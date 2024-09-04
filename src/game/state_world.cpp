#include "gameutils/tiled_map.h"
#include "geometry/edges_to_polygons.h"
#include "geometry/tiles_to_edges.h"
#include "main.h"

#include "entities.h"
#include "tile_grids/core.h"
#include "tile_grids/debug_rendering.h"
#include "tile_grids/entities.h"
#include "tile_grids/high_level.h"
#include "utils/json.h"
#include "utils/ring_multiarray.h"
#include "box2d_physics/math_adapters.h"
#include "utils/bit_manip.h"

#include <box2cpp/box2c.hpp>
#include <box2cpp/debug_imgui_renderer.hpp>

struct PhysicsWorld : Tickable
{
    IMP_STANDALONE_COMPONENT(Game)

    b2::World w;

    b2::DebugImguiRenderer renderer;

    PhysicsWorld()
    {
        w = adjust(b2::World::Params{}, .gravity.y *= -1);
        renderer.camera_pos.x = 12;
        renderer.camera_pos.y = 8;
        renderer.camera_scale = 28;

        renderer.callbacks.drawMass = false;
    }

    void Tick() override
    {
        w.Step(1/60.f, 4);

        renderer.DrawShapes(w);
        renderer.DrawModeToggles();
        renderer.MouseDrag(w);
    }
};

namespace Tiles
{
    struct SystemTraits : TileGrids::DefaultSystemTraits
    {
        // Allow multi-component tiles. 4 bits per tile.
        using TileComponentDesc = std::uint8_t;
        static constexpr int num_tile_component_bits = 4;
    };

    struct System : TileGrids::System<SystemTraits> {};

    static constexpr int chunk_size = 8;

    enum class Tile
    {
        empty,
        wall,
        wall_slope_a, // |/
        wall_slope_b, // \|
        wall_slope_c, // /|
        wall_slope_d, // |\ .
        wall_slope_ac, // |/ /|
        _count,
    };

    // A chunk stores a grid of those.
    struct Cell
    {
        Tile tiles[4]{};

        std::uint8_t bits{};

        // Whether this tile connects to the one to right or below it.
        [[nodiscard]] auto ConnectsToAdjacentTile(this auto &&self, bool vertical) -> BitManip::BitProperty<decltype(self), std::uint8_t, BitManip::any>
        {
            return {self.bits, std::uint8_t(0b100 << vertical)};
        }
        // Whether the internal quaters of this tile connect with each other.
        // The index is: 0 - +X+Y, 1 = -X+Y, 2 = -X-Y, 3 = +X-Y.
        [[nodiscard]] auto ConnectsInternally(this auto &&self, int i) -> BitManip::BitProperty<decltype(self), std::uint8_t, BitManip::any>
        {
            ASSERT(i >= 0 && i < 4);
            return {self.bits, std::uint8_t(1 << i)};
        }
        // Same bits as in `ConnectsInternally()`, but all 4 in a single number.
        [[nodiscard]] std::uint8_t ConnectsInternallyMask() const
        {
            return bits & 0b1111;
        }

        // Returns a bit mask of emptiness, consisting of 4 bits. (0th bit = +X quarter, 1th bit = +Y quarter, then -X, then -Y).
        [[nodiscard]] std::uint8_t ComputeColliderNonEmptinessMask() const
        {
            return
                (tiles[0] != Tile::empty) * 0b0001 |
                (tiles[1] != Tile::empty) * 0b0010 |
                (tiles[2] != Tile::empty) * 0b0100 |
                (tiles[3] != Tile::empty) * 0b1000;
        }

        // Returns which bits in `bits` are allowed to be set with this `tiles[]` configuration.
        // Note that bits 0 and 1 also require a matching adjacent tile in addition to being set in this mask.
        [[nodiscard]] std::uint8_t ComputeAllowedBitsMask() const
        {
            std::uint8_t coll = ComputeColliderNonEmptinessMask();
            std::uint8_t ret = coll & 0b11;

            coll &= coll >> 1 | coll << 3;
            ret |= coll << 2;
            return ret;
        }
    };

    static Geom::TilesToEdges::BakedTileset baked_tileset = []{
        Geom::TilesToEdges::Tileset tileset;
        tileset.tile_size = ivec2(2,2);
        tileset.vertices = {ivec2(1,1), ivec2(0,0), ivec2(2,0), ivec2(2,2), ivec2(0,2)};
        std::vector<unsigned int> sides[4] = {{0,2,3}, {0,3,4}, {0,4,1}, {0,1,2}};
        tileset.tiles = {
            // Padding those with `{}`s to have nice loop indices. After baking they disappear anyway.
            {                                      },
            {sides[0]                              },
            {{},       sides[1]                    },
            {sides[0], sides[1]                    },
            {{},       {},       sides[2]          },
            {sides[0], {},       sides[2]          },
            {{},       sides[1], sides[2]          },
            {sides[0], sides[1], sides[2]          },
            {{},       {},       {},       sides[3]},
            {sides[0], {},       {},       sides[3]},
            {{},       sides[1], {},       sides[3]},
            {sides[0], sides[1], {},       sides[3]},
            {{},       {},       sides[2], sides[3]},
            {sides[0], {},       sides[2], sides[3]},
            {{},       sides[1], sides[2], sides[3]},
            {sides[0], sides[1], sides[2], sides[3]},
        };
        return Geom::TilesToEdges::BakedTileset(std::move(tileset));
    }();

    struct ChunkComponentCollider
    {
        struct Polygon
        {
            b2Polygon shape_data{};

            b2::Shape shape;
        };
        std::vector<Polygon> polygons;
    };

    // The entity for tile grids.
    struct GridEntity : Meta::with_virtual_destructor<GridEntity>
    {
        IMP_STANDALONE_COMPONENT(Game)

        TileGrids::ChunkGrid<System, chunk_size, Cell, ChunkComponentCollider> grid;
        b2::Body body;

        GridEntity(PhysicsWorld &w, fvec2 pos, float angle)
        {
            b2::Body::Params body_params;
            body_params.type = b2_dynamicBody;
            body_params.position = pos;
            body_params.angle = angle;
            body_params.automaticMass = false; // We update it in `FinalizeGridAfterChange()`.

            body = w.w.CreateBody(b2::OwningHandle, body_params);
        }

        void LoadTiles(Stream::ReadOnlyData data);
    };

    [[nodiscard]] b2::Shape::Params GetShapeParams()
    {
        return {};
    }

    // Traits class for `tile_grids/high_level.h`. We wrap it in `EntityHighLevelTraits` below instead of using directly.
    struct BasicHighLevelTraits
    {
        // This part is for `tile_grids/entities.h`, for the `EntityHighLevelTraits` wrapper.

        // The tag for the entity system.
        using EntityTag = Game;
        // The entity that will store the grid.
        using GridEntity = Tiles::GridEntity;

        // This is called when splitting a grid to copy the basic parameters.
        static GridEntity *CreateSplitGrid(typename EntityTag::Controller& controller, const GridEntity *source_grid)
        {
            return &controller.create<GridEntity>(*controller.get<PhysicsWorld>(), fvec2(source_grid->body.GetPosition()), source_grid->body.GetAngle());
        }
        static void FinishSplitGridInit(EntityTag::Controller& controller, const GridEntity *from, GridEntity *to)
        {
            (void)controller;
            (void)from;
            (void)to;
        }

        // This part is directly for `tile_grids/high_level.h`:

        static void OnUpdateGridChunkContents(EntityTag::Controller& controller, GridEntity *grid, vec2<System::WholeChunkCoord> chunk_coord, const System::TileComponentIndices<chunk_size> &comps_per_tile)
        {
            (void)controller;

            static b2::Shape::Params shape_params = GetShapeParams();

            fvec2 chunk_base_offset = fvec2(chunk_coord) * chunk_size;

            auto *chunk = grid->grid.GetChunk(chunk_coord);

            struct CompEntry
            {
                Geom::EdgesToPolygons::TriangulationInput<int> tri_input;
                decltype(tri_input.InsertionCallback()) callback = tri_input.InsertionCallback();

                CompEntry()
                {
                    tri_input.Reserve(100);
                }
                CompEntry(const CompEntry &) = delete;
                CompEntry &operator=(const CompEntry &) = delete;
            };

            // One triangulation input per component.
            std::vector<CompEntry> tri_inputs(chunk->GetComponents().components.size());

            System::ComponentIndex cur_comp = System::ComponentIndex::invalid;
            System::ComponentIndex next_comp = System::ComponentIndex::invalid;

            Geom::TilesToEdges::ConvertTilesToEdges(
                baked_tileset, Geom::TilesToEdges::Mode::closed, ivec2(chunk_size),
                // Tile input.
                [&](ivec2 pos)
                {
                    // The return values are indices into `baked_tileset`.
                    return chunk->GetChunk().at(pos).ComputeColliderNonEmptinessMask();
                },
                // Point output.
                [&](ivec2 pos, Geom::PointInfo info)
                {
                    tri_inputs[std::to_underlying(cur_comp)].callback(pos, info);

                    if (info.type == Geom::PointType::last)
                    {
                        cur_comp = next_comp;
                        next_comp = System::ComponentIndex::invalid;
                    }
                },
                [&](ivec2 pos, int loop_index, ivec2 offset, int other_loop_index) -> bool
                {
                    return
                        comps_per_tile.at(System::CoordInsideChunkWithTileComp{.pos = pos, .comp = DiscoverCellComponent(chunk->GetChunk().at(pos), loop_index)}) ==
                        comps_per_tile.at(System::CoordInsideChunkWithTileComp{.pos = pos + offset, .comp = DiscoverCellComponent(chunk->GetChunk().at(pos), other_loop_index)});
                },
                // Starting tile output for each edge loop.
                [&](ivec2 pos, int loop_index)
                {
                    System::ComponentIndex comp = comps_per_tile.at(System::CoordInsideChunkWithTileComp(pos, std::uint8_t(1) << loop_index));
                    ASSERT(comp != System::ComponentIndex::invalid, "Why does an edge loop start from an empty cell?");
                    ASSERT(next_comp == System::ComponentIndex::invalid, "Too many 'edge loop begins' events?");
                    if (cur_comp == System::ComponentIndex::invalid)
                        cur_comp = comp; // First edge loop.
                    else
                        next_comp = comp; // Subsequent loops.
                }
            );

            for (std::size_t i = 0; i < chunk->GetComponents().components.size(); i++)
            {
                ChunkComponentCollider &coll = chunk->GetMutComponentExtraData(System::ComponentIndex(i));
                coll.polygons.clear(); // This destroys the existing shapes too.

                b2Hull new_hull{};

                Geom::EdgesToPolygons::ConvertToPolygons<float>(
                    tri_inputs[i].tri_input, b2_maxPolygonVertices,
                    [&](ivec2 pos, Geom::PointInfo info)
                    {
                        if (info.type == Geom::PointType::last)
                        {
                            ChunkComponentCollider::Polygon &new_poly = coll.polygons.emplace_back();
                            new_poly.shape_data = b2MakePolygon(&new_hull, 0);

                            new_poly.shape = grid->body.CreateShape(b2::OwningHandle, shape_params, b2MakePolygon(&new_hull, 0));
                            new_hull.count = 0;
                        }
                        else
                        {
                            assert(new_hull.count < b2_maxPolygonVertices);
                            new_hull.points[new_hull.count++] = b2Vec2(fvec2(pos) / baked_tileset.tile_size + chunk_base_offset);
                        }
                    }
                );
            }
        }

        static void OnPreMoveComponentBetweenChunks(typename EntityTag::Controller& controller, GridEntity *source_grid, System::ComponentCoords coords, GridEntity *target_grid)
        {
            (void)controller;
            ChunkComponentCollider &source_collider = source_grid->grid.GetChunk(coords.chunk_coord)->GetMutComponentExtraData(coords.in_chunk_component);
            b2::Shape::Params params = GetShapeParams();
            for (ChunkComponentCollider::Polygon &polygon : source_collider.polygons)
                polygon.shape = target_grid->body.CreateShape(b2::OwningHandle, params, polygon.shape_data);
        }

        static void FinalizeGridAfterChange(typename EntityTag::Controller& controller, GridEntity *grid)
        {
            (void)controller;
            grid->body.ApplyMassFromShapes();
        }

        // Each tile of a chunk stores this.
        using CellType = Cell;

        // Returns the connectivity mask of a cell in the specified direction, for the purposes of splitting unconnected grids.
        // The bit order should NOT be reversed when flipping direction, it's always the same.
        [[nodiscard]] static System::TileEdgeConnectivity CellConnectivity(const CellType &cell, System::NonzeroTileComponentDesc comp, int dir)
        {
            ASSERT(dir >= 0 && dir < 4);

            if (!(comp & 1 << dir))
                return {}; // Not in this component.

            // Respect connectivity bits for +X,+Y directions.
            if (dir == 0 && !cell.ConnectsToAdjacentTile(0))
                return {};
            if (dir == 1 && !cell.ConnectsToAdjacentTile(1))
                return {};

            return cell.tiles[dir] == Tile::empty ? 1 : 0;
        }

        // Calls `func` for every component in the cell.
        static void ForEachCellComponent(const CellType &cell, auto &&func)
        {
            std::uint8_t sides = cell.ComputeColliderNonEmptinessMask();
            if (sides == 0)
                return; // Not strictly necessary, hopefully an optimization.

            std::uint8_t edges = cell.ConnectsInternallyMask();

            // If `sides` ends with 1, move that block on ones 4 bits to the left.
            // This ensures that all blocks of 1s are continuous.
            std::uint8_t t = (std::uint8_t(1) << std::countr_one(sides)) - 1;
            sides |= t << 4;
            sides &= ~t;
            // Duplicate edges 4 bits to the left, to match the extended `sides`.
            edges |= edges << 4;

            std::uint8_t comp = 0;
            for (t = 1; sides || comp; t <<= 1)
            {
                bool has_side = sides & t;
                if (comp && (!has_side || !(edges & (t >> 1))))
                {
                    comp |= comp >> 4;
                    comp &= 0b1111;
                    func(std::as_const(comp));
                    comp = 0;
                }

                if (has_side)
                {
                    sides &= ~t;
                    comp |= t;
                }
            }
        }

        // Given a direction, returns a cell component touching that direction.
        static System::TileComponentDesc DiscoverCellComponent(const CellType &cell, int dir)
        {
            ASSERT(dir >= 0 && dir < 4);

            System::TileComponentDesc bit = System::TileComponentDesc(0b10001) << dir;

            std::uint8_t sides = cell.ComputeColliderNonEmptinessMask();
            if (!(sides & bit))
                return {};

            std::uint8_t edges = cell.ConnectsInternallyMask();

            sides |= sides << 4;
            edges |= edges << 4;

            std::uint8_t comp = bit;

            if ((sides & bit << 1) && (edges & bit))
            {
                comp |= bit << 1;

                if ((sides & bit << 2) && (edges & bit << 1))
                {
                    comp |= bit << 2;

                    if ((sides & bit << 3) && (edges & bit << 2))
                        comp |= bit << 3;
                }
            }

            if ((sides & bit >> 1) && (edges & bit >> 1))
            {
                comp |= bit >> 1;

                if ((sides & bit >> 2) && (edges & bit >> 2))
                {
                    comp |= bit >> 2;

                    if ((sides & bit >> 3) && (edges & bit >> 3))
                        comp |= bit >> 3;
                }
            }

            return comp & 0b1111;
        }

        // Moves parts of a cell listed in `comp` from `from` to `to`.
        static void MoveCellComponent(System::NonzeroTileComponentDesc comp, CellType &from, CellType &to)
        {
            for (int i = 0; i < 4; i++)
            {
                if (comp & (1 << i))
                {
                    to.tiles[i] = from.tiles[i];
                    from.tiles[i] = {};
                }
            }

            from.bits &= from.ComputeAllowedBitsMask();
            to.bits &= to.ComputeAllowedBitsMask();
        }

        // Returns our data from a grid.
        [[nodiscard]] static TileGrids::ChunkGrid<System, chunk_size, CellType, ChunkComponentCollider> &GridToData(GridEntity *grid)
        {
            return grid->grid;
        }
    };
    using HighLevelTraits = TileGrids::EntityHighLevelTraits<BasicHighLevelTraits>;

    // Stores lists of dirty chunks across entities.
    // We have a separate entity that stores those, see `DirtyListsEntity` below.
    using DirtyLists = TileGrids::DirtyChunkLists<System, HighLevelTraits>;

    // This entity stores dirty chunk lists across other entities.
    struct DirtyListsEntity
    {
        IMP_STANDALONE_COMPONENT(Game)

        // The dirty chunk lists.
        DirtyLists dirty;

        // This is used when updating chunks. We're reusing this across frames to reduce heap allocations/deallocations,
        // but nothing is stopping you from creating it per frame as needed.
        DirtyLists::ReusedUpdateData<chunk_size> reused_update_data;
    };

    void GridEntity::LoadTiles(Stream::ReadOnlyData data)
    {
        Json json(data.string(), 32);
        auto tiles = Tiled::LoadTileLayer(Tiled::FindLayer(json, "mid"));

        grid.LoadFromArray(
            game,
            &game.get<DirtyListsEntity>()->dirty,
            dynamic_cast<Game::Entity &>(*this).id(),
            tiles.size().to<System::GlobalTileCoord>(),
            [&](vec2<System::GlobalTileCoord> pos, Cell &cell)
            {
                int tile = tiles.at(pos);
                if (tile < 0 || tile >= int(Tile::_count))
                    throw std::runtime_error(FMT("Unknown tile at {}", pos));
                for (Tile &part : cell.tiles)
                    part = Tile(tile);
                cell.ConnectsInternally(0) = true;
                cell.ConnectsInternally(1) = true;
                cell.ConnectsInternally(2) = true;
                cell.ConnectsInternally(3) = true;
                cell.ConnectsToAdjacentTile(false) = true;
                // cell.ConnectsToAdjacentTile(true) = true;
            }
        );
    }
}

struct MouseCamera : Camera, Tickable
{
    void Tick() override
    {
        pos = mouse.pos();
    }
};

struct TestEntity : Tickable, Renderable
{
    IMP_STANDALONE_COMPONENT(Game)

    int active_entity_index = -2;

    TestEntity()
    {
        auto &w = game.create<PhysicsWorld>();
        auto &d = game.create<Tiles::DirtyListsEntity>();

        auto &e = game.create<Tiles::GridEntity>(w, fvec2(0,0), 0);
        e.LoadTiles(Program::ExeDir() + "assets/map.json");

        d.dirty.Update(game, d.reused_update_data);

        b2::BodyRef platform = w.w.CreateBody(b2::DestroyWithParent, adjust(b2::Body::Params{}, .position = fvec2(10,20)));
        platform.CreateShape(b2::DestroyWithParent, b2::Shape::Params{}, b2MakeBox(20, 1));
    }

    void Tick() override
    {
        auto &renderer = game.get<PhysicsWorld>()->renderer;

        auto &list = game.get<Game::Category<Ent::OrderedList, Tiles::GridEntity>>();
        if (Input::Button(Input::space).pressed())
        {
            ++active_entity_index;
            if (active_entity_index >= list.size())
                active_entity_index = -2;
        }

        for (int index = 0; const auto &e : list)
        {
            if (index++ != active_entity_index && active_entity_index != -1)
                continue;

            const auto &grid = e.get<Tiles::GridEntity>();
            TileGrids::ImguiDebugDraw(grid.grid, [&](fvec2 pos){return fvec2(renderer.Box2dToImguiPoint(grid.body.GetWorldPoint(pos)));}, *ImGui::GetBackgroundDrawList(), TileGrids::DebugDrawFlags::all);
        }
    }

    void Render() const override
    {

    }
};

namespace States
{
    STRUCT( World EXTENDS StateBase )
    {
        MEMBERS()

        void Init() override
        {
            // Configure the audio.
            float audio_distance = screen_size.x * 3;
            Audio::ListenerPosition(fvec3(0, 0, -audio_distance));
            Audio::ListenerOrientation(fvec3(0,0,1), fvec3(0,-1,0));
            Audio::Source::DefaultRefDistance(audio_distance);

            // Entities.
            game = nullptr;

            game.create<MouseCamera>();
            game.create<TestEntity>();
        }

        void Tick(std::string &next_state) override
        {
            (void)next_state;

            for (auto &e : game.get<AllTickable>())
                e.get<Tickable>().Tick();
        }

        void Render() const override
        {
            Graphics::SetClearColor(fvec3(0));
            Graphics::Clear();

            r.BindShader();

            for (auto &e : game.get<AllRenderable>())
                e.get<Renderable>().Render();

            r.Finish();
        }
    };
}
