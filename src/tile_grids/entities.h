#pragma once

namespace TileGrids
{
    // Those traits extend `core.h` traits with some members needed by `high_level.h`, but not all of them.
    // `BaseTraits` must contain following:
    //
    //     // The tag for the entity system.
    //     using EntityTag = ...;
    //     // The entity that will store the grid.
    //     using GridEntity = ...;
    //
    // Also `BaseTraits` (or the class you end up deriving from this one) must implement the rest of funcions needed for `HighLevelTraits`.
    template <typename BaseTraits>
    struct EntityHighLevelTraits : BaseTraits
    {
        using WorldRef = BaseTraits::EntityTag::Controller &;
        using GridHandle = BaseTraits::EntityTag::Id;
        using GridRef = BaseTraits::GridEntity *; // This must be nullable.

        [[nodiscard]] static GridRef HandleToGrid(WorldRef world, GridHandle grid)
        {
            auto ptr = world.get_opt(grid);
            return ptr ? &ptr->template get<typename BaseTraits::GridEntity>() : nullptr;
        }

        // Removes a grid from the world, when it was shrinked into nothing. Having both `grid_handle` and `grid_ref` is redundant.
        static void DestroyGrid(WorldRef world, GridHandle grid_handle, GridRef grid_ref)
        {
            (void)grid_handle;
            world.destroy(*grid_ref);
        }
    };
}
