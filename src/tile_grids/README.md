The tile system consists of multiple layers. You can use more or less layers depending on what you need.

* **`tile_grids/core.h`** is the lowest level.

  It contains `TileGrids::System<>` with all the raw algorithms.

  You can customize some types they use by passing your class as template parameter, inheriting from `TileGrids::DefaultSystemTraits`.

* **`tile_grids/high_level.h`** is based on `tile_grids/core.h`.

  This contains `TileGrids::ChunkGrid`: an implementation of a chunk grid, with its own `...::Chunk` type. The grid contains a 2d ring array of `std::unique_ptr`s to chunks, where each chunk is optional and is automatically removed when it becomes empty, and the whole array of chunks is automatically grown and shrinked as needed.

  It also contains `TileGrids::DirtyChunkLists`, which stores lists of dirty chunks across multiple grids.

  Classes in this file can and must be configured by providing a `HighLevelTraits` template parameter, with must include a bunch of custom functions.

* **`tile_grids/entities.h`** is based on `tile_grids/high_level.h`.

  It includes `TileGrids::Entities`: entities for chunk grids and for a list of dirty chunks across entities.

  This file provides `TileGrids::EntityHighLevelTraits<...>`, which is a wrapper implementing some of the functions of `HighLevelTraits` in terms of some other functions you must provide.

Other files:

* **`tile_grids/debug_rendering.h`** implements debug rendering for chunk grids from `tile_grids/high_level.h`, using ImGui.
