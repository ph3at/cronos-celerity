#pragma once

#include <array>

#include <celerity.h>
#include <sycl/sycl.hpp>

#include "../data-types/faces.h"
#include "../data-types/phys-fields.h"
#include "../grid/padded-grid.h"
#include "boundary-types.h"

namespace Extrapolate {
template <unsigned padding>
void apply(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const unsigned face);

template <unsigned padding>
void applyX(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft);
template <unsigned padding>
void applyY(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft);
template <unsigned padding>
void applyZ(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft);
}; // namespace Extrapolate

template <unsigned padding>
void Extrapolate::applyX(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft) {
    int xInner, xOuter;
    int direction;
    if (isLeft) {
        xInner = grid.xStart() - 1;
        xOuter = -1;
        direction = -1;
    } else {
        xInner = grid.xEnd();
        xOuter = grid.xDim();
        direction = 1;
    }
    for (int x = xInner; x != xOuter; x += direction) {
        for (unsigned y = 0; y < grid.yDim(); y++) {
            for (unsigned z = 0; z < grid.zDim(); z++) {
                grid(x, y, z)[field] = grid(x - direction, y, z)[field];
            }
        }
    }
}

template <unsigned padding>
void Extrapolate::applyY(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft) {
    int yInner, yOuter;
    int direction;
    if (isLeft) {
        yInner = grid.yStart() - 1;
        yOuter = -1;
        direction = -1;
    } else {
        yInner = grid.yEnd();
        yOuter = grid.yDim();
        direction = 1;
    }
    for (unsigned x = 0; x < grid.xDim(); x++) {
        for (int y = yInner; y != yOuter; y += direction) {
            for (unsigned z = 0; z < grid.zDim(); z++) {
                grid(x, y, z)[field] = grid(x, y - direction, z)[field];
            }
        }
    }
}

template <unsigned padding>
void Extrapolate::applyZ(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft) {
    int zInner, zOuter;
    int direction;
    if (isLeft) {
        zInner = grid.zStart() - 1;
        zOuter = -1;
        direction = -1;
    } else {
        zInner = grid.zEnd();
        zOuter = grid.zDim();
        direction = 1;
    }
    for (unsigned x = 0; x < grid.xDim(); x++) {
        for (unsigned y = 0; y < grid.yDim(); y++) {
            for (int z = zInner; z != zOuter; z += direction) {
                grid(x, y, z)[field] = grid(x, y, z - direction)[field];
            }
        }
    }
}

template <unsigned padding>
void Extrapolate::apply(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const unsigned face) {
    if (face == Faces::FaceWest) {
        applyX(grid, field, true);
    } else if (face == Faces::FaceEast) {
        applyX(grid, field, false);
    } else if (face == Faces::FaceSouth) {
        applyY(grid, field, true);
    } else if (face == Faces::FaceNorth) {
        applyY(grid, field, false);
    } else if (face == Faces::FaceBottom) {
        applyZ(grid, field, true);
    } else {
        applyZ(grid, field, false);
    }
}

namespace ExtrapolateSycl {

using Boundary::Axis;
using Boundary::Dir;

template <Axis axis, Dir dir, unsigned padding>
void apply(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid, const unsigned field) {
    constexpr auto idxMap = []() {
        if constexpr (axis == Axis::X) {
            return std::array{ 0, 1, 2 };
        }
        if constexpr (axis == Axis::Y) {
            return std::array{ 1, 0, 2 };
        }
        if constexpr (axis == Axis::Z) {
            return std::array{ 2, 0, 1 };
        }
    }();

    auto inner = 0;
    auto outer = 0;
    auto direction = 0;

    if constexpr (dir == Dir::LEFT) {
        inner = padding - 1;
        outer = -1;
        direction = -1;
    } else {
        inner = grid.get_range()[idxMap[0]] - padding;
        outer = grid.get_range()[idxMap[0]];
        direction = 1;
    }

    queue.submit([&](sycl::handler& cgh) {
        auto gridAccessor = grid.template get_access<sycl::access::mode::read_write>(cgh);

        const auto range = sycl::range<2>(grid.get_range()[idxMap[1]], grid.get_range()[idxMap[2]]);

        cgh.parallel_for(range, [=](const sycl::id<2> id) {
            for (auto d = inner; d != outer; d += direction) {
                auto dstIdx = sycl::id<3>(0, 0, 0);
                dstIdx[idxMap[0]] = d;
                dstIdx[idxMap[1]] = id[0];
                dstIdx[idxMap[2]] = id[1];

                auto srcIdx = dstIdx;
                srcIdx[idxMap[0]] -= direction;

                gridAccessor[dstIdx][field] = gridAccessor[srcIdx][field];
            }
        });
    });
}

template <unsigned padding>
void apply(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid, const unsigned field, const unsigned face) {
    if (face == Faces::FaceWest) {
        apply<Axis::X, Dir::LEFT, padding>(queue, grid, field);
    } else if (face == Faces::FaceEast) {
        apply<Axis::X, Dir::RIGHT, padding>(queue, grid, field);
    } else if (face == Faces::FaceSouth) {
        apply<Axis::Y, Dir::LEFT, padding>(queue, grid, field);
    } else if (face == Faces::FaceNorth) {
        apply<Axis::Y, Dir::RIGHT, padding>(queue, grid, field);
    } else if (face == Faces::FaceBottom) {
        apply<Axis::Z, Dir::LEFT, padding>(queue, grid, field);
    } else {
        apply<Axis::Z, Dir::RIGHT, padding>(queue, grid, field);
    }
}

} // namespace ExtrapolateSycl

namespace ExtrapolateCelerity {

using Boundary::Axis;
using Boundary::Dir;

template <Axis axis, Dir dir, unsigned padding>
void apply(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid, const unsigned field) {
    constexpr auto idxMap = []() {
        if constexpr (axis == Axis::X) {
            return std::array{ 0, 1, 2 };
        }
        if constexpr (axis == Axis::Y) {
            return std::array{ 1, 0, 2 };
        }
        if constexpr (axis == Axis::Z) {
            return std::array{ 2, 0, 1 };
        }
    }();

    auto inner = 0;
    auto outer = 0;
    auto direction = 0;

    if constexpr (dir == Dir::LEFT) {
        inner = padding - 1;
        outer = -1;
        direction = -1;
    } else {
        inner = grid.get_range()[idxMap[0]] - padding;
        outer = grid.get_range()[idxMap[0]];
        direction = 1;
    }

    queue.submit([&grid, &field, &inner, &outer, &direction, &idxMap](celerity::handler& cgh) {
        const auto boundaryMapper = [=](const celerity::chunk<2> chunk) {
            auto offset = celerity::id<3>{ 0, 0, 0 };
            offset[idxMap[0]] = inner - 1;
            offset[idxMap[1]] = chunk.offset[0];
            offset[idxMap[2]] = chunk.offset[1];
            auto range = celerity::range<3>{ 0, 0, 0 };
            range[idxMap[0]] = padding + 1;
            range[idxMap[1]] = chunk.range[0];
            range[idxMap[2]] = chunk.range[1];
            const auto subrange = celerity::subrange<3>{ offset, range };
            return subrange;
        };
        auto gridAccessor = celerity::accessor{ grid, cgh, boundaryMapper, celerity::read_write };

        const auto range = celerity::range<2>{ grid.get_range()[idxMap[1]], grid.get_range()[idxMap[2]] };

        celerity::debug::set_task_name(cgh, "extrapolate");

        cgh.parallel_for(range, [=](const celerity::id<2> id) {
            for (auto d = inner; d != outer; d += direction) {
                auto dstIdx = celerity::id<3>{ 0, 0, 0 };
                dstIdx[idxMap[0]] = d;
                dstIdx[idxMap[1]] = id[0];
                dstIdx[idxMap[2]] = id[1];

                auto srcIdx = dstIdx;
                srcIdx[idxMap[0]] -= direction;

                gridAccessor[dstIdx][field] = gridAccessor[srcIdx][field];
            }
        });
    });
}

template <unsigned padding>
void apply(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid, const unsigned field,
           const unsigned face) {
    if (face == Faces::FaceWest) {
        apply<Axis::X, Dir::LEFT, padding>(queue, grid, field);
    } else if (face == Faces::FaceEast) {
        apply<Axis::X, Dir::RIGHT, padding>(queue, grid, field);
    } else if (face == Faces::FaceSouth) {
        apply<Axis::Y, Dir::LEFT, padding>(queue, grid, field);
    } else if (face == Faces::FaceNorth) {
        apply<Axis::Y, Dir::RIGHT, padding>(queue, grid, field);
    } else if (face == Faces::FaceBottom) {
        apply<Axis::Z, Dir::LEFT, padding>(queue, grid, field);
    } else {
        apply<Axis::Z, Dir::RIGHT, padding>(queue, grid, field);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <Axis axis, Dir dir, unsigned padding>
void apply3D(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid) {
    constexpr auto idxMap = []() {
        if constexpr (axis == Axis::X) {
            return std::array{ 0, 1, 2 };
        }
        if constexpr (axis == Axis::Y) {
            return std::array{ 1, 0, 2 };
        }
        if constexpr (axis == Axis::Z) {
            return std::array{ 2, 0, 1 };
        }
    }();

    auto inner = 0;
    auto outer = 0;
    auto direction = 0;

    if constexpr (dir == Dir::LEFT) {
        inner = padding - 1;
        outer = -1;
        direction = -1;
    } else {
        inner = grid.get_range()[idxMap[0]] - padding;
        outer = grid.get_range()[idxMap[0]];
        direction = 1;
    }

    queue.submit([&grid, &inner, &outer, &direction, &idxMap](celerity::handler& cgh) {
        const auto first = inner;
        const auto last = outer - direction;

        const auto readMapper = [=](const celerity::chunk<3> chunk) {
            if (static_cast<size_t>(inner) >= chunk.offset[idxMap[0]] &&
                static_cast<size_t>(inner) < chunk.offset[idxMap[0]] + chunk.range[idxMap[0]]) {

                const auto firstAccess = first - direction;
                const auto lastAccess = last - direction;

                const auto min = std::min({ firstAccess, lastAccess });
                const auto max = std::max({ firstAccess, lastAccess });

                const auto offset = min;
                const auto range = max - min + 1;

                auto subrange = celerity::subrange<3>{ chunk.offset, chunk.range };
                subrange.offset[idxMap[0]] = offset;
                subrange.range[idxMap[0]] = range;
                return subrange;
            }
            // Working plane of work items is not part of this chunk.
            return celerity::subrange<3>{};
        };

        auto gridRead = celerity::accessor{ grid, cgh, readMapper, celerity::read_only };

        const auto writeMapper = [=](const celerity::chunk<3> chunk) {
            if (static_cast<size_t>(inner) >= chunk.offset[idxMap[0]] &&
                static_cast<size_t>(inner) < chunk.offset[idxMap[0]] + chunk.range[idxMap[0]]) {

                const auto min = std::min({ first, last });
                const auto max = std::max({ first, last });

                const auto offset = min;
                const auto range = max - min + 1;

                auto subrange = celerity::subrange<3>{ chunk.offset, chunk.range };
                subrange.offset[idxMap[0]] = offset;
                subrange.range[idxMap[0]] = range;
                return subrange;
            }
            // Working plane of work items is not part of this chunk.
            return celerity::subrange<3>{};
        };

        auto gridWrite = celerity::accessor{ grid, cgh, writeMapper, celerity::write_only };

        celerity::debug::set_task_name(cgh, "extrapolate3D");

        cgh.parallel_for(grid.get_range(), [=](const celerity::id<3> id) {
            if (id[idxMap[0]] == static_cast<std::size_t>(inner)) {
                for (auto d = inner; d != outer; d += direction) {
                    auto dstIdx = celerity::id<3>{ 0, 0, 0 };
                    dstIdx[idxMap[0]] = d;
                    dstIdx[idxMap[1]] = id[idxMap[1]];
                    dstIdx[idxMap[2]] = id[idxMap[2]];

                    auto srcIdx = dstIdx;
                    srcIdx[idxMap[0]] -= direction;

                    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                        gridWrite[dstIdx][field] = gridRead[srcIdx][field];
                    }
                }
            }
        });
    });
}

template <unsigned padding>
void apply3D(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid, const unsigned face) {
    if (face == Faces::FaceWest) {
        apply3D<Axis::X, Dir::LEFT, padding>(queue, grid);
    } else if (face == Faces::FaceEast) {
        apply3D<Axis::X, Dir::RIGHT, padding>(queue, grid);
    } else if (face == Faces::FaceSouth) {
        apply3D<Axis::Y, Dir::LEFT, padding>(queue, grid);
    } else if (face == Faces::FaceNorth) {
        apply3D<Axis::Y, Dir::RIGHT, padding>(queue, grid);
    } else if (face == Faces::FaceBottom) {
        apply3D<Axis::Z, Dir::LEFT, padding>(queue, grid);
    } else {
        apply3D<Axis::Z, Dir::RIGHT, padding>(queue, grid);
    }
}

} // namespace ExtrapolateCelerity
