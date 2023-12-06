#pragma once

#include <array>

#include <celerity.h>

enum BoundaryType { EMPTY, EXTRAPOLATE, OUTFLOW, USER };

namespace Boundary {

enum class Axis : int {
    X = 0,
    Y = 1,
    Z = 2,
};

enum class Dir {
    LEFT,
    RIGHT,
};

template <Axis axis> constexpr std::array<int, 3> generateIndexMap() {
    if constexpr (axis == Axis::X) {
        return std::array{ 0, 1, 2 };
    }
    if constexpr (axis == Axis::Y) {
        return std::array{ 1, 0, 2 };
    }
    if constexpr (axis == Axis::Z) {
        return std::array{ 2, 0, 1 };
    }
}

inline auto generateBoundaryRangeMapper(const std::array<int, 3>& idxMap, const int inner,
                                        const std::vector<int> accesses) {
    const auto rangeMapper = [idxMap, inner, accesses](const celerity::chunk<3> chunk) {
        if (static_cast<size_t>(inner) >= chunk.offset[idxMap[0]] &&
            static_cast<size_t>(inner) < chunk.offset[idxMap[0]] + chunk.range[idxMap[0]]) {

            const auto [min, max] = std::minmax_element(accesses.begin(), accesses.end());

            const auto offset = *min;
            const auto range = *max - *min + 1;

            auto subrange = celerity::subrange<3>{ chunk.offset, chunk.range };
            subrange.offset[idxMap[0]] = offset;
            subrange.range[idxMap[0]] = range;
            return subrange;
        }
        // Working plane of work items is not part of this chunk.
        return celerity::subrange<3>{};
    };

    return rangeMapper;
}

} // namespace Boundary
