#pragma once

#include <array>
#include <cstddef>

namespace grid::utils {

using dimensions = std::array<std::size_t, 3>;

static inline std::size_t idx3d(const std::size_t x, const std::size_t y, const std::size_t z,
                                const dimensions& dims) {
    return x * dims[1] * dims[2] + y * dims[2] + z;
}

} // namespace grid::utils
