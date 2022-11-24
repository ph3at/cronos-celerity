#pragma once

#include <array>

#include <CL/sycl.hpp>

#include "../configuration/constants.h"
#include "../data-types/face-values.h"
#include "../data-types/phys-values.h"
#include "../grid/padded-grid.h"
#include "../grid/utils.h"

constexpr double DEFAULT_GAMMA = 1.4;

namespace Transformation {
void reconstToConservatives(PhysValues& output, const FieldStruct& reconstructions, const bool isThermal,
                            const double gamma = DEFAULT_GAMMA);

double computeThermalPressure(const FieldStruct& fields, const bool isThermal, const double gamma = DEFAULT_GAMMA);

void primitiveToConservative(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const bool isThermal,
                             const double gamma = DEFAULT_GAMMA);

void conservativeToPrimitive(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const bool isThermal,
                             const double gamma = DEFAULT_GAMMA);
}; // namespace Transformation

namespace TransformationSycl {

void primitiveToConservative(cl::sycl::queue& queue, cl::sycl::buffer<FieldStruct, 3>& grid, const bool isThermal,
                             const double gamma = DEFAULT_GAMMA);

void conservativeToPrimitive(cl::sycl::queue& queue, cl::sycl::buffer<FieldStruct, 3>& grid, const bool isThermal,
                             const double gamma = DEFAULT_GAMMA);

} // namespace TransformationSycl
