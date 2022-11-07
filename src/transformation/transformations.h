#pragma once

#include <array>

#include "../configuration/constants.h"
#include "../data-types/face-values.h"
#include "../data-types/phys-values.h"
#include "../grid/padded-grid.h"
#include "../grid/utils.h"

constexpr double DEFAULT_GAMMA = 1.4;

namespace Transformation {
void reconstToConservatives(PhysValues& output, const FieldStruct& reconstructions,
                            const bool isThermal, const double gamma = DEFAULT_GAMMA);

double computeThermalPressure(const FieldStruct& fields, const bool isThermal,
                              const double gamma = DEFAULT_GAMMA);

void primitiveToConservative(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const bool isThermal,
                             const double gamma = DEFAULT_GAMMA);

void conservativeToPrimitive(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const bool isThermal,
                             const double gamma = DEFAULT_GAMMA);
}; // namespace Transformation

namespace TransformationSycl {

void primitiveToConservative(std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims,
                             const bool isThermal, const double gamma = DEFAULT_GAMMA);

void conservativeToPrimitive(std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims,
                             const bool isThermal, const double gamma = DEFAULT_GAMMA);

} // namespace TransformationSycl
