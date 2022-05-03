#pragma once

#include <array>

#include "../configuration/constants.h"
#include "../field-wrapper/face-values.h"
#include "../field-wrapper/phys-values.h"
#include "../grid/padded-grid.h"

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
