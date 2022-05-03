#pragma once

#include <array>

#include "../field-wrapper/face-values.h"
#include "../field-wrapper/phys-values.h"
#include "../grid/padded-grid.h"
#include "../parameters/constants.h"
#include "../parameters/problem.h"

namespace Transformation {
void reconstToConservatives(PhysValues& output, const FieldStruct& reconstructions,
                            const Problem& problem);

double computeThermalPressure(const FieldStruct& fields, const Problem& problem);

void primitiveToConservative(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const Problem& problem);

void conservativeToPrimitive(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const Problem& problem);
}; // namespace Transformation
