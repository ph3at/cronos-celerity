#pragma once

#include "../field-wrapper/phys-fields.h"
#include "../grid/padded-grid.h"
#include "../parameters/constants.h"
#include "../parameters/problem.h"

namespace Boundary {
void applyAll(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const Problem& problem);
void applyField(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const Problem& problem,
                const unsigned field);
}; // namespace Boundary
