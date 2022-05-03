#include "padded-grid.h"

#include "../field-wrapper/phys-fields.h"
#include "../parameters/constants.h"
#include "../parameters/problem.h"

namespace GridFunctions {
bool checkNaN(PaddedGrid<FieldStruct, GHOST_CELLS>& grid);
}; // namespace GridFunctions
