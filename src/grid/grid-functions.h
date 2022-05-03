#include "padded-grid.h"

#include "../configuration/constants.h"
#include "../field-wrapper/phys-fields.h"

namespace GridFunctions {
bool checkNaN(PaddedGrid<FieldStruct, GHOST_CELLS>& grid);
}; // namespace GridFunctions
