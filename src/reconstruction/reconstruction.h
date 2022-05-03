#pragma once

#include "../configuration/constants.h"
#include "../field-wrapper/fields.h"
#include "../grid/padded-grid.h"
#include "../misc/faces.h"

namespace Reconstruction {
PerFaceValues reconstruct(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned x,
                          const unsigned y, const unsigned z);
};
