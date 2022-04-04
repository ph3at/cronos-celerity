#pragma once

#include "../field-wrapper/fields.h"
#include "../grid/padded-grid.h"
#include "../misc/faces.h"
#include "../parameters/constants.h"

namespace Reconstruction {
PerFaceValues reconstruct(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned x,
                          const unsigned y, const unsigned z);
};
