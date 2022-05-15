#pragma once

#include "../configuration/constants.h"
#include "../data-types/faces.h"
#include "../data-types/fields.h"
#include "../grid/padded-grid.h"

namespace Reconstruction {
PerFaceValues reconstruct(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned x,
                          const unsigned y, const unsigned z);
};
