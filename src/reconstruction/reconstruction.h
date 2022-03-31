#pragma once

#include "../field-wrapper/fields.h"
#include "../grid/padded-grid.h"
#include "../misc/constants.h"
#include "../misc/faces.h"

namespace Reconstruction {
typedef std::array<FieldStruct, Faces::FaceMax> ReconstValues;

ReconstValues reconstruct(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned x,
                          const unsigned y, const unsigned z);
}; // namespace Reconstruction
