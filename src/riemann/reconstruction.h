#pragma once

#include <vector>

#include "../configuration/constants.h"
#include "../data-types/faces.h"
#include "../data-types/fields.h"
#include "../grid/padded-grid.h"
#include "../grid/utils.h"

namespace Reconstruction {
PerFaceValues reconstruct(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned x, const unsigned y,
                          const unsigned z);
};

namespace ReconstructionSycl {
PerFaceValues reconstruct(const std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims, const unsigned x,
                          const unsigned y, const unsigned z);
} // namespace ReconstructionSycl
