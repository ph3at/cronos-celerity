#include <iostream>

#include "boundary.h"

#include "../misc/faces.h"
#include "boundary-types.h"
#include "constant.h"
#include "extrapolate.h"
#include "outflow.h"

namespace Boundary {
void applyAll(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const Problem& problem) {
    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
        applyField(grid, problem, field);
    }
}

void applyField(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const Problem& problem,
                const unsigned field) {
    for (unsigned face = 0; face < Faces::FaceMax; face++) {
        /* GPUs only support compile time polymorphism, so we are using this Frankenstein here.
         * Sorry. */
        if (problem.boundaryTypes[face] == EXTRAPOLATE) {
            Extrapolate::apply(grid, field, face);
        } else if (problem.boundaryTypes[face] == OUTFLOW) {
            Outflow::apply(grid, field, face);
        } else if (problem.boundaryTypes[face] == CONSTANT) {
            Constant::apply(grid, problem, field, face);
        } else {
            std::cerr << "Unknown boundary condition type encountered : "
                      << problem.boundaryTypes[face] << std::endl;
        }
    }
}
}; // namespace Boundary
