#pragma once

#include <iostream>

#include "../configuration/constants.h"
#include "../configuration/problem.h"
#include "../configuration/shock-tube.h"
#include "../data-types/faces.h"
#include "../data-types/phys-fields.h"
#include "../grid/padded-grid.h"
#include "boundary-types.h"
#include "extrapolate.h"
#include "outflow.h"

namespace Boundary {
template <class T>
void applyAll(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const Problem<T>& problem);

template <class T>
void applyField(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const Problem<T>& problem,
                const unsigned field);
}; // namespace Boundary

template <class T>
void Boundary::applyAll(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const Problem<T>& problem) {
    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
        applyField(grid, problem, field);
    }
}

template <class T>
void Boundary::applyField(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const Problem<T>& problem,
                          const unsigned field) {
    for (unsigned face = 0; face < Faces::FaceMax; face++) {
        /* GPUs only support compile time polymorphism, so we are using this Frankenstein here.
         * Sorry. */
        if (grid.boundaryTypes[face] == EMPTY) {
            /* Used for applying boundaries at a different point (e.g. AMR) or for keeping the
             * boundary at its initial values. */
        } else if (grid.boundaryTypes[face] == EXTRAPOLATE) {
            Extrapolate::apply(grid, field, face);
        } else if (grid.boundaryTypes[face] == OUTFLOW) {
            Outflow::apply(grid, field, face);
        } else if (grid.boundaryTypes[face] == USER) {
            problem.applyBoundary(grid, field, face);
        } else {
            std::cerr << "Unknown boundary condition type encountered : "
                      << grid.boundaryTypes[face] << std::endl;
        }
    }
}
