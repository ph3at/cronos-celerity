#pragma once

#include <iostream>
#include <vector>

#include <celerity.h>
#include <sycl/sycl.hpp>

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
template <class T, unsigned padding>
void applyAll(PaddedGrid<FieldStruct, padding>& grid, const Problem<T, FieldStruct, padding>& problem);

template <class T, unsigned padding>
void applyField(PaddedGrid<FieldStruct, padding>& grid, const Problem<T, FieldStruct, padding>& problem,
                const unsigned field);
}; // namespace Boundary

template <class T, unsigned padding>
void Boundary::applyAll(PaddedGrid<FieldStruct, padding>& grid, const Problem<T, FieldStruct, padding>& problem) {
    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
        applyField(grid, problem, field);
    }
}

template <class T, unsigned padding>
void Boundary::applyField(PaddedGrid<FieldStruct, padding>& grid, const Problem<T, FieldStruct, padding>& problem,
                          const unsigned field) {
    for (unsigned face = 0; face < Faces::FaceMax; face++) {
        /* GPUs only support compile time polymorphism, so we are using this Frankenstein here.
         * Sorry. */
        if (problem.boundaryTypes[face] == EMPTY) {
            /* Used for applying boundaries at a different point (e.g. AMR) or for keeping the
             * boundary at its initial values. */
        } else if (problem.boundaryTypes[face] == EXTRAPOLATE) {
            Extrapolate::apply(grid, field, face);
        } else if (problem.boundaryTypes[face] == OUTFLOW) {
            Outflow::apply(grid, field, face);
        } else if (problem.boundaryTypes[face] == USER) {
            problem.applyBoundary(grid, field, face);
        } else {
            std::cerr << "Unknown boundary condition type encountered : " << problem.boundaryTypes[face] << std::endl;
        }
    }
}

namespace BoundarySycl {

template <class T, unsigned padding>
void applyField(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid, const Problem<T, FieldStruct, padding>& problem,
                const unsigned field) {
    for (unsigned face = 0; face < Faces::FaceMax; face++) {
        /* GPUs only support compile time polymorphism, so we are using this Frankenstein here.
         * Sorry. */
        if (problem.boundaryTypes[face] == EMPTY) {
            /* Used for applying boundaries at a different point (e.g. AMR) or for keeping the
             * boundary at its initial values. */
        } else if (problem.boundaryTypes[face] == EXTRAPOLATE) {
            ExtrapolateSycl::apply<padding>(queue, grid, field, face);
        } else if (problem.boundaryTypes[face] == OUTFLOW) {
            OutflowSycl::apply<padding>(queue, grid, field, face);
        } else if (problem.boundaryTypes[face] == USER) {
            problem.applyBoundarySycl(queue, grid, field, face);
        } else {
            std::cerr << "Unknown boundary condition type encountered : " << problem.boundaryTypes[face] << std::endl;
        }
    }
}

template <class T, unsigned padding>
void applyAll(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid, const Problem<T, FieldStruct, padding>& problem) {
    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
        applyField<T, padding>(queue, grid, problem, field);
    }
}

} // namespace BoundarySycl

namespace BoundaryCelerity {

template <class T, unsigned padding>
void applyField(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid,
                const Problem<T, FieldStruct, padding>& problem, const unsigned field) {
    for (unsigned face = 0; face < Faces::FaceMax; face++) {
        /* GPUs only support compile time polymorphism, so we are using this Frankenstein here.
         * Sorry. */
        if (problem.boundaryTypes[face] == EMPTY) {
            /* Used for applying boundaries at a different point (e.g. AMR) or for keeping the
             * boundary at its initial values. */
        } else if (problem.boundaryTypes[face] == EXTRAPOLATE) {
            ExtrapolateCelerity::apply<padding>(queue, grid, field, face);
        } else if (problem.boundaryTypes[face] == OUTFLOW) {
            OutflowCelerity::apply<padding>(queue, grid, field, face);
        } else if (problem.boundaryTypes[face] == USER) {
            problem.applyBoundaryCelerity(queue, grid, field, face);
        } else {
            std::cerr << "Unknown boundary condition type encountered : " << problem.boundaryTypes[face] << std::endl;
        }
    }
}

template <class T, unsigned padding>
void applyAll(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid,
              const Problem<T, FieldStruct, padding>& problem) {
    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
        applyField<T, padding>(queue, grid, problem, field);
    }
}

template <class T, unsigned padding>
void applyAll3D(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid,
                const Problem<T, FieldStruct, padding>& problem) {
    for (unsigned face = 0; face < Faces::FaceMax; ++face) {
        if (problem.boundaryTypes[face] == EMPTY) {
            // Used for applying boundaries at a different point (e.g. AMR) or for keeping the
            // boundary at its initial values.
        } else if (problem.boundaryTypes[face] == EXTRAPOLATE) {
            ExtrapolateCelerity::apply3D<padding>(queue, grid, face);
        } else if (problem.boundaryTypes[face] == OUTFLOW) {
            OutflowCelerity::apply3D<padding>(queue, grid, face);
        } else if (problem.boundaryTypes[face] == USER) {
            ProblemCelerity::apply3D<padding>(queue, grid, face, problem);
        }
    }
}

} // namespace BoundaryCelerity
