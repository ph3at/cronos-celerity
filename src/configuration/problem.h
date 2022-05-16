#pragma once

#include <array>

#include "../boundary/boundary-types.h"
#include "../data-types/direction.h"
#include "../data-types/faces.h"
#include "../data-types/phys-fields.h"
#include "../grid/padded-grid.h"
#include "constants.h"

template <class Specific> class Problem {
  public:
    Problem(const double cflThreshold, const bool thermal, const double timeDelta,
            const double timeStart, const double timeEnd, const double gamma,
            const std::array<BoundaryType, Faces::FaceMax> boundaryTypes);

    const double cflThreshold;
    const bool thermal;
    const double timeDelta;
    const double timeStart;
    const double timeEnd;
    const double gamma;

    const std::array<BoundaryType, Faces::FaceMax> boundaryTypes;

    void initialiseGrid(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const;
    void applyBoundary(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned field,
                       const unsigned face) const;
    void applySource(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const;
};

template <class Specific>
Problem<Specific>::Problem(const double cflThreshold, const bool thermal, const double timeDelta,
                           const double timeStart, const double timeEnd, const double gamma,
                           const std::array<BoundaryType, Faces::FaceMax> boundaryTypes)
    : cflThreshold(cflThreshold), thermal(thermal), timeDelta(timeDelta), timeStart(timeStart),
      timeEnd(timeEnd), gamma(gamma), boundaryTypes(boundaryTypes) {}

template <class Specific>
void Problem<Specific>::initialiseGrid(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const {
    static_cast<const Specific*>(this)->initialiseGrid(grid);
}

template <class Specific>
void Problem<Specific>::applyBoundary(PaddedGrid<FieldStruct, GHOST_CELLS>& grid,
                                      const unsigned field, const unsigned face) const {
    static_cast<const Specific*>(this)->applyBoundary(grid, field, face);
}

template <class Specific>
void Problem<Specific>::applySource(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const {
    static_cast<const Specific*>(this)->applySource(grid);
}
