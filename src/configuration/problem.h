#pragma once

#include <array>

#include "../boundary/boundary-types.h"
#include "../data-types/direction.h"
#include "../data-types/faces.h"
#include "../data-types/phys-fields.h"
#include "../grid/padded-grid.h"
#include "constants.h"

template <class Specific, class Fields, unsigned padding> class Problem {
  public:
    Problem(const double cflThreshold, const bool thermal, const double timeDelta,
            const double timeStart, const double timeEnd, const double gamma);

    const double cflThreshold;
    const bool thermal;
    const double timeDelta;
    const double timeStart;
    const double timeEnd;
    const double gamma;

    void initialiseGrid(PaddedGrid<Fields, padding>& grid) const;
    void applyBoundary(PaddedGrid<Fields, padding>& grid, const unsigned field,
                       const unsigned face) const;
    void applySource(PaddedGrid<Fields, padding>& grid) const;
};

template <class Specific, class Fields, unsigned padding>
Problem<Specific, Fields, padding>::Problem(const double cflThreshold, const bool thermal,
                                            const double timeDelta, const double timeStart,
                                            const double timeEnd, const double gamma)
    : cflThreshold(cflThreshold), thermal(thermal), timeDelta(timeDelta), timeStart(timeStart),
      timeEnd(timeEnd), gamma(gamma) {}

template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::initialiseGrid(PaddedGrid<Fields, padding>& grid) const {
    static_cast<const Specific*>(this)->initialiseGrid(grid);
}

template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::applyBoundary(PaddedGrid<Fields, padding>& grid,
                                                       const unsigned field,
                                                       const unsigned face) const {
    static_cast<const Specific*>(this)->applyBoundary(grid, field, face);
}

template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::applySource(PaddedGrid<Fields, padding>& grid) const {
    static_cast<const Specific*>(this)->applySource(grid);
}
