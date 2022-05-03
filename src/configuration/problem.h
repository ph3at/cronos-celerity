#pragma once

#include <array>

#include "../boundary/boundary-types.h"
#include "../field-wrapper/phys-fields.h"
#include "../grid/padded-grid.h"
#include "../misc/direction.h"
#include "../misc/faces.h"
#include "constants.h"

template <class Specific> class Problem {
  public:
    Problem(const bool thermal, const double timeDelta, const double timeStart,
            const double timeEnd, const double gamma,
            const std::array<std::size_t, Direction::DirMax> numberCells,
            const std::array<double, Direction::DirMax> cellSize,
            const std::array<BoundaryType, Faces::FaceMax> boundaryTypes);

    const bool thermal;
    const double timeDelta;
    const double timeStart;
    const double timeEnd;
    const double gamma;
    const std::array<std::size_t, Direction::DirMax> numberCells;
    const std::array<double, Direction::DirMax> cellSize;
    const std::array<double, Direction::DirMax> inverseCellSize;

    const std::array<BoundaryType, Faces::FaceMax> boundaryTypes;

    void initialiseGrid(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const;
    void applyBoundary(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned field,
                       const unsigned face) const;
    void applySource(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const;
};

template <class Specific>
Problem<Specific>::Problem(const bool thermal, const double timeDelta, const double timeStart,
                           const double timeEnd, const double gamma,
                           const std::array<std::size_t, Direction::DirMax> numberCells,
                           const std::array<double, Direction::DirMax> cellSize,
                           const std::array<BoundaryType, Faces::FaceMax> boundaryTypes)
    : thermal(thermal), timeDelta(timeDelta), timeStart(timeStart), timeEnd(timeEnd), gamma(gamma),
      numberCells(numberCells), cellSize(cellSize),
      inverseCellSize({ 1.0 / cellSize[0], 1.0 / cellSize[1], 1.0 / cellSize[2] }),
      boundaryTypes(boundaryTypes) {}

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
