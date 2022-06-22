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
            const double timeStart, const double timeEnd, const bool preciseEnd, const double gamma,
            const std::array<double, Direction::DirMax> posLeft,
            const std::array<double, Direction::DirMax> posRight,
            const std::array<std::size_t, Direction::DirMax> numberCells,
            const std::array<BoundaryType, Faces::FaceMax> boundaryTypes);
    Problem(const Problem<Specific, Fields, padding>& blueprint);

    double cflThreshold;
    bool thermal;
    double timeDelta;
    double timeStart;
    double timeEnd;
    bool preciseEnd;
    double gamma;

    std::array<double, Direction::DirMax> posLeft;
    std::array<double, Direction::DirMax> posRight;
    std::array<std::size_t, Direction::DirMax> numberCells;
    std::array<double, Direction::DirMax> cellSize;
    std::array<double, Direction::DirMax> inverseCellSize;

    std::array<BoundaryType, Faces::FaceMax> boundaryTypes;

    void initialiseGrid(PaddedGrid<Fields, padding>& grid) const;
    void applyBoundary(PaddedGrid<Fields, padding>& grid, const unsigned field,
                       const unsigned face) const;
    void applySource(PaddedGrid<Fields, padding>& grid) const;
};

template <class Specific, class Fields, unsigned padding>
Problem<Specific, Fields, padding>::Problem(
    const double cflThreshold, const bool thermal, const double timeDelta, const double timeStart,
    const double timeEnd, const bool preciseEnd, const double gamma,
    const std::array<double, Direction::DirMax> posLeft,
    const std::array<double, Direction::DirMax> posRight,
    const std::array<std::size_t, Direction::DirMax> numberCells,
    const std::array<BoundaryType, Faces::FaceMax> boundaryTypes)
    : cflThreshold(cflThreshold), thermal(thermal), timeDelta(timeDelta), timeStart(timeStart),
      timeEnd(timeEnd), preciseEnd(preciseEnd), gamma(gamma), posLeft(posLeft), posRight(posRight),
      numberCells(numberCells), cellSize({ (posRight[Direction::DirX] - posLeft[Direction::DirX]) /
                                               static_cast<double>(numberCells[Direction::DirX]),
                                           (posRight[Direction::DirY] - posLeft[Direction::DirY]) /
                                               static_cast<double>(numberCells[Direction::DirY]),
                                           (posRight[Direction::DirZ] - posLeft[Direction::DirZ]) /
                                               static_cast<double>(numberCells[Direction::DirZ]) }),
      inverseCellSize({ 1.0 / cellSize[Direction::DirX], 1.0 / cellSize[Direction::DirY],
                        1.0 / cellSize[Direction::DirZ] }),
      boundaryTypes(boundaryTypes) {}

template <class Specific, class Fields, unsigned padding>
Problem<Specific, Fields, padding>::Problem(const Problem<Specific, Fields, padding>& blueprint)
    : cflThreshold(blueprint.cflThreshold), thermal(blueprint.thermal),
      timeDelta(blueprint.timeDelta), timeStart(blueprint.timeStart), timeEnd(blueprint.timeEnd),
      preciseEnd(blueprint.preciseEnd), gamma(blueprint.gamma), posLeft(blueprint.posLeft),
      posRight(blueprint.posRight), numberCells(blueprint.numberCells),
      cellSize(blueprint.cellSize), inverseCellSize(blueprint.inverseCellSize),
      boundaryTypes(blueprint.boundaryTypes) {}

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
