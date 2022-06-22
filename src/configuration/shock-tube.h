#pragma once

#include "problem.h"

class ShockTube : public Problem<ShockTube, FieldStruct, GHOST_CELLS> {
  public:
    ShockTube(const double cflThreshold, const bool thermal, const double timeDelta,
              const double timeStart, const double timeEnd, const bool preciseEnd,
              const double gamma, const std::array<double, Direction::DirMax> posLeft,
              const std::array<double, Direction::DirMax> posRight,
              const std::array<std::size_t, Direction::DirMax> numberCells,
              const std::array<BoundaryType, Faces::FaceMax> boundaryTypes,
              const Direction shockDir, const double shockPos, const double densityLeftInit,
              const double densityRightInit, const double velocityLeftInit,
              const double velocityRightInit, const double pressureLeftInit,
              const double pressureRightInit);
    ShockTube(const ShockTube& blueprint);

    void initialiseGrid(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const;
    void applyBoundary(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned field,
                       const unsigned face) const;
    void applySource(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const;

    const Direction shockDir;
    const double shockPos;

    const double densityLeftInit, densityRightInit;
    const double velocityXLeftInit, velocityXRightInit;
    const double velocityYLeftInit, velocityYRightInit;
    const double velocityZLeftInit, velocityZRightInit;
    const double pressureLeftInit, pressureRightInit;

    static ShockTube initialiseTestProblem();
};

