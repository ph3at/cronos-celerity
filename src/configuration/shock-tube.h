#pragma once

#include "problem.h"

#include <vector>

#include <celerity.h>
#include <sycl/sycl.hpp>

class ShockTube : public Problem<ShockTube, FieldStruct, GHOST_CELLS> {
  public:
    ShockTube(const double cflThreshold, const bool thermal, const double timeDelta, const double timeStart,
              const double timeEnd, const bool preciseEnd, const double gamma,
              const std::array<double, Direction::DirMax> posLeft, const std::array<double, Direction::DirMax> posRight,
              const std::array<std::size_t, Direction::DirMax> numberCells,
              const std::array<BoundaryType, Faces::FaceMax> boundaryTypes, const Direction shockDir,
              const double shockPos, const double densityLeftInit, const double densityRightInit,
              const double velocityLeftInit, const double velocityRightInit, const double pressureLeftInit,
              const double pressureRightInit);
    ShockTube(const ShockTube&) = default;
    ShockTube(const toml::table& config);

    void initialiseGrid(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const;
    void initialiseGridSycl(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid) const;
    void initialiseGridCelerity(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid) const;
    void applyBoundary(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned field, const unsigned face) const;
    void applyBoundarySycl(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid, const unsigned field,
                           const unsigned face) const;
    void applyBoundaryCelerity(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid,
                               const unsigned field, const unsigned face) const;
    void applyBoundaryCelerity3D(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid,
                                 const unsigned face) const;
    void applySource(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const;
    void applySourceSycl(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid) const;
    void applySourceCelerity(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid) const;

    Direction shockDir;
    double shockPos;

    double densityLeftInit, densityRightInit;
    double velocityXLeftInit, velocityXRightInit;
    double velocityYLeftInit, velocityYRightInit;
    double velocityZLeftInit, velocityZRightInit;
    double pressureLeftInit, pressureRightInit;
};
