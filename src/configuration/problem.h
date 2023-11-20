#pragma once

#include <array>

#include <toml++/toml.h>

#include <celerity.h>
#include <sycl/sycl.hpp>

#include "../boundary/boundary-types.h"
#include "../data-types/direction.h"
#include "../data-types/faces.h"
#include "../data-types/phys-fields.h"
#include "../grid/padded-grid.h"
#include "../user-interface/value-parser.h"
#include "constants.h"

template <class Specific, class Fields, unsigned padding> class Problem {
  public:
    Problem(const double cflThreshold, const bool thermal, const double timeDelta, const double timeStart,
            const double timeEnd, const bool preciseEnd, const double gamma,
            const std::array<double, Direction::DirMax> posLeft, const std::array<double, Direction::DirMax> posRight,
            const std::array<std::size_t, Direction::DirMax> numberCells,
            const std::array<BoundaryType, Faces::FaceMax> boundaryTypes);
    Problem(const Problem&) = default;
    Problem(const toml::table& config);

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
    void initialiseGridSycl(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid) const;
    void initialiseGridCelerity(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid) const;
    void applyBoundary(PaddedGrid<Fields, padding>& grid, const unsigned field, const unsigned face) const;
    void applyBoundarySycl(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid, const unsigned field,
                           const unsigned face) const;
    void applyBoundaryCelerity(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid,
                               const unsigned field, const unsigned face) const;
    void applyBoundaryCelerity3D(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid,
                                 const unsigned face) const;
    void applySource(PaddedGrid<Fields, padding>& grid) const;
    void applySourceSycl(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid) const;
    void applySourceCelerity(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid) const;
};

template <class Specific, class Fields, unsigned padding>
Problem<Specific, Fields, padding>::Problem(const double cflThreshold, const bool thermal, const double timeDelta,
                                            const double timeStart, const double timeEnd, const bool preciseEnd,
                                            const double gamma, const std::array<double, Direction::DirMax> posLeft,
                                            const std::array<double, Direction::DirMax> posRight,
                                            const std::array<std::size_t, Direction::DirMax> numberCells,
                                            const std::array<BoundaryType, Faces::FaceMax> boundaryTypes)
    : cflThreshold(cflThreshold), thermal(thermal), timeDelta(timeDelta), timeStart(timeStart), timeEnd(timeEnd),
      preciseEnd(preciseEnd), gamma(gamma), posLeft(posLeft), posRight(posRight), numberCells(numberCells),
      cellSize(
          { (posRight[Direction::DirX] - posLeft[Direction::DirX]) / static_cast<double>(numberCells[Direction::DirX]),
            (posRight[Direction::DirY] - posLeft[Direction::DirY]) / static_cast<double>(numberCells[Direction::DirY]),
            (posRight[Direction::DirZ] - posLeft[Direction::DirZ]) /
                static_cast<double>(numberCells[Direction::DirZ]) }),
      inverseCellSize(
          { 1.0 / cellSize[Direction::DirX], 1.0 / cellSize[Direction::DirY], 1.0 / cellSize[Direction::DirZ] }),
      boundaryTypes(boundaryTypes) {}

template <class Specific, class Fields, unsigned padding>
Problem<Specific, Fields, padding>::Problem(const toml::table& config) {
    std::cout << std::endl << "General problem parameters:" << std::endl;
    const toml::node_view<const toml::node>& general = config["general"];
    this->cflThreshold = parseValue<double>(general, "cfl_threshold", 0.4);
    this->thermal = parseValue<bool>(general, "thermal", true);
    this->timeDelta = parseValue<double>(general, "time_delta", 2e-6);
    this->timeStart = parseValue<double>(general, "time_start", 0.0);
    this->timeEnd = parseValue<double>(general, "time_end", 5e-4);
    this->preciseEnd = parseValue<bool>(general, "precise_end", true);
    this->gamma = parseValue<double>(general, "gamma", 1.4);

    std::cout << std::endl << "Grid domain:" << std::endl;
    const toml::node_view<const toml::node>& grid = config["grid"];
    this->posLeft[Direction::DirX] = parseValue<double>(grid, "x_start", 0.45);
    this->posRight[Direction::DirX] = parseValue<double>(grid, "x_end", 0.55);
    this->numberCells[Direction::DirX] = parseValue<std::size_t>(grid, "number_cells_x", 20);
    this->posLeft[Direction::DirY] = parseValue<double>(grid, "y_start", 0.45);
    this->posRight[Direction::DirY] = parseValue<double>(grid, "y_end", 0.55);
    this->numberCells[Direction::DirY] = parseValue<std::size_t>(grid, "number_cells_y", 5);
    this->posLeft[Direction::DirZ] = parseValue<double>(grid, "z_start", 0.45);
    this->posRight[Direction::DirZ] = parseValue<double>(grid, "z_end", 0.55);
    this->numberCells[Direction::DirZ] = parseValue<std::size_t>(grid, "number_cells_z", 5);
    for (unsigned dir = 0; dir < Direction::DirMax; dir++) {
        this->cellSize[dir] = (this->posRight[dir] - this->posLeft[dir]) / static_cast<double>(this->numberCells[dir]);
        this->inverseCellSize[dir] = 1.0 / this->cellSize[dir];
    }

    std::cout << std::endl << "Boundary types:" << std::endl;
    const toml::node_view<const toml::node>& boundary = config["boundary"];
    this->boundaryTypes[Faces::FaceWest] = parseValue<BoundaryType>(boundary, "west", BoundaryType::EMPTY);
    this->boundaryTypes[Faces::FaceEast] = parseValue<BoundaryType>(boundary, "east", BoundaryType::EMPTY);
    this->boundaryTypes[Faces::FaceSouth] = parseValue<BoundaryType>(boundary, "south", BoundaryType::EMPTY);
    this->boundaryTypes[Faces::FaceNorth] = parseValue<BoundaryType>(boundary, "north", BoundaryType::EMPTY);
    this->boundaryTypes[Faces::FaceBottom] = parseValue<BoundaryType>(boundary, "bottom", BoundaryType::EMPTY);
    this->boundaryTypes[Faces::FaceTop] = parseValue<BoundaryType>(boundary, "top", BoundaryType::EMPTY);
}

template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::initialiseGrid(PaddedGrid<Fields, padding>& grid) const {
    static_cast<const Specific*>(this)->initialiseGrid(grid);
}
template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::initialiseGridSycl(sycl::queue& queue,
                                                            sycl::buffer<FieldStruct, 3>& grid) const {
    static_cast<const Specific*>(this)->initialiseGridSycl(queue, grid);
}
template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::initialiseGridCelerity(celerity::distr_queue& queue,
                                                                celerity::buffer<FieldStruct, 3>& grid) const {
    static_cast<const Specific*>(this)->initialiseGridCelerity(queue, grid);
}

template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::applyBoundary(PaddedGrid<Fields, padding>& grid, const unsigned field,
                                                       const unsigned face) const {
    static_cast<const Specific*>(this)->applyBoundary(grid, field, face);
}
template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::applyBoundarySycl(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid,
                                                           const unsigned field, const unsigned face) const {
    static_cast<const Specific*>(this)->applyBoundarySycl(queue, grid, field, face);
}
template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::applyBoundaryCelerity(celerity::distr_queue& queue,
                                                               celerity::buffer<FieldStruct, 3>& grid,
                                                               const unsigned field, const unsigned face) const {
    static_cast<const Specific*>(this)->applyBoundaryCelerity(queue, grid, field, face);
}
template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::applyBoundaryCelerity3D(celerity::distr_queue& queue,
                                                                 celerity::buffer<FieldStruct, 3>& grid,
                                                                 const unsigned face) const {
    static_cast<const Specific*>(this)->applyBoundaryCelerity3D(queue, grid, face);
}

template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::applySource(PaddedGrid<Fields, padding>& grid) const {
    static_cast<const Specific*>(this)->applySource(grid);
}
template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::applySourceSycl(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid) const {
    static_cast<const Specific*>(this)->applySourceSycl(queue, grid);
}
template <class Specific, class Fields, unsigned padding>
void Problem<Specific, Fields, padding>::applySourceCelerity(celerity::distr_queue& queue,
                                                             celerity::buffer<FieldStruct, 3>& grid) const {
    static_cast<const Specific*>(this)->applySourceCelerity(queue, grid);
}

namespace ProblemCelerity {

template <unsigned padding, typename Problem>
void apply3D(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid, const unsigned face,
             const Problem& problem) {
    problem.applyBoundaryCelerity3D(queue, grid, face);
}

} // namespace ProblemCelerity
