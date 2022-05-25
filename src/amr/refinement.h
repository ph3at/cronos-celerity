#pragma once

#include <cmath>
#include <limits>

#include "amr-node.h"
#include "flags.h"

typedef std::array<std::pair<unsigned, unsigned>, 3> GridBoundary;

template <class SolverType, class ProblemType, class Fields, unsigned padding> class Refinery {
  public:
    Refinery(const Problem<ProblemType>& problem, const AMRParameters& configuration);

    void refine();
    AMRNode<SolverType, ProblemType, Fields, padding>&
    initialRefine(PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem);

  private:
    const Problem<ProblemType>& problem;
    const AMRParameters& configuration;

    std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>> amrNodes;

    std::vector<CellFlags::Flags> flagCells() const;
    double efficiency(const GridBoundary& grid, const CellFlags::Flags& flags) const;
    GridBoundary fullGrid(const CellFlags::Flags& flags) const;
    std::vector<GridBoundary> divideGrid(const GridBoundary& grid,
                                         const CellFlags::Flags& flags) const;
    std::vector<GridBoundary> mergeGrids(const std::vector<GridBoundary>& grids,
                                         const CellFlags::Flags& flags) const;
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>
    createGrids(std::vector<GridBoundary>& grids, const unsigned level);
    void initialiseGrid(PaddedGrid<Fields, padding>& grid, const GridBoundary& gridBoundary,
                        const unsigned level);
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
Refinery<SolverType, ProblemType, Fields, padding>::Refinery(const Problem<ProblemType>& problem,
                                                             const AMRParameters& configuration)
    : problem(problem), configuration(configuration) {}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::refine() {
    std::vector<CellFlags::Flags> flags = this->flagCells();
    for (unsigned level = 0; level < flags.size(); level++) {
        CellFlags::Flags& levelFlags = flags[level];
        std::vector<GridBoundary> minimalGrids =
            this->divideGrid(this->fullGrid(levelFlags), levelFlags);
        std::vector<GridBoundary> newGrids = this->mergeGrids(minimalGrids, levelFlags);
        createGrids(newGrids, level + 1);
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<CellFlags::Flags>
Refinery<SolverType, ProblemType, Fields, padding>::flagCells() const {
    const int maxLevel = static_cast<int>(this->amrNodes.size());
    std::vector<CellFlags::Flags> levelFlags;
    levelFlags.reserve(maxLevel);
    for (int level = maxLevel - 1; level >= 0; level--) {
        for (AMRNode<SolverType, ProblemType, Fields, padding> node : this->amrNodes[level]) {
            CellFlags::Flags nodeFlags = {}; // truncation error stuff
            CellFlags::addFlags(levelFlags[level], nodeFlags);
        }
        if (level < maxLevel - 1) {
            CellFlags::Flags higherLevelFlags = CellFlags::translateFlags(
                levelFlags[level + 1], this->configuration.refinementFactor);
            CellFlags::addFlags(levelFlags[level], higherLevelFlags);
        }
        CellFlags::addBuffer(levelFlags[level], this->configuration.bufferSize);
    }
    return levelFlags;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
double Refinery<SolverType, ProblemType, Fields, padding>::efficiency(
    const GridBoundary& grid, const CellFlags::Flags& flags) const {
    unsigned flaggedCells = 0;
    for (unsigned x = std::max(grid[0].first, flags.first);
         x < std::min(grid[0].second, flags.first + static_cast<unsigned>(flags.second.size()));
         x++) {
        const unsigned xx = x - flags.first;
        const std::pair<unsigned, CellFlags::FlagsY>& flagsY = flags.second[xx];
        for (unsigned y = std::max(grid[1].first, flagsY.first);
             y <
             std::min(grid[1].second, flagsY.first + static_cast<unsigned>(flagsY.second.size()));
             y++) {
            const unsigned yy = y - flagsY.first;
            for (std::pair<unsigned, unsigned> pair : flagsY.second[yy]) {
                const int left = std::max(pair.first, grid[3].first);
                const int right = std::min(pair.second + 1, grid[2].second);
                flaggedCells += std::max(0, right - left);
            }
        }
    }
    const double gridSize =
        static_cast<double>((grid[0].second - grid[0].first) * (grid[1].second - grid[1].first) *
                            (grid[2].second - grid[2].first));
    return static_cast<double>(flaggedCells) / gridSize;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
GridBoundary
Refinery<SolverType, ProblemType, Fields, padding>::fullGrid(const CellFlags::Flags& flags) const {
    const unsigned xStart = flags.first;
    const unsigned xEnd = xStart + flags.second.size();
    unsigned yStart = std::numeric_limits<unsigned>::max();
    unsigned yEnd = 0;
    unsigned zStart = std::numeric_limits<unsigned>::max();
    unsigned zEnd = 0;
    for (std::pair<unsigned, CellFlags::FlagsY> flagsY : flags.second) {
        yStart = std::min(yStart, flagsY.first);
        yEnd = std::max(yEnd, flagsY.first + static_cast<unsigned>(flagsY.second.size()));
        for (CellFlags::FlagsZ flagsZ : flagsY.second) {
            const size_t nPairs = flagsZ.size();
            if (nPairs > 0) {
                zStart = std::min(zStart, flagsZ[0].first);
                zEnd = std::max(zEnd, flagsZ[nPairs - 1].second);
            }
        }
    }
    return { std::make_pair(xStart, xEnd), std::make_pair(yStart, yEnd),
             std::make_pair(zStart, zEnd) };
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<GridBoundary> Refinery<SolverType, ProblemType, Fields, padding>::divideGrid(
    const GridBoundary& grid, const CellFlags::Flags& flags) const {
    if (this->efficiency(grid, flags) < this->configuration.efficiencyThreshold) {
        GridBoundary subGrid1;
        GridBoundary subGrid2;
        unsigned longest =
            std::max(grid[0].second - grid[0].first,
                     std::max(grid[1].second - grid[1].first, grid[2].second - grid[2].first));
        if (grid[0].second - grid[0].first == longest) {
            subGrid1 = { std::make_pair(grid[0].first, grid[0].first + grid[0].second / 2), grid[1],
                         grid[2] };
            subGrid2 = { std::make_pair(grid[0].first + grid[0].second / 2 + 1, grid[0].second),
                         grid[1], grid[2] };
        } else if (grid[1].second - grid[1].first == longest) {
            subGrid1 = { grid[0], std::make_pair(grid[1].first, grid[1].first + grid[1].second / 2),
                         grid[2] };
            subGrid2 = { grid[0],
                         std::make_pair(grid[1].first + grid[1].second / 2 + 1, grid[1].second),
                         grid[2] };
        } else {
            subGrid1 = { grid[0], grid[1],
                         std::make_pair(grid[2].first, grid[2].first + grid[2].second / 2) };
            subGrid2 = { grid[0], grid[1],
                         std::make_pair(grid[2].first + grid[2].second / 2 + 1, grid[2].second) };
        }
        std::vector<GridBoundary> grids1 = this->divideGrid(subGrid1, flags);
        std::vector<GridBoundary> grids2 = this->divideGrid(subGrid2, flags);
        grids1.insert(grids1.end(), grids2.begin(), grids2.end());
        return grids1;
    } else {
        return std::vector<GridBoundary>({ grid });
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<GridBoundary> Refinery<SolverType, ProblemType, Fields, padding>::mergeGrids(
    const std::vector<GridBoundary>& grids, const CellFlags::Flags& flags) const {}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>
Refinery<SolverType, ProblemType, Fields, padding>::createGrids(std::vector<GridBoundary>& grids,
                                                                const unsigned level) {

    AMRNode<SolverType, ProblemType, Fields, padding>& root = this->amrNodes[0][0];

    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>> newGrids;
    newGrids.reserve(grids.size());
    const double factor =
        1.0 / static_cast<double>(std::pow(this->configuration.refinementFactor, level));
    const std::array<double, 3> cellSize = { root.grid.cellSize[0] * factor,
                                             root.grid.cellSize[1] * factor,
                                             root.grid.cellSize[2] * factor };
    const double timeDelta = root.solver.timeDelta * factor;
    const double timeCurrent = root.solver.timeCurrent;
    for (GridBoundary grid : grids) {
        std::array<double, 3> dim;
        std::array<double, 3> posLeft;
        std::array<double, 3> posRight;
        for (unsigned dir = 0; dir < 3; dir++) {
            posLeft[dir] = static_cast<double>(grid[dir].first) * cellSize[dir];
            if (posLeft[dir] < (root.grid.posLeft[dir])) {
                posLeft[dir] = root.grid.posLeft[dir];
                grid[dir].first = static_cast<unsigned>(posLeft[dir] / cellSize[dir]);
            }
            posRight[dir] = static_cast<double>(grid[dir].second) * cellSize[dir];
            if (posRight[dir] < (root.grid.posRight[dir])) {
                posRight[dir] = root.grid.posRight[dir];
                grid[dir].second = static_cast<unsigned>(posRight[dir] / cellSize[dir]);
            }
            dim[dir] = grid[dir].second - grid[dir].first;
        }
        PaddedGrid<Fields, padding> newGrid(this->amrNodes[0][0].grid.defaultValue, dim[0], dim[1],
                                            dim[2], posLeft, posRight);
        this->initialiseGrid(newGrid, grid, level);
        newGrids.push_back(AMRNode<SolverType, ProblemType, Fields, padding>(
            newGrid, this->problem, this->configuration, timeDelta, timeCurrent));
    }
    return newGrids;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::initialiseGrid(
    PaddedGrid<Fields, padding>& grid, const GridBoundary& gridBoundary, const unsigned level) {
    // Could be optimised with smaller search space for siblings and parents
    for (unsigned x = 0; x < grid.xDim(); x++) {
        const unsigned xGlobal = x + gridBoundary[0].first;
        for (unsigned y = 0; y < grid.yDim(); y++) {
            const unsigned yGlobal = y + gridBoundary[1].first;
            for (unsigned z = 0; z < grid.zDim(); z++) {
                const unsigned zGlobal = z + gridBoundary[2].first;
                bool foundCell = false;
                if (level < this->amrNodes.size()) {
                    for (AMRNode<SolverType, ProblemType, Fields, padding>& sibling :
                         this->amrNodes[level]) {
                        if (xGlobal >= sibling.gridBoundary[0].first &&
                            xGlobal < sibling.gridBoundary[0].second &&
                            yGlobal >= sibling.gridBoundary[1].first &&
                            yGlobal < sibling.gridBoundary[1].second &&
                            zGlobal >= sibling.gridBoundary[2].first &&
                            zGlobal < sibling.gridBoundary[2].second) {
                            const unsigned xSibling = xGlobal - sibling.gridBoundary[0].first;
                            const unsigned ySibling = yGlobal - sibling.gridBoundary[1].first;
                            unsigned zSibling = zGlobal - sibling.gridBoundary[2].first;
                            while (z < grid.zDim() && zSibling < sibling.grid.zDim()) {
                                grid(x, y, z) = sibling.grid(xSibling, ySibling, zSibling);
                                z++;
                                zSibling++;
                            }
                            foundCell = true;
                            break;
                        }
                    }
                }
                if (!foundCell) {
                    const unsigned xGlobalUp = xGlobal / this->configuration.refinementFactor;
                    const unsigned yGlobalUp = yGlobal / this->configuration.refinementFactor;
                    const unsigned zGlobalUp = zGlobal / this->configuration.refinementFactor;
                    for (AMRNode<SolverType, ProblemType, Fields, padding>& parent :
                         this->amrNodes[level - 1]) {
                        if (xGlobalUp >= parent.gridBoundary[0].first &&
                            xGlobalUp < parent.gridBoundary[0].second &&
                            yGlobalUp >= parent.gridBoundary[1].first &&
                            yGlobalUp < parent.gridBoundary[1].second &&
                            zGlobalUp >= parent.gridBoundary[2].first &&
                            zGlobalUp < parent.gridBoundary[2].second) {
                            const unsigned xParent = xGlobal - parent.gridBoundary[0].first;
                            const unsigned yParent = yGlobal - parent.gridBoundary[1].first;
                            const unsigned zParent = zGlobal - parent.gridBoundary[2].first;
                            // TODO: Interpolation
                            grid(x, y, z) = parent.grid(xParent, yParent, zParent);
                        }
                    }
                }
            }
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRNode<SolverType, ProblemType, Fields, padding>&
Refinery<SolverType, ProblemType, Fields, padding>::initialRefine(
    PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem) {
    AMRNode<SolverType, ProblemType, Fields, padding> root(grid, problem, this->configuration,
                                                           problem.timeDelta, 0.0);
    this->amrNodes.push_back(
        std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>(1, root));
    return this->amrNodes[0][0];
}
