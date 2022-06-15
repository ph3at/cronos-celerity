#pragma once

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>

#include "amr-node.h"
#include "error-estimation.h"
#include "flags.h"
#include "grid-boundary.h"

template <class SolverType, class ProblemType, class Fields, unsigned padding> class Refinery {
  public:
    Refinery(const Problem<ProblemType, Fields, padding>& problem,
             const AMRParameters& configuration);

    std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>>& getNodes();

    void refine();
    void initialRefine(const PaddedGrid<Fields, padding>& grid,
                       const Problem<ProblemType, Fields, padding>& problem);

  private:
    const Problem<ProblemType, Fields, padding>& problem;
    const AMRParameters& configuration;

    std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>> amrNodes;
    std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>> oldNodes;

    std::vector<CellFlags::Flags> flagCells() const;
    GridBoundary::Boundary fullGrid(const CellFlags::Flags& flags) const;
    std::vector<GridBoundary::Boundary> divideGrid(const GridBoundary::Boundary& grid,
                                                   const CellFlags::Flags& flags) const;
    std::vector<GridBoundary::Boundary> mergeGrids(const std::vector<GridBoundary::Boundary>& grids,
                                                   const CellFlags::Flags& flags) const;
    void ensureNesting(std::vector<GridBoundary::Boundary>& grids,
                       const std::vector<GridBoundary::Boundary>& parentGrids);
    std::vector<GridBoundary::Boundary>
    translateGrids(const std::vector<GridBoundary::Boundary>& grids);
    void createGrids(const std::vector<GridBoundary::Boundary>& grids, const unsigned level,
                     const bool initialiseFromProblem);
    void initialiseGrid(const unsigned level, const unsigned nodeId);
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
Refinery<SolverType, ProblemType, Fields, padding>::Refinery(
    const Problem<ProblemType, Fields, padding>& problem, const AMRParameters& configuration)
    : problem(problem), configuration(configuration) {}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>>&
Refinery<SolverType, ProblemType, Fields, padding>::getNodes() {
    return this->amrNodes;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::initialRefine(
    const PaddedGrid<Fields, padding>& grid, const Problem<ProblemType, Fields, padding>& problem) {
    this->amrNodes.clear();
    this->oldNodes.clear();

    GridBoundary::Boundary rootBoundary = { std::make_pair(0, grid.xEnd() - padding - 1),
                                            std::make_pair(0, grid.yEnd() - padding - 1),
                                            std::make_pair(0, grid.zEnd() - padding - 1) };
    this->amrNodes.push_back({});
    this->amrNodes[0].emplace_back(0, grid, rootBoundary, problem, this->configuration,
                                   problem.timeDelta, problem.timeStart);
    AMRNode<SolverType, ProblemType, Fields, padding>& root = this->amrNodes[0][0];

    const bool initialiseFromProblem = true;
    unsigned level = 0;
    std::vector<GridBoundary::Boundary> parents(1, root.gridBoundary);
    while (true) {
        std::cout << "Refining level " << level << std::endl;
        CellFlags::Flags levelFlags;
        bool noFlags = true;
        for (AMRNode<SolverType, ProblemType, Fields, padding>& node : this->amrNodes[level]) {
            std::optional<CellFlags::Flags> gridFlags =
                ErrorEstimation::getFlags<ProblemType, Fields, padding>(
                    node.grid,
                    { node.gridBoundary[0].first, node.gridBoundary[1].first,
                      node.gridBoundary[2].first },
                    node.solver.timeDelta, this->configuration.truncationErrorTreshold,
                    this->problem, root.grid.boundaryTypes);
            if (gridFlags.has_value()) {
                if (noFlags) {
                    noFlags = false;
                    levelFlags.swap(gridFlags.value());
                } else {
                    CellFlags::addFlags(levelFlags, gridFlags.value());
                }
            }
        }
        if (noFlags) {
            break;
        } else {
            CellFlags::addBuffer(levelFlags, this->configuration.bufferSize);
            std::vector<GridBoundary::Boundary> minimalGrids =
                this->divideGrid(this->fullGrid(levelFlags), levelFlags);
            std::vector<GridBoundary::Boundary> newGrids =
                this->mergeGrids(minimalGrids, levelFlags);
            this->ensureNesting(newGrids, parents);
            std::vector<GridBoundary::Boundary> finalNewGrids = this->translateGrids(newGrids);
            this->createGrids(finalNewGrids, level + 1, initialiseFromProblem);
            parents.swap(finalNewGrids);
            level++;
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::refine() {
    std::swap(this->amrNodes, this->oldNodes);
    this->amrNodes.clear();

    std::vector<CellFlags::Flags> flags = this->flagCells();
    this->amrNodes.reserve(flags.size() + 1);
    this->amrNodes.push_back({});
    this->amrNodes[0].swap(this->oldNodes[0]);
    AMRNode<SolverType, ProblemType, Fields, padding>& root = this->amrNodes[0][0];
    root.removeChildren();

    const bool initialiseFromProblem = false;
    std::vector<GridBoundary::Boundary> parents({ root.gridBoundary });
    for (int level = flags.size() - 1; level >= 0; level--) {
        std::cout << "Rerefining level " << flags.size() - 1 - level << ":" << std::endl;
        CellFlags::Flags& levelFlags = flags[level];
        std::vector<GridBoundary::Boundary> minimalGrids =
            this->divideGrid(this->fullGrid(levelFlags), levelFlags);
        std::vector<GridBoundary::Boundary> newGrids = this->mergeGrids(minimalGrids, levelFlags);
        this->ensureNesting(newGrids, parents);
        std::vector<GridBoundary::Boundary> finalNewGrids = this->translateGrids(newGrids);
        this->createGrids(finalNewGrids, flags.size() - level, initialiseFromProblem);
        parents.swap(finalNewGrids);
    }

    this->oldNodes.clear();
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<CellFlags::Flags>
Refinery<SolverType, ProblemType, Fields, padding>::flagCells() const {
    const int maxLevel = static_cast<int>(this->oldNodes.size());
    std::vector<CellFlags::Flags> flags;
    bool refine = false;
    for (int level = maxLevel - 1; level >= 0; level--) {
        bool first = true;
        CellFlags::Flags levelFlags;
        for (const AMRNode<SolverType, ProblemType, Fields, padding>& node :
             this->oldNodes[level]) {
            std::optional<CellFlags::Flags> nodeFlags =
                ErrorEstimation::getFlags<ProblemType, Fields, padding>(
                    node.grid,
                    { node.gridBoundary[0].first, node.gridBoundary[1].first,
                      node.gridBoundary[2].first },
                    node.solver.timeDelta, this->configuration.truncationErrorTreshold,
                    this->problem, this->oldNodes[0][0].grid.boundaryTypes);
            if (nodeFlags.has_value()) {
                if (first) {
                    refine = true;
                    first = false;
                    levelFlags.swap(nodeFlags.value());
                } else {
                    CellFlags::addFlags(levelFlags, nodeFlags.value());
                }
            }
        }
        if (refine) {
            if (flags.size() > 0) {
                CellFlags::Flags higherLevelFlags =
                    CellFlags::translateFlags(flags.back(), this->configuration.refinementFactor);
                CellFlags::addFlags(levelFlags, higherLevelFlags);
            }
            CellFlags::addBuffer(levelFlags, this->configuration.bufferSize);
            flags.push_back(std::make_pair(levelFlags.first, CellFlags::FlagsX({})));
            flags.back().second.swap(levelFlags.second);
        }
    }
    return flags;
}

inline double gridEfficiency(const GridBoundary::Boundary& grid, const CellFlags::Flags& flags) {
    int flaggedCells = 0;
    for (unsigned x = std::max(grid[0].first, flags.first);
         x < std::min(grid[0].second + 1, flags.first + static_cast<unsigned>(flags.second.size()));
         x++) {
        const unsigned xx = x - flags.first;
        const std::pair<unsigned, CellFlags::FlagsY>& flagsY = flags.second[xx];
        for (unsigned y = std::max(grid[1].first, flagsY.first);
             y < std::min(grid[1].second + 1,
                          flagsY.first + static_cast<unsigned>(flagsY.second.size()));
             y++) {
            const unsigned yy = y - flagsY.first;
            for (std::pair<unsigned, unsigned> pair : flagsY.second[yy]) {
                const int left = static_cast<int>(std::max(pair.first, grid[2].first));
                const int right = static_cast<int>(std::min(pair.second, grid[2].second));
                flaggedCells += std::max(0, right - left + 1);
            }
        }
    }
    const double gridSize = static_cast<double>((grid[0].second - grid[0].first + 1) *
                                                (grid[1].second - grid[1].first + 1) *
                                                (grid[2].second - grid[2].first + 1));
    return static_cast<double>(flaggedCells) / gridSize;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
GridBoundary::Boundary
Refinery<SolverType, ProblemType, Fields, padding>::fullGrid(const CellFlags::Flags& flags) const {
    if (flags.second.size() == 0) {
        return { std::make_pair(0, 0), std::make_pair(0, 0), std::make_pair(0, 0) };
    } else {
        const unsigned xStart = flags.first;
        const unsigned xEnd = xStart + flags.second.size() - 1;
        unsigned yStart = std::numeric_limits<unsigned>::max();
        unsigned yEnd = 0;
        unsigned zStart = std::numeric_limits<unsigned>::max();
        unsigned zEnd = 0;
        for (std::pair<unsigned, CellFlags::FlagsY> flagsY : flags.second) {
            yStart = std::min(yStart, flagsY.first);
            const unsigned localYEnd = static_cast<unsigned>(std::max(
                0,
                static_cast<int>(flagsY.first + static_cast<unsigned>(flagsY.second.size())) - 1));
            yEnd = std::max(yEnd, localYEnd);
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
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<GridBoundary::Boundary> Refinery<SolverType, ProblemType, Fields, padding>::divideGrid(
    const GridBoundary::Boundary& grid, const CellFlags::Flags& flags) const {
    const double efficiency = gridEfficiency(grid, flags);
    if (efficiency == 0.0) {
        return std::vector<GridBoundary::Boundary>(0);
    } else if (efficiency < this->configuration.efficiencyThreshold) {
        GridBoundary::Boundary subGrid1;
        GridBoundary::Boundary subGrid2;
        // TODO: What of very small grids?
        if (grid[0].second - grid[0].first > grid[1].second - grid[1].first &&
            grid[0].second - grid[0].first > grid[2].second - grid[2].first) {
            subGrid1 = { std::make_pair(grid[0].first, (grid[0].first + grid[0].second) / 2),
                         grid[1], grid[2] };
            subGrid2 = { std::make_pair((grid[0].first + grid[0].second) / 2 + 1, grid[0].second),
                         grid[1], grid[2] };
        } else if (grid[1].second - grid[1].first > grid[2].second - grid[2].first) {
            subGrid1 = { grid[0],
                         std::make_pair(grid[1].first, (grid[1].first + grid[1].second) / 2),
                         grid[2] };
            subGrid2 = { grid[0],
                         std::make_pair((grid[1].first + grid[1].second) / 2 + 1, grid[1].second),
                         grid[2] };
        } else {
            subGrid1 = { grid[0], grid[1],
                         std::make_pair(grid[2].first, (grid[2].first + grid[2].second) / 2) };
            subGrid2 = { grid[0], grid[1],
                         std::make_pair((grid[2].first + grid[2].second) / 2 + 1, grid[2].second) };
        }
        std::vector<GridBoundary::Boundary> grids1 = this->divideGrid(subGrid1, flags);
        std::vector<GridBoundary::Boundary> grids2 = this->divideGrid(subGrid2, flags);
        grids1.insert(grids1.end(), grids2.begin(), grids2.end());
        return grids1;
    } else {
        return std::vector<GridBoundary::Boundary>({ grid });
    }
}

inline double gridCost(const GridBoundary::Boundary& grid, const CellFlags::Flags& flags) {
    const double x = static_cast<double>(grid[0].second - grid[0].first);
    const double y = static_cast<double>(grid[1].second - grid[1].first);
    const double z = static_cast<double>(grid[2].second - grid[2].first);
    const double efficiency = gridEfficiency(grid, flags);
    return (0.0 - efficiency) * (x * y * z + x * y + x * z + y * z + x + y + z);
}

inline std::optional<GridBoundary::Boundary> tryMerge(const GridBoundary::Boundary& grid1,
                                                      const GridBoundary::Boundary& grid2,
                                                      const CellFlags::Flags& flags) {
    if (GridBoundary::isOverlappingOrAdjaecent(grid1, grid2)) {
        const GridBoundary::Boundary mergedGrid = {
            std::make_pair(std::min(grid1[0].first, grid2[0].first),
                           std::max(grid1[0].second, grid2[0].second)),
            std::make_pair(std::min(grid1[1].first, grid2[1].first),
                           std::max(grid1[1].second, grid2[1].second)),
            std::make_pair(std::min(grid1[2].first, grid2[2].first),
                           std::max(grid1[2].second, grid2[2].second))
        };
        if (gridCost(mergedGrid, flags) < gridCost(grid1, flags) + gridCost(grid2, flags)) {
            return std::make_optional(mergedGrid);
        }
    }
    return std::nullopt;
}

inline std::vector<GridBoundary::Boundary>
mergeGridsAux(const std::vector<GridBoundary::Boundary>& gridsKnown,
              const std::vector<GridBoundary::Boundary>& gridsNew, const CellFlags::Flags& flags) {
    std::vector<bool> keepGrid(gridsNew.size(), true);
    std::vector<GridBoundary::Boundary> newGrids;
    std::vector<GridBoundary::Boundary> knownGrids;
    for (const GridBoundary::Boundary& grid1 : gridsKnown) {
        bool keep = true;
        for (unsigned grid2 = 0; grid2 < newGrids.size(); grid2++) {
            if (keepGrid[grid2]) {
                std::optional<GridBoundary::Boundary> mergedGrid =
                    tryMerge(grid1, newGrids[grid2], flags);
                if (mergedGrid.has_value()) {
                    keep = false;
                    keepGrid[grid2] = false;
                    newGrids.push_back(mergedGrid.value());
                    break;
                }
            }
        }
        if (keep) {
            knownGrids.push_back(grid1);
        } else {
            newGrids.push_back(grid1);
        }
    }
    if (newGrids.size() == 0) {
        knownGrids.insert(knownGrids.end(), gridsNew.begin(), gridsNew.end());
        return knownGrids;
    } else {
        for (unsigned grid = 0; grid < newGrids.size(); grid++) {
            if (keepGrid[grid]) {
                knownGrids.push_back(gridsNew[grid]);
            }
        }
        return mergeGridsAux(knownGrids, newGrids, flags);
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<GridBoundary::Boundary> Refinery<SolverType, ProblemType, Fields, padding>::mergeGrids(
    const std::vector<GridBoundary::Boundary>& grids, const CellFlags::Flags& flags) const {
    std::vector<GridBoundary::Boundary> newGrids;
    std::vector<bool> keepGrid(grids.size(), true);
    for (unsigned grid1 = 0; grid1 + 1 < grids.size(); grid1++) {
        if (keepGrid[grid1]) {
            for (unsigned grid2 = grid1 + 1; grid2 < grids.size(); grid2++) {
                if (keepGrid[grid2]) {
                    std::optional<GridBoundary::Boundary> mergedGrid =
                        tryMerge(grids[grid1], grids[grid2], flags);
                    if (mergedGrid.has_value()) {
                        keepGrid[grid1] = false;
                        keepGrid[grid2] = false;
                        newGrids.push_back(mergedGrid.value());
                        break;
                    }
                }
            }
        }
    }
    if (newGrids.size() == 0) {
        return grids;
    } else {
        std::vector<GridBoundary::Boundary> knownGrids;
        for (unsigned grid = 0; grid < grids.size(); grid++) {
            if (keepGrid[grid]) {
                knownGrids.push_back(grids[grid]);
            }
        }
        return mergeGridsAux(knownGrids, newGrids, flags);
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::ensureNesting(
    std::vector<GridBoundary::Boundary>& grids,
    const std::vector<GridBoundary::Boundary>& parentGrids) {
    std::vector<GridBoundary::Boundary> newGrids;
    for (GridBoundary::Boundary& grid : grids) {
        std::vector<GridBoundary::Boundary> notCovered;
        notCovered.push_back(grid);
        for (const GridBoundary::Boundary& parent : parentGrids) {
            std::optional<GridBoundary::Boundary> overlapping =
                GridBoundary::overlappingArea(grid, parent);
            if (overlapping.has_value()) {
                GridBoundary::Boundary overlap = overlapping.value();
                std::vector<GridBoundary::Boundary> newParts;
                for (GridBoundary::Boundary& part : notCovered) {
                    std::vector<GridBoundary::Boundary> leftOver =
                        GridBoundary::subtract(part, overlap);
                    newParts.insert(newParts.end(), leftOver.begin(), leftOver.end());
                }
                notCovered.swap(newParts);
                if (notCovered.size() == 0) {
                    break;
                }
            }
        }
        if (notCovered.size() > 0) {
            std::vector<GridBoundary::Boundary> remaining =
                GridBoundary::remainingGrids(grid, notCovered);
            grid = remaining[0];
            newGrids.insert(newGrids.end(), remaining.begin() + 1, remaining.end());
        }
    }
    grids.insert(grids.end(), newGrids.begin(), newGrids.end());
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<GridBoundary::Boundary>
Refinery<SolverType, ProblemType, Fields, padding>::translateGrids(
    const std::vector<GridBoundary::Boundary>& grids) {
    std::vector<GridBoundary::Boundary> translatedGrids;
    translatedGrids.reserve(grids.size());
    std::transform(
        grids.begin(), grids.end(), std::back_inserter(translatedGrids),
        [factor = this->configuration.refinementFactor](const GridBoundary::Boundary& parent) {
            return GridBoundary::translateUp(parent, factor);
        });
    return translatedGrids;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::createGrids(
    const std::vector<GridBoundary::Boundary>& grids, const unsigned level,
    const bool initialiseFromProblem) {
    const AMRNode<SolverType, ProblemType, Fields, padding>& root = this->amrNodes[0][0];
    this->amrNodes.emplace_back();
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>& newGrids =
        this->amrNodes[level];
    newGrids.reserve(grids.size());
    const double factor =
        1.0 / static_cast<double>(std::pow(this->configuration.refinementFactor, level));
    const std::array<double, 3> cellSize = { root.grid.cellSize[0] * factor,
                                             root.grid.cellSize[1] * factor,
                                             root.grid.cellSize[2] * factor };
    const double timeDelta = root.solver.timeDelta * factor;
    const double timeCurrent = root.solver.timeCurrent;
    unsigned nodeId = 0;
    std::cout << "Level: " << level << std::endl;
    for (const GridBoundary::Boundary& grid : grids) {
        std::cout << "([" << grid[0].first << "," << grid[0].second << "],[" << grid[1].first << ","
                  << grid[1].second << "],[" << grid[2].first << "," << grid[2].second << "])"
                  << std::endl;
        std::array<unsigned, 3> dim;
        std::array<double, 3> posLeft;
        std::array<double, 3> posRight;
        for (unsigned dir = 0; dir < 3; dir++) {
            posLeft[dir] = static_cast<double>(grid[dir].first) * cellSize[dir];
            posRight[dir] = static_cast<double>(grid[dir].second) * cellSize[dir];
            dim[dir] = grid[dir].second - grid[dir].first;
        }
        newGrids.emplace_back(nodeId, root.grid.defaultValue, dim, posLeft, posRight, grid,
                              this->problem, this->configuration, timeDelta, timeCurrent);
        if (initialiseFromProblem) {
            this->problem.initialiseGrid(newGrids[nodeId].grid);
        } else {
            this->initialiseGrid(level, nodeId);
        }
        nodeId++;
    }

    for (unsigned node1 = 0; node1 < newGrids.size(); node1++) {
        AMRNode<SolverType, ProblemType, Fields, padding>& grid1 = newGrids[node1];
        for (unsigned node2 = node1 + 1; node2 < newGrids.size(); node2++) {
            AMRNode<SolverType, ProblemType, Fields, padding>& grid2 = newGrids[node2];
            if (GridBoundary::isOverlappingOrAdjaecent(grid1.gridBoundary, grid2.gridBoundary)) {
                grid1.addSibling(&grid2);
                grid2.addSibling(&grid1);
            }
        }
        for (AMRNode<SolverType, ProblemType, Fields, padding>& parent :
             this->amrNodes[level - 1]) {
            if (GridBoundary::isOverlappingOrAdjaecent(
                    grid1.gridBoundary,
                    GridBoundary::translateUp(parent.gridBoundary,
                                              this->configuration.refinementFactor))) {
                grid1.addParent(&parent);
                parent.addChild(&grid1);
            }
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::initialiseGrid(const unsigned level,
                                                                        const unsigned nodeId) {
    // Could be optimised with smaller search space for siblings and parents
    PaddedGrid<Fields, padding>& grid = this->amrNodes[level][nodeId].grid;
    const GridBoundary::Boundary& gridBoundary = this->amrNodes[level][nodeId].gridBoundary;
    const double refinementFactor = static_cast<double>(this->configuration.refinementFactor);
    for (unsigned x = 0; x < grid.xDim(); x++) {
        const int xGlobal = static_cast<int>(x + gridBoundary[0].first) - padding;
        const double xGlobalUp = static_cast<double>(xGlobal) / refinementFactor;
        for (unsigned y = 0; y < grid.yDim(); y++) {
            const int yGlobal = static_cast<int>(y + gridBoundary[1].first) - padding;
            const double yGlobalUp = static_cast<double>(yGlobal) / refinementFactor;
            for (unsigned z = 0; z < grid.zDim(); z++) {
                const int zGlobal = static_cast<int>(z + gridBoundary[2].first) - padding;
                bool foundCell = false;
                if (level < this->oldNodes.size()) {
                    for (AMRNode<SolverType, ProblemType, Fields, padding>& sibling :
                         this->oldNodes[level]) {
                        if (xGlobal >= static_cast<int>(sibling.gridBoundary[0].first) &&
                            xGlobal <= static_cast<int>(sibling.gridBoundary[0].second) &&
                            yGlobal >= static_cast<int>(sibling.gridBoundary[1].first) &&
                            yGlobal <= static_cast<int>(sibling.gridBoundary[1].second) &&
                            zGlobal >= static_cast<int>(sibling.gridBoundary[2].first) &&
                            zGlobal <= static_cast<int>(sibling.gridBoundary[2].second)) {
                            const unsigned xSibling =
                                xGlobal - sibling.gridBoundary[0].first + padding;
                            const unsigned ySibling =
                                yGlobal - sibling.gridBoundary[1].first + padding;
                            unsigned zSibling = zGlobal - sibling.gridBoundary[2].first + padding;
                            while (z < grid.zDim() && zSibling < sibling.grid.zDim()) {
                                grid(x, y, z) = sibling.grid(xSibling, ySibling, zSibling);
                                z++;
                                zSibling++;
                            }
                            z--;
                            foundCell = true;
                            break;
                        }
                    }
                }
                if (!foundCell) {
                    const double zGlobalUp = static_cast<double>(zGlobal) / refinementFactor;
                    for (const AMRNode<SolverType, ProblemType, Fields, padding>& parent :
                         this->amrNodes[level - 1]) {
                        std::optional<Fields> fields =
                            parent.valueAtUp(xGlobalUp, yGlobalUp, zGlobalUp);
                        if (fields.has_value()) {
                            grid(x, y, z) = fields.value();
                            break;
                        }
                    }
                }
            }
        }
    }
}
