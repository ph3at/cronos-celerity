#pragma once

#include <algorithm>
#include <cmath>
#include <iterator>
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
    GridBoundary fullGrid(const CellFlags::Flags& flags) const;
    std::vector<GridBoundary> divideGrid(const GridBoundary& grid,
                                         const CellFlags::Flags& flags) const;
    std::vector<GridBoundary> mergeGrids(const std::vector<GridBoundary>& grids,
                                         const CellFlags::Flags& flags) const;
    void ensureNesting(std::vector<GridBoundary>& grids,
                       const std::vector<GridBoundary>& parentGrids);
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>
    createGrids(std::vector<GridBoundary>& grids, const unsigned level,
                std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>& parents);
    void
    initialiseGrid(PaddedGrid<Fields, padding>& grid, const GridBoundary& gridBoundary,
                   const unsigned level,
                   const std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>& parents);
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
Refinery<SolverType, ProblemType, Fields, padding>::Refinery(const Problem<ProblemType>& problem,
                                                             const AMRParameters& configuration)
    : problem(problem), configuration(configuration) {}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRNode<SolverType, ProblemType, Fields, padding>&
Refinery<SolverType, ProblemType, Fields, padding>::initialRefine(
    PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem) {
    AMRNode<SolverType, ProblemType, Fields, padding> root(grid, problem, this->configuration,
                                                           problem.timeDelta, problem.timeStart);
    root.gridBoundary = { std::make_pair(0, grid.xEnd() - padding - 1),
                          std::make_pair(0, grid.yEnd() - padding - 1),
                          std::make_pair(0, grid.zEnd() - padding - 1) };
    this->amrNodes.push_back(
        std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>(1, root));

    unsigned level = 0;
    std::vector<GridBoundary> parents(1, root.gridBoundary);
    while (true) {
        CellFlags::Flags levelFlags;
        bool noFlags = true;
        for (AMRNode<SolverType, ProblemType, Fields, padding>& node : this->amrNodes[level]) {
            std::optional<CellFlags::Flags> gridFlags = {};
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
            std::vector<GridBoundary> minimalGrids =
                this->divideGrid(this->fullGrid(levelFlags), levelFlags);
            std::vector<GridBoundary> newGrids = this->mergeGrids(minimalGrids, levelFlags);
            this->ensureNesting(newGrids, parents);
            std::vector<AMRNode<SolverType, ProblemType, Fields, padding>> levelNodes =
                this->createGrids(newGrids, level + 1, this->amrNodes[level]);
            this->amrNodes.push_back(levelNodes);
            parents.swap(newGrids);
            level++;
        }
    }

    return this->amrNodes[0][0];
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::refine() {
    std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>> newNodes;
    newNodes.push_back(std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>());
    newNodes[0].swap(this->amrNodes[0]);
    newNodes[0][0].removeChildren();

    std::vector<CellFlags::Flags> flags = this->flagCells();
    newNodes.reserve(flags.size() + 1);

    std::vector<GridBoundary> parents({ this->amrNodes[0][0].gridBoundary });
    for (int level = flags.size() - 1; level >= 0; level--) {
        CellFlags::Flags& levelFlags = flags[level];
        std::vector<GridBoundary> minimalGrids =
            this->divideGrid(this->fullGrid(levelFlags), levelFlags);
        std::vector<GridBoundary> newGrids = this->mergeGrids(minimalGrids, levelFlags);
        this->ensureNesting(newGrids, parents);
        std::vector<AMRNode<SolverType, ProblemType, Fields, padding>> levelNodes =
            this->createGrids(newGrids, flags.size() - level, newNodes[flags.size() - level - 1]);
        newNodes.push_back(levelNodes);
        parents.swap(newGrids);
    }

    this->amrNodes.swap(newNodes);
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<CellFlags::Flags>
Refinery<SolverType, ProblemType, Fields, padding>::flagCells() const {
    const int maxLevel = static_cast<int>(this->amrNodes.size());
    std::vector<CellFlags::Flags> flags;
    bool refine = false;
    for (int level = maxLevel - 1; level >= 0; level--) {
        bool first = true;
        CellFlags::Flags levelFlags;
        for (const AMRNode<SolverType, ProblemType, Fields, padding>& node :
             this->amrNodes[level]) {
            std::optional<CellFlags::Flags> nodeFlags = {}; // truncation error stuff
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
            if (level < maxLevel - 1) {
                CellFlags::Flags higherLevelFlags =
                    CellFlags::translateFlags(flags.back(), this->configuration.refinementFactor);
                CellFlags::addFlags(levelFlags, higherLevelFlags);
            }
            CellFlags::addBuffer(levelFlags, this->configuration.bufferSize);
            flags.push_back(levelFlags);
        }
    }
    return flags;
}

inline double gridEfficiency(const GridBoundary& grid, const CellFlags::Flags& flags) {
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
    const double efficiency = gridEfficiency(grid, flags);
    if (efficiency == 0.0) {
        return std::vector<GridBoundary>(0);
    } else if (efficiency < this->configuration.efficiencyThreshold) {
        GridBoundary subGrid1;
        GridBoundary subGrid2;
        unsigned longest =
            std::max(grid[0].second - grid[0].first,
                     std::max(grid[1].second - grid[1].first, grid[2].second - grid[2].first));
        // TODO: What of very small grids?
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

inline double gridCost(const GridBoundary& grid, const CellFlags::Flags& flags) {
    const double x = static_cast<double>(grid[0].second - grid[0].first);
    const double y = static_cast<double>(grid[1].second - grid[1].first);
    const double z = static_cast<double>(grid[2].second - grid[2].first);
    const double efficiency = gridEfficiency(grid, flags);
    return (0.0 - efficiency) * (x * y * z + x * y + x * z + y * z + x + y + z);
}

inline bool isOverlappingOrAdjaecent(const GridBoundary& grid1, const GridBoundary& grid2) {
    return grid1[0].second + 1 >= grid2[0].first && grid1[0].first <= grid2[0].second + 1 &&
           grid1[1].second + 1 >= grid2[1].first && grid1[1].first <= grid2[1].second + 1 &&
           grid1[2].second + 1 >= grid2[2].first && grid1[2].first <= grid2[2].second + 1;
}

inline std::optional<GridBoundary> overlappingArea(const GridBoundary& grid1,
                                                   const GridBoundary& grid2) {
    if (grid1[0].second >= grid2[0].first && grid1[0].first <= grid2[0].second &&
        grid1[1].second >= grid2[1].first && grid1[1].first <= grid2[1].second &&
        grid1[2].second >= grid2[2].first && grid1[2].first <= grid2[2].second) {
        unsigned xStart = std::max(grid1[0].first, grid2[0].first);
        unsigned xEnd = std::min(grid1[0].second, grid2[0].second);
        unsigned yStart = std::max(grid1[1].first, grid2[1].first);
        unsigned yEnd = std::min(grid1[1].second, grid2[1].second);
        unsigned zStart = std::max(grid1[2].first, grid2[2].first);
        unsigned zEnd = std::min(grid1[2].second, grid2[2].second);
        const GridBoundary overlapping = { std::make_pair(xStart, xEnd),
                                           std::make_pair(yStart, yEnd),
                                           std::make_pair(zStart, zEnd) };
        return std::make_optional(overlapping);
    } else {
        return std::nullopt;
    }
}

inline std::optional<GridBoundary> tryMerge(const GridBoundary& grid1, const GridBoundary& grid2,
                                            const CellFlags::Flags& flags) {
    if (isOverlappingOrAdjaecent(grid1, grid2)) {
        const GridBoundary mergedGrid = {
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

inline std::vector<GridBoundary> mergeGridsAux(const std::vector<GridBoundary>& gridsKnown,
                                               const std::vector<GridBoundary>& gridsNew,
                                               const CellFlags::Flags& flags) {
    std::vector<bool> keepGrid(gridsNew.size(), true);
    std::vector<GridBoundary> newGrids;
    std::vector<GridBoundary> knownGrids;
    for (const GridBoundary& grid1 : gridsKnown) {
        bool keep = true;
        for (unsigned grid2 = 0; grid2 < newGrids.size(); grid2++) {
            if (keepGrid[grid2]) {
                std::optional<GridBoundary> mergedGrid = tryMerge(grid1, newGrids[grid2], flags);
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
std::vector<GridBoundary> Refinery<SolverType, ProblemType, Fields, padding>::mergeGrids(
    const std::vector<GridBoundary>& grids, const CellFlags::Flags& flags) const {
    std::vector<GridBoundary> newGrids;
    std::vector<bool> keepGrid(grids.size(), true);
    for (unsigned grid1 = 0; grid1 < grids.size() - 1; grid1++) {
        if (keepGrid[grid1]) {
            for (unsigned grid2 = grid1 + 1; grid2 < grids.size(); grid2++) {
                if (keepGrid[grid2]) {
                    std::optional<GridBoundary> mergedGrid =
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
        std::vector<GridBoundary> knownGrids;
        for (unsigned grid = 0; grid < grids.size(); grid++) {
            if (keepGrid[grid]) {
                knownGrids.push_back(grids[grid]);
            }
        }
        return mergeGridsAux(knownGrids, newGrids, flags);
    }
}

inline GridBoundary translateUp(const GridBoundary& grid, const unsigned factor) {
    return { std::make_pair(grid[0].first * factor, grid[0].second * factor + factor - 1),
             std::make_pair(grid[1].first * factor, grid[1].second * factor + factor - 1),
             std::make_pair(grid[2].first * factor, grid[2].second * factor + factor - 1) };
}

inline std::vector<GridBoundary> subtract(const GridBoundary& grid,
                                          const GridBoundary& subtractor) {
    if (grid[0].first > subtractor[0].second || grid[0].second < subtractor[0].first ||
        grid[1].first > subtractor[1].second || grid[1].second < subtractor[1].first ||
        grid[2].first > subtractor[2].second || grid[2].second < subtractor[2].first) {
        return std::vector<GridBoundary>(1, grid);
    }
    GridBoundary remaining = grid;
    std::vector<GridBoundary> result;
    if (remaining[0].first < subtractor[0].first) {
        result.push_back({ std::make_pair(remaining[0].first, subtractor[0].first - 1),
                           remaining[1], remaining[2] });
        remaining[0].first = subtractor[0].first;
    }
    if (remaining[0].second > subtractor[0].second) {
        result.push_back({ std::make_pair(subtractor[0].second + 1, remaining[0].second),
                           remaining[1], remaining[2] });
        remaining[0].second = subtractor[0].second;
    }
    if (remaining[1].first < subtractor[1].first) {
        result.push_back({ remaining[0],
                           std::make_pair(remaining[1].first, subtractor[1].first - 1),
                           remaining[2] });
        remaining[1].first = subtractor[1].first;
    }
    if (remaining[1].second > subtractor[1].second) {
        result.push_back({ remaining[0],
                           std::make_pair(subtractor[1].second + 1, remaining[1].second),
                           remaining[2] });
        remaining[1].second = subtractor[1].second;
    }
    if (remaining[2].first < subtractor[2].first) {
        result.push_back({ remaining[0], remaining[1],
                           std::make_pair(remaining[2].first, subtractor[2].first - 1) });
        remaining[2].first = subtractor[2].first;
    }
    if (remaining[2].second > subtractor[2].second) {
        result.push_back({ remaining[0], remaining[1],
                           std::make_pair(subtractor[2].second + 1, remaining[2].second) });
        remaining[2].second = subtractor[1].second;
    }
    return result;
}

inline void strictMerge(std::vector<GridBoundary>& grids) {
    std::vector<bool> keepGrid(grids.size(), true);
    std::vector<GridBoundary> newGrids;
    for (unsigned grid1 = 0; grid1 + 1 < grids.size(); grid1++) {
        if (keepGrid[grid1]) {
            const GridBoundary& g1 = grids[grid1];
            for (unsigned grid2 = grid1 + 1; grid2 < grids.size(); grid2++) {
                if (keepGrid[grid2]) {
                    const GridBoundary& g2 = grids[grid2];
                    if (std::max(g1[0].first, g2[0].first) ==
                            std::min(g1[0].second, g2[0].second) + 1 &&
                        g1[1].first == g2[1].first && g1[1].second == g2[1].second &&
                        g1[2].first == g2[2].first && g1[2].second == g2[2].second) {
                        keepGrid[grid1] = false;
                        keepGrid[grid2] = false;
                        newGrids.push_back({ std::make_pair(std::min(g1[0].first, g2[0].first),
                                                            std::max(g2[0].second, g2[0].second)),
                                             g1[1], g1[2] });
                    } else if (std::max(g1[1].first, g2[1].first) ==
                                   std::min(g1[1].second, g2[1].second) + 1 &&
                               g1[0].first == g2[0].first && g1[0].second == g2[0].second &&
                               g1[2].first == g2[2].first && g1[2].second == g2[2].second) {
                        keepGrid[grid1] = false;
                        keepGrid[grid2] = false;
                        newGrids.push_back({ g1[0],
                                             std::make_pair(std::min(g1[1].first, g2[1].first),
                                                            std::max(g2[1].second, g2[1].second)),
                                             g1[2] });
                    } else if (std::max(g1[2].first, g2[2].first) ==
                                   std::min(g1[2].second, g2[2].second) + 1 &&
                               g1[0].first == g2[0].first && g1[0].second == g2[0].second &&
                               g1[1].first == g2[1].first && g1[1].second == g2[1].second) {
                        keepGrid[grid1] = false;
                        keepGrid[grid2] = false;
                        newGrids.push_back({
                            g1[0],
                            g1[1],
                            std::make_pair(std::min(g1[2].first, g2[2].first),
                                           std::max(g2[2].second, g2[2].second)),
                        });
                    }
                }
            }
        }
    }
    if (newGrids.size() > 0) {
        for (unsigned grid = 0; grid < grids.size(); grid++) {
            if (keepGrid[grid]) {
                newGrids.push_back(grids[grid]);
            }
        }
        strictMerge(newGrids);
        grids.swap(newGrids);
    }
}

inline std::vector<GridBoundary> remainingGrids(const GridBoundary& grid,
                                                const std::vector<GridBoundary>& removed) {
    std::vector<GridBoundary> remaining;
    remaining.push_back(grid);
    for (const GridBoundary& remove : removed) {
        std::vector<GridBoundary> stillRemaining;
        for (const GridBoundary& left : remaining) {
            std::vector<GridBoundary> keep = subtract(left, remove);
            stillRemaining.insert(stillRemaining.end(), keep.begin(), keep.end());
        }
        remaining.swap(stillRemaining);
    }
    strictMerge(remaining);
    return remaining;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::ensureNesting(
    std::vector<GridBoundary>& grids, const std::vector<GridBoundary>& parentGrids) {
    std::vector<GridBoundary> translatedParents;
    translatedParents.reserve(parentGrids.size());
    std::transform(parentGrids.begin(), parentGrids.end(), std::back_inserter(translatedParents),
                   [factor = this->configuration.refinementFactor](const GridBoundary& parent) {
                       return translateUp(parent, factor);
                   });
    std::vector<GridBoundary> newGrids;
    for (GridBoundary& grid : grids) {
        std::vector<GridBoundary> notCovered;
        notCovered.push_back(grid);
        for (const GridBoundary& parent : translatedParents) {
            std::optional<GridBoundary> overlapping = overlappingArea(grid, parent);
            if (overlapping.has_value()) {
                GridBoundary overlap = overlapping.value();
                std::vector<GridBoundary> newParts;
                for (GridBoundary& part : notCovered) {
                    std::vector<GridBoundary> leftOver = subtract(part, overlap);
                    newParts.insert(newParts.end(), leftOver.begin(), leftOver.end());
                }
                notCovered.swap(newParts);
                if (notCovered.size() == 0) {
                    break;
                }
            }
        }
        if (notCovered.size() > 0) {
            std::vector<GridBoundary> remaining = remainingGrids(grid, notCovered);
            grid = remaining[0];
            newGrids.insert(newGrids.end(), remaining.begin() + 1, remaining.end());
        }
    }
    grids.insert(grids.end(), newGrids.begin(), newGrids.end());
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>
Refinery<SolverType, ProblemType, Fields, padding>::createGrids(
    std::vector<GridBoundary>& grids, const unsigned level,
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>& parents) {

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
        this->initialiseGrid(newGrid, grid, level, parents);
        AMRNode<SolverType, ProblemType, Fields, padding> newNode(
            newGrid, this->problem, this->configuration, timeDelta, timeCurrent);
        newNode.gridBoundary = grid;
        newGrids.push_back(newNode);
    }

    for (unsigned node1 = 0; node1 < newGrids.size(); node1++) {
        AMRNode<SolverType, ProblemType, Fields, padding>& grid1 = newGrids[node1];
        for (unsigned node2 = node1 + 1; node2 < newGrids.size(); node2++) {
            AMRNode<SolverType, ProblemType, Fields, padding>& grid2 = newGrids[node1];
            if (isOverlappingOrAdjaecent(grid1.gridBoundary, grid2.gridBoundary)) {
                grid1.addSibling(&grid2);
                grid2.addSibling(&grid1);
            }
        }
        for (AMRNode<SolverType, ProblemType, Fields, padding>& parent : parents) {
            if (isOverlappingOrAdjaecent(
                    grid1.gridBoundary,
                    translateUp(parent.gridBoundary, this->configuration.refinementFactor))) {
                grid1.addParent(&parent);
                parent.addChild(&grid1);
            }
        }
    }

    return newGrids;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::initialiseGrid(
    PaddedGrid<Fields, padding>& grid, const GridBoundary& gridBoundary, const unsigned level,
    const std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>& parents) {
    // Could be optimised with smaller search space for siblings and parents
    const double refinementFactor = static_cast<double>(this->configuration.refinementFactor);
    for (unsigned x = 0; x < grid.xDim(); x++) {
        const int xGlobal = x + gridBoundary[0].first - padding;
        const double xGlobalUp = static_cast<double>(xGlobal) / refinementFactor;
        for (unsigned y = 0; y < grid.yDim(); y++) {
            const int yGlobal = y + gridBoundary[1].first - padding;
            const double yGlobalUp = static_cast<double>(yGlobal) / refinementFactor;
            for (unsigned z = 0; z < grid.zDim(); z++) {
                const int zGlobal = z + gridBoundary[2].first - padding;
                bool foundCell = false;
                if (level < this->amrNodes.size()) {
                    for (AMRNode<SolverType, ProblemType, Fields, padding>& sibling :
                         this->amrNodes[level]) {
                        if (xGlobal >= static_cast<int>(sibling.gridBoundary[0].first) &&
                            xGlobal < static_cast<int>(sibling.gridBoundary[0].second) &&
                            yGlobal >= static_cast<int>(sibling.gridBoundary[1].first) &&
                            yGlobal < static_cast<int>(sibling.gridBoundary[1].second) &&
                            zGlobal >= static_cast<int>(sibling.gridBoundary[2].first) &&
                            zGlobal < static_cast<int>(sibling.gridBoundary[2].second)) {
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
                    const double zGlobalUp = static_cast<double>(zGlobal) / refinementFactor;
                    for (const AMRNode<SolverType, ProblemType, Fields, padding>& parent :
                         parents) {
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
