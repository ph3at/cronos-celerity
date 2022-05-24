#pragma once

#include <limits>

#include "amr-node.h"
#include "flags.h"

typedef std::array<std::pair<unsigned, unsigned>, 3> GridBoundary;

template <class SolverType, class ProblemType, class Fields, unsigned padding> class Refinery {
  public:
    Refinery(const AMRParameters& configuration);

    void refine();
    AMRNode<SolverType, ProblemType, Fields, padding>&
    initialRefine(PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem);

  private:
    const AMRParameters& configuration;

    std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>> amrNodes;

    std::vector<CellFlags::Flags> flagCells() const;
    double efficiency(const GridBoundary& grid, const CellFlags::Flags& flags) const;
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
Refinery<SolverType, ProblemType, Fields, padding>::Refinery(const AMRParameters& configuration)
    : configuration(configuration) {}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::refine() {
    std::vector<CellFlags::Flags> flags = this->flagCells();
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
AMRNode<SolverType, ProblemType, Fields, padding>&
Refinery<SolverType, ProblemType, Fields, padding>::initialRefine(
    PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem) {
    AMRNode<SolverType, ProblemType, Fields, padding> root(grid, problem, this->configuration,
                                                           problem.timeDelta, 0.0);
    this->amrNodes.push_back(
        std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>(1, root));
    return this->amrNodes[0][0];
}
