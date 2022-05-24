#pragma once

#include <limits>

#include "amr-node.h"
#include "flags.h"

template <class SolverType, class ProblemType, class Fields, unsigned padding> class Refinery {
  public:
    Refinery(const AMRParameters& configuration);

    void refine();
    AMRNode<SolverType, ProblemType, Fields, padding>&
    initialRefine(PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem);

  private:
    const AMRParameters& configuration;

    std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>> amrNodes;

    std::vector<CellFlags::Flags> flagCells();
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
Refinery<SolverType, ProblemType, Fields, padding>::Refinery(const AMRParameters& configuration)
    : configuration(configuration) {}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::refine() {
    std::vector<CellFlags::Flags> flags = this->flagCells();
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<CellFlags::Flags> Refinery<SolverType, ProblemType, Fields, padding>::flagCells() {
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
AMRNode<SolverType, ProblemType, Fields, padding>&
Refinery<SolverType, ProblemType, Fields, padding>::initialRefine(
    PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem) {
    AMRNode<SolverType, ProblemType, Fields, padding> root(grid, problem, this->configuration,
                                                           problem.timeDelta, 0.0);
    this->amrNodes.push_back(
        std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>(1, root));
    return this->amrNodes[0][0];
}
