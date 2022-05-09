#pragma once

#include <optional>

#include "../configuration/problem.h"
#include "../grid/padded-grid.h"

template <class Specific, class Fields, class ProblemType, unsigned ghostCells> class Solver {
  public:
    void solve(const std::optional<double> untilTime);

  protected:
    Solver(PaddedGrid<Fields, ghostCells>& grid, const Problem<ProblemType>& problem);

    PaddedGrid<Fields, ghostCells>& grid;
    const Problem<ProblemType>& problem;
    double timeDelta;
    double timeCurrent;
    double timeEnd;
    unsigned timeStep;
};

template <class Specific, class Fields, class ProblemType, unsigned ghostCells>
Solver<Specific, Fields, ProblemType, ghostCells>::Solver(PaddedGrid<Fields, ghostCells>& grid,
                                                          const Problem<ProblemType>& problem)
    : grid(grid), problem(problem) {
    this->timeDelta = problem.timeDelta;
    this->timeCurrent = problem.timeStart;
    this->timeEnd = problem.timeEnd;
    this->timeStep = 0;
}

template <class Specific, class Fields, class ProblemType, unsigned ghostCells>
void Solver<Specific, Fields, ProblemType, ghostCells>::solve(
    const std::optional<double> untilTime) {
    this->timeEnd = untilTime.value_or(this->problem.timeEnd);
    static_cast<Specific*>(this)->integrate();
}
