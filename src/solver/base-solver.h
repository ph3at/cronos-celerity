#pragma once

#include <iostream>

#include "../boundary/boundary.h"
#include "../configuration/problem.h"
#include "../grid/padded-grid.h"

template <class Specific, class Fields, class ProblemType, unsigned ghostCells> class Solver {
  public:
    void solve();
    void solve(const double untilTime);
    void step();
    void initialise();
    void report() const;

  protected:
    Solver(PaddedGrid<Fields, ghostCells>& grid, const Problem<ProblemType>& problem);

    PaddedGrid<Fields, ghostCells>& grid;
    const Problem<ProblemType>& problem;
    double timeDelta;
    double timeCurrent;
    double timeEnd;
    unsigned timeStep;

  private:
    void doSolve();
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
void Solver<Specific, Fields, ProblemType, ghostCells>::solve() {
    this->timeEnd = this->problem.timeEnd;
    this->doSolve();
}

template <class Specific, class Fields, class ProblemType, unsigned ghostCells>
void Solver<Specific, Fields, ProblemType, ghostCells>::solve(const double untilTime) {
    this->timeEnd = untilTime;
    this->doSolve();
}

template <class Specific, class Fields, class ProblemType, unsigned ghostCells>
void Solver<Specific, Fields, ProblemType, ghostCells>::doSolve() {
    while (!static_cast<Specific*>(this)->isFinished()) {
        static_cast<Specific*>(this)->singleStep();
    }
}

template <class Specific, class Fields, class ProblemType, unsigned ghostCells>
void Solver<Specific, Fields, ProblemType, ghostCells>::step() {
    static_cast<Specific*>(this)->singleStep();
}

template <class Specific, class Fields, class ProblemType, unsigned ghostCells>
void Solver<Specific, Fields, ProblemType, ghostCells>::initialise() {
    this->problem.initialiseGrid(this->grid);
    Boundary::applyAll(this->grid, this->problem);
}

template <class Specific, class Fields, class ProblemType, unsigned ghostCells>
void Solver<Specific, Fields, ProblemType, ghostCells>::report() const {
    std::cout << "Stopped at time " << this->timeCurrent << ", after " << this->timeStep
              << " steps." << std::endl;
}
