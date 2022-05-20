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
    bool isFinished() const;
    void report() const;

    double timeDelta;
    double timeCurrent;
    unsigned timeStep;
    double timeEnd;

  protected:
    Solver(PaddedGrid<Fields, ghostCells>& grid, const Problem<ProblemType>& problem);

    PaddedGrid<Fields, ghostCells>& grid;
    const Problem<ProblemType>& problem;

  private:
    void doSolve();
};

template <class Specific, class Fields, class ProblemType, unsigned ghostCells>
Solver<Specific, Fields, ProblemType, ghostCells>::Solver(PaddedGrid<Fields, ghostCells>& grid,
                                                          const Problem<ProblemType>& problem)
    : grid(grid), problem(problem) {
    this->timeDelta = problem.timeDelta;
    this->timeCurrent = problem.timeStart;
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
        static_cast<Specific*>(this)->adjust();
        this->timeStep++;
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
    static_cast<Specific*>(this)->init();
}

template <class Specific, class ProblemType, class Fields, unsigned padding>
bool Solver<Specific, ProblemType, Fields, padding>::isFinished() const {
    return this->timeCurrent >= this->timeEnd;
}

template <class Specific, class Fields, class ProblemType, unsigned ghostCells>
void Solver<Specific, Fields, ProblemType, ghostCells>::report() const {
    std::cout << "Stopped at time " << this->timeCurrent << ", after " << this->timeStep
              << " steps." << std::endl;
}
