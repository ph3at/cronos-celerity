#pragma once

#include <iostream>

#include "../boundary/boundary.h"
#include "../configuration/problem.h"
#include "../grid/padded-grid.h"

/* Specifics must implement "init()" for initialising themselves, "singleStep()" which computes one
 * full time step, "adjustConfig()", which adjusts the state of the solver between time steps, and
 * "finaliseResult()", which does post-processing. */
template <class Specific, class Fields, class ProblemType, unsigned padding> class Solver {
  public:
    void solve();
    void solve(const double untilTime);
    void step();
    void initialise();
    void adjust();
    bool isFinished() const;
    void report() const;
    void finalise();

    PaddedGrid<Fields, padding> grid;

    double timeDelta;
    double timeCurrent;
    unsigned timeStep;
    double timeEnd;

    const bool doOutput;

  protected:
    Solver(const ProblemType& problem, const bool doOutput = false);
    Solver(const PaddedGrid<Fields, padding>& grid, const ProblemType& problem,
           const bool doOutput = false);
    const ProblemType& problem;

  private:
    void doSolve();
};

template <class Specific, class Fields, class ProblemType, unsigned padding>
Solver<Specific, Fields, ProblemType, padding>::Solver(const ProblemType& problem,
                                                       const bool doOutput)
    : grid({}, problem.numberCells[Direction::DirX], problem.numberCells[Direction::DirY],
           problem.numberCells[Direction::DirZ]),
      doOutput(doOutput), problem(problem) {
    this->timeEnd = problem.timeEnd;
    this->timeDelta = problem.timeDelta;
    this->timeCurrent = problem.timeStart;
    this->timeStep = 0;
}

template <class Specific, class Fields, class ProblemType, unsigned padding>
Solver<Specific, Fields, ProblemType, padding>::Solver(const PaddedGrid<Fields, padding>& grid,
                                                       const ProblemType& problem,
                                                       const bool doOutput)
    : grid(grid), doOutput(doOutput), problem(problem) {
    this->timeEnd = problem.timeEnd;
    this->timeDelta = problem.timeDelta;
    this->timeCurrent = problem.timeStart;
    this->timeStep = 0;
}

template <class Specific, class Fields, class ProblemType, unsigned padding>
void Solver<Specific, Fields, ProblemType, padding>::solve() {
    this->timeEnd = this->problem.timeEnd;
    this->doSolve();
}

template <class Specific, class Fields, class ProblemType, unsigned padding>
void Solver<Specific, Fields, ProblemType, padding>::solve(const double untilTime) {
    this->timeEnd = untilTime;
    this->doSolve();
}

template <class Specific, class Fields, class ProblemType, unsigned padding>
void Solver<Specific, Fields, ProblemType, padding>::doSolve() {
    while (!static_cast<Specific*>(this)->isFinished()) {
        this->step();
        this->adjust();
        this->timeStep++;
    }
    this->finalise();
}

template <class Specific, class Fields, class ProblemType, unsigned padding>
void Solver<Specific, Fields, ProblemType, padding>::step() {
    if (doOutput) {
        std::cout << "--- Starting time step " << this->timeStep + 1 << " at time "
                  << this->timeCurrent << " with time delta " << this->timeDelta << " ---"
                  << std::endl;
    }
    static_cast<Specific*>(this)->singleStep();
    this->timeCurrent += this->timeDelta;
}

template <class Specific, class Fields, class ProblemType, unsigned padding>
void Solver<Specific, Fields, ProblemType, padding>::initialise() {
    static_cast<Specific*>(this)->init();
}

template <class Specific, class Fields, class ProblemType, unsigned padding>
void Solver<Specific, Fields, ProblemType, padding>::adjust() {
    static_cast<Specific*>(this)->adjustConfig();
}

template <class Specific, class ProblemType, class Fields, unsigned padding>
bool Solver<Specific, ProblemType, Fields, padding>::isFinished() const {
    return this->timeCurrent >= this->timeEnd;
}

template <class Specific, class Fields, class ProblemType, unsigned padding>
void Solver<Specific, Fields, ProblemType, padding>::report() const {
    std::cout << "Stopped at time " << this->timeCurrent << ", after " << this->timeStep
              << " steps." << std::endl;
}

template <class Specific, class Fields, class ProblemType, unsigned padding>
void Solver<Specific, Fields, ProblemType, padding>::finalise() {
    static_cast<Specific*>(this)->finaliseResult();
}
