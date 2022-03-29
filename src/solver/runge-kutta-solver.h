#pragma once

#include "base-solver.h"
#include <algorithm>

template <class T> class RungeKuttaSolver : public Solver<T> {
  public:
    RungeKuttaSolver(Grid<T>& grid, const double timeDelta, const double timeStart,
                     const double timeEnd, unsigned rungeKuttaSteps);

    void solve();

  private:
    double timeCurrent;
    unsigned timeStep;
    double cfl;
    const unsigned rungeKuttaSteps = 2;

    void computeStep();
    double computeSubstep(unsigned substep);
    void advanceTime();
    void adjustTimeDelta();
    bool isFinished() const;
};

template <class T>
RungeKuttaSolver<T>::RungeKuttaSolver(Grid<T>& grid, const double timeDelta, const double timeStart,
                                      const double timeEnd, const unsigned rungeKuttaSteps)
    : Solver<T>(grid, timeDelta, timeStart, timeEnd), rungeKuttaSteps(rungeKuttaSteps) {
    this->timeCurrent = timeStart;
    this->timeStep = 0;
};

template <class T> void RungeKuttaSolver<T>::solve() {
    while (!this->isFinished()) {
        this->computeStep();
        this->advanceTime();
        this->adjustTimeDelta();
    }
}

template <class T> void RungeKuttaSolver<T>::computeStep() {
    this->cfl = 0.0;
    for (unsigned substep = 0; substep < this->rungeKuttaSteps; substep++) {
        double substepCFL = this->computeSubstep(substep);
        this->cfl = std::max(this->cfl, substepCFL);
    }
}

template <class T> double RungeKuttaSolver<T>::computeSubstep(const unsigned substep) {

    return 0.0;
}

template <class T> void RungeKuttaSolver<T>::advanceTime() {
    this->timeCurrent += this->timeDelta;
    this->timeStep++;
}

template <class T> void RungeKuttaSolver<T>::adjustTimeDelta() {}

template <class T> bool RungeKuttaSolver<T>::isFinished() const {
    return this->timeCurrent >= this->timeEnd;
}
