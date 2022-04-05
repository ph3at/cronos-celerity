#pragma once

#include "../field-wrapper/fields.h"
#include "../parameters/problem.h"

namespace RiemannSolver {
void computeFluxes(PhysValues& fields, const FieldStruct& reconstruction, const unsigned face);

std::pair<double, double> characteristicVelocity(const PhysValues& physValsLeft,
                                                 const PhysValues& physValsRight,
                                                 const FieldStruct& reconstLeft,
                                                 const FieldStruct& reconstRight,
                                                 const Problem& problem, const unsigned dir);

NumValues numericalFlux(const std::pair<double, double>& characVelocities,
                        const PhysValues& physValsLeft, const PhysValues& physValsRight,
                        const unsigned dir);
}; // namespace RiemannSolver
