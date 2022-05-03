#pragma once

#include "../field-wrapper/fields.h"

namespace RiemannSolver {
void computeFluxes(PhysValues& fields, const FieldStruct& reconstruction, const unsigned face);

std::pair<double, double> characteristicVelocity(const PhysValues& physValsLeft,
                                                 const PhysValues& physValsRight,
                                                 const FieldStruct& reconstLeft,
                                                 const FieldStruct& reconstRight,
                                                 const double gamma, const unsigned dir);

FieldStruct numericalFlux(const std::pair<double, double>& characteristicVelocities,
                          const PhysValues& physValsLeft, const PhysValues& physValsRight,
                          const FieldStruct& reconstructionLeft,
                          const FieldStruct& reconstructionRight, const unsigned dir);
}; // namespace RiemannSolver
