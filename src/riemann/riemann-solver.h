#pragma once

#include "../field-wrapper/face-values.h"
#include "../field-wrapper/phys-values.h"
#include "../parameters/problem.h"

namespace RiemannSolver {
void computeFluxes(PhysValues& fields, const FieldStruct& reconstruction, const unsigned face);

double computeThermalPressure(const FieldStruct& fields, const Problem& problem);
}; // namespace RiemannSolver
