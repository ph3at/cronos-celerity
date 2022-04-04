#pragma once

#include <array>

#include "../field-wrapper/face-values.h"
#include "../field-wrapper/phys-values.h"
#include "../parameters/problem.h"

namespace Transformation {
void reconstToConservatives(PhysValues& output, const FieldStruct& reconstructions,
                            const Problem& problem);

double computeThermalPressure(const FieldStruct& fields, const Problem& problem);
}; // namespace Transformation
