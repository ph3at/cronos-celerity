#pragma once

#include <array>

#include "../field-wrapper/face-values.h"
#include "../parameters/problem.h"

namespace Transformation {
PerFaceValues reconstToConservatives(const PerFaceValues& reconstructions, const Problem& problem);

PerFaceSingleValue computeThermalPressure(const PerFaceValues& fields, const Problem& problem);
}; // namespace Transformation
