#pragma once

#include <array>

#include "../field-wrapper/face-values.h"
#include "../parameters/problem.h"

namespace Transformation {
PerFaceValues reconstToConservatives(const PerFaceValues& reconstructions, const Problem& problem);

double thermalEnergyToEnergy(const FieldStruct& fields);

double temperatureToEnergy(const FieldStruct& fields, const double gamma);
}; // namespace Transformation
