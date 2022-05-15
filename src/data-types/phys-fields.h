#pragma once

#include <array>
#include <cstddef>

constexpr unsigned NUM_PHYSICAL_FIELDS = 5;
typedef std::array<double, NUM_PHYSICAL_FIELDS> FieldStruct;

enum FieldNames { DENSITY = 0, VELOCITY_X = 1, VELOCITY_Y = 2, VELOCITY_Z = 3, THERMAL_ENERGY = 4 };
