#pragma once

#include "phys-fields.h"

typedef struct physVals {
    double pressureTotal;
    double pressureThermal;
    std::array<double, NUM_PHYSICAL_FIELDS> primitives;
    std::array<double, NUM_PHYSICAL_FIELDS> conservatives;
    std::array<double, NUM_PHYSICAL_FIELDS> fluxes;
} PhysValues;
