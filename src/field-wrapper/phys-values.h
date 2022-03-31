#pragma once

#include "phys-fields.h"

typedef struct physVals {
    // TODO: Better names
    double pTotal;
    double pTherm;
    std::array<double, PHYSICAL_FIELDS> primitives;
    std::array<double, PHYSICAL_FIELDS> conservatives;
    std::array<double, PHYSICAL_FIELDS> fluxes;
} PhysValues;
