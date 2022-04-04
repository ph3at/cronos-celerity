#pragma once

#include "phys-fields.h"

typedef struct numValues {
    double pressureTotal;
    std::array<double, NUM_PHYSICAL_FIELDS> fluxes;
} NumValues;
