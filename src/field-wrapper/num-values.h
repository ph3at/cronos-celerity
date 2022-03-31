#pragma once

#include "phys-fields.h"

typedef struct numVals {
    // TODO: better names
    double pTotal;
    double v_ch_p;
    double v_ch_m;
    std::array<double, PHYSICAL_FIELDS> fluxes;
} NumValues;
