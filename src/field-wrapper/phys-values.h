#pragma once

#include "phys-fields.h"

typedef struct physVals {
    double thermalPressure;
    FieldStruct conservatives;
    FieldStruct fluxes;
} PhysValues;
