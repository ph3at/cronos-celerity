#include "riemann-solver.h"

namespace RiemannSolver {
void computeFluxes(PhysValues& fields, const FieldStruct& reconstruction, const unsigned int face) {
    unsigned dir = face / 2;

    fields.fluxes[FieldNames::DENSITY] =
        reconstruction[FieldNames::VELOCITY_X + dir] * reconstruction[FieldNames::DENSITY];
    for (unsigned fluxDir = FieldNames::VELOCITY_X; fluxDir <= FieldNames::VELOCITY_Z; fluxDir++) {
        if (fluxDir == dir) {
            fields.fluxes[fluxDir] =
                fields.conservatives[fluxDir] * reconstruction[fluxDir] + fields.thermalPressure;
        } else {
            fields.fluxes[fluxDir] =
                fields.conservatives[fluxDir] * reconstruction[FieldNames::VELOCITY_X + dir];
        }
    }
    fields.fluxes[FieldNames::THERMAL_ENERGY] =
        (fields.conservatives[FieldNames::THERMAL_ENERGY] + fields.thermalPressure) *
        reconstruction[FieldNames::VELOCITY_X + dir];
}
}; // namespace RiemannSolver
