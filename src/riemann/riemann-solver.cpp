#include <cmath>

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

std::pair<double, double> characteristicVelocity(const PhysValues& physValsLeft,
                                                 const PhysValues& physValsRight,
                                                 const FieldStruct& reconstLeft,
                                                 const FieldStruct& reconstRight,
                                                 const Problem& problem, const unsigned dir) {
    double flowVelocityLeft = reconstLeft[FieldNames::VELOCITY_X + dir];
    double flowVelocityRight = reconstRight[FieldNames::VELOCITY_X + dir];

    double soundVelocityLeft = std::sqrt(problem.gamma * physValsLeft.thermalPressure /
                                         physValsLeft.conservatives[FieldNames::DENSITY]);
    double soundVelocityRight = std::sqrt(problem.gamma * physValsRight.thermalPressure /
                                          physValsRight.conservatives[FieldNames::DENSITY]);

    double characVelLeft = std::max(
        std::max(soundVelocityLeft + flowVelocityLeft, soundVelocityRight + flowVelocityRight),
        0.0);
    double characVelRight = std::max(
        std::max(soundVelocityLeft - flowVelocityLeft, soundVelocityRight - flowVelocityRight),
        0.0);

    return std::make_pair(characVelLeft, characVelRight);
}

NumValues numericalFlux(const std::pair<double, double>& characVelocities,
                        const PhysValues& physValsLeft, const PhysValues& physValsRight,
                        const unsigned dir) {
    return {};
}
}; // namespace RiemannSolver
