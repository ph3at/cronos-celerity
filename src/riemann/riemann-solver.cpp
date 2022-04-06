#include <cmath>
#include <tuple>

#include "../misc/direction.h"
#include "riemann-solver.h"

constexpr double HLLCSOLVER_HYDRO_VEPS = 1.e-120;

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

FieldStruct numericalFluxOutsideRiemannFan(const PhysValues& physValues) {
    FieldStruct numericalFluxes = {};
    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
        numericalFluxes[field] = physValues.fluxes[field];
    }
    return numericalFluxes;
}

std::tuple<unsigned, unsigned, unsigned> directionalIndices(unsigned dir) {
    switch (dir) {
        case Direction::DirX: {
            return std::make_tuple(FieldNames::VELOCITY_X, FieldNames::VELOCITY_Y,
                                   FieldNames::VELOCITY_Z);
        }
        case Direction::DirY: {
            return std::make_tuple(FieldNames::VELOCITY_Y, FieldNames::VELOCITY_Z,
                                   FieldNames::VELOCITY_X);
        }
        case Direction::DirZ: {
            return std::make_tuple(FieldNames::VELOCITY_Z, FieldNames::VELOCITY_X,
                                   FieldNames::VELOCITY_Y);
        }
        default: {
            /* Impossible */ return std::make_tuple(-1, -1, -1);
        }
    }
}

FieldStruct computeStarValues(const PhysValues& physValues, const FieldStruct& reconstruction,
                              const double characteristicVelocity,
                              const double relativeSignalVelocity,
                              const double intermediateWaveSpeed, const unsigned velocityParallel,
                              const unsigned velocityPerpendicular1,
                              const unsigned velocityPerpendicular2) {
    FieldStruct starValues = {};
    starValues[FieldNames::DENSITY] = physValues.conservatives[FieldNames::DENSITY] *
                                      relativeSignalVelocity /
                                      (characteristicVelocity - intermediateWaveSpeed);

    starValues[velocityParallel] = starValues[FieldNames::DENSITY] * intermediateWaveSpeed;
    starValues[velocityPerpendicular1] =
        starValues[FieldNames::DENSITY] * reconstruction[velocityPerpendicular1];
    starValues[velocityPerpendicular2] =
        starValues[FieldNames::DENSITY] * reconstruction[velocityPerpendicular2];

    starValues[FieldNames::THERMAL_ENERGY] =
        starValues[FieldNames::DENSITY] *
        (physValues.conservatives[FieldNames::THERMAL_ENERGY] /
             physValues.conservatives[FieldNames::DENSITY] +
         (intermediateWaveSpeed - reconstruction[velocityParallel]) *
             (intermediateWaveSpeed +
              physValues.thermalPressure /
                  (physValues.conservatives[FieldNames::DENSITY] * relativeSignalVelocity)));

    return starValues;
}

FieldStruct
computeNumericalFlux(const PhysValues& physValues, const FieldStruct& reconstruction,
                     const double characteristicVelocity, const double relativeSignalVelocity,
                     const double intermediateWaveSpeed, const unsigned velocityParallel,
                     const unsigned velocityPerpendicular1, const unsigned velocityPerpendicular2) {
    FieldStruct numericalFluxes = {};
    FieldStruct starValues = computeStarValues(
        physValues, reconstruction, characteristicVelocity, relativeSignalVelocity,
        intermediateWaveSpeed, velocityParallel, velocityPerpendicular1, velocityPerpendicular2);
    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
        numericalFluxes[field] =
            physValues.fluxes[field] +
            characteristicVelocity * (starValues[field] - physValues.conservatives[field]);
    }
    return numericalFluxes;
}

FieldStruct numericalFluxInsideRiemannFan(const std::pair<double, double>& characteristicVelocities,
                                          const PhysValues& physValsLeft,
                                          const PhysValues physValsRight,
                                          const FieldStruct& reconstructionLeft,
                                          const FieldStruct& reconstructionRight,
                                          const unsigned dir) {
    unsigned velocityParallel, velocityPerpendicular1, velocityPerpendicular2;
    std::tie(velocityParallel, velocityPerpendicular1, velocityPerpendicular2) =
        directionalIndices(dir);

    // Left and right have different meanings for reconstruction and hllc (this).
    // So here we use Minus and Plus instead of right and left.
    double densityMinus = physValsRight.conservatives[FieldNames::DENSITY];
    double densityPlus = physValsLeft.conservatives[FieldNames::DENSITY];

    double thermalPressureMinus = physValsRight.thermalPressure;
    double thermalPressurePlus = physValsLeft.thermalPressure;

    double velocityParallelMinus = reconstructionRight[velocityParallel];
    double velocityParallelPlus = reconstructionLeft[velocityParallel];

    double massFlowRateMinus = physValsRight.conservatives[velocityParallel];
    double massFlowRatePlus = physValsLeft.conservatives[velocityParallel];

    double characteristicVelocityMinus = -characteristicVelocities.second;
    double characteristicVelocityPlus = characteristicVelocities.first;

    double relativeSignalVelocityMinus = characteristicVelocityMinus - velocityParallelMinus;
    double relativeSignalVelocityPlus = characteristicVelocityPlus - velocityParallelPlus;

    double intermediateWaveSpeed =
        (thermalPressurePlus - thermalPressureMinus +
         massFlowRateMinus * relativeSignalVelocityMinus -
         massFlowRatePlus * relativeSignalVelocityPlus) /
        (densityMinus * relativeSignalVelocityMinus - densityPlus * relativeSignalVelocityPlus +
         HLLCSOLVER_HYDRO_VEPS);

    if (intermediateWaveSpeed >= 0.0) {
        return computeNumericalFlux(physValsRight, reconstructionRight, characteristicVelocityMinus,
                                    relativeSignalVelocityMinus, intermediateWaveSpeed,
                                    velocityParallel, velocityPerpendicular1,
                                    velocityPerpendicular2);
    } else {
        return computeNumericalFlux(physValsLeft, reconstructionLeft, characteristicVelocityPlus,
                                    relativeSignalVelocityPlus, intermediateWaveSpeed,
                                    velocityParallel, velocityPerpendicular1,
                                    velocityPerpendicular2);
    }
}

FieldStruct numericalFlux(const std::pair<double, double>& characteristicVelocities,
                          const PhysValues& physValsLeft, const PhysValues& physValsRight,
                          const FieldStruct& reconstructionLeft,
                          const FieldStruct& reconstructionRight, const unsigned dir) {
    // Reference: Toro E.F. - Riemann Solvers and Numerical Methods for Fluid Dynamics, Chapter 10.6
    if (characteristicVelocities.second <= 0) {
        return numericalFluxOutsideRiemannFan(physValsRight);
    } else if (characteristicVelocities.first <= 0) {
        return numericalFluxOutsideRiemannFan(physValsLeft);
    } else /* TODO: 4th case for carbuncle flag (not used in tests) */ {
        return numericalFluxInsideRiemannFan(characteristicVelocities, physValsLeft, physValsRight,
                                             reconstructionLeft, reconstructionRight, dir);
    }
}
}; // namespace RiemannSolver
