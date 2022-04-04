#include "transformations.h"

namespace Transformation {
double velocitySquared(const FieldStruct& fields) {
    return fields[FieldNames::VELOCITY_X] * fields[FieldNames::VELOCITY_X] +
           fields[FieldNames::VELOCITY_Y] * fields[FieldNames::VELOCITY_Y] +
           fields[FieldNames::VELOCITY_Z] * fields[FieldNames::VELOCITY_Z];
}

double thermalEnergyToEnergy(const FieldStruct& fields) {
    return fields[FieldNames::THERMAL_ENERGY] +
           0.5 * velocitySquared(fields) / fields[FieldNames::DENSITY];
}

double temperatureToEnergy(const FieldStruct& fields, const double gamma) {
    double thermalEnergy =
        fields[FieldNames::THERMAL_ENERGY] * fields[FieldNames::DENSITY] / (gamma - 1.0);
    return thermalEnergy + 0.5 * velocitySquared(fields) / fields[FieldNames::DENSITY];
}

void reconstToConservatives(PhysValues& output, const FieldStruct& reconstructions,
                            const Problem& problem) {

    output.conservatives[FieldNames::DENSITY] = reconstructions[FieldNames::DENSITY];
    output.conservatives[FieldNames::VELOCITY_X] =
        reconstructions[FieldNames::VELOCITY_X] * reconstructions[FieldNames::DENSITY];
    output.conservatives[FieldNames::VELOCITY_Y] =
        reconstructions[FieldNames::VELOCITY_Y] * reconstructions[FieldNames::DENSITY];
    output.conservatives[FieldNames::VELOCITY_Z] =
        reconstructions[FieldNames::VELOCITY_Z] * reconstructions[FieldNames::DENSITY];
    if (problem.thermal) {
        output.conservatives[FieldNames::THERMAL_ENERGY] = thermalEnergyToEnergy(reconstructions);
    } else {
        output.conservatives[FieldNames::THERMAL_ENERGY] =
            temperatureToEnergy(reconstructions, problem.gamma);
    }
}

double thermalEnergyToThermalPressure(const FieldStruct& fields, const double gamma) {
    return (gamma - 1.0) * fields[FieldNames::THERMAL_ENERGY];
}

double temperatureToThermalPressure(const FieldStruct& fields) {
    return fields[FieldNames::DENSITY] * fields[FieldNames::THERMAL_ENERGY];
}

double computeThermalPressure(const FieldStruct& fields, const Problem& problem) {
    if (problem.thermal) {
        return thermalEnergyToThermalPressure(fields, problem.gamma);
    } else {
        return temperatureToThermalPressure(fields);
    }
}
}; // namespace Transformation
