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

PerFaceValues reconstToConservatives(const PerFaceValues& reconstructions, const Problem& problem) {
    PerFaceValues conservatives = {};
    for (unsigned face = 0; face < Faces::FaceMax; face++) {
        conservatives[face][FieldNames::DENSITY] = reconstructions[face][FieldNames::DENSITY];
        conservatives[face][FieldNames::VELOCITY_X] =
            reconstructions[face][FieldNames::VELOCITY_X] *
            reconstructions[face][FieldNames::DENSITY];
        conservatives[face][FieldNames::VELOCITY_Y] =
            reconstructions[face][FieldNames::VELOCITY_Y] *
            reconstructions[face][FieldNames::DENSITY];
        conservatives[face][FieldNames::VELOCITY_Z] =
            reconstructions[face][FieldNames::VELOCITY_Z] *
            reconstructions[face][FieldNames::DENSITY];
        if (problem.thermal) {
            conservatives[face][FieldNames::THERMAL_ENERGY] =
                thermalEnergyToEnergy(reconstructions[face]);
        } else {
            conservatives[face][FieldNames::THERMAL_ENERGY] =
                temperatureToEnergy(reconstructions[face], problem.gamma);
        }
    }
    return conservatives;
}

double thermalEnergyToThermalPressure(const FieldStruct& fields, const double gamma) {
    return (gamma - 1.0) * fields[FieldNames::THERMAL_ENERGY];
}

double temperatureToThermalPressure(const FieldStruct& fields) {
    return fields[FieldNames::DENSITY] * fields[FieldNames::THERMAL_ENERGY];
}

PerFaceSingleValue computeThermalPressure(const PerFaceValues& fields, const Problem& problem) {
    PerFaceSingleValue thermalPressures = {};
    for (unsigned face = 0; face < Faces::FaceMax; face++) {
        if (problem.thermal) {
            thermalPressures[face] = thermalEnergyToThermalPressure(fields[face], problem.gamma);
        } else {
            thermalPressures[face] = temperatureToThermalPressure(fields[face]);
        }
    }
    return thermalPressures;
}
}; // namespace Transformation
