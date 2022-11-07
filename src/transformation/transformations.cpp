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

double energyToThermalEnergy(const FieldStruct& fields) {
    return fields[FieldNames::THERMAL_ENERGY] -
           0.5 * velocitySquared(fields) * fields[FieldNames::DENSITY];
}

double temperatureToEnergy(const FieldStruct& fields, const double gamma) {
    double thermalEnergy =
        fields[FieldNames::THERMAL_ENERGY] * fields[FieldNames::DENSITY] / (gamma - 1.0);
    return thermalEnergy + 0.5 * velocitySquared(fields) / fields[FieldNames::DENSITY];
}

double energyToTemperature(const FieldStruct& fields, const double gamma) {
    double thermalEnergy = fields[FieldNames::THERMAL_ENERGY] -
                           0.5 * velocitySquared(fields) * fields[FieldNames::DENSITY];
    return thermalEnergy * (gamma - 1.0) / fields[FieldNames::DENSITY];
}

void reconstToConservatives(PhysValues& output, const FieldStruct& reconstructions,
                            const bool isThermal, const double gamma) {

    output.conservatives[FieldNames::DENSITY] = reconstructions[FieldNames::DENSITY];
    output.conservatives[FieldNames::VELOCITY_X] =
        reconstructions[FieldNames::VELOCITY_X] * reconstructions[FieldNames::DENSITY];
    output.conservatives[FieldNames::VELOCITY_Y] =
        reconstructions[FieldNames::VELOCITY_Y] * reconstructions[FieldNames::DENSITY];
    output.conservatives[FieldNames::VELOCITY_Z] =
        reconstructions[FieldNames::VELOCITY_Z] * reconstructions[FieldNames::DENSITY];
    if (isThermal) {
        output.conservatives[FieldNames::THERMAL_ENERGY] =
            reconstructions[FieldNames::THERMAL_ENERGY];
        output.conservatives[FieldNames::THERMAL_ENERGY] =
            thermalEnergyToEnergy(output.conservatives);
    } else {
        output.conservatives[FieldNames::THERMAL_ENERGY] =
            reconstructions[FieldNames::THERMAL_ENERGY];
        output.conservatives[FieldNames::THERMAL_ENERGY] =
            temperatureToEnergy(output.conservatives, gamma);
    }
}

double thermalEnergyToThermalPressure(const FieldStruct& fields, const double gamma) {
    return (gamma - 1.0) * fields[FieldNames::THERMAL_ENERGY];
}

double temperatureToThermalPressure(const FieldStruct& fields) {
    return fields[FieldNames::DENSITY] * fields[FieldNames::THERMAL_ENERGY];
}

double computeThermalPressure(const FieldStruct& fields, const bool isThermal, const double gamma) {
    if (isThermal) {
        return thermalEnergyToThermalPressure(fields, gamma);
    } else {
        return temperatureToThermalPressure(fields);
    }
}

void velocityToMomentum(FieldStruct& fields) {
    fields[FieldNames::VELOCITY_X] *= fields[FieldNames::DENSITY];
    fields[FieldNames::VELOCITY_Y] *= fields[FieldNames::DENSITY];
    fields[FieldNames::VELOCITY_Z] *= fields[FieldNames::DENSITY];
}

void momentumToVelocity(FieldStruct& fields) {
    fields[FieldNames::VELOCITY_X] /= fields[FieldNames::DENSITY];
    fields[FieldNames::VELOCITY_Y] /= fields[FieldNames::DENSITY];
    fields[FieldNames::VELOCITY_Z] /= fields[FieldNames::DENSITY];
}

void primitiveToConservative(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const bool isThermal,
                             const double gamma) {
    for (unsigned x = 0; x < grid.xDim(); x++) {
        for (unsigned y = 0; y < grid.yDim(); y++) {
            for (unsigned z = 0; z < grid.zDim(); z++) {
                velocityToMomentum(grid(x, y, z));
                if (isThermal) { // Compiler should pull this out of the loop
                    grid(x, y, z)[FieldNames::THERMAL_ENERGY] =
                        thermalEnergyToEnergy(grid(x, y, z));
                } else {
                    grid(x, y, z)[FieldNames::THERMAL_ENERGY] =
                        temperatureToEnergy(grid(x, y, z), gamma);
                }
            }
        }
    }
}

void conservativeToPrimitive(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const bool isThermal,
                             const double gamma) {
    for (unsigned x = 0; x < grid.xDim(); x++) {
        for (unsigned y = 0; y < grid.yDim(); y++) {
            for (unsigned z = 0; z < grid.zDim(); z++) {
                momentumToVelocity(grid(x, y, z));
                if (isThermal) { // Compiler should pull this out of the loop
                    grid(x, y, z)[FieldNames::THERMAL_ENERGY] =
                        energyToThermalEnergy(grid(x, y, z));
                } else {
                    grid(x, y, z)[FieldNames::THERMAL_ENERGY] =
                        energyToTemperature(grid(x, y, z), gamma);
                }
            }
        }
    }
}
}; // namespace Transformation

namespace TransformationSycl {

using grid::utils::idx3d;

void primitiveToConservative(std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims,
                             const bool isThermal, const double gamma) {
    for (std::size_t x = 0; x < dims[0]; x++) {
        for (std::size_t y = 0; y < dims[1]; y++) {
            for (std::size_t z = 0; z < dims[2]; z++) {
                Transformation::velocityToMomentum(grid[idx3d(x, y, z, dims)]);
                if (isThermal) { // Compiler should pull this out of the loop
                    grid[idx3d(x, y, z, dims)][FieldNames::THERMAL_ENERGY] =
                        Transformation::thermalEnergyToEnergy(grid[idx3d(x, y, z, dims)]);
                } else {
                    grid[idx3d(x, y, z, dims)][FieldNames::THERMAL_ENERGY] =
                        Transformation::temperatureToEnergy(grid[idx3d(x, y, z, dims)], gamma);
                }
            }
        }
    }
}

void conservativeToPrimitive(std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims,
                             const bool isThermal, const double gamma) {
    for (std::size_t x = 0; x < dims[0]; ++x) {
        for (std::size_t y = 0; y < dims[1]; ++y) {
            for (std::size_t z = 0; z < dims[2]; ++z) {
                Transformation::momentumToVelocity(grid[idx3d(x, y, z, dims)]);
                if (isThermal) { // Compiler should pull this out of the loop
                    grid[idx3d(x, y, z, dims)][FieldNames::THERMAL_ENERGY] =
                        Transformation::energyToThermalEnergy(grid[idx3d(x, y, z, dims)]);
                } else {
                    grid[idx3d(x, y, z, dims)][FieldNames::THERMAL_ENERGY] =
                        Transformation::energyToTemperature(grid[idx3d(x, y, z, dims)], gamma);
                }
            }
        }
    }
}

} // namespace TransformationSycl
