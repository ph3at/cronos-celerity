#include "transformations.h"

namespace Transformation {
double velocitySquared(const FieldStruct& fields) {
    return fields[FieldNames::VELOCITY_X] * fields[FieldNames::VELOCITY_X] +
           fields[FieldNames::VELOCITY_Y] * fields[FieldNames::VELOCITY_Y] +
           fields[FieldNames::VELOCITY_Z] * fields[FieldNames::VELOCITY_Z];
}

double thermalEnergyToEnergy(const FieldStruct& fields) {
    return fields[FieldNames::THERMAL_ENERGY] + 0.5 * velocitySquared(fields) / fields[FieldNames::DENSITY];
}

double energyToThermalEnergy(const FieldStruct& fields) {
    return fields[FieldNames::THERMAL_ENERGY] - 0.5 * velocitySquared(fields) * fields[FieldNames::DENSITY];
}

double temperatureToEnergy(const FieldStruct& fields, const double gamma) {
    double thermalEnergy = fields[FieldNames::THERMAL_ENERGY] * fields[FieldNames::DENSITY] / (gamma - 1.0);
    return thermalEnergy + 0.5 * velocitySquared(fields) / fields[FieldNames::DENSITY];
}

double energyToTemperature(const FieldStruct& fields, const double gamma) {
    double thermalEnergy =
        fields[FieldNames::THERMAL_ENERGY] - 0.5 * velocitySquared(fields) * fields[FieldNames::DENSITY];
    return thermalEnergy * (gamma - 1.0) / fields[FieldNames::DENSITY];
}

void reconstToConservatives(PhysValues& output, const FieldStruct& reconstructions, const bool isThermal,
                            const double gamma) {

    output.conservatives[FieldNames::DENSITY] = reconstructions[FieldNames::DENSITY];
    output.conservatives[FieldNames::VELOCITY_X] =
        reconstructions[FieldNames::VELOCITY_X] * reconstructions[FieldNames::DENSITY];
    output.conservatives[FieldNames::VELOCITY_Y] =
        reconstructions[FieldNames::VELOCITY_Y] * reconstructions[FieldNames::DENSITY];
    output.conservatives[FieldNames::VELOCITY_Z] =
        reconstructions[FieldNames::VELOCITY_Z] * reconstructions[FieldNames::DENSITY];
    if (isThermal) {
        output.conservatives[FieldNames::THERMAL_ENERGY] = reconstructions[FieldNames::THERMAL_ENERGY];
        output.conservatives[FieldNames::THERMAL_ENERGY] = thermalEnergyToEnergy(output.conservatives);
    } else {
        output.conservatives[FieldNames::THERMAL_ENERGY] = reconstructions[FieldNames::THERMAL_ENERGY];
        output.conservatives[FieldNames::THERMAL_ENERGY] = temperatureToEnergy(output.conservatives, gamma);
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

void primitiveToConservative(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const bool isThermal, const double gamma) {
    for (unsigned x = 0; x < grid.xDim(); x++) {
        for (unsigned y = 0; y < grid.yDim(); y++) {
            for (unsigned z = 0; z < grid.zDim(); z++) {
                velocityToMomentum(grid(x, y, z));
                if (isThermal) { // Compiler should pull this out of the loop
                    grid(x, y, z)[FieldNames::THERMAL_ENERGY] = thermalEnergyToEnergy(grid(x, y, z));
                } else {
                    grid(x, y, z)[FieldNames::THERMAL_ENERGY] = temperatureToEnergy(grid(x, y, z), gamma);
                }
            }
        }
    }
}

void conservativeToPrimitive(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const bool isThermal, const double gamma) {
    for (unsigned x = 0; x < grid.xDim(); x++) {
        for (unsigned y = 0; y < grid.yDim(); y++) {
            for (unsigned z = 0; z < grid.zDim(); z++) {
                momentumToVelocity(grid(x, y, z));
                if (isThermal) { // Compiler should pull this out of the loop
                    grid(x, y, z)[FieldNames::THERMAL_ENERGY] = energyToThermalEnergy(grid(x, y, z));
                } else {
                    grid(x, y, z)[FieldNames::THERMAL_ENERGY] = energyToTemperature(grid(x, y, z), gamma);
                }
            }
        }
    }
}
}; // namespace Transformation

namespace TransformationSycl {

void primitiveToConservative(cl::sycl::queue& queue, cl::sycl::buffer<FieldStruct, 3>& grid, const bool isThermal,
                             const double gamma) {
    queue.submit([&](cl::sycl::handler& cgh) {
        auto gridAccessor = grid.template get_access<cl::sycl::access::mode::read_write>(cgh);

        cgh.parallel_for(grid.get_range(), [=](const cl::sycl::id<3> id) {
            Transformation::velocityToMomentum(gridAccessor[id]);
            if (isThermal) {
                gridAccessor[id][FieldNames::THERMAL_ENERGY] = Transformation::thermalEnergyToEnergy(gridAccessor[id]);
            } else {
                gridAccessor[id][FieldNames::THERMAL_ENERGY] =
                    Transformation::temperatureToEnergy(gridAccessor[id], gamma);
            }
        });
    });
}

void conservativeToPrimitive(cl::sycl::queue& queue, cl::sycl::buffer<FieldStruct, 3>& grid, const bool isThermal,
                             const double gamma) {
    queue.submit([&](cl::sycl::handler& cgh) {
        auto gridAccessor = grid.template get_access<cl::sycl::access::mode::read_write>(cgh);

        cgh.parallel_for(grid.get_range(), [=](const cl::sycl::id<3> id) {
            Transformation::momentumToVelocity(gridAccessor[id]);
            if (isThermal) {
                gridAccessor[id][FieldNames::THERMAL_ENERGY] = Transformation::energyToThermalEnergy(gridAccessor[id]);
            } else {
                gridAccessor[id][FieldNames::THERMAL_ENERGY] =
                    Transformation::energyToTemperature(gridAccessor[id], gamma);
            }
        });
    });
}

} // namespace TransformationSycl
