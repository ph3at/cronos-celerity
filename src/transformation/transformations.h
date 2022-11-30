#pragma once

#include <array>

#include <CL/sycl.hpp>

#include "../configuration/constants.h"
#include "../data-types/face-values.h"
#include "../data-types/phys-values.h"
#include "../grid/padded-grid.h"
#include "../grid/utils.h"

namespace Transformation {

constexpr double DEFAULT_GAMMA = 1.4;

static inline double velocitySquared(const FieldStruct& fields) {
    return fields[FieldNames::VELOCITY_X] * fields[FieldNames::VELOCITY_X] +
           fields[FieldNames::VELOCITY_Y] * fields[FieldNames::VELOCITY_Y] +
           fields[FieldNames::VELOCITY_Z] * fields[FieldNames::VELOCITY_Z];
}

static inline double thermalEnergyToEnergy(const FieldStruct& fields) {
    return fields[FieldNames::THERMAL_ENERGY] + 0.5 * velocitySquared(fields) / fields[FieldNames::DENSITY];
}

static inline double temperatureToEnergy(const FieldStruct& fields, const double gamma) {
    double thermalEnergy = fields[FieldNames::THERMAL_ENERGY] * fields[FieldNames::DENSITY] / (gamma - 1.0);
    return thermalEnergy + 0.5 * velocitySquared(fields) / fields[FieldNames::DENSITY];
}

static inline void reconstToConservatives(PhysValues& output, const FieldStruct& reconstructions, const bool isThermal,
                                          const double gamma = DEFAULT_GAMMA) {

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

static inline double thermalEnergyToThermalPressure(const FieldStruct& fields, const double gamma) {
    return (gamma - 1.0) * fields[FieldNames::THERMAL_ENERGY];
}

static inline double temperatureToThermalPressure(const FieldStruct& fields) {
    return fields[FieldNames::DENSITY] * fields[FieldNames::THERMAL_ENERGY];
}

static inline double computeThermalPressure(const FieldStruct& fields, const bool isThermal,
                                            const double gamma = DEFAULT_GAMMA) {
    if (isThermal) {
        return thermalEnergyToThermalPressure(fields, gamma);
    } else {
        return temperatureToThermalPressure(fields);
    }
}

void primitiveToConservative(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const bool isThermal,
                             const double gamma = DEFAULT_GAMMA);

void conservativeToPrimitive(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const bool isThermal,
                             const double gamma = DEFAULT_GAMMA);

} // namespace Transformation

namespace TransformationSycl {

void primitiveToConservative(cl::sycl::queue& queue, cl::sycl::buffer<FieldStruct, 3>& grid, const bool isThermal,
                             const double gamma = Transformation::DEFAULT_GAMMA);

void conservativeToPrimitive(cl::sycl::queue& queue, cl::sycl::buffer<FieldStruct, 3>& grid, const bool isThermal,
                             const double gamma = Transformation::DEFAULT_GAMMA);

} // namespace TransformationSycl
