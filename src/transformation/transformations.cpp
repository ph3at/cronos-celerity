#include "transformations.h"

namespace Transformation {

double energyToThermalEnergy(const FieldStruct& fields) {
    return fields[FieldNames::THERMAL_ENERGY] - 0.5 * velocitySquared(fields) * fields[FieldNames::DENSITY];
}

double energyToTemperature(const FieldStruct& fields, const double gamma) {
    double thermalEnergy =
        fields[FieldNames::THERMAL_ENERGY] - 0.5 * velocitySquared(fields) * fields[FieldNames::DENSITY];
    return thermalEnergy * (gamma - 1.0) / fields[FieldNames::DENSITY];
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

void primitiveToConservative(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid, const bool isThermal,
                             const double gamma) {
    queue.submit([&](sycl::handler& cgh) {
        auto gridAccessor = grid.template get_access<sycl::access::mode::read_write>(cgh);

        cgh.parallel_for(grid.get_range(), [=](const sycl::id<3> id) {
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

void conservativeToPrimitive(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid, const bool isThermal,
                             const double gamma) {
    queue.submit([&](sycl::handler& cgh) {
        auto gridAccessor = grid.template get_access<sycl::access::mode::read_write>(cgh);

        cgh.parallel_for(grid.get_range(), [=](const sycl::id<3> id) {
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

namespace TransformationCelerity {

void primitiveToConservative(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid, const bool isThermal,
                             const double gamma) {
    queue.submit([=](celerity::handler& cgh) {
        auto gridAccessor = celerity::accessor{ grid, cgh, celerity::access::one_to_one{}, celerity::read_write };

        cgh.parallel_for(grid.get_range(), [=](const celerity::id<3> id) {
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

void conservativeToPrimitive(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid, const bool isThermal,
                             const double gamma) {
    queue.submit([=](celerity::handler& cgh) {
        auto gridAccessor = celerity::accessor{ grid, cgh, celerity::access::one_to_one{}, celerity::read_write };

        cgh.parallel_for(grid.get_range(), [=](const celerity::id<3> id) {
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

} // namespace TransformationCelerity