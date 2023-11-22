#include "grid-functions.h"

#include <cmath>
#include <fstream>

#include "../data-types/faces.h"
#include "../sycl/reduction.h"

namespace GridFunctions {
bool checkNaN(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid) {
    for (unsigned x = 0; x < grid.xDim(); x++) {
        for (unsigned y = 0; y < grid.yDim(); y++) {
            for (unsigned z = 0; z < grid.zDim(); z++) {
                for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                    if (std::isnan(grid(x, y, z)[field])) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

SimpleGrid<FieldStruct> readFromFile(const std::string filename) {
    std::ifstream inputFile(filename, std::ios::in);
    if (!inputFile.is_open()) {
        std::cerr << "File \"" << filename << "\" could not be opened." << std::endl;
        return SimpleGrid<FieldStruct>({}, 0, 0, 0);
    } else {
        unsigned xDim, yDim, zDim;
        inputFile >> xDim >> yDim >> zDim;
        SimpleGrid<FieldStruct> result = SimpleGrid<FieldStruct>({}, xDim, yDim, zDim);
        double value;
        for (unsigned x = 0; x < xDim; x++) {
            for (unsigned y = 0; y < yDim; y++) {
                for (unsigned z = 0; z < zDim; z++) {
                    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                        inputFile >> value;
                        result(x, y, z)[field] = value;
                    }
                }
            }
        }
        return result;
    }
}
}; // namespace GridFunctions

namespace GridFunctionsSycl {

bool checkNaN(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid) {
    const auto range = grid.get_range();
    auto isNanBuffer = sycl::buffer<bool, 1>(range.size());
    queue.submit([&](sycl::handler& cgh) {
        auto gridAccessor = grid.template get_access<sycl::access::mode::read>(cgh);
        auto isNanBufferAccessor = isNanBuffer.template get_access<sycl::access::mode::discard_write>(cgh);
        cgh.parallel_for(range, [=](const sycl::item<3> item) {
            isNanBufferAccessor[item.get_linear_id()] = false;
            for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                isNanBufferAccessor[item.get_linear_id()] |= std::isnan(gridAccessor[item][field]);
            }
        });
    });
    return sycl_utils::reduce(queue, isNanBuffer, std::bit_or<bool>{});
}

} // namespace GridFunctionsSycl

namespace GridFunctionsCelerity {

bool checkNaN(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid) {
    auto resultBuffer = celerity::buffer<bool, 1>{ celerity::range<1>{ 1 } };

    queue.submit([&resultBuffer, &grid](celerity::handler& cgh) {
        auto gridAccessor = celerity::accessor{ grid, cgh, celerity::access::one_to_one{}, celerity::read_only };
        auto reduction = celerity::reduction(resultBuffer, cgh, sycl::bit_or<bool>{},
                                             celerity::property::reduction::initialize_to_identity{});

        const auto range = grid.get_range();
        cgh.parallel_for(range, reduction, [=](const celerity::item<3> item, auto& val) {
            auto hasNan = false;
            for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; ++field) {
                hasNan |= std::isnan(gridAccessor[item][field]);
            }
            val.combine(hasNan);
        });
    });

    auto nanOccurred = false;

    queue.submit([&resultBuffer, &nanOccurred](celerity::handler& cgh) {
        celerity::accessor bufferAccessor{ resultBuffer, cgh, celerity::access::all{}, celerity::read_only_host_task };
        cgh.host_task(
            celerity::experimental::collective,
            [=, &nanOccurred](celerity::experimental::collective_partition) { nanOccurred |= bufferAccessor[0]; });
    });

    queue.slow_full_sync();

    return nanOccurred;
}

} // namespace GridFunctionsCelerity
