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

bool checkNaN(cl::sycl::queue& queue, cl::sycl::buffer<FieldStruct, 3>& grid) {
    const auto range = grid.get_range();
    auto isNanBuffer = cl::sycl::buffer<bool, 1>(range.size());
    queue.submit([&](cl::sycl::handler& cgh) {
        auto gridAccessor = grid.template get_access<cl::sycl::access::mode::read>(cgh);
        auto isNanBufferAccessor = isNanBuffer.template get_access<cl::sycl::access::mode::discard_write>(cgh);
        cgh.parallel_for(range, [=](const cl::sycl::item<3> item) {
            isNanBufferAccessor[item.get_linear_id()] = false;
            for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                isNanBufferAccessor[item.get_linear_id()] |= std::isnan(gridAccessor[item][field]);
            }
        });
    });
    queue.wait_and_throw();
    return utils::sycl::reduce(queue, isNanBuffer, std::bit_or<bool>{});
}

} // namespace GridFunctionsSycl
