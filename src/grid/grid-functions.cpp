#include <cmath>

#include "../misc/faces.h"
#include "grid-functions.h"

namespace GridFunctions {
bool checkNaN(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) {
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
}; // namespace GridFunctions
