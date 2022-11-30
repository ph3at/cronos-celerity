#pragma once

#include <vector>

#include <CL/sycl.hpp>

#include "../configuration/constants.h"
#include "../data-types/faces.h"
#include "../data-types/fields.h"
#include "../grid/padded-grid.h"
#include "../grid/utils.h"

namespace Reconstruction {

PerFaceValues reconstruct(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned x, const unsigned y,
                          const unsigned z);

} // namespace Reconstruction

namespace ReconstructionSycl {

static inline double limitMinmod(const double deltaLeft, const double deltaRight) {
    if ((deltaLeft < 0 && deltaRight < 0) || (deltaLeft > 0 && deltaRight > 0)) {
        if (deltaRight > 0) {
            return std::min(deltaRight, deltaLeft);
        } else {
            return std::max(deltaRight, deltaLeft);
        }
    } else {
        return 0.0;
    }
}

static inline double derivX(const cl::sycl::accessor<FieldStruct, 3, cl::sycl::access::mode::read,
                                                     cl::sycl::access::target::global_buffer>& gridAccessor,
                            const cl::sycl::id<3>& idx, const unsigned field) {
    const auto prevX = cl::sycl::id<3>(idx[0] - 1, idx[1], idx[2]);
    const auto nextX = cl::sycl::id<3>(idx[0] + 1, idx[1], idx[2]);
    double deltaLeft = gridAccessor[idx][field] - gridAccessor[prevX][field];
    double deltaRight = gridAccessor[nextX][field] - gridAccessor[idx][field];

    return limitMinmod(deltaLeft, deltaRight);
}

static inline double derivY(const cl::sycl::accessor<FieldStruct, 3, cl::sycl::access::mode::read,
                                                     cl::sycl::access::target::global_buffer>& gridAccessor,
                            const cl::sycl::id<3>& idx, const unsigned field) {
    const auto prevY = cl::sycl::id<3>(idx[0], idx[1] - 1, idx[2]);
    const auto nextY = cl::sycl::id<3>(idx[0], idx[1] + 1, idx[2]);
    double deltaLeft = gridAccessor[idx][field] - gridAccessor[prevY][field];
    double deltaRight = gridAccessor[nextY][field] - gridAccessor[idx][field];

    return limitMinmod(deltaLeft, deltaRight);
}

static inline double derivZ(const cl::sycl::accessor<FieldStruct, 3, cl::sycl::access::mode::read,
                                                     cl::sycl::access::target::global_buffer>& gridAccessor,
                            const cl::sycl::id<3>& idx, const unsigned field) {
    const auto prevZ = cl::sycl::id<3>(idx[0], idx[1], idx[2] - 1);
    const auto nextZ = cl::sycl::id<3>(idx[0], idx[1], idx[2] + 1);
    double deltaLeft = gridAccessor[idx][field] - gridAccessor[prevZ][field];
    double deltaRight = gridAccessor[nextZ][field] - gridAccessor[idx][field];

    return limitMinmod(deltaLeft, deltaRight);
}

static inline PerFaceValues reconstruct(const cl::sycl::accessor<FieldStruct, 3, cl::sycl::access::mode::read,
                                                                 cl::sycl::access::target::global_buffer>& gridAccessor,
                                        const cl::sycl::id<3>& idx) {
    const auto prevX = cl::sycl::id<3>(idx[0] - 1, idx[1], idx[2]);
    const auto prevY = cl::sycl::id<3>(idx[0], idx[1] - 1, idx[2]);
    const auto prevZ = cl::sycl::id<3>(idx[0], idx[1], idx[2] - 1);

    PerFaceValues reconst = {};
    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
        reconst[Faces::FaceWest][field] = gridAccessor[idx][field] - 0.5 * derivX(gridAccessor, idx, field);
        reconst[Faces::FaceEast][field] = gridAccessor[prevX][field] + 0.5 * derivX(gridAccessor, prevX, field);

        reconst[Faces::FaceSouth][field] = gridAccessor[idx][field] - 0.5 * derivY(gridAccessor, idx, field);
        reconst[Faces::FaceNorth][field] = gridAccessor[prevY][field] + 0.5 * derivY(gridAccessor, prevY, field);

        reconst[Faces::FaceBottom][field] = gridAccessor[idx][field] - 0.5 * derivZ(gridAccessor, idx, field);
        reconst[Faces::FaceTop][field] = gridAccessor[prevZ][field] + 0.5 * derivZ(gridAccessor, prevZ, field);
    }

    return reconst;
}

} // namespace ReconstructionSycl
