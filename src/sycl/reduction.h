#pragma once

#include <algorithm>
#include <bit>
#include <functional>

#include <CL/sycl.hpp>

namespace utils::sycl {

template <typename T, typename BinOp = std::plus<T>>
T reduce(cl::sycl::queue& queue, cl::sycl::buffer<T, 1>& data, BinOp op = std::plus<T>{},
         const T neutralElement = T()) {
    auto device = queue.get_device();

    const auto actualSize = data.get_count();
    const auto ceiledSize = std::bit_ceil(actualSize);
    const size_t workgroupSize = std::min(ceiledSize, device.get_info<cl::sycl::info::device::max_work_group_size>());
    size_t remainingSize = ceiledSize;

    data.set_final_data(nullptr);

    {
        do {
            const auto reductionKernel = [actualSize, remainingSize, workgroupSize, &data, op,
                                          neutralElement](cl::sycl::handler& cgh) mutable {
                const auto range =
                    cl::sycl::nd_range<1>{ cl::sycl::range<1>{ std::max(remainingSize, workgroupSize) },
                                           cl::sycl::range<1>{ std::min(remainingSize, workgroupSize) } };

                auto dataAccessor = data.template get_access<cl::sycl::access::mode::read_write>(cgh);
                auto scratchBuffer =
                    cl::sycl::accessor<T, 1, cl::sycl::access::mode::read_write, cl::sycl::access::target::local>(
                        cl::sycl::range<1>(workgroupSize), cgh);

                cgh.parallel_for(range, [dataAccessor, scratchBuffer, actualSize, workgroupSize, remainingSize, op,
                                         neutralElement](cl::sycl::nd_item<1> id) {
                    const auto gid = id.get_global_id(0);
                    const auto lid = id.get_local_id(0);

                    if (gid < actualSize) {
                        scratchBuffer[lid] = dataAccessor[gid];
                    } else {
                        scratchBuffer[lid] = neutralElement;
                    }
                    id.barrier(cl::sycl::access::fence_space::local_space);

                    if (gid < remainingSize) {
                        const auto minSize = std::min(remainingSize, workgroupSize);
                        for (size_t offset = minSize / 2; offset > 0; offset /= 2) {
                            if (lid < offset) {
                                scratchBuffer[lid] = op(scratchBuffer[lid], scratchBuffer[lid + offset]);
                            }
                            id.barrier(cl::sycl::access::fence_space::local_space);
                        }

                        if (lid == 0) {
                            dataAccessor[id.get_group(0)] = scratchBuffer[lid];
                        }
                    }
                });
            };
            queue.submit(reductionKernel);

            remainingSize = remainingSize / workgroupSize;
        } while (remainingSize > 1);
    }

    const auto dataAccessor = data.template get_access<cl::sycl::access::mode::read>();
    return dataAccessor[0];
}

} // namespace utils::sycl
