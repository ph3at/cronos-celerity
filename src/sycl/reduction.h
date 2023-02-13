#pragma once

#include <algorithm>
#include <climits>
#include <cstddef>
#include <functional>

#include <sycl/sycl.hpp>

namespace sycl_utils {

template <typename T> constexpr T bit_ceil(const T num) {
    static_assert(std::is_unsigned_v<T> && std::is_integral_v<T>,
                  "bit_ceil is only defined for unsigned integral types");

    if (num == 0) {
        return 1;
    }

    auto ceiled = num - 1;
    for (std::size_t i = 1; i < sizeof(T) * CHAR_BIT; i *= 2) {
        ceiled |= ceiled >> i;
    }
    ++ceiled;
    return ceiled;
}

template <typename T, typename BinOp = std::plus<T>>
T reduce(sycl::queue& queue, sycl::buffer<T, 1>& data, BinOp op = std::plus<T>{}, const T neutralElement = T()) {
    auto device = queue.get_device();

    const auto actualSize = data.get_count();
    const auto ceiledSize = bit_ceil(actualSize);
    const size_t workgroupSize = std::min(ceiledSize, device.get_info<sycl::info::device::max_work_group_size>());
    size_t remainingSize = ceiledSize;

    data.set_final_data(nullptr);

    {
        do {
            const auto reductionKernel = [actualSize, remainingSize, workgroupSize, &data, op,
                                          neutralElement](sycl::handler& cgh) mutable {
                const auto range = sycl::nd_range<1>{ sycl::range<1>{ std::max(remainingSize, workgroupSize) },
                                                      sycl::range<1>{ std::min(remainingSize, workgroupSize) } };

                auto dataAccessor = data.template get_access<sycl::access::mode::read_write>(cgh);
                auto scratchBuffer = sycl::accessor<T, 1, sycl::access::mode::read_write, sycl::access::target::local>(
                    sycl::range<1>(workgroupSize), cgh);

                cgh.parallel_for(range, [dataAccessor, scratchBuffer, actualSize, workgroupSize, remainingSize, op,
                                         neutralElement](sycl::nd_item<1> id) {
                    const auto gid = id.get_global_id(0);
                    const auto lid = id.get_local_id(0);

                    if (gid < actualSize) {
                        scratchBuffer[lid] = dataAccessor[gid];
                    } else {
                        scratchBuffer[lid] = neutralElement;
                    }
                    id.barrier(sycl::access::fence_space::local_space);

                    if (gid < remainingSize) {
                        const auto minSize = std::min(remainingSize, workgroupSize);
                        for (size_t offset = minSize / 2; offset > 0; offset /= 2) {
                            if (lid < offset) {
                                scratchBuffer[lid] = op(scratchBuffer[lid], scratchBuffer[lid + offset]);
                            }
                            id.barrier(sycl::access::fence_space::local_space);
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

    const auto dataAccessor = data.template get_access<sycl::access::mode::read>();
    return dataAccessor[0];
}

} // namespace sycl_utils
