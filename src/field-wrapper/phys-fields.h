#pragma once

#include <array>
#include <cstddef>

constexpr unsigned PHYSICAL_FIELDS = 5;
typedef std::array<double, PHYSICAL_FIELDS> FieldStruct;

// TODO: Find out meaning of fields
enum FieldNames { FIELD1 = 0, FIELD2 = 1, FIELD3 = 2, FIELD4 = 3, FIELD5 = 4 };
