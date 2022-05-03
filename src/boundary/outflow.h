#pragma once

#include "../field-wrapper/phys-fields.h"
#include "../grid/padded-grid.h"

namespace Outflow {
template <unsigned padding>
void apply(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const unsigned face);
}; // namespace Outflow
