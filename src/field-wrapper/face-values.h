#pragma once

#include <array>

#include "../misc/faces.h"
#include "phys-fields.h"

typedef std::array<FieldStruct, Faces::FaceMax> PerFaceValues;
typedef std::array<double, Faces::FaceMax> PerFaceSingleValue;
