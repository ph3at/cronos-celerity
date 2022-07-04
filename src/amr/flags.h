#pragma once

#include <vector>

namespace CellFlags {
/* Each vector represents one of the dimensions. The unsigned is the offset to the first index in
 * the component with a flagged cell. The innermost pairs tell from where to where cells are
 * flagged, i.e. (4, 9) means cells with index 4,5,...,9 are flagged. */
typedef std::vector<std::pair<unsigned, unsigned>> FlagsZ;
typedef std::vector<FlagsZ> FlagsY;
typedef std::vector<std::pair<unsigned, FlagsY>> FlagsX;
typedef std::pair<unsigned, FlagsX> Flags;

/* Flags in new flags are added into knownFlags, inplace. */
void addFlags(Flags& knownFlags, const Flags& newFlags);

/* Flags are translated to a grid with resolution lower by resolution factor. */
Flags translateFlags(const Flags& higherFlags, const unsigned resolutionFactor);

/* Flags all cells at most bufferSize cells away from already flagged cells */
void addBuffer(Flags& flags, const unsigned bufferSize);
}; // namespace CellFlags
