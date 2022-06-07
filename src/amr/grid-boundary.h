#pragma once

#include <array>
#include <optional>
#include <vector>

namespace GridBoundary {
typedef std::array<std::pair<unsigned, unsigned>, 3> Boundary;

bool isOverlappingOrAdjaecent(const Boundary& grid1, const Boundary& grid2);
std::optional<GridBoundary::Boundary> overlappingArea(const Boundary& grid1, const Boundary& grid2);
Boundary translateUp(const Boundary& grid, const unsigned factor);
Boundary translateDown(const Boundary& grid, const unsigned factor);
std::vector<Boundary> subtract(const Boundary& grid, const Boundary& subtractor);
void strictMerge(std::vector<Boundary>& grids);
std::vector<Boundary> remainingGrids(const Boundary& grid, const std::vector<Boundary>& removed);
}; // namespace GridBoundary
