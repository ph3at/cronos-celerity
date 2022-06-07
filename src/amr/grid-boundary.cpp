#include "grid-boundary.h"

namespace GridBoundary {
bool isOverlappingOrAdjaecent(const Boundary& grid1, const Boundary& grid2) {
    return grid1[0].second + 1 >= grid2[0].first && grid1[0].first <= grid2[0].second + 1 &&
           grid1[1].second + 1 >= grid2[1].first && grid1[1].first <= grid2[1].second + 1 &&
           grid1[2].second + 1 >= grid2[2].first && grid1[2].first <= grid2[2].second + 1;
}

std::optional<Boundary> overlappingArea(const Boundary& grid1, const Boundary& grid2) {
    if (grid1[0].second >= grid2[0].first && grid1[0].first <= grid2[0].second &&
        grid1[1].second >= grid2[1].first && grid1[1].first <= grid2[1].second &&
        grid1[2].second >= grid2[2].first && grid1[2].first <= grid2[2].second) {
        unsigned xStart = std::max(grid1[0].first, grid2[0].first);
        unsigned xEnd = std::min(grid1[0].second, grid2[0].second);
        unsigned yStart = std::max(grid1[1].first, grid2[1].first);
        unsigned yEnd = std::min(grid1[1].second, grid2[1].second);
        unsigned zStart = std::max(grid1[2].first, grid2[2].first);
        unsigned zEnd = std::min(grid1[2].second, grid2[2].second);
        const GridBoundary::Boundary overlapping = { std::make_pair(xStart, xEnd),
                                                     std::make_pair(yStart, yEnd),
                                                     std::make_pair(zStart, zEnd) };
        return std::make_optional(overlapping);
    } else {
        return std::nullopt;
    }
}

Boundary translateUp(const Boundary& grid, const unsigned factor) {
    return { std::make_pair(grid[0].first * factor, grid[0].second * factor + factor - 1),
             std::make_pair(grid[1].first * factor, grid[1].second * factor + factor - 1),
             std::make_pair(grid[2].first * factor, grid[2].second * factor + factor - 1) };
}

Boundary translateDown(const Boundary& grid, const unsigned factor) {
    return { std::make_pair(grid[0].first / factor, grid[0].second / factor),
             std::make_pair(grid[1].first / factor, grid[1].second / factor),
             std::make_pair(grid[2].first / factor, grid[2].second / factor) };
}

std::vector<Boundary> subtract(const Boundary& grid, const Boundary& subtractor) {
    if (grid[0].first > subtractor[0].second || grid[0].second < subtractor[0].first ||
        grid[1].first > subtractor[1].second || grid[1].second < subtractor[1].first ||
        grid[2].first > subtractor[2].second || grid[2].second < subtractor[2].first) {
        return std::vector<Boundary>(1, grid);
    }
    Boundary remaining = grid;
    std::vector<Boundary> result;
    if (remaining[0].first < subtractor[0].first) {
        result.push_back({ std::make_pair(remaining[0].first, subtractor[0].first - 1),
                           remaining[1], remaining[2] });
        remaining[0].first = subtractor[0].first;
    }
    if (remaining[0].second > subtractor[0].second) {
        result.push_back({ std::make_pair(subtractor[0].second + 1, remaining[0].second),
                           remaining[1], remaining[2] });
        remaining[0].second = subtractor[0].second;
    }
    if (remaining[1].first < subtractor[1].first) {
        result.push_back({ remaining[0],
                           std::make_pair(remaining[1].first, subtractor[1].first - 1),
                           remaining[2] });
        remaining[1].first = subtractor[1].first;
    }
    if (remaining[1].second > subtractor[1].second) {
        result.push_back({ remaining[0],
                           std::make_pair(subtractor[1].second + 1, remaining[1].second),
                           remaining[2] });
        remaining[1].second = subtractor[1].second;
    }
    if (remaining[2].first < subtractor[2].first) {
        result.push_back({ remaining[0], remaining[1],
                           std::make_pair(remaining[2].first, subtractor[2].first - 1) });
        remaining[2].first = subtractor[2].first;
    }
    if (remaining[2].second > subtractor[2].second) {
        result.push_back({ remaining[0], remaining[1],
                           std::make_pair(subtractor[2].second + 1, remaining[2].second) });
        remaining[2].second = subtractor[1].second;
    }
    return result;
}

void strictMerge(std::vector<Boundary>& grids) {
    std::vector<bool> keepGrid(grids.size(), true);
    std::vector<Boundary> newGrids;
    for (unsigned grid1 = 0; grid1 + 1 < grids.size(); grid1++) {
        if (keepGrid[grid1]) {
            const Boundary& g1 = grids[grid1];
            for (unsigned grid2 = grid1 + 1; grid2 < grids.size(); grid2++) {
                if (keepGrid[grid2]) {
                    const Boundary& g2 = grids[grid2];
                    if (std::max(g1[0].first, g2[0].first) ==
                            std::min(g1[0].second, g2[0].second) + 1 &&
                        g1[1].first == g2[1].first && g1[1].second == g2[1].second &&
                        g1[2].first == g2[2].first && g1[2].second == g2[2].second) {
                        keepGrid[grid1] = false;
                        keepGrid[grid2] = false;
                        newGrids.push_back({ std::make_pair(std::min(g1[0].first, g2[0].first),
                                                            std::max(g2[0].second, g2[0].second)),
                                             g1[1], g1[2] });
                    } else if (std::max(g1[1].first, g2[1].first) ==
                                   std::min(g1[1].second, g2[1].second) + 1 &&
                               g1[0].first == g2[0].first && g1[0].second == g2[0].second &&
                               g1[2].first == g2[2].first && g1[2].second == g2[2].second) {
                        keepGrid[grid1] = false;
                        keepGrid[grid2] = false;
                        newGrids.push_back({ g1[0],
                                             std::make_pair(std::min(g1[1].first, g2[1].first),
                                                            std::max(g2[1].second, g2[1].second)),
                                             g1[2] });
                    } else if (std::max(g1[2].first, g2[2].first) ==
                                   std::min(g1[2].second, g2[2].second) + 1 &&
                               g1[0].first == g2[0].first && g1[0].second == g2[0].second &&
                               g1[1].first == g2[1].first && g1[1].second == g2[1].second) {
                        keepGrid[grid1] = false;
                        keepGrid[grid2] = false;
                        newGrids.push_back({
                            g1[0],
                            g1[1],
                            std::make_pair(std::min(g1[2].first, g2[2].first),
                                           std::max(g2[2].second, g2[2].second)),
                        });
                    }
                }
            }
        }
    }
    if (newGrids.size() > 0) {
        for (unsigned grid = 0; grid < grids.size(); grid++) {
            if (keepGrid[grid]) {
                newGrids.push_back(grids[grid]);
            }
        }
        strictMerge(newGrids);
        grids.swap(newGrids);
    }
}

std::vector<Boundary> remainingGrids(const Boundary& grid, const std::vector<Boundary>& removed) {
    std::vector<Boundary> remaining;
    remaining.push_back(grid);
    for (const Boundary& remove : removed) {
        std::vector<Boundary> stillRemaining;
        for (const Boundary& left : remaining) {
            std::vector<Boundary> keep = subtract(left, remove);
            stillRemaining.insert(stillRemaining.end(), keep.begin(), keep.end());
        }
        remaining.swap(stillRemaining);
    }
    strictMerge(remaining);
    return remaining;
}
} // namespace GridBoundary
