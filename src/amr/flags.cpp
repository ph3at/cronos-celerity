#include <limits>

#include "flags.h"

namespace CellFlags {
void addFlags(Flags& knownFlags, const Flags& newFlags) {
    const int knownXStart = knownFlags.first;
    const int newXStart = newFlags.first;
    const int knownXSize = knownFlags.second.size();
    const int newXSize = newFlags.second.size();
    const int knownXOffset = std::max(0, newXStart - knownXStart);
    const int newXOffset = std::max(0, knownXStart - newXStart);
    const int xUntil = std::min(knownXSize - knownXOffset, newXSize - newXOffset);
    for (int x = 0; x < xUntil; x++) {
        std::pair<unsigned, FlagsY>& knownFlagsX = knownFlags.second[x + knownXOffset];
        const std::pair<unsigned, FlagsY>& newFlagsX = newFlags.second[x + newXOffset];
        const int knownYStart = knownFlagsX.first;
        const int newYStart = newFlagsX.first;
        const int knownYSize = knownFlagsX.second.size();
        const int newYSize = newFlagsX.second.size();
        const int knownYOffset = std::max(0, newYStart - knownYStart);
        const int newYOffset = std::max(0, knownYStart - newYStart);
        const int yUntil = std::min(knownYSize - knownYOffset, newYSize - newYOffset);
        for (int y = 0; y < yUntil; y++) {
            unsigned knownZ = 0;
            unsigned newZ = 0;
            const FlagsZ& knownPairs = knownFlagsX.second[y + knownYOffset];
            const FlagsZ& newPairs = newFlagsX.second[y + newYOffset];
            const unsigned knownZEnd = knownPairs.size();
            const unsigned newZEnd = newPairs.size();
            FlagsZ pairs;
            if (knownZ < knownZEnd && newZ < newZEnd) {
                unsigned left = std::min(knownPairs[0].first, newPairs[0].first);
                unsigned right = left;
                while (true) {
                    if (knownZ < knownZEnd && right >= knownPairs[knownZ].first - 1) {
                        right = std::max(right, knownPairs[knownZ].second);
                        knownZ++;
                    } else if (newZ < newZEnd && right >= newPairs[newZ].first - 1) {
                        right = std::max(right, newPairs[newZ].second);
                        newZ++;
                    } else {
                        pairs.push_back(std::make_pair(left, right));
                        if (knownZ >= knownZEnd) {
                            pairs.insert(pairs.end(), newPairs.begin() + newZ, newPairs.end());
                            break;
                        } else if (newZ >= newZEnd) {
                            pairs.insert(pairs.end(), knownPairs.begin() + knownZ,
                                         knownPairs.end());
                            break;
                        } else {
                            left = std::min(knownPairs[knownZ].first, newPairs[newZ].first);
                            right = left;
                        }
                    }
                }
            }
            knownFlagsX.second[y + knownYOffset] = pairs;
        }

        if (knownYStart > newYStart) {
            knownFlagsX.second.insert(knownFlagsX.second.begin(), newFlagsX.second.begin(),
                                      newFlagsX.second.begin() + knownYStart - newYStart);
            knownFlagsX.first = newFlagsX.first;
        }
        const int yExtra = newYStart + newYSize - (knownYStart + knownYSize);
        if (yExtra > 0) {
            knownFlagsX.second.insert(knownFlagsX.second.end(), newFlagsX.second.end() - yExtra,
                                      newFlagsX.second.end());
        }
    }
    if (knownXStart > newXStart) {
        knownFlags.second.insert(
            knownFlags.second.begin(), newFlags.second.begin(),
            newFlags.second.begin() +
                std::min(knownXStart - newXStart, static_cast<int>(newFlags.second.size())));
        knownFlags.first = newFlags.first;
    }
    const int xExtra = newXStart + newXSize - (knownXStart + knownXSize);
    if (xExtra > 0) {
        knownFlags.second.insert(knownFlags.second.end(), newFlags.second.end() - xExtra,
                                 newFlags.second.end());
    }
}

unsigned minLeft(const std::vector<unsigned>& indices, const FlagsY& candidates) {
    unsigned minValue = std::numeric_limits<unsigned>::max();
    for (unsigned candidate = 0; candidate < indices.size(); candidate++) {
        if (indices[candidate] < candidates[candidate].size()) {
            minValue = std::min(minValue, candidates[candidate][indices[candidate]].first);
        }
    }
    return minValue;
}

Flags translateFlags(const Flags& flags, const unsigned resolutionFactor) {
    const unsigned xFirst = flags.first;
    const unsigned xLast = xFirst + flags.second.size();
    const unsigned xStart = flags.first / resolutionFactor;
    const unsigned xEnd = (flags.second.size() + flags.first - 1) / resolutionFactor;
    FlagsX flagsX;
    flagsX.reserve(xEnd - xStart + 1);
    for (unsigned x = xStart; x <= xEnd; x++) {
        FlagsY flagsY;
        unsigned yFirst = std::numeric_limits<unsigned>::max();
        unsigned yLast = 0;
        for (unsigned xOffset = 0; xOffset < resolutionFactor; xOffset++) {
            unsigned xx = x * resolutionFactor + xOffset;
            if (xx >= xFirst && xx < xLast) {
                const std::pair<unsigned, FlagsY>& flagsHere = flags.second[xx - flags.first];
                yFirst = std::min(yFirst, flagsHere.first);
                yLast = std::max(yLast,
                                 static_cast<unsigned>(flagsHere.second.size()) + flagsHere.first);
            }
        }
        unsigned yStart = yFirst / resolutionFactor;
        unsigned yEnd = (yLast - 1) / resolutionFactor;
        for (unsigned y = yStart; y <= yEnd; y++) {
            FlagsY flagsZCandidates;
            for (unsigned xOffset = 0; xOffset < resolutionFactor; xOffset++) {
                const unsigned xx = x * resolutionFactor + xOffset;
                if (xx >= xFirst && xx < xLast) {
                    const std::pair<unsigned, FlagsY>& flagsHere = flags.second[xx - flags.first];
                    for (unsigned yOffset = 0; yOffset < resolutionFactor; yOffset++) {
                        unsigned yy = y * resolutionFactor + yOffset;
                        if (yy >= flagsHere.first &&
                            yy < flagsHere.first + flagsHere.second.size()) {
                            flagsZCandidates.push_back(flagsHere.second[yy - flagsHere.first]);
                        }
                    }
                }
            }
            FlagsZ flagsZ;
            std::vector<unsigned> indices(flagsZCandidates.size(), 0);
            unsigned counter = static_cast<unsigned>(flagsZCandidates.size());
            unsigned left = minLeft(indices, flagsZCandidates) / resolutionFactor;
            unsigned right = left;
            while (counter > 0) {
                for (unsigned candidate = 0; candidate <= indices.size(); candidate++) {
                    if (candidate == indices.size()) {
                        flagsZ.push_back(std::make_pair(left, right));
                        left = minLeft(indices, flagsZCandidates) / resolutionFactor;
                        right = left;
                    } else if (indices[candidate] == flagsZCandidates[candidate].size()) {
                        counter--;
                        indices[candidate]++;
                    } else if (indices[candidate] < flagsZCandidates[candidate].size()) {
                        std::pair<unsigned, unsigned>& pair =
                            flagsZCandidates[candidate][indices[candidate]];
                        if (right + 1 >= pair.first / resolutionFactor) {
                            right = std::max(right, pair.second / resolutionFactor);
                            indices[candidate]++;
                            break;
                        }
                    }
                }
            }
            flagsY.push_back(flagsZ);
        }
        flagsX.push_back(std::make_pair(yStart, flagsY));
    }
    return std::make_pair(xStart, flagsX);
}

void addBuffer(Flags& flags, const unsigned bufferSize) {}
}; // namespace CellFlags
