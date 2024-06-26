#include <limits>

#include "flags.h"

namespace CellFlags {
void addFlags(Flags& knownFlags, const Flags& newFlags) {
    const int knownXStart = static_cast<int>(knownFlags.first);
    const int newXStart = static_cast<int>(newFlags.first);
    const int knownXSize = static_cast<int>(knownFlags.second.size());
    const int newXSize = static_cast<int>(newFlags.second.size());
    const int knownXOffset = std::max(0, newXStart - knownXStart);
    const int newXOffset = std::max(0, knownXStart - newXStart);
    const int xUntil = std::min(knownXSize - knownXOffset, newXSize - newXOffset);
    for (int x = 0; x < xUntil; x++) {
        std::pair<unsigned, FlagsY>& knownFlagsX = knownFlags.second[x + knownXOffset];
        const std::pair<unsigned, FlagsY>& newFlagsX = newFlags.second[x + newXOffset];
        const int knownYStart = static_cast<int>(knownFlagsX.first);
        const int newYStart = static_cast<int>(newFlagsX.first);
        const int knownYSize = static_cast<int>(knownFlagsX.second.size());
        const int newYSize = static_cast<int>(newFlagsX.second.size());
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
                    if (knownZ < knownZEnd && right + 1 >= knownPairs[knownZ].first) {
                        right = std::max(right, knownPairs[knownZ].second);
                        knownZ++;
                    } else if (newZ < newZEnd && right + 1 >= newPairs[newZ].first) {
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

FlagsZ mergeFlagsZ(const FlagsY& flags) {
    bool allEmpty = true;
    for (FlagsZ candidate : flags) {
        if (candidate.size() > 0) {
            allEmpty = false;
            break;
        }
    }
    if (allEmpty) {
        return FlagsZ({});
    } else {
        FlagsZ newFlags;
        std::vector<unsigned> indices(flags.size(), 0);
        unsigned counter = static_cast<unsigned>(flags.size());
        unsigned left = minLeft(indices, flags);
        unsigned right = left;
        while (counter > 0) {
            for (unsigned candidate = 0; candidate <= indices.size(); candidate++) {
                if (candidate == indices.size()) {
                    newFlags.push_back(std::make_pair(left, right));
                    left = minLeft(indices, flags);
                    right = left;
                } else if (indices[candidate] == flags[candidate].size()) {
                    counter--;
                    indices[candidate]++;
                } else if (indices[candidate] < flags[candidate].size()) {
                    const std::pair<unsigned, unsigned>& pair =
                        flags[candidate][indices[candidate]];
                    if (right + 1 >= pair.first) {
                        right = std::max(right, pair.second);
                        indices[candidate]++;
                        break;
                    }
                }
            }
        }
        return newFlags;
    }
}

std::pair<unsigned, FlagsY> mergeFlagsY(const FlagsX& flags) {
    FlagsY newFlags;
    unsigned yMin = std::numeric_limits<unsigned>::max();
    unsigned yEnd = 0;
    for (std::pair<unsigned, FlagsY> flagsY : flags) {
        yMin = std::min(yMin, flagsY.first);
        yEnd = std::max(yEnd, flagsY.first + static_cast<unsigned>(flagsY.second.size()));
    }
    newFlags.reserve(yEnd - yMin);
    for (unsigned y = yMin; y < yEnd; y++) {
        FlagsY flagCandidates;
        for (std::pair<unsigned, FlagsY> flagsY : flags) {
            if (y >= flagsY.first &&
                y < flagsY.first + static_cast<unsigned>(flagsY.second.size())) {
                flagCandidates.push_back(flagsY.second[y - flagsY.first]);
            }
        }
        newFlags.push_back(mergeFlagsZ(flagCandidates));
    }
    return std::make_pair(yMin, newFlags);
}

void addBuffer(Flags& flags, const unsigned bufferSize) {
    for (std::pair<unsigned, FlagsY>& flagsY : flags.second) {
        FlagsY newFlagsY;
        for (FlagsZ& flagsZ : flagsY.second) {
            FlagsZ newPairs;
            for (unsigned pair = 0; pair < flagsZ.size(); pair++) {
                unsigned left = std::max(flagsZ[pair].first, bufferSize) - bufferSize;
                unsigned right = flagsZ[pair].second + bufferSize;
                while (pair + 1 < flagsZ.size() &&
                       right + bufferSize + 1 >= flagsZ[pair + 1].first) {
                    pair++;
                    right = flagsZ[pair].second + bufferSize;
                }
                newPairs.push_back(std::make_pair(left, right));
            }
            newFlagsY.push_back(newPairs);
        }
        FlagsY newNewFlagsY;
        newNewFlagsY.reserve(newFlagsY.size() + 4);
        for (unsigned y = std::max(flagsY.first, bufferSize) - flagsY.first; y < bufferSize; y++) {
            const unsigned offset = std::min(y + 1, static_cast<unsigned>(flagsY.second.size()));
            newNewFlagsY.push_back(
                mergeFlagsZ(FlagsY(newFlagsY.begin(), newFlagsY.begin() + offset)));
        }
        for (unsigned y = 0; y < flagsY.second.size(); y++) {
            const unsigned beginOffset = std::max(bufferSize, y) - bufferSize;
            const unsigned endOffset =
                std::min(y + bufferSize + 1, static_cast<unsigned>(flagsY.second.size()));
            newNewFlagsY.push_back(mergeFlagsZ(
                FlagsY(newFlagsY.begin() + beginOffset, newFlagsY.begin() + endOffset)));
        }
        for (unsigned y = bufferSize; y > 0; y--) {
            const unsigned offset = std::min(static_cast<unsigned>(flagsY.second.size()), y);
            FlagsZ newFlagsZ = mergeFlagsZ(FlagsY(newFlagsY.end() - offset, newFlagsY.end()));
            if (newFlagsZ.size() == 0) {
                break;
            } else {
                newNewFlagsY.push_back(newFlagsZ);
            }
        }
        flagsY.second = newNewFlagsY;
        flagsY.first = std::max(flagsY.first, bufferSize) - bufferSize;
    }
    FlagsX newFlagsX;
    newFlagsX.reserve(flags.second.size() + 4);
    for (unsigned x = std::max(flags.first, bufferSize) - flags.first; x < bufferSize; x++) {
        const unsigned offset = std::min(x + 1, static_cast<unsigned>(flags.second.size()));
        newFlagsX.push_back(
            mergeFlagsY(FlagsX(flags.second.begin(), flags.second.begin() + offset)));
    }
    for (unsigned x = 0; x < flags.second.size(); x++) {
        const unsigned beginOffset = std::max(bufferSize, x) - bufferSize;
        const unsigned endOffset =
            std::min(x + bufferSize + 1, static_cast<unsigned>(flags.second.size()));
        newFlagsX.push_back(mergeFlagsY(
            FlagsX(flags.second.begin() + beginOffset, flags.second.begin() + endOffset)));
    }
    for (unsigned x = bufferSize; x > 0; x--) {
        const unsigned offset = std::min(static_cast<unsigned>(flags.second.size()), x);
        newFlagsX.push_back(mergeFlagsY(FlagsX(flags.second.end() - offset, flags.second.end())));
    }
    flags.second = newFlagsX;
    flags.first = std::max(flags.first, bufferSize) - bufferSize;
}
}; // namespace CellFlags
