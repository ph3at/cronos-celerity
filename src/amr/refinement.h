#pragma once

#include <limits>

#include "amr-node.h"

/* Each vector represents one of the dimensions. The unsigned is the offset to the first index in
 * the component with a flagged cell. The innermost pairs tell from where to where cells are
 * flagged, i.e. (4, 9) means cells with index 4,5,...,9 are flagged. */
typedef std::pair<
    unsigned,
    std::vector<std::pair<unsigned, std::vector<std::vector<std::pair<unsigned, unsigned>>>>>>
    Flags;

template <class SolverType, class ProblemType, class Fields, unsigned padding> class Refinery {
  public:
    Refinery(const AMRParameters& configuration);

    void refine();
    AMRNode<SolverType, ProblemType, Fields, padding>&
    initialRefine(PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem);

  private:
    const AMRParameters& configuration;

    std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>> amrNodes;

    std::vector<Flags> flagCells();
    void addFlags(Flags& levelFlags, const Flags& nodeFlags);
    Flags translateFlags(const Flags& subLevelFlags);
    void addBuffer(Flags& flags);
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
Refinery<SolverType, ProblemType, Fields, padding>::Refinery(const AMRParameters& configuration)
    : configuration(configuration) {}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::refine() {
    std::vector<Flags> flags = this->flagCells();
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<Flags> Refinery<SolverType, ProblemType, Fields, padding>::flagCells() {
    const int maxLevel = static_cast<int>(this->amrNodes.size());
    std::vector<Flags> levelFlags;
    levelFlags.reserve(maxLevel);
    for (int level = maxLevel - 1; level >= 0; level--) {
        for (AMRNode<SolverType, ProblemType, Fields, padding> node : this->amrNodes[level]) {
            Flags nodeFlags = {}; // truncation error stuff
            this->addFlags(levelFlags[level], nodeFlags);
        }
        if (level < maxLevel - 1) {
            Flags higherLevelFlags = this->translateFlags(levelFlags[level + 1]);
            this->addFlags(levelFlags[level], higherLevelFlags);
        }
        this->addBuffer(levelFlags[level]);
    }
    return levelFlags;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::addFlags(Flags& levelFlags,
                                                                  const Flags& nodeFlags) {
    const int levelXStart = levelFlags.first;
    const int nodeXStart = nodeFlags.first;
    const int levelXSize = levelFlags.second.size();
    const int nodeXSize = nodeFlags.second.size();
    const int levelXOffset = std::max(0, nodeXStart - levelXStart);
    const int nodeXOffset = std::max(0, levelXStart - nodeXStart);
    const int xUntil = std::min(levelXSize - levelXOffset, nodeXSize - nodeXOffset);
    for (int x = 0; x < xUntil; x++) {
        const int levelYStart = levelFlags.second[x + levelXOffset].first;
        const int nodeYStart = nodeFlags.second[x + nodeXOffset].first;
        const int levelYSize = levelFlags.second[x + levelXOffset].second.size();
        const int nodeYSize = nodeFlags.second[x + nodeXOffset].second.size();
        const int levelYOffset = std::max(0, nodeYStart - levelYStart);
        const int nodeYOffset = std::max(0, levelYStart - nodeYStart);
        const int yUntil = std::min(levelYSize - levelYOffset, nodeYSize - nodeYOffset);

        for (int y = 0; y < yUntil; y++) {
            unsigned levelZ = 0;
            unsigned nodeZ = 0;
            const std::vector<std::pair<unsigned, unsigned>>& levelPairs =
                levelFlags.second[x + levelXOffset].second[y + levelYOffset];
            const std::vector<std::pair<unsigned, unsigned>>& nodePairs =
                nodeFlags.second[x + nodeXOffset].second[y + nodeYOffset];
            const unsigned levelZEnd = levelPairs.size();
            const unsigned nodeZEnd = nodePairs.size();
            std::vector<std::pair<unsigned, unsigned>> pairs;
            if (levelZ < levelZEnd && nodeZ < nodeZEnd) {
                unsigned left = std::min(levelPairs[0].first, nodePairs[0].first);
                unsigned right = left;
                while (true) {
                    std::cout << left << " " << right << " | " << levelZ << " " << nodeZ
                              << std::endl;
                    if (levelZ < levelZEnd && right >= levelPairs[levelZ].first - 1) {
                        right = std::max(right, levelPairs[levelZ].second);
                        levelZ++;
                    } else if (nodeZ < nodeZEnd && right >= nodePairs[nodeZ].first - 1) {
                        right = std::max(right, nodePairs[nodeZ].second);
                        nodeZ++;
                    } else {
                        pairs.push_back(std::make_pair(left, right));
                        if (levelZ >= levelZEnd) {
                            pairs.insert(pairs.end(), nodePairs.begin() + nodeZ, nodePairs.end());
                            break;
                        } else if (nodeZ >= nodeZEnd) {
                            pairs.insert(pairs.end(), levelPairs.begin() + levelZ,
                                         levelPairs.end());
                            break;
                        } else {
                            left = std::min(levelPairs[levelZ].first, nodePairs[nodeZ].first);
                            right = left;
                        }
                    }
                }
            }
            levelFlags.second[x + levelXOffset].second[y + levelYOffset] = pairs;
        }

        if (levelYStart > nodeYStart) {
            levelFlags.second[x + levelXOffset].second.insert(
                levelFlags.second[x + levelXOffset].second.begin(),
                nodeFlags.second[x + nodeXOffset].second.begin(),
                nodeFlags.second[x + nodeXOffset].second.begin() + levelYStart - nodeYStart);
            levelFlags.second[x + levelXOffset].first = nodeFlags.second[x + nodeXOffset].first;
        }
        const int yExtra = nodeYStart + nodeYSize - (levelYStart + levelYSize);
        if (yExtra > 0) {
            levelFlags.second[x + levelXOffset].second.insert(
                levelFlags.second[x + levelXOffset].second.end(),
                nodeFlags.second[x + nodeXOffset].second.end() - yExtra,
                nodeFlags.second[x + nodeXOffset].second.end());
        }
    }
    if (levelXStart > nodeXStart) {
        levelFlags.second.insert(
            levelFlags.second.begin(), nodeFlags.second.begin(),
            nodeFlags.second.begin() +
                std::min(levelXStart - nodeXStart, static_cast<int>(nodeFlags.second.size())));
        levelFlags.first = nodeFlags.first;
    }
    const int xExtra = nodeXStart + nodeXSize - (levelXStart + levelXSize);
    if (xExtra > 0) {
        levelFlags.second.insert(levelFlags.second.end(), nodeFlags.second.end() - xExtra,
                                 nodeFlags.second.end());
    }
}

inline unsigned minLeft(const std::vector<unsigned>& indices,
                        const std::vector<std::vector<std::pair<unsigned, unsigned>>>& candidates) {
    unsigned minValue = std::numeric_limits<unsigned>::max();
    for (unsigned candidate = 0; candidate < indices.size(); candidate++) {
        if (indices[candidate] < candidates[candidate].size()) {
            minValue = std::min(minValue, candidates[candidate][indices[candidate]].first);
        }
    }
    return minValue;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
Flags Refinery<SolverType, ProblemType, Fields, padding>::translateFlags(const Flags& flags) {
    const unsigned xStart = flags.first / this->configuration.refinementFactor;
    const unsigned xSize = flags.second.size() / this->configuration.refinementFactor - xStart + 1;
    std::vector<std::pair<unsigned, std::vector<std::vector<std::pair<unsigned, unsigned>>>>>
        flagsX;
    flagsX.reserve(xSize);
    for (unsigned x = 0; x < xSize; x++) {
        std::vector<std::vector<std::pair<unsigned, unsigned>>> flagsY;

        unsigned yStart = std::numeric_limits<unsigned>::max();
        unsigned yEnd = 0;
        for (unsigned xOffset = 0; xOffset < this->configuration.refinementFactor; xOffset++) {
            if (x * this->configuration.refinementFactor + xOffset >= flags.first) {
                const std::pair<
                    unsigned, std::vector<std::vector<std::pair<unsigned, unsigned>>>>& flagsHere =
                    flags.second[x * this->configuration.refinementFactor + xOffset - flags.first];
                yStart = std::min(yStart, flagsHere.first);
                yEnd = std::max(yEnd,
                                static_cast<unsigned>(flagsHere.second.size()) + flagsHere.first);
            }
        }
        yStart /= this->configuration.refinementFactor;
        yEnd = yEnd / this->configuration.refinementFactor - yStart + 1;
        for (unsigned y = 0; y < yEnd; y++) {
            std::vector<std::vector<std::pair<unsigned, unsigned>>> flagsZCandidates;
            for (unsigned xOffset = 0; xOffset < this->configuration.refinementFactor; xOffset++) {
                if (x * this->configuration.refinementFactor + xOffset >= flags.first) {
                    const std::pair<unsigned,
                                    std::vector<std::vector<std::pair<unsigned, unsigned>>>>&
                        flagsHere = flags.second[x * this->configuration.refinementFactor +
                                                 xOffset - flags.first];
                    for (unsigned yOffset = 0; yOffset < this->configuration.refinementFactor;
                         yOffset++) {
                        if (y * this->configuration.refinementFactor + yOffset >= flags.first) {
                            flagsZCandidates.push_back(
                                flagsHere.second[y * this->configuration.refinementFactor +
                                                 yOffset - flagsHere.first]);
                        }
                    }
                }
            }
            std::vector<std::pair<unsigned, unsigned>> flagsZ;
            std::vector<unsigned> indices(flagsZCandidates.size(), 0);
            unsigned counter = static_cast<unsigned>(flagsZCandidates.size());
            unsigned left = minLeft(indices, flagsZCandidates);
            unsigned right = left;
            while (counter > 0) {
                for (unsigned candidate = 0; candidate <= indices.size(); candidate++) {
                    if (candidate == indices.size()) {
                        flagsZ.push_back(std::make_pair(left, right));
                        left = minLeft(indices, flagsZCandidates);
                        right = left;
                    } else if (indices[candidate] == flagsZCandidates[candidate].size()) {
                        counter--;
                        indices[candidate]++;
                    } else if (indices[candidate] < flagsZCandidates[candidate].size()) {
                        std::pair<unsigned, unsigned>& pair =
                            flagsZCandidates[candidate][indices[candidate]];
                        if (right >= pair.first / this->configuration.refinementFactor) {
                            right =
                                std::max(right, pair.second / this->configuration.refinementFactor);
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

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::addBuffer(Flags& flags) {}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRNode<SolverType, ProblemType, Fields, padding>&
Refinery<SolverType, ProblemType, Fields, padding>::initialRefine(
    PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem) {
    AMRNode<SolverType, ProblemType, Fields, padding> root(grid, problem, this->configuration,
                                                           problem.timeDelta, 0.0);
    this->amrNodes.push_back(
        std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>(1, root));
    return this->amrNodes[0][0];
}
