#pragma once

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
    void refine();
    AMRNode<SolverType, ProblemType, Fields, padding>&
    initialRefine(PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem,
                  const AMRParameters& configuration);

  private:
    std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>> amrNodes;

    std::vector<Flags> flagCells();
    void addFlags(Flags& levelFlags, const Flags& nodeFlags);
    void translateAndAddFlags(Flags& levelFlags, const Flags& subLevelFlags);
    void addBuffer(Flags& flags);
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::refine() {
    std::vector<Flags> flags = this->flagCells();
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::vector<Flags> Refinery<SolverType, ProblemType, Fields, padding>::flagCells() {
    const unsigned maxLevel = this->amrNodes.size();
    std::vector<Flags> levelFlags;
    levelFlags.reserve(maxLevel);
    for (unsigned level = maxLevel - 1; level >= 0; level--) {
        for (AMRNode<SolverType, ProblemType, Fields, padding> node : this->amrNodes[level]) {
            Flags nodeFlags = {}; // truncation error stuff
            this->addFlags(levelFlags[level], nodeFlags);
        }
        if (level < maxLevel) {
            this->translateAndAddFlags(levelFlags[level], levelFlags[level + 1]);
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
                while (levelZ < levelZEnd && nodeZ < nodeZEnd) {
                    if (right >= levelPairs[levelZ].first - 1) {
                        right = std::max(right, levelPairs[levelZ].second);
                        levelZ++;
                    } else if (right >= nodePairs[nodeZ].first - 1) {
                        right = std::max(right, nodePairs[nodeZ].second);
                        nodeZ++;
                    } else {
                        pairs.push_back(std::make_pair(left, right));
                        left = std::min(levelPairs[levelZ].first, nodePairs[nodeZ].first);
                    }
                }
            }
            if (nodeZ < nodeZEnd) {
                pairs.insert(pairs.end(), nodePairs.begin() + nodeZ, nodePairs.end());
            }
            levelFlags.second[x + levelXOffset].second[y + levelYOffset] = pairs;
        }
        if (levelYStart > nodeYStart) {
            levelFlags.second[x + levelXOffset].second.insert(
                levelFlags.second[x + levelXOffset].second.begin(),
                nodeFlags.second[x + nodeXOffset].second.begin(),
                nodeFlags.second[x + nodeXOffset].second.begin() + nodeYStart - levelYStart);
            levelFlags.second[x + levelXOffset].first = nodeFlags.second[x + nodeXOffset].first;
        }
        const unsigned yExtra = nodeYStart + nodeYSize - (levelYStart + levelYSize);
        if (yExtra > 0) {
            levelFlags.second[x + levelXOffset].second.insert(
                levelFlags.second[x + levelXOffset].second.end(),
                nodeFlags.second[x + nodeXOffset].second.end() - yExtra,
                nodeFlags.second[x + nodeXOffset].second.end());
        }
    }
    if (levelXStart > nodeXStart) {
        levelFlags.second.insert(levelFlags.second.begin(), nodeFlags.second.begin(),
                                 nodeFlags.second.begin() + nodeXStart - levelXStart);
        levelFlags.first = nodeFlags.first;
    }
    const unsigned xExtra = nodeXStart + nodeXSize - (levelXStart + levelXSize);
    if (xExtra > 0) {
        levelFlags.second.insert(levelFlags.second.end(), nodeFlags.second.end() - xExtra,
                                 nodeFlags.second.end());
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::translateAndAddFlags(
    Flags& levelFlags, const Flags& subLevelFlags) {}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void Refinery<SolverType, ProblemType, Fields, padding>::addBuffer(Flags& flags) {}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRNode<SolverType, ProblemType, Fields, padding>&
Refinery<SolverType, ProblemType, Fields, padding>::initialRefine(
    PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem,
    const AMRParameters& configuration) {}
