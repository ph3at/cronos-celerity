#include "shock-tube.h"

ShockTube::ShockTube(const double cflThreshold, const bool thermal, const double timeDelta,
                     const double timeStart, const double timeEnd, const double gamma,
                     const std::array<std::size_t, Direction::DirMax> numberCells,
                     const std::array<double, Direction::DirMax> cellSize,
                     const std::array<BoundaryType, Faces::FaceMax> boundaryTypes,
                     const Direction shockDir, const double shockPos, const double densityLeftInit,
                     const double densityRightInit, const double velocityLeftInit,
                     const double velocityRightInit, const double pressureLeftInit,
                     const double pressureRightInit)
    : Problem<ShockTube>(cflThreshold, thermal, timeDelta, timeStart, timeEnd, gamma, numberCells,
                         cellSize, boundaryTypes),
      shockDir(shockDir), shockPos(shockPos), densityLeftInit(densityLeftInit),
      densityRightInit(densityRightInit),
      velocityXLeftInit(shockDir == Direction::DirX ? velocityLeftInit : 0.0),
      velocityXRightInit(shockDir == Direction::DirX ? velocityRightInit : 0.0),
      velocityYLeftInit(shockDir == Direction::DirY ? velocityLeftInit : 0.0),
      velocityYRightInit(shockDir == Direction::DirY ? velocityRightInit : 0.0),
      velocityZLeftInit(shockDir == Direction::DirZ ? velocityLeftInit : 0.0),
      velocityZRightInit(shockDir == Direction::DirZ ? velocityRightInit : 0.0),
      pressureLeftInit(pressureLeftInit), pressureRightInit(pressureRightInit) {}

void ShockTube::initialiseGrid(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const {
    double posParallel;
    for (unsigned x = 0; x < grid.xDim(); x++) {
        for (unsigned y = 0; y < grid.yDim(); y++) {
            for (unsigned z = 0; z < grid.zDim(); z++) {
                if (this->shockDir == Direction::DirX) {
                    posParallel = this->cellSize[Direction::DirX] *
                                  (static_cast<double>(x) - GHOST_CELLS + 0.5);
                } else if (this->shockDir == Direction::DirY) {
                    posParallel = this->cellSize[Direction::DirY] *
                                  (static_cast<double>(y) - GHOST_CELLS + 0.5);
                } else {
                    posParallel = this->cellSize[Direction::DirZ] *
                                  (static_cast<double>(z) - GHOST_CELLS + 0.5);
                }
                if (posParallel < this->shockPos) {
                    grid(x, y, z)[FieldNames::DENSITY] = this->densityLeftInit;
                    grid(x, y, z)[FieldNames::VELOCITY_X] = this->velocityXLeftInit;
                    grid(x, y, z)[FieldNames::VELOCITY_Y] = this->velocityYLeftInit;
                    grid(x, y, z)[FieldNames::VELOCITY_Z] = this->velocityZLeftInit;
                    if (this->thermal) {
                        grid(x, y, z)[FieldNames::THERMAL_ENERGY] =
                            this->pressureLeftInit / (this->gamma - 1.0);
                    } else {
                        grid(x, y, z)[FieldNames::THERMAL_ENERGY] =
                            this->pressureLeftInit / (this->densityLeftInit);
                    }
                } else {
                    grid(x, y, z)[FieldNames::DENSITY] = this->densityRightInit;
                    grid(x, y, z)[FieldNames::VELOCITY_X] = this->velocityXRightInit;
                    grid(x, y, z)[FieldNames::VELOCITY_Y] = this->velocityYRightInit;
                    grid(x, y, z)[FieldNames::VELOCITY_Z] = this->velocityZRightInit;
                    if (this->thermal) {
                        grid(x, y, z)[FieldNames::THERMAL_ENERGY] =
                            this->pressureRightInit / (this->gamma - 1.0);
                    } else {
                        grid(x, y, z)[FieldNames::THERMAL_ENERGY] =
                            this->pressureRightInit / (this->densityRightInit);
                    }
                }
            }
        }
    }
}

void ShockTube::applyBoundary(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned field,
                              const unsigned face) const {
    if (face == Faces::FaceEast) {

        double initValue;
        if (field == FieldNames::DENSITY) {
            initValue = this->densityRightInit;
        } else if (field == FieldNames::VELOCITY_X) {
            initValue = this->velocityXRightInit;
        } else if (field == FieldNames::VELOCITY_Y) {
            initValue = this->velocityYRightInit;
        } else if (field == FieldNames::VELOCITY_Z) {
            initValue = this->velocityZRightInit;
        } else if (field == FieldNames::THERMAL_ENERGY) {
            initValue = this->pressureRightInit / this->densityRightInit / (this->gamma - 1.0);
        } else {
            initValue = 0.0;
        }

        for (unsigned x = grid.xEnd(); x < grid.xDim(); x++) {
            for (unsigned y = grid.yStart(); y < grid.yEnd(); y++) {
                for (unsigned z = grid.zStart(); z < grid.zEnd(); z++) {
                    grid(x, y, z)[field] = initValue;
                }
            }
        }
    }
}

void ShockTube::applySource([[maybe_unused]] PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const {
    // No source in shock tube
}

constexpr double CFL_THRESHOLD = 0.4;
constexpr bool THERMAL = true;
constexpr double TIME_DELTA = 0.000002;
constexpr double TIME_START = 0.0;
constexpr double TIME_END = 0.0005;
constexpr double GAMMA = 1.4;
constexpr std::size_t NUMBER_CELLS_X = 200;
constexpr double X_START = 0.0;
constexpr double X_END = 1.0;
constexpr std::size_t NUMBER_CELLS_Y = 5;
constexpr double Y_START = 0.0;
constexpr double Y_END = 1.0;
constexpr std::size_t NUMBER_CELLS_Z = 5;
constexpr double Z_START = 0.0;
constexpr double Z_END = 1.0;
constexpr Direction SHOCK_DIR = Direction::DirX;
constexpr double SHOCK_POS = 0.5;
constexpr double DENSITY_LEFT_INIT = 1.0;
constexpr double DENSITY_RIGHT_INIT = 1.0;
constexpr double VELOCITY_LEFT_INIT = -19.59745;
constexpr double VELOCITY_RIGHT_INIT = -19.59745;
constexpr double PRESSURE_LEFT_INIT = 1000.0;
constexpr double PRESSURE_RIGHT_INIT = 0.01;

ShockTube ShockTube::initialiseTestProblem() {
    const std::array<std::size_t, Direction::DirMax> numberCells = { NUMBER_CELLS_X, NUMBER_CELLS_Y,
                                                                     NUMBER_CELLS_Z };

    const std::array<double, Direction::DirMax> cellSize = { (X_END - X_START) / NUMBER_CELLS_X,
                                                             (Y_END - Y_START) / NUMBER_CELLS_Y,
                                                             (Z_END - Z_START) / NUMBER_CELLS_Z };

    const std::array<BoundaryType, Faces::FaceMax> boundaryTypes = {
        BoundaryType::OUTFLOW,     BoundaryType::USER,        BoundaryType::EXTRAPOLATE,
        BoundaryType::EXTRAPOLATE, BoundaryType::EXTRAPOLATE, BoundaryType::EXTRAPOLATE
    };

    return ShockTube(CFL_THRESHOLD, THERMAL, TIME_DELTA, TIME_START, TIME_END, GAMMA, numberCells,
                     cellSize, boundaryTypes, SHOCK_DIR, SHOCK_POS, DENSITY_LEFT_INIT,
                     DENSITY_RIGHT_INIT, VELOCITY_LEFT_INIT, VELOCITY_RIGHT_INIT,
                     PRESSURE_LEFT_INIT, PRESSURE_RIGHT_INIT);
}
