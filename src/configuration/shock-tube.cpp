#include "shock-tube.h"

ShockTube::ShockTube(const double cflThreshold, const bool thermal, const double timeDelta,
                     const double timeStart, const double timeEnd, const bool preciseEnd,
                     const double gamma, const std::array<double, Direction::DirMax> posLeft,
                     const std::array<double, Direction::DirMax> posRight,
                     const std::array<std::size_t, Direction::DirMax> numberCells,
                     const std::array<BoundaryType, Faces::FaceMax> boundaryTypes,
                     const Direction shockDir, const double shockPos, const double densityLeftInit,
                     const double densityRightInit, const double velocityLeftInit,
                     const double velocityRightInit, const double pressureLeftInit,
                     const double pressureRightInit)
    : Problem<ShockTube, FieldStruct, GHOST_CELLS>(cflThreshold, thermal, timeDelta, timeStart,
                                                   timeEnd, preciseEnd, gamma, posLeft, posRight,
                                                   numberCells, boundaryTypes),
      shockDir(shockDir), shockPos(shockPos), densityLeftInit(densityLeftInit),
      densityRightInit(densityRightInit),
      velocityXLeftInit(shockDir == Direction::DirX ? velocityLeftInit : 0.0),
      velocityXRightInit(shockDir == Direction::DirX ? velocityRightInit : 0.0),
      velocityYLeftInit(shockDir == Direction::DirY ? velocityLeftInit : 0.0),
      velocityYRightInit(shockDir == Direction::DirY ? velocityRightInit : 0.0),
      velocityZLeftInit(shockDir == Direction::DirZ ? velocityLeftInit : 0.0),
      velocityZRightInit(shockDir == Direction::DirZ ? velocityRightInit : 0.0),
      pressureLeftInit(pressureLeftInit), pressureRightInit(pressureRightInit) {}

ShockTube::ShockTube(const ShockTube& blueprint)
    : Problem<ShockTube, FieldStruct, GHOST_CELLS>(blueprint), shockDir(blueprint.shockDir),
      shockPos(blueprint.shockPos), densityLeftInit(blueprint.densityLeftInit),
      densityRightInit(blueprint.densityRightInit), velocityXLeftInit(blueprint.velocityXLeftInit),
      velocityXRightInit(blueprint.velocityXRightInit),
      velocityYLeftInit(blueprint.velocityYLeftInit),
      velocityYRightInit(blueprint.velocityYRightInit),
      velocityZLeftInit(blueprint.velocityZLeftInit),
      velocityZRightInit(blueprint.velocityZRightInit),
      pressureLeftInit(blueprint.pressureLeftInit), pressureRightInit(blueprint.pressureRightInit) {
}

ShockTube::ShockTube(const toml::table& config)
    : Problem<ShockTube, FieldStruct, GHOST_CELLS>(config) {
    std::cout << std::endl << "Shock tube specific parameters:" << std::endl;
    const toml::node_view<const toml::node>& specific = config["specific"];
    this->shockDir = parseValue<Direction>(specific, "shock_dir", Direction::DirX);
    this->shockPos = parseValue<double>(specific, "shock_pos", 0.5);
    this->densityLeftInit = parseValue<double>(specific, "density_left_init", 1.0);
    this->densityRightInit = parseValue<double>(specific, "density_right_init", 1.0);
    this->velocityXLeftInit = this->shockDir == Direction::DirX
                                  ? parseValue<double>(specific, "velocity_left_init", -20.0)
                                  : 0.0;
    this->velocityXRightInit = this->shockDir == Direction::DirX
                                   ? parseValue<double>(specific, "velocity_right_init", -20.0)
                                   : 0.0;
    this->velocityYLeftInit = this->shockDir == Direction::DirY
                                  ? parseValue<double>(specific, "velocity_left_init", -20.0)
                                  : 0.0;
    this->velocityYRightInit = this->shockDir == Direction::DirY
                                   ? parseValue<double>(specific, "velocity_right_init", -20.0)
                                   : 0.0;
    this->velocityZLeftInit = this->shockDir == Direction::DirZ
                                  ? parseValue<double>(specific, "velocity_left_init", -20.0)
                                  : 0.0;
    this->velocityZRightInit = this->shockDir == Direction::DirZ
                                   ? parseValue<double>(specific, "velocity_right_init", -20.0)
                                   : 0.0;
    this->pressureLeftInit = parseValue<double>(specific, "pressure_left_init", 1000.0);
    this->pressureRightInit = parseValue<double>(specific, "pressure_right_init", 0.01);
}

void ShockTube::initialiseGrid(PaddedGrid<FieldStruct, GHOST_CELLS>& grid) const {
    double posParallel;
    for (unsigned x = 0; x < grid.xDim(); x++) {
        for (unsigned y = 0; y < grid.yDim(); y++) {
            for (unsigned z = 0; z < grid.zDim(); z++) {
                if (this->shockDir == Direction::DirX) {
                    posParallel = this->posLeft[Direction::DirX] +
                                  (static_cast<double>(x) - GHOST_CELLS + 0.5) *
                                      this->cellSize[Direction::DirX];
                } else if (this->shockDir == Direction::DirY) {
                    posParallel = this->posLeft[Direction::DirY] +
                                  (static_cast<double>(y) - GHOST_CELLS + 0.5) *
                                      this->cellSize[Direction::DirY];
                } else {
                    posParallel = this->posLeft[Direction::DirZ] +
                                  (static_cast<double>(z) - GHOST_CELLS + 0.5) *
                                      this->cellSize[Direction::DirZ];
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
