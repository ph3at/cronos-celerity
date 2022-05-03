#include "problem.h"

constexpr bool THERMAL = true;
constexpr double TIME_DELTA = 0.001;
constexpr double TIME_START = 0.0;
constexpr double TIME_END = 1.0;
constexpr double GAMMA = 1.4;
constexpr std::size_t NUMBER_CELLS_X = 10;
constexpr double X_START = 0.0;
constexpr double X_END = 1.0;
constexpr std::size_t NUMBER_CELLS_Y = 10;
constexpr double Y_START = 0.0;
constexpr double Y_END = 1.0;
constexpr std::size_t NUMBER_CELLS_Z = 10;
constexpr double Z_START = 0.0;
constexpr double Z_END = 1.0;

Problem Problem::initialiseTestProblem() {
    const std::array<std::size_t, Direction::DirMax> numberCells = { NUMBER_CELLS_X, NUMBER_CELLS_Y,
                                                                     NUMBER_CELLS_Z };
    const std::array<double, Direction::DirMax> cellSize = { (X_END - X_START) / NUMBER_CELLS_X,
                                                             (Y_END - Y_START) / NUMBER_CELLS_Y,
                                                             (Z_END - Z_START) / NUMBER_CELLS_Z };
    const std::array<double, Direction::DirMax> inverseCellSize = { 1.0 / cellSize[0],
                                                                    1.0 / cellSize[1],
                                                                    1.0 / cellSize[2] };

    const std::array<BoundaryType, Faces::FaceMax> boundaryTypes = {
        BoundaryType::OUTFLOW,     BoundaryType::CONSTANT,    BoundaryType::EXTRAPOLATE,
        BoundaryType::EXTRAPOLATE, BoundaryType::EXTRAPOLATE, BoundaryType::EXTRAPOLATE
    };

    return { .thermal = THERMAL,
             .timeDelta = TIME_DELTA,
             .timeStart = TIME_START,
             .timeEnd = TIME_END,
             .gamma = GAMMA,
             .numberCells = numberCells,
             .cellSize = cellSize,
             .inverseCellSize = inverseCellSize,
             .boundaryTypes = boundaryTypes };
}
