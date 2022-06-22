#include "value-parser.h"

template <>
BoundaryType parseValue<BoundaryType>(const toml::node_view<const toml::node>& config,
                                      const std::string& name, const BoundaryType& defaultValue) {
    const std::optional<std::string> maybeValue = config[name].value<std::string>();
    if (maybeValue.has_value()) {
        const std::string& boundaryName = maybeValue.value();
        if (boundaryName.compare("empty") == 0) {
            std::cout << name << " = Empty" << std::endl;
            return BoundaryType::EMPTY;
        } else if (boundaryName.compare("extrapolate") == 0) {
            std::cout << name << " = Extrapolate" << std::endl;
            return BoundaryType::EXTRAPOLATE;
        } else if (boundaryName.compare("outflow") == 0) {
            std::cout << name << " = Outflow" << std::endl;
            return BoundaryType::OUTFLOW;
        } else if (boundaryName.compare("user") == 0) {
            std::cout << name << " = User" << std::endl;
            return BoundaryType::USER;
        } else {
            std::cerr << "Unknown boundary type " << boundaryName << std::endl;
            exit(5);
        }
    } else {
        std::cout << name << " = " << defaultValue << " (default)" << std::endl;
        return defaultValue;
    }
}

template <>
Direction parseValue<Direction>(const toml::node_view<const toml::node>& config,
                                const std::string& name, const Direction& defaultValue) {
    const std::optional<std::string> maybeValue = config[name].value<std::string>();
    if (maybeValue.has_value()) {
        const std::string& directionName = maybeValue.value();
        if (directionName.compare("x") == 0) {
            std::cout << name << " = x" << std::endl;
            return Direction::DirX;
        } else if (directionName.compare("y") == 0) {
            std::cout << name << " = y" << std::endl;
            return Direction::DirY;
        } else if (directionName.compare("z") == 0) {
            std::cout << name << " = z" << std::endl;
            return Direction::DirZ;
        } else {
            std::cerr << "Unknown direction name " << directionName << std::endl;
            exit(5);
        }
    } else {
        std::cout << name << " = " << defaultValue << " (default)" << std::endl;
        return defaultValue;
    }
}
