#pragma once

#include <iostream>
#include <toml++/toml.h>

#include "../boundary/boundary-types.h"
#include "../data-types/direction.h"

template <class T>
T parseValue(const toml::node_view<const toml::node>& config, const std::string& name,
             const T& defaultValue) {
    const std::optional<T> maybeValue = config[name].value<T>();
    if (maybeValue.has_value()) {
        const T& value = maybeValue.value();
        std::cout << name << " = " << value << std::endl;
        return value;
    } else {
        std::cout << name << " = " << defaultValue << " (default)" << std::endl;
        return defaultValue;
    }
}

template <>
BoundaryType parseValue<BoundaryType>(const toml::node_view<const toml::node>& config,
                                      const std::string& name, const BoundaryType& defaultValue);

template <>
Direction parseValue<Direction>(const toml::node_view<const toml::node>& config,
                                const std::string& name, const Direction& defaultValue);
