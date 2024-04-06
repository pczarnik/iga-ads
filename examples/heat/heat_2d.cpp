// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "heat_2d.hpp"
#include <string>

int main(int argc, char* argv[]) {
    double dt = 1e-5;
    if (argc > 1) {
        dt = std::stod(argv[1]);
    }

    ads::dim_config dim{2, 40};
    ads::timesteps_config steps{10000, dt};
    int ders = 1;

    ads::config_2d c{dim, dim, steps, ders};
    ads::problems::heat_2d sim{c};
    sim.run();
}
