#pragma once

#include <vector>

#include "config.hpp"

void get_tracers_grid(size_t num_tracers, vect_t L, double dx, std::vector<vect_t>& pos);
void get_tracers_poisson_disk(size_t num_tracers, vect_t L, sub_t N, double dx, std::vector<RawType>& cell_type, std::vector<vect_t>& pos);