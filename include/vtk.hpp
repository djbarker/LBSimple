#pragma once

#include <memory>
#include <string>

#include <include/config.hpp>

using grid_vect_t = std::unique_ptr<vect_t[]>;
using grid_double_t = std::unique_ptr<double[]>;
using grid_cell_t = std::unique_ptr<RawType[]>;

void write_grid(std::string fname, double dx, sub_t N, grid_cell_t& cell_type, grid_vect_t& v, grid_double_t& rho);
void write_tracers(std::string fname, size_t num_tracers, grid_vect_t& pos, grid_vect_t& vel, grid_double_t& id);