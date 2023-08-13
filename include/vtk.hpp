#pragma once

#include <memory>
#include <string>
#include <vector>

#include <include/config.hpp>

using grid_vect_t = std::vector<vect_t>;
using grid_double_t = std::vector<double>;
using grid_cell_t = std::vector<RawType>;

void write_grid(std::string fname, double dx, sub_t N, grid_cell_t &cell_type, grid_vect_t &v, grid_double_t &rho);
void write_tracers(std::string fname, size_t num_tracers, std::vector<vect_t> &pos, std::vector<vect_t> &vel, std::vector<double> &id, std::vector<vect_t> &pos_init, std::vector<double> &colour);