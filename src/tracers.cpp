#include <cmath>
#include <vector>
#include <random>
#include <iostream>

#include "config.hpp"
#include "utils.hpp"

using namespace std;

void get_tracers_grid(size_t num_tracers, vect_t L, double dx, vector<vect_t>& pos) {
    for (int i=0; i<num_tracers; ++i) {
        if (Dims == 2)
        {
            int M = (int)sqrt(num_tracers*(L[0] / L[1]));
            double dx_ = L[0] / M;
            pos.push_back(vect_t {
                (i - (i / M)*M) * dx_,
                (i / M) * dx_,
            });
        }
        else if (Dims == 3)
        {
            // TODO: implement
        }
    }
}

//TODO: bundle L, N, dx into a "domain" struct or something
void get_tracers_poisson_disk(size_t num_tracers, vect_t L, sub_t N, double dx, vector<RawType>& cell_type, vector<vect_t>& pos) {
    
	// (Lx/(f*r)) * (Ly/(f*r)) = N  => (f*r)^2 = Lx*Ly/N  => r = sqrt(Lx*Ly/N) / f 
    double r_pois = sqrt(L[0] * L[1] / num_tracers) / 1.3;
    double dx_pois = r_pois / sqrt(2.0);
	sub_t N_pois = (L / dx_pois).as<int>();
	int k_pois = 50;

	auto rng = mt19937_64(1);
	auto uniform = uniform_real_distribution<double>(0, 1.0);
    double eps = 0.0;
	
    // initialize data structures
	auto active_idx = vector<int>();
	auto grid_idx = make_unique<int[]>(trace(N_pois));
	for (int idx = 0; idx < trace(N_pois); ++idx) {
		grid_idx[idx] = -1; 
	}   

    // func for adding points to the data structures
    auto add_point = [&](vect_t x) {
        active_idx.push_back(pos.size());
        grid_idx[sub2idx((x/dx_pois).as<int>(), N_pois)] = pos.size();
        pos.push_back(x);
    };

    // uniform seed point
    vect_t seed = vect_t { uniform(rng),  uniform(rng) };
    seed = seed * (L + dx * vect_t::ones());
    add_point(seed);

    while (active_idx.size() >  0) {

        int aidx = (int) (uniform(rng) * active_idx.size());
        int idx = active_idx[aidx];
        sub_t sub = (pos[idx] / dx_pois).as<int>();

        bool valid = false;
        for (int j=0; j < k_pois; ++j) {
            double radius = r_pois * (1 + uniform(rng));
            double theta = 2 * M_PI * uniform(rng);

            vect_t pos2 = pos[idx] + vect_t {
                cos(theta),
                sin(theta)
            } * radius;

            pos2 = periodic(pos2, L);
            sub_t sub2 = (pos2/dx_pois).as<int>();
            int idx2 = sub2idx(sub2, N_pois);

            // NOTE: If Lx and Ly were exactly divisible by r_pois, then this would not happen.
            if (grid_idx[idx2] >= 0) {
                continue;
            }

            // check valid point - cell type
            auto type = cell_type[sub2idx((pos2/dx).as<int>(), N)];
            if (type == Empty || type == Wall) {
                continue;
            }

            // check valid point - distance
            bool valid2 = true;
            for (int ix=-1; ix <= 1; ++ix) {
                for (int iy=-1; iy <= 1; ++iy) {
                    sub_t sub3 = sub2 + sub_t { ix, iy };
                    sub3 = periodic(sub3, N_pois);

                    int idx3 = sub2idx(sub3, N_pois);
                    idx3 = grid_idx[idx3];

                    if (idx3 >= 0) {
                        vect_t diff = periodic_delta(pos[idx3] - pos2, L);
                        if (dot(diff, diff) < (r_pois * r_pois) - eps) {
                            valid2 = false;
                        }
                    }
                    
                    if (!valid2) break;
                }

                if (!valid2) break;
            }

            if (valid2) {
                add_point(pos2);
                valid = true;
                break;  // loop over j
            }
        }

        if (!valid) {
            active_idx.erase(active_idx.begin() + aidx);
        }
    };

}