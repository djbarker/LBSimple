#pragma once

#include <cmath>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <unordered_map>
#include <sys/stat.h>
#include <vector>

#include "Vect.hpp"

using namespace std;

enum CellType
{
	Fluid = 0,
	Wall = 1,
	Empty = 2,
};

// compile time factorial
constexpr size_t factorial(size_t n)
{
	return n == 0 ? 1 : n * factorial(n - 1);
}

// compile time integer power
constexpr size_t pow_int(size_t exp, size_t pow)
{
	return pow == 0 ? 1 : pow_int(exp, pow - 1)*exp;
}

// number of elements with dimension < n on the boundary of a hypercube of dimension n
constexpr size_t hc_elements(size_t n, size_t m = 0)
{
	return m == n ? 0 : hc_elements(n, m + 1) + pow_int(2, n - m)*factorial(n) / (factorial(m)*factorial(n - m));
}

template<class T>
T periodic(T idx, T period)
{
	//static_assert(std::is_unsigned<T>::value, "Attempting to use periodic clamp with unsigned type!");
	if (idx < (T)0) idx += period;
	else if (idx >= period) idx -= period;
	return idx;
}

template<class T, size_t D>
Vect<T, D> periodic(Vect<T, D> idx, Vect<T, D> period)
{
	Vect<T, D> out;
	for (int d = 0; d < D; ++d)
		out[d] = periodic(idx[d], period[d]);
	return out;
}

template<class T>
T periodic_delta(T delta, T period) {
	if (delta < -0.5*period) delta += period;
	else if (delta >= 0.5 * period) delta -= period;
	return delta;
}

template<class T, size_t D>
Vect<T, D> periodic_delta(Vect<T, D> delta, Vect<T, D> period) {
	Vect<T, D> out;
	for (int d = 0; d < D; ++d) {
		out[d] = periodic_delta(delta[d], period[d]);
	}
	return out;
}

template<size_t D>
void raster(Vect<int, D>& sub, const Vect<int, D>& extent)
{
	sub[0] += 1;
	for (int i = 0; i < D-1;++i)
		if (sub[i] == extent[i])
		{
			sub[i] = 0;
			sub[i + 1] += 1;
		}
}

template<size_t D>
Vect<int, D> raster_end(const Vect<int, D>& extent)
{
	Vect<int, D> end = extent - Vect<int,D>::ones();
	raster(end, extent);
	return end;
}

template<size_t D> int sub2idx(const Vect<int, D>& sub, const Vect<int, D>& extent);
template<size_t D> Vect<int, D> idx2sub(size_t idx, const Vect<int, D>& extent);

// In situations where we know if we have 2 or 3 dims we can use these to avoid constructing the sub_t (although maybe that doesn't help)
int sub2idx_2d(int subx, int suby, const Vect<int, 2>& extent);
int sub2idx_3d(int subx, int suby, int subz, const Vect<int, 3>& extent);

template<size_t D> Vect<int, D> calc_num_domains(size_t procs);

template<class T>
T bilinear_interp(Vect<double, 2> x, Vect<double, 2> x0, double dx, T f00, T f10, T f01, T f11)
{
	double x_ = (x[0] - x0[0]) / dx;
	double y_ = (x[1] - x0[1]) / dx;
	T b1 = f00;
	T b2 = f10 - f00;
	T b3 = f01 - f00;
	T b4 = f00 - f10 - f01 + f11;
	T out = b1 + b2*x_ + b3*y_ + b4*x_*y_;

	return out;
}

template<class T>
T trilinear_interp(Vect<double, 3> x, Vect<double, 3> x0, double dx, T f000, T f100, T f010, T f001, T f110, T f101, T f011, T f111)
{
	auto x_ = (x - x0) / dx;
	T c00 = f000*(1 - x_[0]) + f100*x_[0];
	T c10 = f010*(1 - x_[0]) + f110*x_[0];
	T c01 = f001*(1 - x_[0]) + f101*x_[0];
	T c11 = f011*(1 - x_[0]) + f111*x_[0];
	T c0 = c00*(1 - x_[1]) + c10*x_[1];
	T c1 = c01*(1 - x_[1]) + c11*x_[1];
	return c0*(1 - x_[2]) + c1*x_[1];
}

template<class T>
std::vector<T> read_raw(const std::string& fname, size_t header, int Nx, int Ny)
{
	// check size of file matches
	struct stat filestatus;
	stat(fname.c_str(), &filestatus);
	
	if ((filestatus.st_size - header) != Nx*Ny*sizeof(T))
	{
		std::stringstream msg;
		msg << "File size does not match specified size:" << endl;
		msg << '\t' << filestatus.st_size << " bytes. Expected " << header << " + " << Nx*Ny*sizeof(T) << " bytes.";
		throw std::runtime_error(msg.str());
	}

	auto out = std::vector<T>(Nx*Ny);

	ifstream fin;
	fin.open(fname.c_str());

	std::unordered_map<T, int> counts;

	// read in raw data
	Vect<int, 2> N = { Nx, Ny };
	T tmp;
	int type;
	for (int j = 0; j < Ny; ++j)
	for (int i = 0; i < Nx; ++i)
	{
		fin.read((char*)&tmp, sizeof(T));
		if (tmp == std::numeric_limits<T>::max()) {
			tmp = Fluid;
		} else if (tmp == 0) {
			tmp = Wall;
		}
		out[sub2idx(Vect<int, 2>{i, Ny - j - 1}, N)] = tmp;

		if (counts.count((int)tmp) == 0) {
			counts[(int)tmp] = 1;
		}
		
		counts[(int)tmp] += 1;
	}

	cout << "Cell summary ... " << endl;
	for (auto pair : counts) {
		cout << " " << (int)pair.first << ": " << pair.second << endl;
	}

	for (int i = 0; i < Nx;++i)
	for (int j = 0; j < Ny; ++j)
	{
		for (int di = -1; di <= 1;++di)
		for (int dj = -1; dj <= 1; ++dj)
		{
			int idx = sub2idx(Vect<int, 2>{i, j}, N);
			int nidx = sub2idx(periodic(Vect<int, 2>{i + di, j + dj}, N), N);
			
			if (out[idx] == Empty && out[nidx] == Fluid)
				out[idx] = Wall;
		}
	}

	return out;
}