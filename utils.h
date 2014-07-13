#pragma once

#include <string>
#include <memory>
#include <iostream>
#include <sys/stat.h>

#include "Vect.hpp"

enum CellType
{
	Fluid,
	Wall,
	Empty
};

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

template<size_t D>
void raster(Vect<int, D>& sub, const Vect<int, D>& extent);

template<size_t D>
Vect<int, 2> raster_end(const Vect<int, D>& extent);

template<size_t D>
int sub2idx(const Vect<int, D>& sub, const Vect<int, D>& extent);

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
void read_raw(const std::string& fname, size_t header, int Nx, int Ny, std::unique_ptr<CellType[]>& ptr)
{
	// check size of file matches
	struct stat filestatus;
	stat(fname.c_str(), &filestatus);
	
	if ((filestatus.st_size - header) != Nx*Ny*sizeof(T))
	{
		stringstream msg;
		msg << "File size does not match specified size:" << endl;
		msg << '\t' << filestatus.st_size << " bytes. Expected " << header << " + " << Nx*Ny*sizeof(T) << " bytes.";
		throw runtime_error(msg.str());
	}

	ptr = make_unique<CellType[]>(Nx*Ny);

	ifstream fin;
	fin.open(fname.c_str());

	// read in raw data
	Vect<int, 2> N = { Nx, Ny };
	T tmp;
	for (int j = 0; j < Ny; ++j)
	for (int i = 0; i < Nx; ++i)
	{
		fin.read((char*)&tmp, sizeof(T));
		ptr[sub2idx(Vect<int, 2>{i, Ny - j - 1}, N)] = (tmp>0 ? Fluid : Empty);
	}

	for (int i = 0; i < Nx;++i)
	for (int j = 0; j < Ny; ++j)
	{
		for (int di = -1; di <= 1;++di)
		for (int dj = -1; dj <= 1; ++dj)
		{
			int idx = sub2idx(Vect<int, 2>{i, j}, N);
			int nidx = sub2idx(periodic(Vect<int, 2>{i + di, j + dj}, N), N);
			
			if (ptr[idx] == Empty && ptr[nidx] == Fluid)
				ptr[idx] = Wall;
		}
	}
}