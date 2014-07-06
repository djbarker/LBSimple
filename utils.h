#pragma once

#include <string>
#include <memory>
#include <sys/stat.h>

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

int sub2idx(int i, int j, int Nx, int Ny);

double bilinear_interp(double x, double y, double x0, double y0, double dx, double f00, double f10, double f01, double f11);

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
	T tmp;
	for (int j = 0; j < Ny; ++j)
	for (int i = 0; i < Nx; ++i)
	{
		fin.read((char*)&tmp, sizeof(T));
		ptr[sub2idx(i,Ny-j-1,Nx,Ny)] = (tmp>0 ? Fluid : Empty);
	}

	for (int i = 0; i < Nx;++i)
	for (int j = 0; j < Ny; ++j)
	{
		for (int di = -1; di <= 1;++di)
		for (int dj = -1; dj <= 1; ++dj)
		{
			int idx = sub2idx(i, j, Nx, Ny);
			int nidx = sub2idx(periodic(i + di, Nx), periodic(j + dj, Ny), Nx, Ny);
			
			if (ptr[idx] == Empty && ptr[nidx] == Fluid)
				ptr[idx] = Wall;
		}
	}
}