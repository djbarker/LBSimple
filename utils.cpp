#include "utils.h"

#include <memory>

#include <vtkPNGReader.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>

using namespace std;

int sub2idx(int i, int j, int Nx, int Ny)
{
	return i + j*Nx;
}

double bilinear_interp(double x, double y, double x0, double y0, double dx, double f00, double f10, double f01, double f11)
{
	double x_ = (x - x0) / dx;
	double y_ = (y - y0) / dx;
	double b1 = f00;
	double b2 = f10 - f00;
	double b3 = f01 - f00;
	double b4 = f00 - f10 - f01 + f11;
	double out = b1 + b2*x_ + b3*y_ + b4*x_*y_;

	return out;
}