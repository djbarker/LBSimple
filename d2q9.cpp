#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <iomanip>
#include <cmath>

#include <vtkFileOutputWindow.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkZLibDataCompressor.h>
#include <vtkXMLImageDataWriter.h>

using namespace std;

// model params
const int qs[9][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }, { 1, 1 }, { -1, 1 }, { -1, -1 }, { 1, -1 } };
const double Ws[9] = { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. };


int periodic(int idx, int period)
{
	if (idx < 0) idx += period;
	else if (idx >= period) idx -= period;
	return idx;
}

int sub2idx(int i, int j,int Nx,int Ny)
{
	return i + j*Nx;
}

double get_feq(double vx, double vy, double rho, double c, int q)
{
	double v_dot_e = vx*qs[q][0] + vy*qs[q][1];
	double v_dot_v = vx*vx + vy*vy;
	double s = (3. / c)*(v_dot_e)+(9. / 2.)*(v_dot_e*v_dot_e / (c*c)) - (3. / 2.)*v_dot_v / (c*c);
	return Ws[q] * rho*(1.0 + s);
}

void get_vortex(double x0, double y0, double x, double y, double vx, double vy, double vmax, double sigma, bool clockwise, double &vxout, double& vyout)
{
	double r = sqrt((x0 - x)*(x0 - x) + (y0 - y)*(y0 - y));
	double v = (2 * r / (sigma*sigma))*exp(-r*r/(sigma*sigma));
	double theta = atan2(y - y0, x - x0);

	vxout = +vmax*v*sin(theta)*(clockwise ? 1 : -1) + 2*exp(-r*r/(sigma*sigma))*vmax*vx;
	vyout = -vmax*v*cos(theta)*(clockwise ? 1 : -1) + 2*exp(-r*r / (sigma*sigma))*vmax*vy;
}

#define Nx 401
#define Ny 401

int main(int argc, char* argv[])
{
	// redirect vtk messages to a file
	vtkSmartPointer<vtkFileOutputWindow> vtkout = vtkSmartPointer<vtkFileOutputWindow>::New();
	vtkout->SetInstance(vtkout);
	vtkout->SetFileName("vtk_output.txt");

	// params - hard coded for now
	double dx = 0.005;
	double dt = 0.02;
	double rho0 = 1000.;
	double nu = 0.001/rho0;

	double c = dx / dt;
	double tau = 3. * nu*(1. /(c* dx)) + 0.5;

	cout << "C = " << c << endl;
	cout << "Tau = " << tau << endl;

	// allocate memory
	array<double, 9>* f   = new array<double,9>[Nx*Ny];
	array<double, 9>* feq = new array<double, 9>[Nx*Ny];
	double* vx = new double[Nx*Ny];
	double* vy = new double[Nx*Ny];
	double* rho = new double[Nx*Ny];

	// initialize f
	for (int i = 0; i < Nx;++i)
	for (int j = 0; j < Ny; ++j)
	{
		double vx_, vy_;
		vy_ = vx_ = 0.0;
		
		double vx__, vy__;
		get_vortex(100 * dx, 210 * dx, i*dx, j*dx, 10,0, 0.001*c, 0.066, false, vx__, vy__);
		vx_ += vx__;
		vy_ += vy__;
		get_vortex(100 * dx, 190 * dx, i*dx, j*dx, 10,0, 0.001*c, 0.066, true, vx__, vy__);
		vx_ += vx__;
		vy_ += vy__;
		get_vortex(300 * dx, 210 * dx, i*dx, j*dx, -10, 0, 0.001*c, 0.066, true, vx__, vy__);
		vx_ += vx__;
		vy_ += vy__;
		get_vortex(300 * dx, 190 * dx, i*dx, j*dx, -10, 0, 0.001*c, 0.066, false, vx__, vy__);
		vx_ += vx__;
		vy_ += vy__;

		int idx = sub2idx(i, j, Nx, Ny);
		for (int q = 0; q < 9; ++q)
		{
			feq[idx][q]=f[idx][q] = get_feq(vx_, vy_, rho0, c, q);
		}

		vx[idx] = vx_;
		vy[idx] = vy_;
		rho[idx] = rho0;
	}

	int fact = 50;
	for (int iteration = 0; iteration < fact*1000; ++iteration)
	{
		// write out velocities
		if (iteration % fact == 0)
		{
			stringstream fname;
			fname << "./out/output_" << setw(6) << setfill('0') << iteration / fact << ".vti";
			
			// setup arrays and vtkImageData objects
			vtkSmartPointer<vtkImageData> data = vtkSmartPointer<vtkImageData>::New();
			data->SetExtent(0, Nx-1, 0, Ny-1, 0, 0);
			data->SetSpacing(dx, dx, 0);
			data->SetOrigin(0, 0, 0);

			vtkSmartPointer<vtkDoubleArray> rho_arr = vtkSmartPointer<vtkDoubleArray>::New();
			rho_arr->SetName("Density");
			rho_arr->SetNumberOfComponents(1);
			rho_arr->SetNumberOfTuples(Nx*Ny);

			vtkSmartPointer<vtkDoubleArray> vel_arr = vtkSmartPointer<vtkDoubleArray>::New();
			vel_arr->SetName("Velocity");
			vel_arr->SetNumberOfComponents(3);
			vel_arr->SetNumberOfTuples(Nx*Ny);

			vtkSmartPointer<vtkDoubleArray> curl_arr = vtkSmartPointer<vtkDoubleArray>::New();
			curl_arr->SetName("Curl");
			curl_arr->SetNumberOfComponents(1);
			curl_arr->SetNumberOfTuples(Nx*Ny);

			vtkPointData* pdata = data->GetPointData(); // belongs to data => no smart pointer necessary
			pdata->AddArray(vel_arr);
			pdata->AddArray(rho_arr);
			pdata->AddArray(curl_arr);

			for (int i = 0; i < Nx;++i)
			for (int j = 0; j < Ny; ++j)
			{
				int idx = sub2idx(i, j, Nx, Ny);

				double curl = (vy[sub2idx(periodic(i + 1, Nx), j, Nx, Ny)] - vy[sub2idx(periodic(i - 1, Nx), j, Nx, Ny)]) - (vx[sub2idx(i, periodic(j + 1, Ny), Nx, Ny)] - vx[sub2idx(i, periodic(j - 1, Ny), Nx, Ny)]);
				curl /= (2 * dx);

				vel_arr->SetTuple3(idx, vx[idx], vy[idx], 0.0);
				rho_arr->SetTuple1(idx, rho[idx]);
				curl_arr->SetTuple1(idx, curl);
			}

			vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();
			vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();

			writer->SetDataModeToBinary();
			writer->SetFileName(fname.str().c_str());
			writer->SetInputData(data);
			writer->SetCompressor(compressor);

			if (!writer->Write())
			{
				throw runtime_error("Error writing imagedata to vtk file!");
			}

			cout << iteration << endl;
		}

		// stream - use feq as temporary storage
		for (int i = 0; i < Nx; ++i)
		for (int j = 0; j < Ny; ++j)
		{
			for (int q = 0; q < 9; ++q)
			{
				int qi = qs[q][0];
				int qj = qs[q][1];

				feq[sub2idx(periodic(i + qi, Nx),
					        periodic(j + qj, Ny),Nx,Ny)][q] = f[sub2idx(i,j,Nx,Ny)][q];
			}
		}
		std::swap(f, feq);

		// calc macroscopic properties
		for (int i = 0; i < Nx; ++i)
		for (int j = 0; j < Ny; ++j)
		{
			int idx = sub2idx(i, j, Nx, Ny);

			rho[idx] = 0.0;
			for (int q = 0; q < 9; ++q)
				rho[idx] += f[idx][q];

			vx[idx] = vy[idx] = 0.0;
			for (int q = 0; q < 9; ++q)
			{
				int qi = qs[q][0];
				int qj = qs[q][1];

				vx[idx] += qi*f[idx][q];
				vy[idx] += qj*f[idx][q];
			}
			vx[idx] *= c / rho[idx];
			vy[idx] *= c / rho[idx];
		}

		// calculate equilibrium distribution & update f
		for (int i = 0; i < Nx; ++i)
		for (int j = 0; j < Ny; ++j)
		{
			for (int q = 0; q < 9; ++q)
			{
				feq[sub2idx(i, j, Nx, Ny)][q] = get_feq(vx[sub2idx(i, j, Nx, Ny)], vy[sub2idx(i, j, Nx, Ny)],rho[sub2idx(i,j,Nx,Ny)],c,q);

				f[sub2idx(i,j,Nx,Ny)][q] = f[sub2idx(i,j,Nx,Ny)][q] - (1. / tau)*(f[sub2idx(i,j,Nx,Ny)][q] - feq[sub2idx(i,j,Nx,Ny)][q]);
			}
		}
	}

	return 0;
}
