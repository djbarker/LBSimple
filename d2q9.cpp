#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <array>
#include <memory>
#include <iomanip>
#include <cmath>

#include <vtkFileOutputWindow.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkZLibDataCompressor.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>

using namespace std;

// model params
const int qs[9][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }, { 1, 1 }, { -1, 1 }, { -1, -1 }, { 1, -1 } };
const double Ws[9] = { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. };

template<class T>
T periodic(T idx, T period)
{
	//static_assert(std::is_unsigned<T>::value, "Attempting to use periodic clamp with unsigned type!");
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

void get_vortex(double x0, double y0, double x, double y, double vx, double vy, double vmax, double sigma, bool clockwise, double Lx, double Ly, double &vxout, double& vyout)
{
	double dxp = periodic(x - x0 + Lx / 2., Lx) - Lx / 2.;
	double dyp = periodic(y - y0 + Ly / 2., Ly) - Ly / 2.;
	double r = sqrt(dxp*dxp + dyp*dyp);
	double v = (2 * r / (sigma*sigma))*exp(-r*r/(sigma*sigma));
	double theta = atan2(dyp, dxp);

	vxout = +vmax*v*sin(theta)*(clockwise ? 1 : -1) + 2*exp(-r*r/(sigma*sigma))*vmax*vx;
	vyout = -vmax*v*cos(theta)*(clockwise ? 1 : -1) + 2*exp(-r*r / (sigma*sigma))*vmax*vy;
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

template<typename T>
std::unique_ptr<T[]> make_unique_arr(size_t N)
{
	return std::unique_ptr<T[]>(new T[N]);
}

int main(int argc, char* argv[])
{
	// redirect vtk messages to a file
	vtkSmartPointer<vtkFileOutputWindow> vtkout = vtkSmartPointer<vtkFileOutputWindow>::New();
	vtkout->SetInstance(vtkout);
	vtkout->SetFileName("vtk_output.txt");

	// params - hard coded for now
	double L = 1.0;
	double dx = 0.004;
	double dt = 0.02;
	double rho0 = 1000.;
	double nu = 0.001/rho0;

	int Nx = L / dx + 1;
	int Ny = L / dx + 1;
	double c = dx / dt;
	double tau = 3. * nu*(1. /(c* dx)) + 0.5;

	cout << "C = " << c << endl;
	cout << "Tau = " << tau << endl;
	cout << "N = (" << Nx << ", " << Ny << ")" << endl;

	// allocate memory
	auto f = make_unique_arr<array<double, 9>>(Nx*Ny);
	auto feq = make_unique_arr<array<double, 9>>(Nx*Ny);
	auto vx = make_unique_arr<double>(Nx*Nx);
	auto vy = make_unique_arr<double>(Nx*Nx);
	auto rho = make_unique_arr<double>(Nx*Nx);

	// initialize tracers
	int num_tracers = 1E4;
	auto posx = make_unique_arr<double>(num_tracers);
	auto posy = make_unique_arr<double>(num_tracers);
	auto velx = make_unique_arr<double>(num_tracers);
	auto vely = make_unique_arr<double>(num_tracers);
	auto ids = make_unique_arr<int>(num_tracers);

	mt19937_64 rng(42);
	uniform_real_distribution<double> dist_x(0.0, L);
	uniform_real_distribution<double> dist_y(0.0, L);

	for (int i = 0; i < num_tracers; ++i)
	{
		posx[i] = dist_x(rng);
		posy[i] = dist_y(rng);
		velx[i] = vely[i] = 0.0;
		ids[i] = sub2idx(posy[i] / dx, posx[i] / dx, Ny, Nx); // i/j deliberately swapped here
	}

	// initialize f
	for (int i = 0; i < Nx;++i)
	for (int j = 0; j < Ny; ++j)
	{
		double vx_, vy_;

		vx_ = 0.05*c*cos((j*dx / L - L*0.75)*M_PI*2.);
		vy_ = 0.05*c*cos((i*dx / L - L*0.25)*M_PI*2.);

		int idx = sub2idx(i, j, Nx, Ny);
		for (int q = 0; q < 9; ++q)
		{
			feq[idx][q]=f[idx][q] = get_feq(vx_, vy_, rho0, c, q);
		}

		vx[idx] = vx_;
		vy[idx] = vy_;
		rho[idx] = rho0;
	}

	int fact = 25;
	for (int iteration = 0; iteration < fact*1000; ++iteration)
	{
		// write out data files
		if (iteration % fact == 0)
		{
			stringstream fname;
			fname << "./out/output_" << setw(6) << setfill('0') << iteration / fact << ".vti";

			vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();

			{
				// setup arrays and vtkImageData objects
				vtkSmartPointer<vtkImageData> data = vtkSmartPointer<vtkImageData>::New();
				data->SetExtent(0, Nx - 1, 0, Ny - 1, 0, 0);
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
				curl_arr->SetName("Vorticity");
				curl_arr->SetNumberOfComponents(1);
				curl_arr->SetNumberOfTuples(Nx*Ny);

				vtkPointData* pdata = data->GetPointData(); // belongs to data => no smart pointer necessary
				pdata->AddArray(vel_arr);
				pdata->AddArray(rho_arr);
				pdata->AddArray(curl_arr);

				for (int i = 0; i < Nx; ++i)
				for (int j = 0; j < Ny; ++j)
				{
					int idx = sub2idx(i, j, Nx, Ny);

					double curl = (vy[sub2idx(periodic(i + 1, Nx), j, Nx, Ny)] - vy[sub2idx(periodic(i - 1, Nx), j, Nx, Ny)]) - (vx[sub2idx(i, periodic(j + 1, Ny), Nx, Ny)] - vx[sub2idx(i, periodic(j - 1, Ny), Nx, Ny)]);
					curl /= (2 * dx);

					vel_arr->SetTuple3(idx, vx[idx], vy[idx], 0.0);
					rho_arr->SetTuple1(idx, rho[idx]);
					curl_arr->SetTuple1(idx, curl);
				}

				vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();

				writer->SetDataModeToBinary();
				writer->SetFileName(fname.str().c_str());
				writer->SetInputData(data);
				writer->SetCompressor(compressor);

				if (!writer->Write())
				{
					throw runtime_error("Error writing imagedata to vtk file!");
				}
			}

			fname.str("");
			fname.clear();
			fname << "./out/tracer_" << setw(6) << setfill('0') << iteration / fact << ".vtp";

			{
				vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
				vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();

				vtkSmartPointer<vtkDoubleArray> vel_arr = vtkSmartPointer<vtkDoubleArray>::New();
				vel_arr->SetName("Velocity");
				vel_arr->SetNumberOfComponents(3);
				vel_arr->SetNumberOfTuples(num_tracers);

				vtkSmartPointer<vtkDoubleArray> id_arr = vtkSmartPointer<vtkDoubleArray>::New();
				id_arr->SetName("Id");
				id_arr->SetNumberOfComponents(1);
				id_arr->SetNumberOfTuples(num_tracers);

				for (int i = 0; i < num_tracers; ++i)
				{
					points->InsertNextPoint(posx[i], posy[i], 0.0);
					vel_arr->SetTuple3(i, velx[i], vely[i], 0.0);
					id_arr->SetTuple1(i, ids[i]);
					vtkIdType id[1] = { i };
					vertices->InsertNextCell(1, id);
				}

				vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
				polydata->Allocate(num_tracers);
				polydata->SetPoints(points);
				polydata->SetVerts(vertices);
				polydata->GetPointData()->AddArray(vel_arr);
				polydata->GetPointData()->AddArray(id_arr);

				vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
				writer->SetFileName(fname.str().c_str());
				writer->SetInputData(polydata);
				writer->SetCompressor(compressor);

				if (!writer->Write())
				{
					throw runtime_error("Error writing tracers to vtk file!");
				}
			}

			cout << "iteration " << iteration << "\tt = " << iteration *dt << " s" << endl;
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

		// advect tracers using improved Euler method
		for (int i = 0; i < num_tracers; ++i)
		{
			int xidx = (int)(posx[i] / dx);
			int yidx = (int)(posy[i] / dx);
			int idx00 = sub2idx(xidx, yidx, Nx, Ny);
			int idx10 = sub2idx(xidx+1, yidx, Nx, Ny);
			int idx01 = sub2idx(xidx, yidx+1, Nx, Ny);
			int idx11 = sub2idx(xidx+1, yidx+1, Nx, Ny);

			double vx1, vx2, vy1, vy2;
			vx1 = bilinear_interp(posx[i], posy[i], xidx*dx, yidx*dx, dx, vx[idx00], vx[idx10], vx[idx01], vx[idx11]);
			vy1 = bilinear_interp(posx[i], posy[i], xidx*dx, yidx*dx, dx, vy[idx00], vy[idx10], vy[idx01], vy[idx11]);
			
			// initial guess
			double posx1 = periodic(posx[i] + dt*vx1, L);
			double posy1 = periodic(posy[i] + dt*vy1, L);

			// could have moved out of the grid cell => recalculate idxs
			xidx = (int)(posx1 / dx);
			yidx = (int)(posy1 / dx);
			idx00 = sub2idx(xidx, yidx, Nx, Ny);
			idx10 = sub2idx(xidx + 1, yidx, Nx, Ny);
			idx01 = sub2idx(xidx, yidx + 1, Nx, Ny);
			idx11 = sub2idx(xidx + 1, yidx + 1, Nx, Ny);

			vx2 = bilinear_interp(posx1, posy1, xidx*dx, yidx*dx, dx, vx[idx00], vx[idx10], vx[idx01], vx[idx11]);
			vy2 = bilinear_interp(posx1, posy1, xidx*dx, yidx*dx, dx, vy[idx00], vy[idx10], vy[idx01], vy[idx11]);

			// corrector step
			posx[i] = periodic(posx[i] + 0.5*dt*(vx1 + vx2), L);
			posy[i] = periodic(posy[i] + 0.5*dt*(vy1 + vy2), L);

			// strictly should be recalculated at the new position but this is for visualization only
			velx[i] = vx2; 
			vely[i] = vy2;
		}
	}

	return 0;
}