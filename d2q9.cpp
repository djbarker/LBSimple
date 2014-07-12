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

#include "utils.h"

using namespace std;

// model params
const int qs[9][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }, { 1, 1 }, { -1, 1 }, { -1, -1 }, { 1, -1 } };
const int qneg[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
const double Ws[9] = { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. };

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


int run_main(int argc, char* argv[])
{
	// redirect vtk messages to a file
	vtkSmartPointer<vtkFileOutputWindow> vtkout = vtkSmartPointer<vtkFileOutputWindow>::New();
	vtkout->SetInstance(vtkout);
	vtkout->SetFileName("vtk_output.txt");

	// params - hard coded for now
	double dx = 0.0001;
	double dt = 0.00005;
	double dtout = 0.001;
	double rho0 = 1000.;
	double nu = 0.01/rho0;
	double gx = 0.001;
	double gy = 0.0;
	double perturbation = 0.5; // in percent

	int Nx = 500;
	int Ny = 300;
	double c = dx / dt;
	double tau = 3. * nu*(1. /(c* dx)) + 0.5;

	string fname = "domain4.raw";
	unique_ptr<CellType[]> cell_type;
	read_raw<unsigned char>(fname, 0, Nx, Ny, cell_type);

	double Lx = Nx*dx;
	double Ly = Ny*dx;

	cout << "C = " << c << endl;
	cout << "Tau = " << tau << endl;
	cout << "N = (" << Nx << ", " << Ny << ")" << endl;
	cout << "L = (" << Lx << ", " << Ly << ")" << endl;

	// allocate memory
	auto f = make_unique<array<double, 9>[]>(Nx*Ny);
	auto feq = make_unique<array<double, 9>[]>(Nx*Ny);
	auto vx = make_unique<double[]>(Nx*Nx);
	auto vy = make_unique<double[]>(Nx*Nx);
	auto rho = make_unique<double[]>(Nx*Nx);

	mt19937_64 rng(42);
	uniform_real_distribution<double> pert_dist(-perturbation / 100., perturbation / 100.);

	// initialize arrays
	for (int i = 0; i < Nx;++i)
	for (int j = 0; j < Ny; ++j)
	{
		double vx_, vy_;

		//vx_ = 0.015*c*cos((j*dx / L - L*0.75)*M_PI*2.) + 0.025*c*cos((j*dx / L - L*0.5)*M_PI * 4) + 0.025*c*cos((j*dx / L - L*0.15)*M_PI * 6) + 0.01*c*cos((j*dx / L - L*0.5)*M_PI * 8);
		//vy_ = 0.015*c*cos((i*dx / L - L*0.25)*M_PI*2.) + 0.025*c*cos((i*dx / L - L*0.125)*M_PI * 4) + 0.025*c*cos((i*dx / L - L*0.2)*M_PI * 6) + 0.01*c*cos((i*dx / L - L*0.7)*M_PI * 8);
		
		vx_ = 0.0;
		vy_ = 0.0;

		int idx = sub2idx(i, j, Nx, Ny);

		if(cell_type[idx] != Fluid) vx_ = vy_ = 0.0;

		for (int q = 0; q < 9; ++q)
		{
			feq[idx][q]=f[idx][q] = get_feq(vx_, vy_, rho0, c, q);
			f[idx][q] *= (1 + pert_dist(rng));
		}

		vx[idx] = vx_;
		vy[idx] = vy_;
		rho[idx] = rho0;
	}

	// initialize tracers
	int num_tracers = 5E4;
	auto posx = make_unique<double[]>(num_tracers);
	auto posy = make_unique<double[]>(num_tracers);
	auto velx = make_unique<double[]>(num_tracers);
	auto vely = make_unique<double[]>(num_tracers);
	auto ids = make_unique<int[]>(num_tracers);

	uniform_real_distribution<double> dist_x(0.0, Lx);
	uniform_real_distribution<double> dist_y(0.0, Ly);
	bool random = true;

	for (int i = 0; i < num_tracers; ++i)
	{
		if (random)
		{
			bool valid = false;
			while (!valid)
			{
				posx[i] = dist_x(rng);
				posy[i] = dist_y(rng);

				valid = cell_type[sub2idx(posx[i] / dx, posy[i] / dx, Nx, Ny)] == Fluid;
			}
		}
		else
		{
			int N = (int)sqrt(num_tracers*(Lx / Ly));
			double dx_ = Lx / N;
			posx[i] = (i - (i / N)*N) * dx_;
			posy[i] = (i / N) * dx_;
		}

		int xidx = posx[i] / dx;
		int yidx = posy[i] / dx;
		int idx = sub2idx(xidx, yidx, Nx, Ny);
		
		velx[i] = vx[idx];
		vely[i] = vy[idx];
		ids[i] = sub2idx(posy[i] / dx, posx[i] / dx, Ny, Nx); // i/j deliberately swapped here
	}

	auto interp_vel = [&](double x, double y, double& vx_, double& vy_){
		int xidx = (int)(x / dx);
		int yidx = (int)(y / dx);
		int idx00 = sub2idx(xidx, yidx, Nx, Ny);
		int idx10 = sub2idx(periodic(xidx + 1, Nx), yidx, Nx, Ny);
		int idx01 = sub2idx(xidx, periodic(yidx + 1, Ny), Nx, Ny);
		int idx11 = sub2idx(periodic(xidx + 1, Nx), periodic(yidx + 1, Ny), Nx, Ny);

		vx_ = bilinear_interp(x, y, xidx*dx, yidx*dx, dx, vx[idx00], vx[idx10], vx[idx01], vx[idx11]);
		vy_ = bilinear_interp(x, y, xidx*dx, yidx*dx, dx, vy[idx00], vy[idx10], vy[idx01], vy[idx11]);
	};

	// run the simulation
	int fact = dtout / dt;
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
				data->SetSpacing(dx*(Nx)/(Nx-1), dx*(Ny)/(Ny-1), 0);
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
				pdata->SetActiveAttribute("Velocity", vtkDataSetAttributes::VECTORS);

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
			int idx = sub2idx(i, j, Nx, Ny);
			if (cell_type[idx] == Empty) continue;

			for (int q = 0; q < 9; ++q)
			{
				int qi = qs[q][0];
				int qj = qs[q][1];

				int neighbour_idx = sub2idx(periodic(i + qi, Nx), periodic(j + qj, Ny), Nx, Ny);
				if (cell_type[neighbour_idx] == Fluid)
				{
					feq[neighbour_idx][q] = f[idx][q];
				}				
			}

			if (cell_type[idx] == Wall)
			{
				for (int q = 0; q < 9; ++q)
				{
					int qi = qs[q][0];
					int qj = qs[q][1];

					int neighbour_idx = sub2idx(periodic(i + qi, Nx), periodic(j + qj, Ny), Nx, Ny);
					if (cell_type[neighbour_idx] != Fluid)
					{
						feq[idx][qneg[q]] = f[idx][q];
					}
				}
			}
		}
		std::swap(f, feq);

		// calc macroscopic properties
		for (int i = 0; i < Nx; ++i)
		for (int j = 0; j < Ny; ++j)
		{
			int idx = sub2idx(i, j, Nx, Ny);
			if (cell_type[idx] == Empty) continue;

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
			int idx = sub2idx(i, j, Nx, Ny);
			if (cell_type[idx] == Empty) continue;

			for (int q = 0; q < 9; ++q)
			{
				if (cell_type[idx] == Fluid)
					feq[idx][q] = get_feq(vx[idx] + gx*tau, vy[idx] + gy*tau, rho[idx], c, q);
				else
					feq[idx][q] = get_feq(vx[idx], vy[idx] , rho[idx], c, q);

				f[idx][q] = f[idx][q] - (1. / tau)*(f[idx][q] - feq[idx][q]);
			}
		}

		// advect tracers using RK4 integration
		for (int i = 0; i < num_tracers; ++i)
		{
			double vx1, vx2, vx3, vx4;
			double vy1, vy2, vy3, vy4;
						
			interp_vel(posx[i], posy[i], vx1, vy1);
			double posx2 = periodic(posx[i] + 0.5*dt*vx1, Lx);
			double posy2 = periodic(posy[i] + 0.5*dt*vy1, Ly);

			interp_vel(posx2, posy2, vx2, vy2);
			double posx3 = periodic(posx[i] + 0.5*dt*vx2, Lx);
			double posy3 = periodic(posy[i] + 0.5*dt*vy2, Ly);

			interp_vel(posx3, posy3, vx3, vy3);
			double posx4 = periodic(posx[i] + dt*vx3, Lx);
			double posy4 = periodic(posy[i] + dt*vy3, Ly);

			interp_vel(posx4, posy4, vx4, vy4);
			posx[i] = periodic(posx[i] + (1. / 6.)*dt*(vx1 + 2.*vx2 + 2.*vx3 + vx4), Lx);
			posy[i] = periodic(posy[i] + (1. / 6.)*dt*(vy1 + 2.*vy2 + 2.*vy3 + vy4), Ly);

			// strictly should be recalculated at the new position but this is for visualization only
			velx[i] = vx2; 
			vely[i] = vy2;
		}
	}

	return 0;
}

int main(int argc, char* argv[])
{
	try
	{
		run_main(argc, argv);
	}
	catch (runtime_error& e)
	{
		cerr << "Runtime error detected in LBSimple:" << endl;
		cerr << "\t" << e.what() << endl;
	}
}
