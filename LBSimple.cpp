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
#include "Model.hpp"

using namespace std;

template<class Model>
double get_feq(const Vect<double, 2>& v, double rho, double c, int q)
{
	double v_dot_e = dot((Model::Es[q]).as<double>(), v);
	double v_dot_v = dot(v, v);
	double s = (3. / c)*(v_dot_e)+(9. / 2.)*(v_dot_e*v_dot_e / (c*c)) - (3. / 2.)*v_dot_v / (c*c);
	return Model::Ws[q] * rho*(1.0 + s);
}

void get_vortex(double x0, double y0, double x, double y, double vx, double vy, double vmax, double sigma, bool clockwise, double Lx, double Ly, double &vxout, double& vyout)
{
	double dxp = periodic(x - x0 + Lx / 2., Lx) - Lx / 2.;
	double dyp = periodic(y - y0 + Ly / 2., Ly) - Ly / 2.;
	double r = sqrt(dxp*dxp + dyp*dyp);
	double v = (2 * r / (sigma*sigma))*exp(-r*r / (sigma*sigma));
	double theta = atan2(dyp, dxp);

	vxout = +vmax*v*sin(theta)*(clockwise ? 1 : -1) + 2 * exp(-r*r / (sigma*sigma))*vmax*vx;
	vyout = -vmax*v*cos(theta)*(clockwise ? 1 : -1) + 2 * exp(-r*r / (sigma*sigma))*vmax*vy;
}

// select the model
#define the_model D2Q9
#define Dims 2
#define Q 9

typedef Vect<double, Dims> vect_t;
typedef Vect<int, Dims> sub_t;

int run_main(int argc, char* argv[])
{
	// redirect vtk messages to a file
	vtkSmartPointer<vtkFileOutputWindow> vtkout = vtkSmartPointer<vtkFileOutputWindow>::New();
	vtkout->SetInstance(vtkout);
	vtkout->SetFileName("vtk_output.txt");

	// params - hard coded for now
	double dx = 0.0001;
	double dt = 0.0001;
	double dtout = 0.005;
	double rho0 = 1000.;
	double nu = 0.001 / rho0;
	vect_t g = { 1E-4, 0.0 };
	double perturbation = 0.5; // in percent

	sub_t N = { 500, 300 };
	double c = dx / dt;
	double tau = 3. * nu*(1. / (c* dx)) + 0.5;

	string fname = "domain4.raw";
	unique_ptr<CellType[]> cell_type;
	read_raw<unsigned char>(fname, 0, N[0], N[1], cell_type);

	vect_t L = N.as<double>()*dx;

	cout << "C = " << c << endl;
	cout << "Tau = " << tau << endl;
	cout << "N = (" << N[0] << ", " << N[1] << ")" << endl;
	cout << "L = (" << L[0] << ", " << L[1] << ")" << endl;

	// allocate memory
	auto f = make_unique<Vect<double, Q>[]>(trace(N));
	auto feq = make_unique<Vect<double, Q>[]>(trace(N));
	auto v = make_unique<vect_t[]>(trace(N));
	auto rho = make_unique<double[]>(trace(N));

	mt19937_64 rng(42);
	uniform_real_distribution<double> pert_dist(-perturbation / 100., perturbation / 100.);

	cout << "Initializing cells..." << endl;

	// initialize arrays
	for (sub_t sub; sub != raster_end(N); raster(sub, N))
	{
		int idx = sub2idx(sub, N);

		vect_t v_ = { 0.0, 0.0 };

		//v_[0] = 0.015*c*cos((j*dx / L - L*0.75)*M_PI*2.) + 0.025*c*cos((j*dx / L - L*0.5)*M_PI * 4) + 0.025*c*cos((j*dx / L - L*0.15)*M_PI * 6) + 0.01*c*cos((j*dx / L - L*0.5)*M_PI * 8);
		//v_[1] = 0.015*c*cos((i*dx / L - L*0.25)*M_PI*2.) + 0.025*c*cos((i*dx / L - L*0.125)*M_PI * 4) + 0.025*c*cos((i*dx / L - L*0.2)*M_PI * 6) + 0.01*c*cos((i*dx / L - L*0.7)*M_PI * 8);

		if (cell_type[idx] != Fluid) v_ = vect_t(); // zero

		for (int q = 0; q < 9; ++q)
		{
			feq[idx][q] = f[idx][q] = get_feq<Model<D2Q9>>(v_, rho0, c, q);
			f[idx][q] *= (1 + pert_dist(rng));
		}

		v[idx] = v_;
		rho[idx] = rho0;
	}

	cout << "Initializing tracers..." << endl;

	// initialize tracers
	int num_tracers = 5E4;
	auto pos = make_unique<vect_t[]>(num_tracers);
	auto vel = make_unique<vect_t[]>(num_tracers);
	auto ids = make_unique<int[]>(num_tracers);

	array<uniform_real_distribution<double>, 2> dist;
	for (int d = 0; d < Dims; ++d)
		dist[d] = uniform_real_distribution<double>(0, L[d]);

	bool random = true;

	for (int i = 0; i < num_tracers; ++i)
	{
		if (random)
		{
			bool valid = false;
			while (!valid)
			{
				for (int d = 0; d < Dims; ++d)
					pos[i][d] = dist[d](rng);

				valid = cell_type[sub2idx((pos[i] / dx).as<int>(), N)] == Fluid;
			}
		}
		else
		{
			if (Dims == 2)
			{
				int M = (int)sqrt(num_tracers*(L[0] / L[1]));
				double dx_ = L[0] / M;
				pos[i][0] = (i - (i / M)*M) * dx_;
				pos[i][1] = (i / M) * dx_;
			}
			else if (Dims == 3)
			{
				// TODO: implement
			}

		}

		int idx = sub2idx((pos[i] / dx).as<int>(), N);

		for (int d = 0; d < Dims; ++d)
			vel[i] = v[idx];

		ids[i] = idx;
	}

	auto interp_vel = [&](vect_t x, vect_t& v_){
		if (Dims == 2)
		{
			auto idx = (x / dx).as<int>();
			int idx00 = sub2idx(idx, N);
			int idx10 = sub2idx(sub_t{ periodic(idx[0] + 1, N[0]), idx[1] }, N);
			int idx01 = sub2idx(sub_t{ idx[0], periodic(idx[1] + 1, N[1]) }, N);
			int idx11 = sub2idx(sub_t{ periodic(idx[0] + 1, N[0]), periodic(idx[1] + 1, N[1]) }, N);

			v_ = bilinear_interp(x, (idx).as<double>()*dx, dx, v[idx00], v[idx10], v[idx01], v[idx11]);
		}
		else
		{
			// TODO: implement
		}
	};

	// run the simulation
	int fact = dtout / dt;
	for (int iteration = 0; iteration < fact * 1000; ++iteration)
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
				data->SetExtent(0, N[0] - 1, 0, N[1] - 1, 0, 0);
				data->SetSpacing(dx*(N[0]) / (N[0] - 1), dx*(N[1]) / (N[1] - 1), 0);
				data->SetOrigin(0, 0, 0);

				vtkSmartPointer<vtkDoubleArray> rho_arr = vtkSmartPointer<vtkDoubleArray>::New();
				rho_arr->SetName("Density");
				rho_arr->SetNumberOfComponents(1);
				rho_arr->SetNumberOfTuples(trace(N));

				vtkSmartPointer<vtkDoubleArray> vel_arr = vtkSmartPointer<vtkDoubleArray>::New();
				vel_arr->SetName("Velocity");
				vel_arr->SetNumberOfComponents(3);
				vel_arr->SetNumberOfTuples(trace(N));

				vtkSmartPointer<vtkDoubleArray> curl_arr = vtkSmartPointer<vtkDoubleArray>::New();
				curl_arr->SetName("Vorticity");
				curl_arr->SetNumberOfComponents(1);
				curl_arr->SetNumberOfTuples(trace(N));

				vtkSmartPointer<vtkDoubleArray> type_arr = vtkSmartPointer<vtkDoubleArray>::New();
				type_arr->SetName("CellType");
				type_arr->SetNumberOfComponents(1);
				type_arr->SetNumberOfTuples(trace(N));

				vtkPointData* pdata = data->GetPointData(); // belongs to data => no smart pointer necessary
				pdata->AddArray(vel_arr);
				pdata->AddArray(rho_arr);
				pdata->AddArray(curl_arr);
				pdata->AddArray(type_arr);
				pdata->SetActiveAttribute("Velocity", vtkDataSetAttributes::VECTORS);

				for (int i = 0; i < N[0]; ++i)
				for (int j = 0; j < N[1]; ++j)
				{
					int idx = sub2idx(sub_t{ i, j }, N);

					double curl = (v[sub2idx(sub_t{ periodic(i + 1, N[0]), j }, N)][1] - v[sub2idx(sub_t{ periodic(i - 1, N[0]), j }, N)][1])
						- (v[sub2idx(sub_t{ i, periodic(j + 1, N[1]) }, N)][0] - v[sub2idx(sub_t{ i, periodic(j - 1, N[1]) }, N)][0]);
					curl /= (2 * dx);

					vel_arr->SetTuple3(idx, v[idx][0], v[idx][1], 0.0);
					rho_arr->SetTuple1(idx, rho[idx]);
					curl_arr->SetTuple1(idx, curl);
					type_arr->SetTuple1(idx, cell_type[idx]);
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
					points->InsertNextPoint(pos[i][0], pos[i][1], 0.0);
					vel_arr->SetTuple3(i, vel[i][0], vel[i][1], 0.0);
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
		for (sub_t sub; sub != raster_end(N); raster(sub, N))
		{
			int idx = sub2idx(sub, N);
			if (cell_type[idx] == Empty) continue;

			for (int q = 0; q < Q; ++q)
			{
				int neighbour_idx = sub2idx(periodic(sub + Model<the_model>::Es[q], N), N);
				if (cell_type[neighbour_idx] != Empty)
				{
					feq[neighbour_idx][q] = f[idx][q];
				}
			}

			if (cell_type[idx] == Wall)
			{
				for (int q = 0; q < Q; ++q)
				{
					int qi = Model<the_model>::Es[q][0];
					int qj = Model<the_model>::Es[q][1];

					int neighbour_idx = sub2idx(periodic(sub + Model<the_model>::Es[q], N), N);
					if (cell_type[neighbour_idx] != Fluid)
					{
						feq[idx][Model<the_model>::Qneg[q]] = f[idx][q];
					}
				}
			}
		}
		std::swap(f, feq);

		// calc macroscopic properties
		for (sub_t sub; sub != raster_end(N); raster(sub, N))
		{
			int idx = sub2idx(sub, N);
			if (cell_type[idx] == Empty) continue;

			rho[idx] = 0.0;
			for (int q = 0; q < Q; ++q)
				rho[idx] += f[idx][q];

			v[idx] = { 0.0, 0.0 };
			for (int q = 0; q < Q; ++q)
			{
				v[idx] += (Model<the_model>::Es[q]).as<double>() * f[idx][q];
			}
			v[idx] *= c / rho[idx];
		}

		// calculate equilibrium distribution & update f
		for (sub_t sub; sub != raster_end(N); raster(sub, N))
		{
			int idx = sub2idx(sub, N);
			if (cell_type[idx] == Empty) continue;

			for (int q = 0; q < Q; ++q)
			{
				if (cell_type[idx] == Fluid)
					feq[idx][q] = get_feq<Model<the_model>>(v[idx] + g*tau, rho[idx], c, q);
				else
					feq[idx][q] = get_feq<Model<the_model>>(v[idx], rho[idx], c, q);

				f[idx][q] = f[idx][q] - (1. / tau)*(f[idx][q] - feq[idx][q]);
			}
		}

		// advect tracers using RK4 integration
		for (int i = 0; i < num_tracers; ++i)
		{
			vect_t v1, v2, v3, v4;

			interp_vel(pos[i], v1);
			auto pos2 = periodic(pos[i] + 0.5*dt*v1, L);

			interp_vel(pos2, v2);
			auto pos3 = periodic(pos[i] + 0.5*dt*v2, L);

			interp_vel(pos3, v3);
			auto pos4 = periodic(pos[i] + dt*v3, L);

			interp_vel(pos4, v4);
			pos[i] = periodic(pos[i] + (1. / 6.)*dt*(v1 + 2.*v2 + 2.*v3 + v4), L);

			// strictly should be recalculated at the new position but this is for visualization only
			vel[i] = v1;
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
