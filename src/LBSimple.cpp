#include <array>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFileOutputWindow.h>
#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkZLibDataCompressor.h>

#include "json.hpp" // nlohmann

#include "utils.hpp"
#include "Model.hpp"

// select the model
#define the_model D2Q9
#define Dims 2
#define Q 9
// #define the_model D3Q19
// #define Dims 3
// #define Q 19

typedef Vect<double, Dims> vect_t;
typedef Vect<int, Dims> sub_t;

using namespace std;

template<ModelType M>
double get_feq(const vect_t& v, double rho, double c, int q)
{
	double v_dot_e = dot((Model<M>::Es[q]).template as<double>(), v);
	double v_dot_v = dot(v, v);
	double s = (3. / c)*(v_dot_e)+(9. / 2.)*(v_dot_e*v_dot_e / (c*c)) - (3. / 2.)*v_dot_v / (c*c);
	return Model<M>::Ws[q] * rho*(1.0 + s);
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

Vect<double, 3> calculate_curl(Vect<int, 2> sub, Vect<int, 2> N, double dx, const unique_ptr<Vect<double, 2>[]>& v)
{
	Vect<double, 3> curl;
	curl[2] = v[sub2idx(periodic(sub + Vect<int, 2>{1, 0}, N), N)][1] - v[sub2idx(periodic(sub + Vect<int, 2>{-1, 0}, N), N)][1]
		    - v[sub2idx(periodic(sub + Vect<int, 2>{0, 1}, N), N)][0] + v[sub2idx(periodic(sub + Vect<int, 2>{0, -1}, N), N)][0];
	curl[2] /= 2 * dx;

	return curl;
}

Vect<double, 3> calculate_curl(Vect<int, 3> sub, Vect<int, 3> N, double dx, const unique_ptr<Vect<double, 3>[]>& v)
{
	Vect<double, 3> curl;
	curl[2] = v[sub2idx(periodic(sub + Vect<int, 3>{1, 0, 0}, N), N)][1] - v[sub2idx(periodic(sub + Vect<int, 3>{-1, 0, 0}, N), N)][1]
		    - v[sub2idx(periodic(sub + Vect<int, 3>{0, 1, 0}, N), N)][0] + v[sub2idx(periodic(sub + Vect<int, 3>{0, -1, 0}, N), N)][0];

	curl[1] = v[sub2idx(periodic(sub + Vect<int, 3>{0, 0, 1}, N), N)][0] - v[sub2idx(periodic(sub + Vect<int, 3>{0, 0, -1}, N), N)][0]
		    - v[sub2idx(periodic(sub + Vect<int, 3>{1, 0, 0}, N), N)][2] + v[sub2idx(periodic(sub + Vect<int, 3>{-1, 0, 0}, N), N)][2];

	curl[0] = v[sub2idx(periodic(sub + Vect<int, 3>{0, 1, 0}, N), N)][2] - v[sub2idx(periodic(sub + Vect<int, 3>{0, -1, 0}, N), N)][2]
		    - v[sub2idx(periodic(sub + Vect<int, 3>{0, 0, 1}, N), N)][1] + v[sub2idx(periodic(sub + Vect<int, 3>{0, 0, -1}, N), N)][1];
	curl /= 2 * dx;

	return curl;
}


#include <json.hpp>
using json = nlohmann::json;
json parse_config(string fname) {
	ifstream i(fname);
	return json::parse(i);
}


typedef unsigned char RawType;


int run_main(int argc, char* argv[])
{
	// redirect vtk messages to a file
	vtkSmartPointer<vtkFileOutputWindow> vtkout = vtkSmartPointer<vtkFileOutputWindow>::New();
	vtkout->SetInstance(vtkout);
	vtkout->SetFileName("vtk_output.txt");

	ifstream ifs = ifstream(argv[1]);
	auto config = json::parse(ifs);
	cout << "Read config" << endl;

	// params - hard coded for now
	double dx = config["dx"];
	double dt = config["dt"];
	double dtout = config["dt_out"]; 
	double rho0 = config["rho0"];
	double nu = 0.001 / rho0;
	// vect_t g = { 0.0001, 0.0 };
	vect_t g = config["gravity"].get<std::array<double, 2>>();
	double perturbation = config["perturbation"]; // in percent

#if Dims == 2
	const int Nx = config["domain"]["size"][0];
	const int Ny = config["domain"]["size"][1];
	sub_t N = { Nx, Ny };
#elif Dims == 3
	const int Nx = config["domain"]["size"][0];
	const int Ny = config["domain"]["size"][1];
	const int Nz = config["domain"]["size"][2];
	sub_t N = { Nx, Ny, Nz };
#endif

	double c = dx / dt;
	double tau = 3. * nu*(1. / (c* dx)) + 0.5;
	
	const bool output_grid    = config["output"]["grid"]; // can eat lots of space
	const bool output_tracers = config["output"]["tracers"];

	unique_ptr<RawType[]> cell_type;
	if (Dims == 2 && config["domain"].count("file") > 0)
	{
		string fname = config["domain"]["file"];
		read_raw<RawType>(fname, 0, N[0], N[1], cell_type);
	}
	else
	{
		cell_type = make_unique<RawType[]>(trace(N));
	}

	vect_t L = N.as<double>()*dx;

	cout << "C = " << c << endl;
	cout << "Tau = " << tau << endl;
	cout << "N = " << N << " [" << trace(N) << "]" << endl;
	cout << "L = " << L << endl;

	// load boundaries
	cout << "Reading boundary config..." << endl;
	auto velocity_boundary = unordered_map<int, vect_t>();
	for (auto& c : config["boundaries"].items()) {
		if (c.value()["type"] == "velocity") {
			vect_t v = c.value()["vel"].get<std::array<double, Dims>>();
			int key = std::atoi(c.key().c_str());
			velocity_boundary[key] = v;

			cout << " " << key  << " -> velocity " << v << endl;
		}
	}

	// allocate memory
	auto f = make_unique<Vect<double, Q>[]>(trace(N));
	auto feq = make_unique<Vect<double, Q>[]>(trace(N));
	auto v = make_unique<vect_t[]>(trace(N));
	auto rho = make_unique<double[]>(trace(N));

	mt19937_64 rng(42);
	uniform_real_distribution<double> pert_dist(-perturbation / 100., perturbation / 100.);

	cout << "Initializing cells... ";
	cout.flush();

	// initialize arrays
	for (sub_t sub; sub != raster_end(N); raster(sub, N))
	{
		int idx = sub2idx(sub, N);

		vect_t v_ = vect_t::zero();
		vect_t x_ = sub.as<double>() / dx;

		//v_[0] = 0.015*c*cos((sub[1] * dx / L[1] - L[1] * 0.75)*M_PI*2.) + 0.025*c*cos((sub[1] * dx / L[1] - L[1] * 0.5)*M_PI * 4) + 0.025*c*cos((sub[1] * dx / L[1] - L[1] * 0.15)*M_PI * 6) + 0.01*c*cos((sub[1] * dx / L[1] - L[1] * 0.5)*M_PI * 8);
		//v_[2] = v_[1] = 0.015*c*cos((sub[0] * dx / L[0] - L[0] * 0.25)*M_PI*2.) + 0.025*c*cos((sub[0] * dx / L[0] - L[0] * 0.125)*M_PI * 4) + 0.025*c*cos((sub[0] * dx / L[0] - L[0] * 0.2)*M_PI * 6) + 0.01*c*cos((sub[0] * dx / L[0] - L[0] * 0.7)*M_PI * 8);
		
		/*int i, j, k;
		i = sub[0];
		j = sub[1];
		k = sub[2];
		if (i > 20 && i <= 30 && j > 45 && j <= 55 && k > 45 && k <= 55)
		{
			v_[0] = 0.2*c;
		}

		if (i > 70 && i <= 80 && j > 45 && j <= 55 && k > 45 && k <= 55)
		{
			v_[0] = -0.2*c;
		}*/

		if (cell_type[idx] != Fluid) v_ = vect_t::zero();

		double rho_ = 1.0;

		// if (
		// 	(sub[0] > Nx * 0.35) && (sub[0] < Nx * 0.625) &&
		//   (sub[1] > Ny * 0.45) && (sub[1] < Ny * 0.55)
		//   ) {
		// 	rho_ = 1.2;
		// 	// rho_ = 1.0;
		// 	v_ = {0.0, 0.4};
		// }

		rho_ *= rho0;

		for (int q = 0; q < Q; ++q)
		{
			feq[idx][q] = f[idx][q] = get_feq<the_model>(v_, rho_, c, q);
			if (cell_type[idx] == Fluid)
				f[idx][q] *= (1 + pert_dist(rng));
		}

		v[idx] = v_;
		rho[idx] = rho_;
	}

	cout << "Done." << endl;
	cout << "Initializing tracers... ";
	cout.flush();

	// initialize tracers
	int num_tracers = 5E4;
	auto pos = make_unique<vect_t[]>(num_tracers);
	auto vel = make_unique<vect_t[]>(num_tracers);
	auto ids = make_unique<int[]>(num_tracers);

	array<uniform_real_distribution<double>, Dims> dist;
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

	cout << "Done." << endl;

	auto interp_vel = [&](vect_t x, vect_t& v_){
#if Dims==2
			auto idx = (x / dx).as<int>();
			int idx00 = sub2idx(idx, N);
			int idx10 = sub2idx(sub_t{ periodic(idx[0] + 1, N[0]), idx[1] }, N);
			int idx01 = sub2idx(sub_t{ idx[0], periodic(idx[1] + 1, N[1]) }, N);
			int idx11 = sub2idx(sub_t{ periodic(idx[0] + 1, N[0]), periodic(idx[1] + 1, N[1]) }, N);

			v_ = bilinear_interp(x, (idx).as<double>()*dx, dx, v[idx00], v[idx10], v[idx01], v[idx11]);
#elif Dims==3
			auto idx = (x / dx).as<int>();
			int idx000 = sub2idx(idx, N);
			int idx100 = sub2idx(periodic(idx + sub_t{ 1, 0, 0 }, N), N);
			int idx010 = sub2idx(periodic(idx + sub_t{ 0, 1, 0 }, N), N);
			int idx001 = sub2idx(periodic(idx + sub_t{ 0, 0, 1 }, N), N);
			int idx110 = sub2idx(periodic(idx + sub_t{ 1, 1, 0 }, N), N);
			int idx101 = sub2idx(periodic(idx + sub_t{ 1, 0, 1 }, N), N);
			int idx011 = sub2idx(periodic(idx + sub_t{ 0, 1, 1 }, N), N);
			int idx111 = sub2idx(periodic(idx + sub_t{ 1, 1, 1 }, N), N);

			v_ = trilinear_interp(x, (idx).as<double>()*dx, dx, v[idx000], v[idx100], v[idx010], v[idx001], v[idx110], v[idx101], v[idx011], v[idx111]);
#endif
	};

	cout << "Running simulation..." << endl;
	cout << endl;

	// run the simulation
	auto then = std::clock();
	int fact = dtout / dt;
	for (int iteration = 0; iteration < fact * 1000; ++iteration)
	{
		// write out data files
		if (iteration % fact == 0)
		{
			stringstream fname;
			fname << "./out/output_" << setw(6) << setfill('0') << iteration / fact << ".vti";

			vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();

			if (output_grid) {
				// setup arrays and vtkImageData objects
				vtkSmartPointer<vtkImageData> data = vtkSmartPointer<vtkImageData>::New();
				if (Dims == 2)
				{
					data->SetExtent(0, N[0] - 1, 0, N[1] - 1, 0, 0);
					data->SetSpacing(dx*(N[0]) / (N[0] - 1), dx*(N[1]) / (N[1] - 1), 0);
				}
				else if (Dims == 3)
				{
					data->SetExtent(0, N[0] - 1, 0, N[1] - 1, 0, N[2] - 1);
					data->SetSpacing(dx*(N[0]) / (N[0] - 1), dx*(N[1]) / (N[1] - 1), dx*(N[2])/(N[2]-1));
				}
				
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
				if (Dims==2)
					curl_arr->SetNumberOfComponents(1);
				else
					curl_arr->SetNumberOfComponents(3);
				curl_arr->SetNumberOfTuples(trace(N));

				vtkSmartPointer<vtkIntArray> type_arr = vtkSmartPointer<vtkIntArray>::New();
				type_arr->SetName("CellType");
				type_arr->SetNumberOfComponents(1);
				type_arr->SetNumberOfTuples(trace(N));

				vtkPointData* pdata = data->GetPointData(); // belongs to data => no smart pointer necessary
				pdata->AddArray(vel_arr);
				pdata->AddArray(rho_arr);
				pdata->AddArray(curl_arr);
				pdata->AddArray(type_arr);
				pdata->SetActiveAttribute("Velocity", vtkDataSetAttributes::VECTORS);

				for (sub_t sub; sub != raster_end(N); raster(sub, N))
				{
					int idx = sub2idx(sub, N);
					Vect<double, 3> curl = calculate_curl(sub, N, dx, v);
					
					type_arr->SetTuple1(idx, cell_type[idx]);
					rho_arr->SetTuple1(idx, rho[idx]);
					if (Dims == 2)
					{
						vel_arr->SetTuple3(idx, v[idx][0], v[idx][1], 0.0);
						curl_arr->SetTuple1(idx, curl[2]);
					}
					else
					{
						vel_arr->SetTuple3(idx, v[idx][0], v[idx][1], v[idx][2]);
						curl_arr->SetTuple3(idx, curl[0], curl[1], curl[2]);
					}
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

			if (output_tracers) {
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
					if (Dims == 2)
					{
						points->InsertNextPoint(pos[i][0], pos[i][1], 0.0);
						vel_arr->SetTuple3(i, vel[i][0], vel[i][1], 0.0);
					}
					else
					{
						points->InsertNextPoint(pos[i][0], pos[i][1], pos[i][2]);
						vel_arr->SetTuple3(i, vel[i][0], vel[i][1], vel[i][2]);
					}
					
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
		}

		auto now = std::clock();
		cout << iteration << " (" << (int)(iteration / fact) << "): " << "<time/itr> = " << (int)((now-then)*1000.0 / (CLOCKS_PER_SEC * (iteration+1.0))) << " ms\tt = " << iteration * dt << " s\r";
		cout.flush();

		// stream - use feq as temporary storage
		for (sub_t sub; sub != raster_end(N); raster(sub, N))
		{
			const int idx = sub2idx(sub, N);
			const RawType type = cell_type[idx];
			if (type == Empty) continue;

			for (int q = 0; q < Q; ++q)
			{
				int neighbour_idx = sub2idx(periodic(sub + Model<the_model>::Es[q], N), N);
				if (cell_type[neighbour_idx] == Fluid)
				{
					feq[neighbour_idx][q] = f[idx][q];
				}
			}

			if (type == Wall)
			{
				for (int q = 0; q < Q; ++q)
				{
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
			const int idx = sub2idx(sub, N);
			const RawType type = cell_type[idx];

			if (type == Empty) continue;

			if (type > Empty) {
				if (velocity_boundary.count(type) != 0) {
					vect_t vel = velocity_boundary[type];
					for (int q=0; q<Q; q++) {
						feq[idx][q] = f[idx][q] = get_feq<the_model>(vel, rho0, c, q);
					}
				}
				continue;
			}

			for (int q = 0; q < Q; ++q)
			{
				if (type == Fluid) 
					feq[idx][q] = get_feq<the_model>(v[idx] + g*tau, rho[idx], c, q);
				else
					feq[idx][q] = get_feq<the_model>(v[idx], rho[idx], c, q);

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

			// TODO: reflect back if end up in wall

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
	catch (exception& e)
	{
		cerr << "Exception detected in LBSimple:" << endl;
		cerr << "\t" << e.what() << endl;
	}
}
