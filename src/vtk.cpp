#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkZLibDataCompressor.h>

#include <include/vtk.hpp>
#include <include/utils.hpp>

using namespace std;


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



void write_grid(string fname, double dx, sub_t N,  grid_cell_t& cell_type, grid_vect_t& v, grid_double_t& rho) {

	vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();

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
	writer->SetFileName(fname.c_str());
	writer->SetInputData(data);
	writer->SetCompressor(compressor);

	if (!writer->Write())
	{
		throw runtime_error("Error writing imagedata to vtk file!");
	}
}

void write_tracers(string fname, size_t num_tracers, vector<vect_t>& pos, vector<vect_t>& vel, vector<double>& ids, vector<vect_t>& pos_init, vector<double>& colour) {

	vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();

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

	vtkSmartPointer<vtkDoubleArray> pos_init_arr = vtkSmartPointer<vtkDoubleArray>::New();
	pos_init_arr->SetName("InitialPos");
	pos_init_arr->SetNumberOfComponents(3);
	pos_init_arr->SetNumberOfTuples(num_tracers);

	vtkSmartPointer<vtkDoubleArray> colour_arr = vtkSmartPointer<vtkDoubleArray>::New();
	colour_arr->SetName("Colour");
	colour_arr->SetNumberOfComponents(1);
	colour_arr->SetNumberOfTuples(num_tracers);

	for (int i = 0; i < num_tracers; ++i)
	{
		if (Dims == 2)
		{
			points->InsertNextPoint(pos[i][0], pos[i][1], 0.0);
			vel_arr->SetTuple3(i, vel[i][0], vel[i][1], 0.0);		
			pos_init_arr->SetTuple3(i, pos_init[i][0], pos_init[i][1], 0.0);
		}
		else
		{
			points->InsertNextPoint(pos[i][0], pos[i][1], pos[i][2]);
			vel_arr->SetTuple3(i, vel[i][0], vel[i][1], vel[i][2]);
			pos_init_arr->SetTuple3(i, pos_init[i][0], pos_init[i][1], pos_init[i][2]);
		}
		
		id_arr->SetTuple1(i, ids[i]);
		colour_arr->SetTuple1(i, colour[i]);
		vtkIdType id[1] = { i };
		vertices->InsertNextCell(1, id);
	}

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->Allocate(num_tracers);
	polydata->SetPoints(points);
	polydata->SetVerts(vertices);
	polydata->GetPointData()->AddArray(vel_arr);
	polydata->GetPointData()->AddArray(id_arr);
	polydata->GetPointData()->AddArray(pos_init_arr);
	polydata->GetPointData()->AddArray(colour_arr);

	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(fname.c_str());
	writer->SetInputData(polydata);
	writer->SetCompressor(compressor);

	if (!writer->Write())
	{
		throw runtime_error("Error writing tracers to vtk file!");
	}
}