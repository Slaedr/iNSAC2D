/** @file aoutput_struct.hpp
 * @brief Provides a class to manage output of mesh data to VTK-type files.
 * @author Aditya Kashi
*/

#include "aoutput_struct.hpp"

namespace acfd {

Structdata2d::Structdata2d(Structmesh2d* mesh, amat::Array2d<double>* _u, amat::Array2d<double>* _res, 
		const std::string _title)
{
	m = mesh;
	u = _u;
	res = _res;
	title = _title;
	upoint = new amat::Array2d<double>[3];
	for(int i = 0; i < 3; i++)
		upoint[i].setup(m->gimx()+1,m->gjmx()+1);
}

Structdata2d::~Structdata2d()
{
	delete [] upoint;
}

void Structdata2d::writevtk(const std::string fname)
{
	std::cout << "Structdata2d: writevtk(): Writing data to file " << fname << std::endl;
	std::ofstream fout(fname);
	fout << "# vtk DataFile Version 2.0\n";
	fout << title << '\n';
	std::cout << title << '\n';
	fout << "ASCII\n";
	fout << "DATASET STRUCTURED_GRID\n";
	fout << "DIMENSIONS " << m->gimx() << " " << m->gjmx() << " 1\n";
	fout << "POINTS " << m->gimx()*m->gjmx() << " float\n";
	for(int j = 1; j <= m->gjmx(); j++)
		for(int i = 1; i <= m->gimx(); i++)
			fout << m->gx(i,j) << " " << m->gy(i,j) << " " << 0.0 << '\n';
	
	// calculate point data
	for(int k = 0; k < 3; k++)
		for(int i = 1; i <= m->gimx(); i++)
			for(int j = 1; j <= m->gjmx(); j++)
			{
				upoint[k](i,j) = (u[k].get(i,j) + u[k].get(i-1,j) 
						+ u[k].get(i,j-1) + u[k].get(i-1,j-1))*0.25;
			}

	// Now output data
	fout << "POINT_DATA " << m->gimx()*m->gjmx() << '\n';

	fout << "SCALARS " << "pressure" << " float 1\n";
	fout << "LOOKUP_TABLE default\n";
	for(int j = 1; j <= m->gjmx(); j++)
		for(int i = 1; i <= m->gimx(); i++)
			fout << upoint[0].get(i,j) << '\n';
	

	fout << "VECTORS " << "velocity" << " float\n";
	for(int j = 1; j <= m->gjmx(); j++)
		for(int i = 1; i <= m->gimx(); i++)
		{
			for(int idim = 1; idim < NDIM+1; idim++)
				fout << upoint[idim].get(i,j) << " ";
			for(int idim = 3-NDIM; idim > 0; idim--)
				fout << "0 ";
			fout << '\n';
		}
	
	fout << "CELL_DATA " << (m->gimx()-1)*(m->gjmx()-1) << '\n';

	fout << "SCALARS " << "pressure_cell" << " float 1\n";
	fout << "LOOKUP_TABLE default\n";
	for(int j = 1; j <= m->gjmx()-1; j++)
		for(int i = 1; i <= m->gimx()-1; i++)
			fout << u[0].get(i,j) << '\n';
	
	fout << "SCALARS " << "mass_res_cell" << " float 1\n";
	fout << "LOOKUP_TABLE default\n";
	for(int j = 1; j <= m->gjmx()-1; j++)
		for(int i = 1; i <= m->gimx()-1; i++)
			fout << res[0].get(i,j) << '\n';
	

	fout << "VECTORS " << "velocity" << "_cell float\n";
	for(int j = 1; j <= m->gjmx()-1; j++)
		for(int i = 1; i <= m->gimx()-1; i++)
		{
			for(int idim = 1; idim < NDIM+1; idim++)
				fout << u[idim].get(i,j) << " ";
			for(int idim = 3-NDIM; idim > 0; idim--)
				fout << "0 ";
			fout << '\n';
		}
	
	fout << "VECTORS " << "momentum_res" << "_cell float\n";
	for(int j = 1; j <= m->gjmx()-1; j++)
		for(int i = 1; i <= m->gimx()-1; i++)
		{
			for(int idim = 1; idim < NDIM+1; idim++)
				fout << res[idim].get(i,j) << " ";
			for(int idim = 3-NDIM; idim > 0; idim--)
				fout << "0 ";
			fout << '\n';
		}
}

}
