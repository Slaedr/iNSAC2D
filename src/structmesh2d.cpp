
#include <structmesh2d.hpp>
#include <exception>

namespace acfd {


Structmesh2d::Structmesh2d()
{
	allocdel = false;
	ndelcomp = 4;
}

Structmesh2d::Structmesh2d(int nxpoin, int nypoin)
{
	imx = nxpoin;
	jmx = nypoin;

	// Allocate space for point coordinates. 
	// Two extra in each direction for ghost cells.
	x.setup(imx+2,jmx+2);
	y.setup(imx+2,jmx+2);

	// Allocate space for cell centers. 
	// Space required = (imx-1) real cells + 2 ghost cells; similarly for j direction.
	xc.setup(imx+1,jmx+1);
	yc.setup(imx+1,jmx+1);

	vol.setup(imx+1,jmx+1);

	ndelcomp = 4;
	del = new amat::Array2d<double>[ndelcomp];
	for(int i = 0; i < ndelcomp; i++)
		del[i].setup(imx+1,jmx+1);
	allocdel = true;
}

Structmesh2d::~Structmesh2d()
{
	if(allocdel)
	{
		delete [] del;
		allocdel = false;
	}
}

void Structmesh2d::setup(int nxpoin, int nypoin)
{
	imx = nxpoin;
	jmx = nypoin;

	// Allocate space for point coordinates. Two extra in each direction for ghost cells.
	x.setup(imx+2,jmx+2);
	y.setup(imx+2,jmx+2);

	// Allocate space for cell centers. 
	// Space required = (imx-1) real cells + 2 ghost cells; similarly for j direction.
	xc.setup(imx+1,jmx+1);
	yc.setup(imx+1,jmx+1);

	vol.setup(imx+1,jmx+1);

	ndelcomp = 4;
	del = new amat::Array2d<double>[ndelcomp];
	for(int i = 0; i < ndelcomp; i++)
		del[i].setup(imx+1,jmx+1);
	allocdel = true;
}

void Structmesh2d::readmesh(const std::string meshname)
{
	std::ifstream fin;
	try {
		fin.open(meshname);
	}
	catch(std::exception& e) {
		throw "Could not open file!";
	}
	
	fin >> imx >> jmx;
	
	std::cout << " Structmesh2d: readmesh(): Reading file " << meshname 
		<< ", size = " << imx << " x " << jmx << std::endl;

	/// Allocate space for point coordinates. Two extra in each direction for ghost cells.
	x.setup(imx+2,jmx+2);
	y.setup(imx+2,jmx+2);

	// Allocate space for cell centers. 
	// Space required = (imx-1) real cells + 2 ghost cells; similarly for j direction.
	xc.setup(imx+1,jmx+1);
	yc.setup(imx+1,jmx+1);
	
	if(allocdel)
		delete [] del;

	del = new amat::Array2d<double>[ndelcomp];
	for(int i = 0; i < ndelcomp; i++)
		del[i].setup(imx+1,jmx+1);
	allocdel = true;
	vol.setup(imx+1,jmx+1);

	for(int j = 1; j <= jmx; j++)
		for(int i = 1; i <= imx; i++)
			fin >> x(i,j) >> y(i,j);
	
	std::cout << "Structmesh2d: readmesh(): Mesh read. Points in i-dir: " << imx 
		<< ", points in j-dir: " << jmx << "." << std::endl;
}

void Structmesh2d::preprocess()
{
	// Add ghost points. 
	// The ghost point is such that the boundary point is the 
	// arithmetic mean of the ghost point and the first interior point.
	for(int j = 1; j <= jmx; j++)
	{
		x(0,j) = 2*x(1,j) - x(2,j);
		y(0,j) = 2*y(1,j) - y(2,j);
		x(imx+1,j) = 2*x(imx,j) - x(imx-1,j);
		y(imx+1,j) = 2*y(imx,j) - y(imx-1,j);
	}
	for(int i = 1; i <= imx; i++)
	{
		x(i,0) = 2*x(i,1) - x(i,2);
		y(i,0) = 2*y(i,1) - y(i,2);
		x(i,jmx+1) = 2*x(i,jmx) - x(i,jmx-1);
		y(i,jmx+1) = 2*y(i,jmx) - y(i,jmx-1);
	}
	// corner points:
	x(0,0) = 2*x(1,0) - x(2,0);
	y(0,0) = 2*y(1,0) - y(2,0);
	x(0,jmx+1) = 2*x(1,jmx+1) - x(2,jmx+1);
	y(0,jmx+1) = 2*y(1,jmx+1) - y(2,jmx+1);
	x(imx+1,0) = 2*x(imx,0) - x(imx-1,0);
	y(imx+1,0) = 2*y(imx,0) - y(imx-1,0);
	x(imx+1,jmx+1) = 2*x(imx,jmx+1) - x(imx-1,jmx+1);
	y(imx+1,jmx+1) = 2*y(imx,jmx+1) - y(imx-1,jmx+1);

	// Calculate cell centers, volumes and face areas. We calculate volumes by Heron's formula.
	std::cout << "Structmesh2d: preprocess(): Calculating cell centers, volume and del" << std::endl;
	double a,b,c,s;
	for(int i = 0; i <= imx; i++)
		for(int j = 0; j <= jmx; j++)
		{
			if((i>0 || j>0) && (i < imx || j < imx))	// ignore corner ghost cells 0,0 and imx,jmx
			{
				// set cell centers
				xc(i,j) = 0.25*(x(i,j) + x(i+1,j) + x(i,j+1) + x(i+1,j+1));
				yc(i,j) = 0.25*(y(i,j) + y(i+1,j) + y(i,j+1) + y(i+1,j+1));
				// next, calculate volumes
				a = sqrt( (x(i,j)-x(i,j+1))*(x(i,j)-x(i,j+1)) + (y(i,j)-y(i,j+1))*(y(i,j)-y(i,j+1)) );
				b = sqrt( (x(i,j+1)-x(i+1,j+1))*(x(i,j+1)-x(i+1,j+1)) + (y(i,j+1)-y(i+1,j+1))*(y(i,j+1)-y(i+1,j+1)) );
				c = sqrt( (x(i,j)-x(i+1,j+1))*(x(i,j)-x(i+1,j+1)) + (y(i,j)-y(i+1,j+1))*(y(i,j)-y(i+1,j+1)) );
				s = (a+b+c)*0.5;
				vol(i,j) = sqrt(s*(s-a)*(s-b)*(s-c));
				
				a = sqrt( (x(i+1,j)-x(i+1,j+1))*(x(i+1,j)-x(i+1,j+1)) + (y(i+1,j)-y(i+1,j+1))*(y(i+1,j)-y(i+1,j+1)) );
				b = sqrt( (x(i,j)-x(i+1,j))*(x(i,j)-x(i+1,j)) + (y(i,j)-y(i+1,j))*(y(i,j)-y(i+1,j)) );
				s = (a+b+c)*0.5;
				vol(i,j) += sqrt(s*(s-a)*(s-b)*(s-c));
			}

			// Compute face area vectors
			// Note that face normals point in the positive i or positive j directions.
			del[0](i,j) = y(i+1,j+1) - y(i+1,j);
			del[1](i,j) = -(x(i+1,j+1) - x(i+1,j));
			del[2](i,j) = -(y(i+1,j+1) - y(i,j+1));
			del[3](i,j) = x(i+1,j+1) - x(i,j+1);
		}
}


}
