/** 
\file structmesh2d.hpp
\brief This file provides a handler class for structured meshes.
\author Aditya Kashi
*/

#ifndef INSAC_STRUCTMESH2D_H
#define INSAC_STRUCTMESH2D_H 1

#include <iostream>
#include <string>

#include "aarray2d.hpp"

namespace acfd {

/**	\brief This class encapsulates the 2D structured mesh.

	It reads from a mesh file and computes coordinates of cell centres, cell volumes and area vectors of each face.
	NOTE that the index of a cell is same as that of its lower-left corner node, as opposed to the convention in class.
*/
class Structmesh2d
{
	amat::Array2d<double> x;
	amat::Array2d<double> y;
	int imx;                                             ///< Number of points along x
	int jmx;                                             ///< Number of points along y

	amat::Array2d<double> xc;                            ///< x-coordinates of cell centres
	amat::Array2d<double> yc;                            ///< y-coordinates of cell centres
	
	/**	\brief Area vectors of faces 1 and 2 for each face.
	 *
	 * del[0](i,j) and del[1](i,j) are x- and y-components of area of face 1; 
	 * del[2](i,j) and del[3](i,j) are x- and y-components of area of face 2.
	 * Note that the vectors point along +i direction for the for the i-faces and 
	 * the +j direction for the j-faces (rather than -i and -j).
	 */
	amat::Array2d<double>* del;

	int ndelcomp;                                    ///< Number of components of del
	bool allocdel;                                   ///< Stores whether or not del has been allocated

	amat::Array2d<double> vol;                       ///< Stores the measure (area in 2D) of each cell

public:
	Structmesh2d();
	Structmesh2d(int nxpoin, int nypoin);

	void setup(int nxpoin, int nypoin);

	~Structmesh2d();

	void readmesh(const std::string fname);

	void writemesh_vtk(std::string fname) const;

	/// Calculates cell centers, volumes and area vectors
	/** Stores in xc, yc, vol  and del respectively. 
	 */
	void preprocess();

	int gimx() const						///< Returns the number of (real) grid points in i-direction.
	{ return imx; }
	int gjmx() const						///< Returns the number of (real) grid points in j-direction.
	{ return jmx; }
	double gx(int i, int j) const			///< Returns x-coordinates of grid points.
	{ return x.get(i,j); }
	double gy(int i, int j) const			///< Returns y-coordinates of grid points.
	{ return y.get(i,j); }
	double gxc(int i, int j) const			///< Returns x-coordinates of cell centers.
	{ return xc.get(i,j); }
	double gyc(int i, int j) const			///< Returns y-coordinates of cell centers.
	{ return yc.get(i,j); }
	double gdel(int i, int j, int idat) const		///< Returns del[idat](i,j); refer to [del](@ref del)
	{ return del[idat].get(i,j); }
	double gvol(int i, int j) const		///< Returns 'volune' (area) of a cell.
	{ return vol.get(i,j); }
};

} // end namespace
#endif
