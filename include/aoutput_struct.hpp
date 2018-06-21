/** @file aoutput_struct.hpp
 * @brief Provides a class to manage output of mesh data to VTK-type files.
 * @author Aditya Kashi
*/

#ifndef INSAC_OUTPUT_STRUCT_H
#define INSAC_OUTPUT_STRUCT_H 1

#include <iostream>
#include <fstream>
#include <string>

#include "aarray2d.hpp"
#include "structmesh2d.hpp"

namespace acfd {

///	Class for managing output of analysis data for simulations on structured 2D grids.
/**	We assume 1-based arrays for all array-quantities.
 * Cell data is taken as input and both cell data and point data are written to the output file.
 */
class Structdata2d
{
	Structmesh2d* m;
	amat::Array2d<double>* u;
	amat::Array2d<double>* res;
	amat::Array2d<double>* upoint;
	std::string title;

public:
	Structdata2d(Structmesh2d* mesh, amat::Array2d<double>* _u, amat::Array2d<double>* _res, 
			const std::string _title);
	
	~Structdata2d();

	/// Writes data to a file in legacy VTK format
	void writevtk(const std::string fname);
};

}
#endif
