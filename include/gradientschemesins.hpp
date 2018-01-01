/** \file gradientschemesins.hpp 
 * \brief Gradient computation schemes for viscous terms
 */

#ifndef INSAC_GRADIENTSCHEMESINS_H
#define INSAC_GRADIENTSCHEMESINS_H 1

#include "aarray2d.hpp"
#include "structmesh2d.hpp"

#include <vector>

namespace acfd {

/// \brief Base class for gradient reconstruction schemes for AC formulation of iNS equations.
/**	\note compute_fluxes() *increments* the residual res. If res contains anything non-zero, 
 * contribution by gradients is *added* to it; it is not replaced.
 */
class GradientSchemeIns
{
protected:
	int nvar;				///< Number of variables in u.
	int ndim;
	double mu;				///< Viscosity (dynamic)
	Structmesh2d* m;		///< Associated mesh
	amat::Array2d<double>* u;		///< The unknown with which to compute the graient flux
	amat::Array2d<double>* res;	///< The residual containing fluxes at each grid cell
	
	/** \brief LHS of AF scheme corresponding to a particular (compact) gradient model.
	
		- a(i,j,0) corresponds to coeff of u(i-1,j)
		- a(i,j,1) corresponds to coeff of u(i,j-1)
		- a(i,j,2) corresponds to coeff of u(i,j)
		- a(i,j,3) corresponds to coeff of u(i,j+1)
		- a(i,j,4) corresponds to coeff of u(i+1,j),
	 in the row of Ax=b corresponding to the (i,j) cell.
	*/
	amat::Array2d<double>* a;

	/// \brief Volumes of the thin-layer CVs for the gradient.
	///
	/// Two components: first component for +i face and the other for +j face. Arranged like Structmesh2d::del.
	std::vector<amat::Array2d<double>> dualvol;

public:
	virtual ~GradientSchemeIns();
	virtual void setup(Structmesh2d* mesh, amat::Array2d<double>* unknown, amat::Array2d<double>* residual, amat::Array2d<double>* lhs, double visc);
	void compute_CV_volumes();		///< Computes [volumes](@ref dualvol) of thin-layer CVs. To be precomputed just once.

	virtual void compute_s();		///< Required for Normal tangent gradient decomposition scheme.
	
	virtual void compute_fluxes() = 0;
	
	/** \brief Computes LHS arrays corresponding to a particular gradient scheme.
	
		Make sure to execcute [calculate_CV_volumes](@ref calculate_CV_volumes) before calling this function.
	*/
	virtual void compute_lhs() = 0;
};

/// Parallel CV model for gradient reconstruction
class ParallelCVGradientIns : public GradientSchemeIns
{
	/// flux across any i-face
	std::vector<double> g;
	/// flux across any j-face
	std::vector<double> h;

	/// CV normal
	std::vector<double> deln;
	/// gradient
	std::vector<double> grad;
	
public:
	void setup(Structmesh2d* mesh, amat::Array2d<double>* unknown, amat::Array2d<double>* residual, amat::Array2d<double>* lhs, double visc);
	
	/// This function is a dummy for this class
	void compute_lhs();

	void compute_fluxes();
};

/**	\brief Thin layer gradient reconstruction scheme. 
*
* We consider grad u at a face to be influenced by only the change in normal component of u at the face.
*/
class ThinLayerGradientIns : public GradientSchemeIns
{
public:
	/// Initializes data and computes LHS
	void setup(Structmesh2d* mesh, amat::Array2d<double>* unknown, amat::Array2d<double>* residual, amat::Array2d<double>* lhs, double visc);
	/// Adds viscous flux contribution to the residual
	void compute_fluxes();
	void compute_lhs();
};

/**	\brief Computes LHS and residual for FVM solution of Poisson equation using 
 * Normal-tangent gradient decomposition model.
 *
 *	Applied properly, the scheme is non-compact and has 13 terms in the LHS for each cell. 
 *	However, we implement only the 5 compact terms for the LHS. But we treat the residual fully. 
 *	This requires specification of more than one layer of ghost cells. 
 *	This is taken care of without introducing another layer;
 *	rather, we compute a `ghost state' on-the-fly using the BCs.
 */
class NormTanGradientIns : public GradientSchemeIns
{
	/** \brief Stores a unit vector in the direction of the line joining the two cells across each face.
	 *
	 *  Organized like Structmesh2d::del.
	 */
	std::vector<amat::Array2d<double>> svect;

	/** \brief Magnitude of the corresponding vectors in [svect](@ref svect).
	 *
	 * dels[0] refers to magnitude of the vector for +i face and 
	 * dels[1] refers to magnitude of the vector for +j face.
	 */
	std::vector<amat::Array2d<double>> dels;

public:
	void compute_s();
	void compute_lhs();
	void compute_fluxes();
};

}
#endif
