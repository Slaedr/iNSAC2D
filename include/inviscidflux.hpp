/** \file inviscidflux.hpp 
 * \brief Inviscid flux schemes for the artificial compressiblity (AC) formulation of the iNS equations
*/

#ifndef INSAC_INVISCIDFLUX_H
#define INSAC_INVISCIDFLUX_H

#include <vector>

#include <aarray2d.hpp>
#include <structmesh2d.hpp>

namespace acfd {

/// Abstract class to compute pressure difference dp across each face as (p_L - p_R)
class PressureReconstruction
{
public:
	PressureReconstruction(const Structmesh2d& mesh);
	virtual ~PressureReconstruction();
	virtual void compute_pressure_difference(const amat::Array2d<double> unknowns[NDIM],
	                                         amat::Array2d<double> delta_p[NDIM]) const = 0;
protected:
	const Structmesh2d& m;
};

/** \brief Implements the basic first-order pressure reconstruction needed for the mass flux
 * in the Rhie-Chow method.
 *
 * Computes pressure difference across each face in the mesh
 * (including faces shared between ghost cells).
 */
class BasicPR : public PressureReconstruction
{
public:
	BasicPR(const Structmesh2d& mesh);
	void compute_pressure_difference(const amat::Array2d<double> unknowns[NDIM],
	                                 amat::Array2d<double> delta_p[NDIM]) const;
};

/** \brief Implements a second-order accurate TVD pressure reconstruction for the mass flux needed
 * for the Rhie-Chow scheme for colocated grid.
 */
class TVDPR : public PressureReconstruction
{
public:
	TVDPR(const Structmesh2d& mesh);
	void compute_pressure_difference(const amat::Array2d<double> unknowns[NDIM],
	                                 amat::Array2d<double> delta_p[NDIM]) const;
};

/// Class to add contributions of inviscid fluxes to the residual.
/** Implements a Rhie-Chow stabilization in the mass flux, for which pressure reconstruction is needed.
* @see PressureReconstruction
*/
class InviscidFlux
{
	Structmesh2d* m;
	amat::Array2d<double>* u;			///< Contains one (imx+1) x (jmx+1) array each for pressure, x-velocity and y-velocity.
	amat::Array2d<double>* res;		///< Contains residuals corresponding to [u](@ref u).
	amat::Array2d<double>* beta;		///< Artificial compressiblity factor for each cell.
	double rho;					///< Density of fluid.
	PressureReconstruction* pr;

	/// Pressure difference across each face.
	/** dp[0] contains pressure difference across i-faces. 
	 * dp[1] contains pressure differences across j-faces. 
	 */
	amat::Array2d<double> dp[NDIM];	
	
	double c;					///< Rhie-Chow constant.
	bool isallocdp;

	/** \brief Integers indicating the type of boundary for each of the 4 boundaries.
	 *
	 *	0. velocity inflow
	 *	1. pressure outflow
	 *	2. no-slip wall
	 *	3. slip wall
	 */
	std::vector<int> bcflags;
	
	/// Values of prescribed quantities corresponding to the flags in bcflags.
	/** The value for each boundary consists of two numbers - this is needed for inflow velocity.
	* For other types of boundary, we just read the first value corresponding to that boundary.
	*/
	std::vector<std::vector<double>> bvalues;

	/// one flag for each of the 4 boundaries, indicating whether or not that boundary is a wall.
	std::vector<int> wall;

public:
	void setup(Structmesh2d* mesh, amat::Array2d<double>* unknown, amat::Array2d<double>* residuals, 
			amat::Array2d<double>* _beta, double _rho, std::string pressurereconstruction, 
			std::vector<int> _bcflag, std::vector<std::vector<double>> _bvalues);

	~InviscidFlux();

	/// Add the inviscid flux contribution to the [residual](@ref res)
	void compute_fluxes();
};

} // end namespace
#endif
