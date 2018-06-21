
#ifndef INSAC_INS_H
#define INSAC_INS_H

#include <array>

#include <inviscidflux.hpp>
#include <gradientschemesins.hpp>

namespace acfd {

/** \brief Implements the main solution loop of stead-state iNS equations in artificial compressibility form.
*/
class Steady_insac
{
	Structmesh2d* m;				///< Pointer to mesh.

	/// Values of the 3 unkowns - index 0 is pressure, 1 is x-velocity and 2 is y-velocity.
	amat::Array2d<double> u[NVARS];				
	amat::Array2d<double> res[NVARS];			///< Residuals for each of the 3 quanitities.
	amat::Array2d<double> visc_lhs[NDIM*2+1];		///< Required for gradient computation.
	amat::Array2d<double> beta;			///< Artificial compressibility in each cell.
	amat::Array2d<double> dt;				///< Local time step for each cell.
	double cfl;						///< CFL number

	/// Flags indicating the type of boundary for each of the 4 boundaries.
	std::vector<int> bcflags;			
	/// Boundary values for each of the 4 boundaries (2 each, at most)
	std::vector<std::vector<double>> bvalues;

	double rho;						///< Density
	double mu;						///< Viscosity

	InviscidFlux* invf;				///< Inviscid flux object
	GradientSchemeIns* grad;		///< Gradient reconstruction object for viscous flux

	std::string pressurescheme;			///< Denotes which pressure reconstruction scheme to use in mass flux.
	double uref;					///< Some reference fluid velocity
	int ndim;						///< No. of spatial dimensions (2)
	double tol;						///< Relative residual tolerance to decide convergence to steady-state
	int maxiter;					///< Max number of time steps

	bool isalloc;
	bool isinviscid;				///< If true, viscous flux calculation is not carried out

	amat::Array2d<double>* vel[2];			///< Required for output of velocity

public:
	Steady_insac();

	/// Sets up the iNS problem.
	/** \param _bcflags contains 4 integers denoting the type of boundary for each boundary.
	 - 0 is a velocity inlet
	 - 1 is a pressure outlet
	 - 2 is a no-slip wall
	 - 3 is a slip wall
	\param bvalues contains the corresponding boundary values - at most two per boundary.
	\param gradscheme is string which is either "thinlayer" or "normtan",
	  describing the gradient reconstruction scheme to use for viscous fluxes.
	\param pressure_scheme is a string (either "basic" or "tvd") describing
	  the pressure reconstruction to use for the Rhie-Chow mass flux.
	\param refvel is some reference fluid velocity.
	\param CFL is the C.F.L. number to use.
	\param tolerance is the relative tolerance
	\param maxiters is the maximum number of time steps
	*/
	void setup(Structmesh2d* mesh,
	           double dens, double visc,
	           std::vector<int> _bcflags, std::vector<std::vector<double>> _bvalues,
	           std::string gradscheme, std::string pressure_scheme,
	           double refvel,
	           double CFL, double tolerance, int maxiters, bool is_inviscid = false);

	~Steady_insac();

	/// Computes artificial conmpressibility ([beta](@ref beta) ) for each cell.
	void compute_beta();
	
	/// Sets quantities to ghost cells.
	void setBCs();
	
	void setInitialConditions();
	
	void compute_timesteps();
	
	/// Contains the main solver time loop.
	void solve();
	
	/// Computes solutions at grid nodes by averaging cell values around each node
	void getPointSolution();

	amat::Array2d<double>* getpressure();
	amat::Array2d<double>** getvelocity();
	
	amat::Array2d<double>* getVariables();
	amat::Array2d<double>* getResiduals();
	
	std::vector<double> reslist;			///< For convergence history
};

inline amat::Array2d<double>* Steady_insac::getVariables()
{
	return u;
}

inline amat::Array2d<double>* Steady_insac::getResiduals()
{
	return res;
}

}
#endif
