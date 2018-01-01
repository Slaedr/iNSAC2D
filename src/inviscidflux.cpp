/** \file inviscidflux.cpp
 * \brief Implements inviscid flux schemes for iNS equations
 */

#include <inviscidflux.hpp>

namespace acfd {

void PressureReconstruction::setup(Structmesh2d* mesh, amat::Array2d<double>* unknowns, amat::Array2d<double>* delp)
{
	std::cout << "PressureReconstruction: Contructor called." << std::endl;
	m = mesh;
	u = unknowns;
	dp = delp;
}

PressureReconstruction::~PressureReconstruction() { }

void BasicPR::compute_pressure_difference()
{
	// we include the i=0 and j=0 ghost cells
	for(int j = 0; j <= m->gjmx()-1; j++)
#pragma omp simd
		for(int i = 0; i <= m->gimx()-1; i++)
		{
			dp[0](i,j) = u[0].get(i,j) - u[0].get(i+1,j);		// p_L - p_R
			dp[1](i,j) = u[0].get(i,j) - u[0].get(i,j+1);		// p_L - p_R in j-direction, so it's more like p_down - p_up
		}
}

void TVDPR::compute_pressure_difference()
{
	int i,j;
	for(j = 1; j <= m->gjmx()-2; j++)
		for(i = 1; i <= m->gimx()-2; i++)
		{
			dp[0](i,j) = u[0].get(i,j) +0.5*minmod_avg(u[0].get(i+1,j)-u[0](i,j), u[0].get(i,j)-u[0].get(i-1,j));				// p_L
			dp[0](i,j) -= u[0].get(i+1,j) - 0.5*minmod_avg(u[0].get(i+1,j)-u[0].get(i,j), u[0].get(i+2,j)-u[0].get(i+1,j));		// p_R
			dp[1](i,j) = u[0].get(i,j) +0.5*minmod_avg(u[0].get(i,j+1)-u[0](i,j), u[0].get(i,j)-u[0].get(i,j-1));				// p_down
			dp[1](i,j) -= u[0].get(i,j+1) - 0.5*minmod_avg(u[0].get(i,j+1)-u[0].get(i,j), u[0].get(i,j+2)-u[0].get(i,j+1));		// p_up
		}

	// Now values for remaining cells (real or ghost) - just copy values from adjacent interior cells.
	// We do not need values for cells i = imx and j = jmx (ghost cells). We also probably don't need dp for corner cells.
	for(j = 1; j <= m->gjmx()-2; j++)
	{
		dp[0](0,j) = dp[0](1,j);
		dp[0](m->gimx()-1,j) = dp[0](m->gimx()-2,j);
		dp[1](0,j) = u[0].get(0,j) +0.5*minmod_avg(u[0].get(0,j+1)-u[0](0,j), u[0].get(0,j)-u[0].get(0,j-1));				// p_down
		dp[1](0,j) -= u[0].get(0,j+1) - 0.5*minmod_avg(u[0].get(0,j+1)-u[0].get(0,j), u[0].get(0,j+2)-u[0].get(0,j+1));		// p_up
		dp[1](m->gimx()-1,j) = u[0].get(m->gimx()-1,j) +0.5*minmod_avg(u[0].get(m->gimx()-1,j+1)-u[0](m->gimx()-1,j), u[0].get(m->gimx()-1,j)-u[0].get(m->gimx()-1,j-1));
		dp[1](m->gimx()-1,j) -= u[0].get(m->gimx()-1,j+1) -0.5*minmod_avg(u[0].get(m->gimx()-1,j+1)-u[0].get(m->gimx()-1,j), u[0].get(m->gimx()-1,j+2)-u[0].get(m->gimx()-1,j+1));
	}
	for(i = 1; i <= m->gimx()-2; i++)
	{
		dp[1](i,0) = dp[1](i,1);
		dp[1](i,m->gjmx()-1) = dp[1](i,m->gjmx()-2);
		dp[0](i,0) = u[0].get(i,0) +0.5*minmod_avg(u[0].get(i+1,0)-u[0](i,0), u[0].get(i,0)-u[0].get(i-1,0));				// p_L
		dp[0](i,0) -= u[0].get(i+1,0) - 0.5*minmod_avg(u[0].get(i+1,0)-u[0].get(i,0), u[0].get(i+2,0)-u[0].get(i+1,0));		// p_R
		dp[0](i,m->gjmx()-1) = u[0].get(i,j) +0.5*minmod_avg(u[0].get(i+1,j)-u[0](i,j), u[0].get(i,j)-u[0].get(i-1,j));
		dp[0](i,m->gjmx()-1) -= u[0].get(i+1,m->gjmx()-1) - 0.5*minmod_avg(u[0].get(i+1,m->gjmx()-1)-u[0].get(i,m->gjmx()-1), u[0].get(i+2,m->gjmx()-1)-u[0].get(i+1,m->gjmx()-1));
	}
}

void InviscidFlux::setup(Structmesh2d* mesh, amat::Array2d<double>* unknown, amat::Array2d<double>* residuals, amat::Array2d<double>* _beta, double _rho, std::string pressurereconstruction, std::vector<int> _bcflag, std::vector<std::vector<double>> _bvalues)
{
	std::cout << "InviscidFlux: setup(): Setting up inviscid flux calculator" << std::endl;
	m = mesh;
	u = unknown;
	res = residuals;
	beta = _beta;
	rho = _rho;
	bvalues = _bvalues;
	bcflags = _bcflag;
	dp = new amat::Array2d<double>[2];
	for(int i = 0; i<2; i++)
		dp[i].setup(m->gimx()+1, m->gjmx()+1);
	isallocdp = true;
	
	/// Sets the pressure reconstruction scheme to be used based on the last argument - "basic" or "TVD".
	if(pressurereconstruction=="basic")
	{
		pr = new BasicPR;
		std::cout << "InviscidFlux: setup(): BasicPR selected.\n";
	}
	else if(pressurereconstruction=="tvd")
		pr = new TVDPR;
	else {
		std::cout << "InviscidFlux: setup(): Pressure reconstruction scheme requested does not exist. Choosing basic first-order scheme." << std::endl;
		pr = new BasicPR;
	}
	pr->setup(m, unknown, dp);
	
	/// Rhie-Chow constant [c](@ref c) is maximum 0.5.
	c = 0.5;
	
	nvar = 3;
	
	/*std::cout << "BC flags ";
	for(int i = 0; i < 4; i++)
		std::cout << bcflags[i] << " ";
	std::cout << "\nB values:\n";
	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 2; j++)
			std::cout << bvalues[i][j] << " ";
		std::cout << std::endl;
	}*/
	
	// account for solid walls; ie bcflag values of 2 or 3
	wall.resize(4);		// one flag for each of the 4 boundaries.
	for(int i = 0; i < 4; i++)
		if(bcflags[i] == 2 || bcflags[i] == 3)
			wall[i] = 0.0;
		else
			wall[i] = 1.0;
	
	/*std::cout << "Wall:";
	for(int i = 0; i < 4; i++)
		std::cout << " " << wall[i];
	std::cout << std::endl;*/
}

InviscidFlux::~InviscidFlux()
{
	if(isallocdp) {
		delete [] dp;
		delete pr;
	}
}

void InviscidFlux::compute_fluxes()
{
	pr->compute_pressure_difference();

	// add inviscid flux contribution to residuals
	//std::cout << "InviscidFlux: compute_flux(): Computing inviscid fluxes now..." << std::endl;
	int i, j, k;
	double area, nx, ny, eigen, vdotn;
	double bhalf2; std::vector<double> uhalf(nvar);			// interface values for each unknown
	amat::Array2d<double> g(m->gimx(),nvar);
	
	for(j = 1; j <= m->gjmx()-1; j++)
	{
		for(i = 0; i <= m->gimx()-1; i++)
		{
			area = sqrt(m->gdel(i,j,0)*m->gdel(i,j,0) + m->gdel(i,j,1)*m->gdel(i,j,1));
			nx = m->gdel(i,j,0)/area;
			ny = m->gdel(i,j,1)/area;
			for(k = 0; k < nvar; k++)
				uhalf[k] = 0.5*(u[k](i,j) + u[k](i+1,j));

			// get average of beta^2
			bhalf2 = 0.5*(beta->get(i,j)*beta->get(i,j) + beta->get(i+1,j)*beta->get(i+1,j));

			vdotn = uhalf[1]*nx + uhalf[2]*ny;
			
			// now get eigenvalue for Rhie-Chow
			eigen = 0.5*( fabs(vdotn) + sqrt(vdotn*vdotn + 4.0*bhalf2) );

			// Boundary faces need to be treated differently to account for wall BCs.
			if(i == 0)
			{
				g(i,0) = area*(rho*vdotn + 0.5*c*eigen*dp[0](i,j)/bhalf2)*wall[2];
				g(i,1) = area*(rho*vdotn*uhalf[1]*wall[2] + uhalf[0]*nx);
				g(i,2) = area*(rho*vdotn*uhalf[2]*wall[2] + uhalf[0]*ny);
			}
			else if(i == m->gimx()-1)
			{
				g(i,0) = area*(rho*vdotn + 0.5*c*eigen*dp[0](i,j)/bhalf2)*wall[0];
				g(i,1) = area*(rho*vdotn*uhalf[1]*wall[0] + uhalf[0]*nx);
				g(i,2) = area*(rho*vdotn*uhalf[2]*wall[0] + uhalf[0]*ny);
			}
			else
			{
				g(i,0) = area*(rho*vdotn + 0.5*c*eigen*dp[0](i,j)/bhalf2);
				g(i,1) = area*(rho*vdotn*uhalf[1] + uhalf[0]*nx);
				g(i,2) = area*(rho*vdotn*uhalf[2] + uhalf[0]*ny);
			}
		}
		for(i = 1; i <= m->gimx()-1; i++)
			for(k = 0; k < nvar; k++)
				res[k](i,j) += g(i,k) - g(i-1,k);
	}

	// now we add contribution of j-fluxes
	amat::Array2d<double> h(m->gjmx(),nvar);

	for(i = 1; i <= m->gimx()-1; i++)
	{
		for(j = 0; j <= m->gjmx()-1; j++)
		{
			area = sqrt(m->gdel(i,j,2)*m->gdel(i,j,2) + m->gdel(i,j,3)*m->gdel(i,j,3));
			nx = m->gdel(i,j,2)/area;
			ny = m->gdel(i,j,3)/area;
			for(k = 0; k < nvar; k++)
				uhalf[k] = 0.5*(u[k](i,j) + u[k](i,j+1));
			bhalf2 = 0.5*(beta->get(i,j)*beta->get(i,j) + beta->get(i,j+1)*beta->get(i,j+1));
			vdotn = uhalf[1]*nx + uhalf[2]*ny;
			eigen = 0.5*( fabs(vdotn) + sqrt(vdotn*vdotn + 4.0*bhalf2) );
			
			if(j == 0)
			{
				h(j,0) = area*(rho*vdotn + 0.5*c*eigen*dp[1](i,j)/bhalf2)*wall[3];
				h(j,1) = area*(rho*vdotn*uhalf[1]*wall[3] + uhalf[0]*nx);
				h(j,2) = area*(rho*vdotn*uhalf[2]*wall[3] + uhalf[0]*ny); 
			}
			else if(j == m->gjmx()-1)
			{
				h(j,0) = area*(rho*vdotn + 0.5*c*eigen*dp[1](i,j)/bhalf2)*wall[1];
				h(j,1) = area*(rho*vdotn*uhalf[1]*wall[1] + uhalf[0]*nx);
				h(j,2) = area*(rho*vdotn*uhalf[2]*wall[1] + uhalf[0]*ny); 
			}
			else
			{
				h(j,0) = area*(rho*vdotn + 0.5*c*eigen*dp[1](i,j)/bhalf2);
				h(j,1) = area*(rho*vdotn*uhalf[1] + uhalf[0]*nx);
				h(j,2) = area*(rho*vdotn*uhalf[2] + uhalf[0]*ny); 
			}
		}
		for(j = 1; j <= m->gjmx()-1; j++)
			for(k = 0; k < nvar; k++)
				res[k](i,j) += h(j,k) - h(j-1,k);
	}
}

}
