
#include <cmath>
#include <limits>
#include "ins.hpp"

namespace acfd {

using namespace amat;

Steady_insac::Steady_insac() {
	isalloc = false;
}

void Steady_insac::setup(Structmesh2d* mesh, double dens, double visc,
                         std::vector<int> _bcflags, std::vector<std::vector<double>> _bvalues,
                         std::string gradscheme, std::string pressure_scheme,
                         double refvel,
                         double CFL, double tolerance, int maxiters, bool is_inviscid)
{
	m = mesh;
	rho = dens;
	mu = visc;
	bcflags = _bcflags;
	bvalues = _bvalues;
	pressurescheme = pressure_scheme;
	uref = refvel;
	if(gradscheme == "normtan")
		grad = new NormTanGradientIns;
	else if(gradscheme == "parallelcv")
		grad = new ParallelCVGradientIns;
	else
		grad = new ThinLayerGradientIns;
	
	invf = new InviscidFlux;

	isalloc = true;

	isinviscid = is_inviscid;

	for(int i = 0; i < NVARS; i++)
	{
		u[i].setup(m->gimx()+1, m->gjmx()+1);
		res[i].setup(m->gimx()+1, m->gjmx()+1);
	}
	for(int i = 0; i < NDIM*2+1; i++)
		visc_lhs[i].setup(m->gimx()+1, m->gjmx()+1);
	beta.setup(m->gimx()+1,m->gjmx()+1);
	dt.setup(m->gimx()+1,m->gjmx()+1);
	cfl = CFL;
	tol = tolerance;
	maxiter = maxiters;

	invf->setup(m, u, res, &beta, rho, pressure_scheme, _bcflags, _bvalues);
	grad->setup(m,u,res,visc_lhs, mu);
	
	std::cout << "Steady_insac: setup():\n";
	std::cout << "BC flags ";
	for(int i = 0; i < 4; i++)
		std::cout << bcflags[i] << " ";
	std::cout << "\nB values:\n";
	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 2; j++)
			std::cout << bvalues[i][j] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

Steady_insac::~Steady_insac()
{
	if(isalloc) {
		delete grad;
		delete invf;
	}
}

/** We currently do not consider viscosity in calculating the artificial compressibility [beta](@ref beta).
*/
void Steady_insac::compute_beta()
{
	for(int i = 0; i <= m->gimx(); i++)
		for(int j = 0; j <= m->gjmx(); j++)
		{
			double vmag = sqrt(u[1].get(i,j)*u[1].get(i,j) + u[2].get(i,j)*u[2].get(i,j));
			if(vmag >= uref)
				beta(i,j) = vmag;
			else
				beta(i,j) = uref;
			
			// account for viscous effect
			/*
			double h = sqrt((m->gx(i+1,j+1)-m->gx(i,j+1))*(m->gx(i+1,j+1)-m->gx(i,j+1)) + (m->gy(i+1,j+1)-m->gy(i,j+1))*(m->gy(i+1,j+1)-m->gy(i,j+1)));
			h += sqrt((m->gx(i+1,j+1)-m->gx(i+1,j))*(m->gx(i+1,j+1)-m->gx(i+1,j)) + (m->gy(i+1,j+1)-m->gy(i+1,j))*(m->gy(i+1,j+1)-m->gy(i+1,j)));
			h /= 2.0;
			double uvisc = (mu/rho)/h;
			if(beta(i,j) < uvisc) beta(i,j) = uvisc;*/
		}
}

/** Note that for an inlet boundary, a parabolic profile is imposed with maximum velocity as that given by the [bvalues](@ref bvalues) entries. We assume that the respective boundaries are parallel to x- or y-axis, so the inlets work only for straight boundaries parallel to the axes.
* \note For the time being, inlet is only allowed for boundary 2, ie, for the i=1 boundary.
* \note The corner ghost cells are not directly given values; they are just set as the average of their neighboring ghost cells.
*/
void Steady_insac::setBCs()
{
	// boundary 0
	{
		int i = m->gimx();
		if(bcflags[0] == 1)		// pressure outlet
			for(int j = 1; j <= m->gjmx()-1; j++)
			{
				u[0](i,j) = bvalues[0][0];
				u[1](i,j) = u[1].get(i-1,j);
				u[2](i,j) = u[2].get(i-1,j);
			}
		else if(bcflags[0] == 2)	// no-slip wall
		{
			for(int j = 1; j <= m->gjmx()-1; j++)
			{
				u[0](i,j) = u[0](i-1,j);
				u[1](i,j) = 2.0*bvalues[0][0] - u[1](i-1,j);
				u[2](i,j) = 2.0*bvalues[0][1] - u[2](i-1,j);
			}
		}
		else 			// slip-wall
			for(int j = 1; j <= m->gjmx()-1; j++)
			{
				u[0](i,j) = u[0](i-1,j);
				double area = sqrt(m->gdel(i-1,j,0)*m->gdel(i-1,j,0) + m->gdel(i-1,j,1)*m->gdel(i-1,j,1));
				double nx = m->gdel(i-1,j,0)/area;
				double ny = m->gdel(i-1,j,1)/area;
				double vdotn = u[1](i-1,j)*nx + u[2](i-1,j)*ny;
				u[1](i,j) = u[1](i-1,j) - 2*vdotn*nx;
				u[2](i,j) = u[2](i-1,j) - 2*vdotn*ny;
			}
	}

	// boundary 2
	{
		int i = 0;
		if(bcflags[2] == 0)			// parabolic velocity inlet
		{
			double rm = (m->gy(1,1) + m->gy(1,m->gjmx()))/2.0;		// mid point
			double a = bvalues[2][0] / ( rm*rm - 2.0*rm*rm + m->gy(1,1)*m->gy(1,m->gjmx()) );
			double b = -a*2.0*rm;
			double c = -a*m->gy(1,1)*m->gy(1,1) - b*m->gy(1,1);
			for(int j = 1; j <= m->gjmx()-1; j++) 
			{
				u[1](i,j) = a*m->gyc(i,j)*m->gyc(i,j) + b*m->gyc(i,j) + c;
				u[2](i,j) = 0;
				u[0](i,j) = u[0](i+1,j);
			}
		}
		else if(bcflags[2] == 1)		// pressure outlet
		{
			for(int j = 1; j <= m->gjmx()-1; j++) 
			{
				u[0](i,j) = bvalues[2][0];
				u[1](i,j) = u[1](i+1,j);
				u[2](i,j) = u[2](i+1,j);
			}
		}
		else if(bcflags[2] == 2)	// no-slip wall
			for(int j = 1; j <= m->gjmx()-1; j++)
			{
				u[0](i,j) = u[0](i+1,j);
				u[1](i,j) = 2.0*bvalues[2][0] - u[1](i+1,j);
				u[2](i,j) = 2.0*bvalues[2][1] - u[2](i+1,j);
			}
		else 			// slip-wall
			for(int j = 1; j <= m->gjmx()-1; j++)
			{
				u[0](i,j) = u[0](i+1,j);
				double area = sqrt(m->gdel(i,j,0)*m->gdel(i,j,0) + m->gdel(i,j,1)*m->gdel(i,j,1));
				// we use negative of the del value as the normal always points in the positive i direction.
				double nx = -m->gdel(i,j,0)/area;
				double ny = -m->gdel(i,j,1)/area;
				double vdotn = u[1](i+1,j)*nx + u[2](i+1,j)*ny;
				u[1](i,j) = u[1](i+1,j) - 2*vdotn*nx;
				u[2](i,j) = u[2](i+1,j) - 2*vdotn*ny;
			}
	}

	// boundary 1
	{
		int j = m->gjmx();
		if(bcflags[1] == 1)		// pressure outlet
			for(int i = 1; i <= m->gimx()-1; i++)
			{
				u[0](i,j) = bvalues[1][0];
				u[1](i,j) = u[1](i,j-1);
				u[2](i,j) = u[2](i,j-1);
			}
		else if(bcflags[1] == 2)	// no-slip wall
		{
			for(int i = 1; i <= m->gimx()-1; i++)
			{
				u[0](i,j) = u[0](i,j-1);
				u[1](i,j) = 2.0*bvalues[1][0] - u[1](i,j-1);
				u[2](i,j) = 2.0*bvalues[1][1] - u[2](i,j-1);
			}
		}
		else if(bcflags[1] == 3)	// slip wall
			for(int i = 1; i <= m->gimx()-1; i++)
			{
				u[0](i,j) = u[0](i,j-1);
				double area = sqrt(m->gdel(i,j-1,2)*m->gdel(i,j-1,2) + m->gdel(i,j-1,3)*m->gdel(i,j-1,3));
				double nx = m->gdel(i,j-1,2)/area;
				double ny = m->gdel(i,j-1,3)/area;
				double vdotn = u[1](i,j-1)*nx + u[2](i,j-1)*ny;
				u[1](i,j) = u[1](i,j-1) - 2.0*vdotn*nx;
				u[2](i,j) = u[2](i,j-1) - 2.0*vdotn*ny;
			}
	}

	// boundary 3
	{
		int j = 0;
		if(bcflags[3] == 1)		// pressure outlet
			for(int i = 1; i <= m->gimx()-1; i++)
			{
				u[0](i,j) = bvalues[3][0];
				u[1](i,j) = u[1](i,j+1);
				u[2](i,j) = u[2](i,j+1);
			}
		else if(bcflags[3] == 2)	// no-slip wall
			for(int i = 1; i <= m->gimx()-1; i++)
			{
				u[0](i,j) = u[0](i,j+1);
				u[1](i,j) = 2.0*bvalues[3][0] - u[1](i,j+1);
				u[2](i,j) = 2.0*bvalues[3][1] - u[2](i,j+1);
			}
		else if(bcflags[3] == 3)	// slip wall
			for(int i = 1; i <= m->gimx()-1; i++)
			{
				u[0](i,j) = u[0](i,j+1);
				double area = sqrt(m->gdel(i,j,2)*m->gdel(i,j,2) + m->gdel(i,j,3)*m->gdel(i,j,3));
				// we use negative of the del values as the normal always points in the positive j-direction (this really does not matter, though)
				double nx = m->gdel(i,j,2)/area;
				double ny = m->gdel(i,j,3)/area;
				double vdotn = u[1](i,j+1)*nx + u[2](i,j+1)*ny;
				u[1](i,j) = u[1](i,j+1) - 2.0*vdotn*nx;
				u[2](i,j) = u[2](i,j+1) - 2.0*vdotn*ny;
			}
	}

	// for corner ghost cells
	for(int k = 0; k < NVARS; k++)
	{
		u[k](0,0) = 0.5*(u[k](1,0)+u[k](0,1));
		u[k](0,m->gjmx()) = 0.5*(u[k](1,m->gjmx())+u[k](0,m->gjmx()-1));
		u[k](m->gimx(),0) = 0.5*(u[k](m->gimx(),1)+u[k](m->gimx()-1,0));
		u[k](m->gimx(),m->gjmx()) = 0.5*(u[k](m->gimx()-1,m->gjmx())+u[k](m->gimx(),m->gjmx()-1));
	}
}

/** Computes local time-step for each cell using characteristics of the system in the normal directions.
 * Face normals 'in a cell' are calculated by averaging those of the two faces bounding the cell
 * in each the i- and j-directions.
 */
void Steady_insac::compute_timesteps()
{
	//std::cout << "Steady_insac: compute_timesteps(): Now computing time steps for next iteration..." << endl;

	for(int i = 1; i <= m->gimx()-1; i++)
		for(int j = 1; j <= m->gjmx()-1; j++)
		{
			// area vectors and unit normal vectors
			std::array<a_real,NDIM> areavi, areavj, inormal, jnormal;

			for(int dim = 0; dim < NDIM; dim++) 
			{
				areavi[dim] = 0.5*(m->gdel(i,j,dim) + m->gdel(i-1,j,dim));
				areavj[dim] = 0.5*(m->gdel(i,j,2+dim) + m->gdel(i,j-1,2+dim));
			}
			a_real areai = sqrt(areavi[0]*areavi[0] + areavi[1]*areavi[1]);
			a_real areaj = sqrt(areavj[0]*areavj[0] + areavj[1]*areavj[1]);
			
			for(int dim = 0; dim < NDIM; dim++)
			{
				inormal[dim] = areavi[dim]/areai;
				jnormal[dim] = areavj[dim]/areaj;
			}
			// we now have face geometrical info "at the cell centers"

			a_real vdotni = u[1](i,j)*inormal[0] * u[2](i,j)*inormal[1];
			a_real vdotnj = u[1](i,j)*jnormal[0] * u[2](i,j)*jnormal[1];

			a_real eigeni = 0.5*( fabs(vdotni) + sqrt(vdotni*vdotni + 4.0*beta(i,j)*beta(i,j)) );
			a_real eigenj = 0.5*( fabs(vdotnj) + sqrt(vdotnj*vdotnj + 4.0*beta(i,j)*beta(i,j)) );
			
			a_real voldt = eigeni*areai + eigenj*areaj;
			dt(i,j) = m->gvol(i,j)/voldt * cfl;
		}
}

/** Sets all quantities in all cells to zero. */
void Steady_insac::setInitialConditions()
{
//#pragma omp parallel for default(shared)
	for(int i = 0; i <= m->gimx(); i++)
		for(int j = 0; j <= m->gjmx(); j++)
		{
			for(int k = 0; k < NVARS; k++)
			{
				u[k](i,j) = 0;
				res[k](i,j) = 0;
			}
		}
}

/** Make sure the [Steady_insac](@ref Steady_insac) object has been [setup](@ref setup)
 * and initialized with some [initial condition](@ref setInitialConditions).
 * Both the momentum-magnitude residual and the mass flux are taken as convergence criteria;
 * the tolerance for the mass flux is the square-root of the tolerance for
 * the relative momentum-magnitude residual. [tol](@ref tol) is the latter.
 */
void Steady_insac::solve()
{
	setInitialConditions();
	setBCs();
	compute_beta();

	double resnorm, resnorm0=1.0, massflux;
	std::cout << "Steady_insac: solve(): Beginning the time-march." << std::endl;
	for(int n = 0; n < maxiter; n++)
	{
		// calculate stuff needed for this iteration
		setBCs();
		
		compute_beta();
		
		compute_timesteps();

		for(int k = 0; k < NVARS; k++)
			res[k].zeros();
		
		// compute fluxes
		invf->compute_fluxes();
		if(!isinviscid)
			grad->compute_fluxes();

		// check convergence
		resnorm = 0; massflux = 0;
//#pragma omp parallel for default(shared) reduction(+:resnorm) reduction(+:massflux)
		for(int i = 1; i <= m->gimx()-1; i++)
			for(int j = 1; j <= m->gjmx()-1; j++)
			{
				resnorm += (res[1].get(i,j)*res[1].get(i,j) + res[2].get(i,j)*res[2].get(i,j))*m->gvol(i,j);
				massflux += res[0].get(i,j);
			}
		
		resnorm = std::sqrt(resnorm);

		// It is important to check both the quantities in the if below!
		if(!std::isfinite(resnorm) || !std::isfinite(resnorm/resnorm0)) {
			throw "Steady_insac: Pseudo time stepper diverged to inf or NaN!";
		}

		reslist.push_back(resnorm);
		if(n == 0) resnorm0 = resnorm;
		if(n == 1 || n%20 == 0) {
			std::cout << "Steady_insac: solve(): Iteration " << n << ": relative L2 norm of residual = "
			          << resnorm/resnorm0 << std::endl;
			//std::cout << "  L2 norm of residual = " << resnorm << endl;
		}

		if(resnorm/resnorm0 < tol /*&& fabs(massflux) < sqrt(tol)*/)
		{
			std::cout << "Steady_insac: solve(): Converged in " << n
			          << " iterations. Norm of final residual = " << resnorm
			          << ", final net mass flux = " << massflux << std::endl;

			if(massflux > std::numeric_limits<a_real>::epsilon())
				throw "Mass flux is greater than machine epsilon!";

			break;
		}

		// update u
		for(int i = 1; i <= m->gimx()-1; i++)
			for(int j = 1; j <= m->gjmx()-1; j++)
			{
				double dtv = dt.get(i,j)/m->gvol(i,j);
				u[0](i,j) = u[0].get(i,j) - res[0].get(i,j)*beta.get(i,j)*beta.get(i,j)*dtv;
				for(int k = 1; k < NVARS; k++)
					u[k](i,j) = u[k].get(i,j) - res[k].get(i,j)*dtv/rho;
			}
	}
}

}
