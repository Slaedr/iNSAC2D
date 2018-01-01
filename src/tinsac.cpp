/** @file tinsac.cpp
 * @brief Driver function for time-dependent iNS AC solution
 */

#include <instime.hpp>
#ifndef __AOUTPUT_STRUCT_H
#include <aoutput_struct.hpp>
#endif

using namespace amat;
using namespace acfd;
using namespace std;

int main(int argc, char* argv[])
{
	if(argc < 2) {
		cout << "Please provide the name of a control file. Quitting.\n";
		return -1;
	}

	string meshfile, outfile, pressurescheme, gradscheme, dum;
	double tol, rho, mu, cfl, refvel, totaltime;
	int maxiter;
	vector<int> bcflags(4);
	vector<vector<double>> bvalues(4);
	for(int i = 0; i < 4; i++)
		bvalues[i].resize(2,-1);

	ifstream conf(argv[1]);

	// read control file
	conf >> dum; conf >> meshfile;
	conf >> dum; conf >> outfile;
	conf >> dum; conf >> rho;
	conf >> dum; conf >> mu;
	conf >> dum; conf >> refvel;
	conf >> dum; conf >> gradscheme;
	conf >> dum; conf >> pressurescheme;
	conf >> dum; conf >> cfl;
	conf >> dum; conf >> tol;
	conf >> dum; conf >> maxiter;
	conf >> dum; conf >> totaltime;
	conf >> dum;
	for(int i = 0; i < 4; i++)
		conf >> bcflags[i];
	conf >> dum;
	for(int i = 0; i < 4; i++)
	{
		// for a np-slip wall, we can take x- and y-velocities
		// for a velocity inlet, we just take the max velocity of a parabolic profile

		if(bcflags[i] == 2)
			conf >> bvalues[i][0] >> bvalues[i][1];
		else
			conf >> bvalues[i][0];
	}

	conf.close();

	cout << "Input data:\n";
	cout << "Pressure scheme = " << pressurescheme << endl;
	cout << "Gradient scheme = " << gradscheme << endl << endl;
	cout << "CFL = " << cfl << endl;
	cout << "Dynamic viscosity = " << mu << endl;
	cout << "Density = " << rho << endl;
	cout << "Reference velocity = " << refvel << endl;
	cout << "Final time = " << totaltime << endl;

	Structmesh2d m;
	m.readmesh(meshfile);
	m.preprocess();
	
	cout << "Length = " << m.gy(1,m.gjmx()) - m.gy(1,1) << endl;

	Time_insac ins;
	ins.setup(&m, rho, mu, bcflags, bvalues, gradscheme, pressurescheme, refvel, cfl, tol, maxiter, totaltime);
	ins.solve();

	Array2d<double>* u = ins.getVariables();
	Array2d<double>* residuals = ins.getResiduals();

	Structdata2d strd(&m, u, residuals, "time_insac");
	strd.writevtk(outfile);

	cout << endl;
	return 0;
}
