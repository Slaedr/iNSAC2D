/** @file create_rect_grid.cpp
 * @brief Creates a rectilinear grid from data passed as command-line arguments.
 * 
 * Requires 5 arguments:
 * - Number of points in i-direction
 * - Number of points in j-direction
 * - Length of domain in i-direction
 * - Length of domain in j-direction
 * - Output mesh file name.
*/ 

#include <vector>
#include <structmesh2d.hpp>

using namespace amat;
using namespace acfd;
using namespace std;

/** Returns a sequence of N+1 points which divide a segment of length L into N intervals
 * according to a geometric progression with a given common ratio.
 */
std::vector<a_real> geometricProgression(const a_real len, const a_real ratio, const int N)
{
	// get first term of the geom prog
	const a_real a0 = (ratio-1.0)*len/(std::pow(ratio,N)-1.0);
	std::vector<a_real> points(N+1);
	for(int i = 0; i <= N; i++)
		points[i] = a0 * ratio^i;
	return points;
}

int main(int argc, char* argv[])
{
	if(argc < 4) {
		cout << "Need 5 arguments:\n"
		     << " No. of points in i-dir, number of points in j-dir, x-length, y-length\n"
		     << " and output file name."
		     << endl;
		return -1;
	}
	int imx = stoi(argv[1]);
	int jmx = stoi(argv[2]);
	double xlength = stod(argv[3]);
	double ylength = stod(argv[4]);

	ofstream fout(argv[5]);
	fout << setprecision(20);

	fout << imx << " " << jmx << '\n';
	for(int j = 0; j < jmx; j++)
		for(int i = 0; i < imx; i++)
			fout << xlength/(imx-1)*i << " " << ylength/(jmx-1)*j << '\n';
	
	fout.close();
	cout << endl;
	return 0;
}
