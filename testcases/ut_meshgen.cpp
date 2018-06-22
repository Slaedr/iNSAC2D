#undef NDEBUG

#include "create_rect_grid.hpp"

int test_geomProg() {
	std::vector<a_real> pvec = geometricProgression(1.0, 2.0, 4);
	assert(pvec.size() == 5);

	constexpr a_real epsilon = std::numeric_limits<a_real>::epsilon();
	assert(std::fabs(pvec[0]) < epsilon);
	assert(std::fabs(pvec[1]-1.0/15) < epsilon);
	assert(std::fabs(pvec[2]-1.0/5) < epsilon);
	assert(std::fabs(pvec[3]-7.0/15) < epsilon);
	assert(std::fabs(pvec[4]-1.0) < epsilon);
	return 0;
}
