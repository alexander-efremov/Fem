#include <gtest/gtest.h>
#include "common.h"
#include "test_utils.h"
#include "LowOrdOper.h"

class cpu : public testing::Test
{
protected:

	double* solve_internal(ComputeParameters& p)
	{
		return compute_density(p.b, p.lb, p.rb, p.bb,
		                       p.ub, p.tau, p.t_count, p.x_size,
		                       p.y_size, p.norm);
	}

	double* get_model_result(ComputeParameters& p, int lvl)
	{
		return solByEqualVolWithVarStepPlusPrint1(p.a, p.b, p.lb, p.rb, p.bb,
		                                          p.ub, p.tau, p.t_count, p.x_size,
		                                          p.y_size, lvl);
	}
};

inline static double sign(const dp_t1& p1, const dp_t1 p2, const dp_t1 p3)
{
	return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

inline static bool is_points_belong_to_one_line(const dp_t1& p1, const dp_t1 p2, const dp_t1 p3)
{
	return sign(p1, p2, p3) == FLT_MIN;
}

inline static bool is_point_in_triangle(dp_t1 pt, dp_t1 v1, dp_t1 v2, dp_t1 v3)
{
	bool b1, b2, b3;
	b1 = sign(pt, v1, v2) < 0.0;
	b2 = sign(pt, v2, v3) < 0.0;
	b3 = sign(pt, v3, v1) < 0.0;
	return b1 == b2 && b2 == b3;
}


TEST_F(cpu, test_tri)
{
	dp_t1 p[4];
	p[0] = dp_t1(1,1);
	p[1] = dp_t1(1,3);
	p[2] = dp_t1(3,2);
	p[3] = dp_t1(1.72,2.01);
	bool b = is_point_in_triangle(p[3], p[0], p[1], p[2]);
	ASSERT_TRUE(b);

}



TEST_F(cpu, DISABLED_test_to_model)
//TEST_F(cpu, test_to_model)
{
	int first = 0, last = 3;
	double norm_test, norm_model;
	ComputeParameters p = ComputeParameters();
	for (int lvl = first; lvl < last; ++lvl)
	{
		p.recompute_params(lvl);
		double* data = solve_internal(p);
		double norm_test = p.norm;
		double* model = get_model_result(p, lvl);
		double norm_model = p.norm;

		if (lvl < 2)
			for (int i = 0; i < p.get_size(); i++)
			{
				ASSERT_NEAR(model[i], data[i], 1e-12);
			}
		printf("model norm = %f\n", norm_model);
		printf("test norm = %f\n", norm_test);
		ASSERT_NEAR(norm_model, norm_test, 1e-12);
		delete[] data;
		delete[] model;
	}
}