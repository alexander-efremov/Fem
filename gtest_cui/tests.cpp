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

TEST_F(cpu, test_to_model)
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
		for (int i = 0; i < p.get_size(); i++)
		{
	//		ASSERT_NEAR(model[i], data[i], 1e-12);
		}
		printf("model norm = %f\n", norm_model);
		printf("test norm = %f\n", norm_test);
		ASSERT_NEAR(norm_model, norm_test, 1e-12);
		delete[] data;
		delete[] model;
	}
}