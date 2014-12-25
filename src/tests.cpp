#include <gtest/gtest.h>
#include "common.h"
#include "test_utils.h"
#include "LowOrdOper.h"
#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

class cpu : public testing::Test
{
public:
	cpu()
	{
	    #ifdef VER
		 printf("%s\n", STRINGIZE_VALUE_OF(VER));
	    #endif
	}
	
protected:
	double* solve_internal(ComputeParameters& p, double& time)
	{
		return compute_density(p.b, p.lb, p.rb, p.bb,
		                       p.ub, p.tau, p.t_count, p.x_size,
		                       p.y_size, p.norm, time);
	}
	
	double* solve_internal_cuda(ComputeParameters& p, float& time)
	{
		return compute_density_cuda(p.b, p.lb, p.rb, p.bb, p.ub, p.tau, p.t_count, p.x_size, p.y_size, p.norm, time);
	}
	

	double* get_model_result(ComputeParameters& p, int lvl)
	{
		return solByEqualVolWithVarStepPlusPrint1(p.a, p.b, p.lb, p.rb, p.bb,
		                                          p.ub, p.tau, p.t_count, p.x_size,
		                                          p.y_size, lvl);
	}
};

//TEST_F(cpu, DISABLED_test_to_model)
TEST_F(cpu, test_to_model)
{
	int first = 0, last = 1;
	double norm_test, norm_model;
	double time;
	ComputeParameters p = ComputeParameters();
	for (int lvl = first; lvl < last; ++lvl)
	{
		p.recompute_params(lvl);
		double* data = solve_internal(p, time);
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
		printf("Time: %lf\n", time);
		ASSERT_NEAR(norm_model, norm_test, 1e-12);
		delete[] data;
		delete[] model;
	}
	
}

//TEST_F(cpu, DISABLED_test_to_model)
TEST_F(cpu, test_to_model_cuda)
{
	int first = 0, last = 1;
	double norm_test, norm_model;
	float time = 0;
	ComputeParameters p = ComputeParameters();
	for (int lvl = first; lvl < last; ++lvl)
	{
		p.recompute_params(lvl);
		double* data = solve_internal_cuda(p, time);
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
		printf("Time: %lf\n", time);
		ASSERT_NEAR(norm_model, norm_test, 1e-12);
		delete[] data;
		delete[] model;
	}
}
