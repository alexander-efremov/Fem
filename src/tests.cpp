#include <gtest/gtest.h>
#include "common.h"
#include "test_utils.h"
#include "LowOrdOper.h"
#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)
#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 1
  #define omp_get_num_threads() 1
#endif

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
	int first = 0, last = 3;
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
		printf("GPU\n");
		_print_matrix(data, 11, 11);
		printf("CPU\n");
		_print_matrix(model, 11, 11);

		double norm_model = p.norm;
//		if (lvl < 2)
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

TEST_F(cpu, china2015_test)
{
	int first = 4, last = 9;
//	int first = 6, last = 7;
	double time = 0;
	ComputeParameters p = ComputeParameters();
	for (int lvl = first; lvl < last; ++lvl)
	{
		printf("Start level %d\n", lvl);
		fflush(stdout);
		time = 0;
		p.recompute_params(lvl);
		double* data = solve_internal(p, time);
		delete[] data;
	}
}

TEST_F(cpu, valgrind_test)
{
//	int first = 4, last = 8;
	int first = 1, last = 2;
	double time = 0;
	ComputeParameters p = ComputeParameters();
	for (int lvl = first; lvl < last; ++lvl)
	{
		printf("Start level %d\n", lvl);
		fflush(stdout);
		time = 0;
		p.recompute_params(lvl);
		double* data = solve_internal(p, time);
		delete[] data;
	}
}
