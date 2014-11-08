#include <gtest/gtest.h>
#include "common.h"
#include "test_utils.h"
#include "LowOrdOper.h"

class cpu : public testing::Test
{
protected:

	double* solve_internal(ComputeParameters* p)
	{
		return compute_density(p->b, p->lb, p->rb, p->bb,
			p->ub, p->tau, p->t_count, p->x_size,
			p->y_size);
	}

	double* get_model_result(ComputeParameters* p, int lvl)
	{
		return solByEqualVolWithVarStepPlusPrint1(p->a, p->b, p->lb, p->rb, p->bb,
			p->ub, p->tau, p->t_count, p->x_size,
			p->y_size, lvl);
	}
};

TEST_F(cpu, test_to_model)
{
	const int first = 0;
	const int last = 2;
	double norm_test = 0;
	double norm_model = 0;
	ComputeParameters* p = new ComputeParameters();

	for (int lvl = first; lvl < last; ++lvl)
	{
		printf("level = %d\n", lvl);
		p->recompute_params(lvl);
		double* data = solve_internal(p);
		norm_test = p->norm;
		double* model = get_model_result(p, lvl);
		norm_model = p->norm;
		// print_matrix(model, p->x_length(), p->y_length());
		for (int i = 0; i < p->get_size(); i++)
		{
			ASSERT_NEAR(model[i], data[i], 1e-12);
		}
//		ASSERT_NEAR(norm_model, norm_test, 1e-12);
//		printf("norm test = %f\n", norm_test);
//		printf("norm model = %f\n", norm_model);

		delete[] data;
		delete[] model;
	}
	delete p;
}