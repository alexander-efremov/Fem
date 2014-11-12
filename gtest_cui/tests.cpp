#include <gtest/gtest.h>
#include "common.h"
#include "test_utils.h"
#include "LowOrdOper.h"

struct dp_t2 {
	double x;
	double y;

	dp_t2(double x, double y)
		: x(x), y(y) {
	}

	dp_t2()
		: x(0), y(0) {
	}

	inline double operator[](int i) {
		return (i == 0) ? x : y;
	}

	inline dp_t2& operator=(const dp_t2& p) {
		if (this != &p) {
			x = p.x;
			y = p.y;
		}

		return *this;
	}

	inline dp_t2 operator+(const dp_t2& p) const {
		return dp_t2(p.x + x, p.y + y);
	}

	inline bool operator==(const dp_t2& p) const {
		return (x == p.x && y == p.y);
	}

	inline int operator<(const dp_t2& p) {
		return ((x < p.x) || ((x == p.x) && (y < p.y)));
	}

	inline int operator>(const dp_t2& p) {
		return ((x > p.x) || ((x == p.x) && (y > p.y)));
	}
};

void sort_by_x_asc3(dp_t2* a)
{
	for (int i = 1; i < 4; i++)
	{
		for (int j = i; j > 0 && a[j - 1].x > a[j].x; j--)
		{
			double t = a[j].x;
			a[j].x = a[j - 1].x;
			a[j - 1].x = t;
			t = a[j].y;
			a[j].y = a[j - 1].y;
			a[j - 1].y = t;
		}
	}
}

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
	dp_t2 a[4];
	a[0] = dp_t2(5,1);
	a[1] = dp_t2(1,5);
	a[2] = dp_t2(-1,5);
	a[3] = dp_t2(-1,-1);
	printf("x = %f y = %f\n", a[0]);
	printf("x = %f y = %f\n", a[1]);
	printf("x = %f y = %f\n", a[2]);
	printf("x = %f y = %f\n", a[3]);
	printf("\n");
	sort_by_x_asc3(a);
	if (a[0].y > a[1].y)
	{
		double t = a[0].x;
		a[0].x = a[1].x;
		a[1].x = t;
		t = a[0].y;
		a[0].y = a[1].y;
		a[1].y = t;
	}
	if (a[2].y < a[3].y)
	{
		double t = a[2].x;
		a[2].x = a[3].x;
		a[3].x = t;
		t = a[2].y;
		a[2].y = a[3].y;
		a[3].y = t;
	}
	printf("x = %f y = %f\n", a[0]);
	printf("x = %f y = %f\n", a[1]);
	printf("x = %f y = %f\n", a[2]);
	printf("x = %f y = %f\n", a[3]);
}

TEST_F(cpu, DISABLED_test_to_model)
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