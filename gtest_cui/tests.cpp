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

struct dp_t1 {
	double x;
	double y;

	dp_t1(double x, double y)
		: x(x), y(y) {
	}

	dp_t1()
		: x(0), y(0) {
	}

	inline double operator[](int i) {
		return (i == 0) ? x : y;
	}

	inline dp_t1& operator=(const dp_t1& p) {
		if (this != &p) {
			x = p.x;
			y = p.y;
		}

		return *this;
	}

	inline dp_t1 operator+(const dp_t1& p) const {
		return dp_t1(p.x + x, p.y + y);
	}

	inline bool operator==(const dp_t1& p) const {
		return (x == p.x && y == p.y);
	}

	inline int operator<(const dp_t1& p) {
		return ((x < p.x) || ((x == p.x) && (y < p.y)));
	}

	inline int operator>(const dp_t1& p) {
		return ((x > p.x) || ((x == p.x) && (y > p.y)));
	}
};

inline void sort_by_y_desc_31(dp_t1* a)
{
	for (int i = 2; i < 4; i++)
	{
		for (int j = i; j > 1 && a[j - 1].y < a[j].y; j--)
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

inline void sort_by_x_asc1(dp_t1* a)
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

TEST_F(cpu, test_to_model)
{
	dp_t1 p[4];
	p[0] = dp_t1(4, 9);
	p[1] = dp_t1(3, 1);
	p[2] = dp_t1(2, 2);
	p[3] = dp_t1(1, 3);
	printf("x = %f y = %f\n", p[0].x,p[0].y);
	printf("x = %f y = %f\n", p[1].x,p[1].y);
	printf("x = %f y = %f\n", p[2].x,p[2].y);
	printf("x = %f y = %f\n", p[3].x,p[3].y);
	printf("================\n");
	sort_by_x_asc1(p);
	printf("x = %f y = %f\n", p[0].x, p[0].y);
	printf("x = %f y = %f\n", p[1].x, p[1].y);
	printf("x = %f y = %f\n", p[2].x, p[2].y);
	printf("x = %f y = %f\n", p[3].x, p[3].y);
	sort_by_y_desc_31(p);
	printf("================\n");
	printf("x = %f y = %f\n", p[0].x, p[0].y);
	printf("x = %f y = %f\n", p[1].x, p[1].y);
	printf("x = %f y = %f\n", p[2].x, p[2].y);
	printf("x = %f y = %f\n", p[3].x, p[3].y);
	
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