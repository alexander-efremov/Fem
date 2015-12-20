#include "common.h"
#include "point.h"
#include "utils.h"
#include <algorithm>
#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num() 1
	#define omp_get_num_threads() 1
#endif

#define sqr(x) ((x)*(x))
#define cub(x) ((x)*(x)*(x))
#define quad(x) ((x)*(x)*(x)*(x))

#ifdef __GNUC__
#define __pure
//#define __pure __attribute__((const))
/*почему то на GCC не работает и ломает результаты расчета*/
#elif __INTEL_COMPILER
#define __pure __declspec(const)
#elif __NVCC__
#define __pure
#else
#define __pure
#endif

static double B; //-V707
static double UB; //-V707
static double BB; //-V707
static double LB; //-V707
static double RB; //-V707
static double TAU;
static int OX_LEN;
static int OX_LEN_1; // OX_LEN_1
static int OY_LEN;
static int XY_LEN;
static int TIME_STEP_CNT;
static double HX; //-V707
static double HY; //-V707
static double* OX; //-V707
static double* OY; //-V707
static double* PREV_DENSITY;
static double TIME;
static double PREV_TIME; // tau * (tl - 1)
static double INVERTED_HX_HY;

__pure inline static double analytical_solution(double t, double x, double y)
{
	return 1.1 + sin(t * x * y);
}

inline static void __print_matrix11(double* a, int n, int m, int precision = 8) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			int k = i * n + j;
			switch (precision) {
			case 1:
				printf("%.1f ", a[k]);
				break;
			case 2:
				printf("%.2f ", a[k]);
				break;
			case 3:
				printf("%.3f ", a[k]);
				break;
			case 4:
				printf("%.4f ", a[k]);
				break;
			case 5:
				printf("%.5f ", a[k]);
				break;
			case 6:
				printf("%.6f ", a[k]);
				break;
			case 7:
				printf("%.7f ", a[k]);
				break;
			case 8:
				printf("%.8f ", a[k]);
				break;
			}
		}
		printf("\n");
	}
}

inline static void __print_vector(double* a, int n, int precision = 8) {
	for (int k = 0; k < n; ++k) 
	{
			switch (precision) 
			{
			case 1:
				printf("%.1f ", a[k]);
				break;
			case 2:
				printf("%.2f ", a[k]);
				break;
			case 3:
				printf("%.3f ", a[k]);
				break;
			case 4:
				printf("%.4f ", a[k]);
				break;
			case 5:
				printf("%.5f ", a[k]);
				break;
			case 6:
				printf("%.6f ", a[k]);
				break;
			case 7:
				printf("%.7f ", a[k]);
				break;
			case 8:
				printf("%.8f ", a[k]);
				break;
			}
	}
}

static int f = 0;

static double func_u()
{
	return 0;
}
static double func_v()
{
	return 0;
}

// Calculate coordintes in tk-1 square by coordintes from 
// ideal square and 4 points of tk-1 square
static dp_t get_point_real_pt_by_ideal_pt(dp_t p1, dp_t p2, dp_t p3, dp_t p4,
	dp_t ideal)
{
	dp_t r;
	double x = p1.x + (p2.x - p1.x)*ideal.x + (p4.x - p1.x)*ideal.y 
			+ (p1.x + p3.x - p2.x - p4.x)*ideal.x*ideal.y;
	double y = p1.y + (p2.y - p1.y)*ideal.x + (p4.y - p1.y)*ideal.y
			+ (p1.y + p3.y - p2.y - p4.y)*ideal.x*ideal.y;
	r.x = x;
	r.y = y;
	return r;
}
// Calculate density in real point by 4 points of tk-1 square
static double get_density_in_real_pt(dp_t p1, dp_t p2, dp_t p3, dp_t p4,
	dp_t real)
{
	// find out in which square real point was placed
	int i = static_cast<int>((real.x - FLT_MIN) / HX);
	if (real.x - FLT_MIN <= 0) i -= 1;
	int j = static_cast<int>((real.y - FLT_MIN) / HY);
	if (real.y - FLT_MIN <= 0) j -= 1;
	// if i,j = square index
	// compute 4 points of this square
	dp_t b_left;
	b_left.x = i*HX;
	b_left.y = j*HY;
	dp_t b_right;
	b_right.x = b_left.x + HX;
	b_right.y = b_left.y;
	dp_t u_left;
	u_left.x = b_left.x;
	u_left.y = b_left.y + HY;
	dp_t u_right;
	u_right.x = u_left.x + HX;
	u_right.y = u_left.y;

	double r = 0.;

	// formula 4
	r = PREV_DENSITY[j * OX_LEN_1 + i]*(1-(real.x - b_left.x)/HX)*(1 - (real.y - b_left.y)/HY)
	+ PREV_DENSITY[(j+1) * OX_LEN_1 + i]*((real.x - b_left.x)/HX)*(1 - (real.y - b_left.y)/HY)
	+ PREV_DENSITY[(j+1) * OX_LEN_1 + i + 1]*((real.x - b_left.x)/HX)*((real.y - b_left.y)/HY)
	+ PREV_DENSITY[j * OX_LEN_1 + i + 1]*(1 - (real.x - b_left.x)/HX)*((real.y - b_left.y)/HY);

	return r;
}
static double jakobian(dp_t p1, dp_t p2, dp_t p3, dp_t p4, dp_t ideal)
{
	double a11 = (p2.x - p1.x) + (p1.x + p3.x - p2.x - p4.x)*ideal.x;
	double a12 = (p4.x - p1.x) + (p1.x + p3.x - p2.x - p4.x)*ideal.y;
	double a21 = (p2.y - p1.y) + (p1.y + p3.y - p2.y - p4.y)*ideal.x;
	double a22 = (p4.y - p1.y) + (p1.y + p3.y - p2.y - p4.y)*ideal.y;
	double r = a11*a22 - a21*a12;
	return r;
}

static double* create_mesh_vector(int n)
{
	double* a = new double[n];
	double h = 1.0/n;
	for(int i = 0; i < n; ++i)
		a[i] = i*h;
	return a;
}

static double get_phi(int i, int j)
{
	dp_t u_left(OX[i-1], OY[j+1]); // left upper
	dp_t u_right(OX[i+1], OY[j+1]); // right upper
	dp_t b_left(OX[i-1], OY[j-1]); // left bottom
	dp_t b_right(OX[i+1], OY[j-1]); // right bottom
	dp_t center(OX[i], OY[j]);

	double u = func_u();
	double v = func_v();
	center.x = center.x - TAU * u;
	center.y = center.y - TAU * v;
	u = func_u();
	v = func_v();
	u_left.x = u_left.x - TAU * u;
	u_left.y = u_left.y - TAU * v;
	u = func_u();
	v = func_v();
	u_right.x = u_right.x - TAU * u;
	u_right.y = u_right.y - TAU * v;
	u = func_u();
	v = func_v();
	b_left.x = b_left.x - TAU * u;
	b_left.y = b_left.y - TAU * v;
	u = func_u();
	v = func_v();
	b_right.x = b_right.x - TAU * u;
	b_right.y = b_right.y - TAU * v;

	int nx = 10;
	int ny = 10;
	double x_step = 1.0/nx;
	double y_step = 1.0/ny;
	double* x_mesh = create_mesh_vector(nx);
	double* y_mesh = create_mesh_vector(ny);

	if (i == 1 && j == 1 && f == 0)
	{
		__print_vector(x_mesh, nx);
		__print_vector(y_mesh, ny);
		f = 1;
	}

	// get right part for jakoby
	double phi = 0;
	for(int i = 0; i < nx; ++i)
	{
		for(int j = 0; j < ny; ++j)
		{
			dp_t ideal;
			ideal.x = x_mesh[i] + x_step/2.;
			ideal.y = y_mesh[j] + y_step/2.;
			dp_t real = get_point_real_pt_by_ideal_pt(b_left, b_right, u_right, u_left, 
				ideal);
			double dens = get_density_in_real_pt(b_left, b_right, u_right, u_left, 
				real);
			double jakob = jakobian(b_left, b_right, u_right, u_left, ideal);
			phi += dens * jakob;
		}
	}
	phi = (16*phi/9*x_step*y_step);
	delete[] x_mesh;
	delete[] y_mesh;
	return phi;
}

inline static double get_norm_of_error(double* density, double ts_count_mul_steps)
{
	double r = 0;
	for (int k = 1; k < OY_LEN; ++k)
		for (int j = 1; j < OX_LEN; ++j)
			r += fabs(analytical_solution(ts_count_mul_steps, OX[j], OY[k])
				- density[(OY_LEN + 1) * k + j]);
	return HX * HY * r;
}

static void solve(double* density, double& time)
{
	PREV_DENSITY = new double[XY_LEN];
	for (int j = 0; j < OY_LEN + 1; j++)
	{
		for (int i = 0; i < OX_LEN_1; i++)
		{
			PREV_DENSITY[OX_LEN_1 * j + i] = analytical_solution(0, OX[i], OY[j]);
		}
	}

	int i = 0, j = 0, tl = 0;
	double timeStart = 0, timeEnd=0;
#ifdef _OPENMP
	// printf("OPENMP THREADS COUNT = %d\n", omp_get_max_threads());
	long count = 0;
	// dummy parallel section to get all threads running
	#pragma omp parallel private(i,j)
	{
		_InterlockedIncrement(&count);
	}
#endif

#ifdef _OPENMP
//	printf("OPENMP timer function is used!\n");
	timeStart = omp_get_wtime();
#else
//	printf("Standart timer function is used!\n");
	StartTimer();
#endif
	fflush(stdout);
	double* phi = new double[XY_LEN];
	for (tl = 1; tl <= TIME_STEP_CNT; tl++)
	{
		for (int k = 0; k <= OX_LEN; k++)
		{
			density[k] = analytical_solution(OX[k], BB, TIME);
			density[OX_LEN_1 * OY_LEN + k] = analytical_solution(OX[k], UB, TIME);
		}
		for (int u = 0; u <= OY_LEN; u++)
		{
			density[OX_LEN_1 * u] = analytical_solution(LB, OY[u], TIME);
			density[OX_LEN_1 * u + OX_LEN] = analytical_solution(RB, OY[u], TIME);
		}
#ifdef _OPENMP
	#pragma omp parallel for collapse(2) private(i, j)
#endif
		for (j = 1; j < OY_LEN; ++j)
			for (i = 1; i < OX_LEN; ++i)
				phi[OX_LEN_1 * j + i] = get_phi(i, j);
		// jakoby
		int iter_count = 20;
		int i = 0;
		while(i < iter_count)
		{
			for (j = 1; j < OY_LEN; ++j)
			{
				for (i = 1; i < OX_LEN; ++i)
				{
					density[OX_LEN_1 * j + i] = -1/9*(
						1.5*(
							PREV_DENSITY[OX_LEN_1 * j + i+1] +
							PREV_DENSITY[OX_LEN_1 *(j+1) + i] +
							PREV_DENSITY[OX_LEN_1 * j + i - 1] +
							PREV_DENSITY[OX_LEN_1 * (j-1) + i]) + 
						0.25*(
						PREV_DENSITY[OX_LEN_1 * (j+1) + i+1] + 
						PREV_DENSITY[OX_LEN_1 * (j+1) + i-1] + 
						PREV_DENSITY[OX_LEN_1 * (j-1) + i-1] + 
						PREV_DENSITY[OX_LEN_1 * (j-1) + i+1])) + 
						phi[OX_LEN_1 * j + i];
				}
			}
			i++;
		}

		memcpy(PREV_DENSITY, density, XY_LEN * sizeof(double));// заменить на быструю версию из agnerasmlib
	}
#ifdef _OPENMP
	timeEnd = omp_get_wtime();
	time = (timeEnd-timeStart);
//	printf("time %f s.\n", time);
#else
	time = GetTimer()/1000;
//	printf("time %f s.\n", time/1000);
#endif
	delete[] PREV_DENSITY;
	delete[] phi;
}

inline static void init(double b, double lb, double rb, double bb, double ub,
						double tau, int time_step_count, int ox_length, int oy_length)
{
	B = b;
	UB = ub;
	BB = bb;
	LB = lb;
	RB = rb;
	TAU = tau;
	TIME_STEP_CNT = time_step_count;
	XY_LEN = (ox_length + 1) * (oy_length + 1);
	OX_LEN = ox_length;
	OX_LEN_1 = ox_length + 1;
	OY_LEN = oy_length;
	OX = new double[OX_LEN_1];
	OY = new double[OY_LEN + 1];
	for (int i = 0; i <= OX_LEN; ++i) OX[i] = lb + i * (rb - lb) / OX_LEN;
	for (int i = 0; i <= OY_LEN; ++i) OY[i] = bb + i * (ub - bb) / OY_LEN;
	HX = OX[1] - OX[0];
	HY = OY[1] - OY[0];
	INVERTED_HX_HY = 1 / HX / HY;
}

inline static void clean()
{
	B = 0;
	UB = 0;
	BB = 0;
	LB = 0;
	RB = 0;
	TAU = 0;
	TIME = 0;
	OX_LEN = 0;
	OY_LEN = 0;
	OX_LEN_1 = 0;
	TIME_STEP_CNT = 0;
	XY_LEN = 0;
	HX = 0;
	HY = 0;
	INVERTED_HX_HY = 0;
	delete[] OX;
	delete[] OY;
}

double* compute_quad2_density(double b, double lb, double rb, double bb, double ub,
								double tau, int time_step_count, int ox_length, int oy_length, double& norm, double& time)
{
	init(b, lb, rb, bb, ub, tau, time_step_count, ox_length, oy_length);
	double* density = new double[XY_LEN];
	solve(density, time);
	norm = get_norm_of_error(density, TIME_STEP_CNT * TAU);
	clean();
	return density;
}
