#include "common.h"
#include "point.h"
#include "utils.h"
#include "compute_density_cuda.cuh"
#include <algorithm>
#include <cuda.h>
#include <hemi.h>
__constant__ double C_B;
__constant__ double C_LB;
__constant__ double C_RB;
__constant__ double C_UB;
__constant__ double C_BB;
__constant__ double C_INVERTED_HX_HY;
__constant__ int C_HX;
__constant__ int C_HY;
__constant__ int C_OY_LEN;
__constant__ int C_OX_LEN;
__constant__ int C_OX_LEN_1;
__constant__ int C_XY_LEN;
__constant__ double C_PREV_TIME; // tau * (tl - 1)
__constant__ double C_TIME;
__constant__ double C_TAU;
__device__ double *OX_DEVICE, *OY_DEVICE;

#define sqr(x) ((x)*(x))
#define cub(x) ((x)*(x)*(x))
#define quad(x) ((x)*(x)*(x)*(x))

#ifdef __NVCC__
#define __pure __device__
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
static double TIME;
static double INVERTED_HX_HY;

/*
__pure inline static void print_params_const(int index, int needed_index,
	double b,
	double lb,
	double rb,
	double bb,
	double ub,
	double tau,
	int tl,
	int tl_count,
	int ox_length,
	int oy_length) {
	if (index == needed_index) {
		printf("index = %d\n", index);
		printf("b = %f\n", b);
		printf("lbDom = %f\n", lb);
		printf("rbDom = %f\n", rb);
		printf("bbDom = %f\n", bb);
		printf("ubDom = %f\n", ub);
		printf("tau = %f\n", tau);
		printf("Time level count = %d\n", tl_count);
		printf("current time level = %d\n", tl);
		printf("ox length = %d\n", ox_length + 1);
		printf("oy length = %d\n", oy_length + 1);
	}
}*/

__pure inline static void sort_by_y_asc(c_dp_t& x, c_dp_t& y, c_dp_t& z)
{	
	double t;
	if (x.y < y.y)
	{
		if (z.y < x.y) 
		{
			//swap(x, z);
			double t = x.x;
			x.x = z.x;
			z.x = t;
			t = x.y;
			x.y = z.y;
			z.y = t;
		}
	}
	else
	{
		if (y.y < z.y) 
		{
			//swap(x, y);
			t = x.x;
			x.x = y.x;
			y.x = t;
			t = x.y;
			x.y = y.y;
			y.y = t;
		}
		else 
		{
			//swap(x, z);
			t = x.x;
			x.x = z.x;
			z.x = t;
			t = x.y;
			x.y = z.y;
			z.y = t;
		}
	}
	if (z.y < y.y) 
	{
		//swap(y, z);
		t = y.x;
		y.x = z.x;
		z.x = t;
		t = y.y;
		y.y = z.y;
		z.y = t;
	}
}

__pure inline void sort_by_y(c_dp_t* a)
{
	for (int i = 1; i < 4; i++)
	{
		for (int j = i; j > 0 && a[j - 1].y > a[j].y; j--)
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

__pure inline void sort_by_y_desc_3(c_dp4_t* a)
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
			t = a[j].x_initial;
			a[j].x_initial = a[j - 1].x_initial;
			a[j - 1].x_initial = t;
			t = a[j].y_initial;
			a[j].y_initial = a[j - 1].y_initial;
			a[j - 1].y_initial = t;
		}
	}
}

__pure inline void sort_by_x_asc(c_dp4_t* a)
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
			t = a[j].x_initial;
			a[j].x_initial = a[j - 1].x_initial;
			a[j - 1].x_initial = t;
			t = a[j].y_initial;
			a[j].y_initial = a[j - 1].y_initial;
			a[j - 1].y_initial = t;
		}
	}

	if (a[0].y > a[1].y)
	{
		double t = a[0].x;
		a[0].x = a[1].x;
		a[1].x = t;
		t = a[0].y;
		a[0].y = a[1].y;
		a[1].y = t;
		t = a[0].x_initial;
		a[0].x_initial = a[1].x_initial;
		a[1].x_initial = t;
		t = a[0].y_initial;
		a[0].y_initial = a[1].y_initial;
		a[1].y_initial = t;
	}
	if (a[2].y < a[3].y)
	{
		double t = a[2].x;
		a[2].x = a[3].x;
		a[3].x = t;
		t = a[2].y;
		a[2].y = a[3].y;
		a[3].y = t;

		t = a[2].x_initial;
		a[2].x_initial = a[3].x_initial;
		a[3].x_initial = t;
		t = a[2].x_initial;
		a[2].x_initial = a[3].x_initial;
		a[3].x_initial = t;
	}
}

// получается порядок
/*

a[1]    a[2]
a[0]   a[3]
*/
__pure inline void sort_by_xy_wall_2(c_dp4_t* a)
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
			t = a[j].x_initial;
			a[j].x_initial = a[j - 1].x_initial;
			a[j - 1].x_initial = t;
			t = a[j].y_initial;
			a[j].y_initial = a[j - 1].y_initial;
			a[j - 1].y_initial = t;
		}
	}
}

__pure inline static bool try_get_slope_ratio(const c_dp_t& bv, const c_dp_t& uv, double& value)
{
	if (fabs(bv.x - uv.x) < 1e-12)
	{
		return false;
	}
	value = fabs((uv.y - bv.y) / (uv.x - bv.x)); // угловой коэффициент прямой
	if (value < 1e-12)
	{
		return false;
	}
	return true;
}


__pure inline static c_dp_t get_intersection_point(const c_dp4_t& alpha, const c_dp4_t& beta, const c_dp4_t& gamma, const c_dp4_t& theta)
{
	double a1 = gamma.y - alpha.y;
	double b1 = alpha.x - gamma.x; //double b1 = -(gamma.x - alpha.x);
	double c1 = a1 * alpha.x + b1 * alpha.y;
	double a2 = theta.y - beta.y;
	double b2 = beta.x - theta.x; //double b2 = -(theta.x - beta.x);
	double c2 = a2 * beta.x + b2 * beta.y;
	return c_dp_t((b1 * c2 - b2 * c1) / (b1 * a2 - b2 * a1), (a1 * c2 - a2 * c1) / (-b1 * a2 + b2 * a1));
}

// __pure inline static double sign(const c_dp_t& p1, const c_dp_t p2, const c_dp_t p3)
// {
// 	return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
// }

__pure inline static double sign(const c_dp4_t& p1, const c_dp4_t p2, const c_dp4_t p3)
{
	return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

__pure inline static bool is_points_belong_to_one_line(const c_dp4_t& p1, const c_dp4_t p2, const c_dp4_t p3)
{
	return sign(p1, p2, p3) == FLT_MIN ;
}

// __pure inline static bool is_point_in_triangle(c_dp_t pt, c_dp_t v1, c_dp_t v2, c_dp_t v3)
// {
// 	bool b1, b2, b3;
// 	b1 = sign(pt, v1, v2) < 0.0;
// 	b2 = sign(pt, v2, v3) < 0.0;
// 	b3 = sign(pt, v3, v1) < 0.0;
// 	return b1 == b2 && b2 == b3;
// }

__host__ __pure inline static double analytical_solution(double t, double x, double y)
{
	return 1.1 + sin(t * x * y);
}

__pure inline static double func_u(double b, double x, double y)
{
	return b * y * (1 - y) * (M_PI_2 + atan(-x));
}

__pure inline static double func_u(double b, const c_dp_t& p)
{
	return func_u(b, p.x, p.y);
}

__pure inline static double func_v(double ub, double bb, double lb, double rb, double time, double x, double y)
{
	return atan(0.1 * (x - lb) * (x - rb) * (1 + time) * (y - ub) * (y - bb));
}

__pure inline static double func_v(double ub, double bb, double lb, double rb, double time, const c_dp_t& p)
{
	return func_v(ub, bb, lb, rb, time, p.x, p.y);
}

__pure inline static double func_f(double b, double time, double ub, double bb, double lb, double rb, double x, double y)
{
	double arg_v = 0.1 * (x - lb) * (x - rb) * (1 + time) * (y - ub) * (y - bb);
	double rho = analytical_solution(time, x, y);
	double drho_dt = x * y * cos(time * x * y);
	double drho_dx = time * y * cos(time * x * y);
	double dtho_dy = time * x * cos(time * x * y);
	double u = func_u(b, x, y);
	double v = func_v(ub, bb, lb, rb, time, x, y);
	double du_dx = -b * y * (1 - y) / (1 + sqr(x));
	double dv_dx = 0.1 * (x - lb) * (x - rb) * (1 + time) * (y - bb + y - ub);
	dv_dx /= (1 + arg_v * arg_v);
	double res = drho_dt + rho * du_dx + u * drho_dx + rho * dv_dx + v * dtho_dy;
	// print_f_params()...
	return res;
}

__pure inline static double integrate_rectangle(double py, double qy, double gx, double hx, double a, double b)
{
	return 0.25 * (sqr(hx - a) - sqr(gx - a)) * (sqr(qy - b) - sqr(py - b));
}

__pure inline static double integrate_triangle(double py, double qy, double alpha, double beta, double a, double b)
{
	return (((qy - alpha) * cub(a * qy + b - beta) - (py - alpha) * cub(a * py + b - beta)) / (6 * a))
		- (quad(a * qy + b - beta) - quad(a * py + b - beta)) / (24 * sqr(a));
}
__device__ int flag;
__pure static double integrate_rectangle_one_cell(double* prev_dens, double py, double qy, double gx, double hx, const c_ip_t& sx, const c_ip_t& sy)
{
	double result, a, b;
	a = sx.y >= 0 && sy.y >= 0 ? OX_DEVICE[sx.y] : C_HX * sx.y;
	b = sx.y >= 0 && sy.y >= 0 ? OY_DEVICE[sy.y] : C_HY * sy.y; // ЭТО ПЛОТНОСТЬ С ПРЕДЫДУЩЕГО СЛОЯ ДЛЯ ДАННОЙ ЯЧЕЙКИ
	result = integrate_rectangle(py, qy, gx, hx, a, b) * (sx.x >= 0 && sy.x >= 0 ? prev_dens[C_OX_LEN_1 * sy.x + sx.x] : analytical_solution(C_PREV_TIME, sx.x * C_HX, sy.x * C_HY));
	a = sx.x >= 0 && sy.y >= 0 ? OX_DEVICE[sx.x] : C_HX * sx.x;
	b = sx.x >= 0 && sy.y >= 0 ? OY_DEVICE[sy.y] : C_HY * sy.y;
	result -= integrate_rectangle(py, qy, gx, hx, a, b) * (sx.x >= 0 && sy.x >= 0 ? prev_dens[C_OX_LEN_1 * sy.x + sx.y] : analytical_solution(C_PREV_TIME, sx.y * C_HX, sy.x * C_HY));
	a = sx.y >= 0 && sy.x >= 0 ? OX_DEVICE[sx.y] : C_HX * sx.y;
	b = sx.y >= 0 && sy.x >= 0 ? OY_DEVICE[sy.x] : C_HY * sy.x;
	result -= integrate_rectangle(py, qy, gx, hx, a, b) * (sx.x >= 0 && sy.x >= 0 ? prev_dens[C_OX_LEN_1 * sy.y + sx.x] : analytical_solution(C_PREV_TIME, sx.x * C_HX, sy.y * C_HY));
	a = sx.x >= 0 && sy.x >= 0 ? OX_DEVICE[sx.x] : C_HX * sx.x;
	b = sx.x >= 0 && sy.x >= 0 ? OY_DEVICE[sy.x] : C_HY * sy.x;
	result += integrate_rectangle(py, qy, gx, hx, a, b) * (sx.x >= 0 && sy.x >= 0 ? prev_dens[C_OX_LEN_1 * sy.y + sx.y] : analytical_solution(C_PREV_TIME, sx.y * C_HX, sy.y * C_HY));
	// if (flag == 1) 
	// {
	// 	printf("%s\n", "integrate_rectangle_one_cell");
	// 	printf("%lf\n", result);
	// }
	return result * C_INVERTED_HX_HY;
}

__pure static double integrate_triangle_left_one_cell(double* prev_dens, const c_dp_t& bv, const c_dp_t& uv, double hx,
                                               const c_ip_t& sx, const c_ip_t& sy)
{
	double a_sl = (bv.x - uv.x) / (bv.y - uv.y); //   Coefficients of slant line: x = a_SL *y  +  b_SL.
	if (fabs(a_sl) <= FLT_MIN) return 0;
	double b_sl = uv.x - a_sl * uv.y;
	double result = 0, tmp, alpha, beta;
	alpha = sx.y >= 0 && sy.y >= 0 ? OY_DEVICE[sy.y] : C_HY * sy.y;
	beta = sx.y >= 0 && sy.y >= 0 ? OX_DEVICE[sx.y] : C_HX * sx.y;
	tmp = 0.25 * (sqr(uv.y - OY_DEVICE[sy.y]) - sqr(bv.y - OY_DEVICE[sy.y])) * sqr(hx - beta) - integrate_triangle(bv.y, uv.y, alpha, beta, a_sl, b_sl);
	
	result += tmp * (sx.x >= 0 && sx.y <= C_OX_LEN && sy.x >= 0 && sy.y <= C_OY_LEN ? prev_dens[C_OX_LEN_1 * sy.x + sx.x] : analytical_solution(C_PREV_TIME, sx.x * C_HX, sy.x * C_HY));
	// if (flag == 1) 
	// {
	// 	printf("%s\n", "integrate_triangle_left_one_cell result 1");
	// 	printf("%lf\n", result);
	// }
	beta = sx.x >= 0 && sy.y >= 0 ? OX_DEVICE[sx.x] : C_HX * sx.x;
	tmp = sqr(uv.y - OY_DEVICE[sy.y]) - sqr(bv.y - OY_DEVICE[sy.y]);
	tmp = -0.25 * tmp * sqr(hx - beta) + integrate_triangle(bv.y, uv.y, alpha, beta, a_sl, b_sl);

	result += tmp * (sx.x >= 0 && sx.y <= C_OX_LEN && sy.x >= 0 && sy.y <= C_OY_LEN ? prev_dens[C_OX_LEN_1 * sy.x + sx.y] : analytical_solution(C_PREV_TIME, sx.y * C_HX, sy.x * C_HY));
	// if (flag == 1) 
	// {
	// 	printf("%s\n", "integrate_triangle_left_one_cell result 2");
	// 	printf("%lf\n", result);
	// }
	alpha = sx.y >= 0 && sy.x >= 0 ? OY_DEVICE[sy.x] : C_HY * sy.x;
	beta = sx.y >= 0 && sy.x >= 0 ? OX_DEVICE[sx.y] : C_HX * sx.y;
	tmp = sqr(uv.y - OY_DEVICE[sy.x]) - sqr(bv.y - OY_DEVICE[sy.x]);
	tmp = -0.25 * tmp * sqr(hx - beta) + integrate_triangle(bv.y, uv.y, alpha, beta, a_sl, b_sl);
	result += tmp * (sx.x >= 0 && sx.y <= C_OX_LEN && sy.x >= 0 && sy.y <= C_OY_LEN ? prev_dens[C_OX_LEN_1 * sy.y + sx.x] : analytical_solution(C_PREV_TIME, sx.x * C_HX, sy.y * C_HY));
		if (flag == 1) 
	{
		printf("%s\n", "integrate_triangle_left_one_cell result 3");
		printf("%lf\n", sy.y);
		printf("%lf\n", sx.x);
		printf("%lf\n", C_OX_LEN_1);
		printf("%lf\n", C_OX_LEN_1 * sy.y + sx.x);
		printf("%lf\n", prev_dens[C_OX_LEN_1 * sy.y + sx.x]);
		printf("%lf\n", result);
	}
	alpha = sx.x >= 0 && sy.x >= 0 ? OY_DEVICE[sy.x] : C_HY * sy.x;
	beta = sx.x >= 0 && sy.x >= 0 ? OX_DEVICE[sx.x] : C_HX * sx.x;
	tmp = sqr(uv.y - OY_DEVICE[sy.x]) - sqr(bv.y - OY_DEVICE[sy.x]);
	tmp = 0.25 * tmp * sqr(hx - beta) - integrate_triangle(bv.y, uv.y, alpha, beta, a_sl, b_sl);
	result += tmp * (sx.x >= 0 && sx.y <= C_OX_LEN && sy.x >= 0 && sy.y <= C_OY_LEN ? prev_dens[C_OX_LEN_1 * sy.y + sx.y] : analytical_solution(C_PREV_TIME, sx.y * C_HX, sy.y * C_HY));	
// if (flag == 1) 
// 	{
// 		printf("%s\n", "integrate_triangle_left_one_cell result 4");
// 		printf("%lf\n", tmp);
// 	}
	return result * C_INVERTED_HX_HY;
}

__pure static double integrate_right_slant_chanel(double* prev_dens, const c_dp_t& bv, const c_dp_t& uv, bool is_rect_truncated, const c_ip_t& sx, double b, const c_ip_t& sb, const c_ip_t& sy)
{
	if (fabs(uv.y - bv.y) <= FLT_MIN) return FLT_MIN ;
	double result = 0, gx = 0;
	double x = uv.x <= bv.x ? uv.x : bv.x;

	//   A. Under rectangle.
	result += -1 * integrate_triangle_left_one_cell(prev_dens, bv, uv, x, sx, sy);

	// case B: неполный прямоугольник    
	if (is_rect_truncated)
	{
		if (sx.x == sb.x) gx = b;
		if (sx.x > sb.x)
		{
			gx = sx.x >= 0 ? OX_DEVICE[sx.x] : C_HX * sx.x;
		}
		result += integrate_rectangle_one_cell(prev_dens, bv.y, uv.y, gx, x, sx, sy);
	}

	//   А теперь прибавим все прямоугольные куски, которые помещаются в ячейку
	c_ip_t ch_pos(sb.x, sb.x + 1);
	for (int j = sb.x; j < sx.x; j++)
	{
		if (j == sb.x) gx = b;
		else gx = ch_pos.x >= 0 ? OX_DEVICE[ch_pos.x] : C_HX * ch_pos.x;
		result += integrate_rectangle_one_cell(prev_dens, bv.y, uv.y, gx, ch_pos.x >= 0 ? OX_DEVICE[ch_pos.y] : C_HX * ch_pos.y, ch_pos, sy);
		ch_pos.x += 1;
		ch_pos.y = ch_pos.x + 1;
	}
	return result;
}

// используется для upper left и для bottom left треугольника
// т.е. случай
// UPPERLEFTTR
//
//                  CENTRE
//
// BOTTOMLEFTTR

__pure static double integrate_left_slant_chanel(double* prev_dens, const c_dp_t& bv, const c_dp_t& uv,
                                          bool is_rect_trunc, const c_ip_t& sx, const c_ip_t& sy,
                                          double b, const c_ip_t& sb)
{
	if (fabs(uv.y - bv.y) <= FLT_MIN) return FLT_MIN;
	double result = 0, hx = 0; //   -  Left and right boundary for each integration.   
	double x = uv.x <= bv.x ? bv.x : uv.x;

	// case A: triangle
	result += integrate_triangle_left_one_cell(prev_dens, bv, uv, x, sx, sy);

	// case B: не полный прямоугольник
	if (is_rect_trunc)
	{ // это значит, что прямоугольник занимает не всю ячейку  
		hx = sx.x == sb.x ? b : (sx.y >= 0 ? OX_DEVICE[sx.y] : C_HX * sx.y);
		result += integrate_rectangle_one_cell(prev_dens, bv.y, uv.y, x, hx, sx, sy);
	}

	//   А теперь прибавим все прямоугольные куски, которые помещаются в ячейку
	c_ip_t ch_pos(sx.x + 1, sx.x + 2); //   - координаты канала
	for (int j = sx.x + 1; j < sb.x + 1; j++)
	{
		hx = ch_pos.y <= 0 ? C_HX * ch_pos.y : hx = OX_DEVICE[ch_pos.y];
		if (j == sb.x) hx = b;
		result += integrate_rectangle_one_cell(prev_dens, bv.y, uv.y, ch_pos.y <= 0 ? C_HX * ch_pos.x : OX_DEVICE[ch_pos.x], hx, ch_pos, sy);
		ch_pos.x += 1;
		ch_pos.y = ch_pos.x + 1;
	}
	return result;
}

// определим целочисленные индексы квадратов в которых лежат верхняя и нижняя точки треугольника
// sx = (x,y) координаты квадрата в которой лежит нижняя точка
// sy = (x,y) координаты квадрата в которой лежит верхняя точка
// в случае успешной проверки, k = будет  угловой коэфициент прямой

__pure static double integrate_right_triangle_bottom_left(double* prev_dens, const c_dp_t& bv, const c_dp_t& uv)
{
	double k = 0;
	if (!try_get_slope_ratio(bv, uv, k)) return k;

	//   -  Index of current square by Ox and Oy axes. 
	c_ip_t sx, sy;
	sx.x = static_cast<int>((bv.x - FLT_MIN) / C_HX);
	if (bv.x - FLT_MIN <= 0) sx.x -= 1;
	sx.y = sx.x + 1;
	sy.x = static_cast<int>((bv.y + FLT_MIN) / C_HY);
	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
	sy.y = sy.x + 1;

	c_ip_t ib(sx.x, sx.x + 1);
	double result = 0;
	int curr_i = 0, next_i;
	c_dp_t curr = bv, next;
	while (true)
	{
		//TODO: sx.x и sx.y должны быть положительными всегда? Кажется для sx.x это всегда верно...
		double slope = sx.y >= 0 ? OY_DEVICE[sy.y] - curr.y : fabs(C_HY * sy.y - curr.y);
		slope /= sx.x >= 0 ? curr.x - OX_DEVICE[sx.x] : fabs(curr.x - C_HX * sx.x);
		if (slope <= k)
		{
			next_i = 1;
			next.y = sy.y >= 0 ? OY_DEVICE[sy.y] : C_HY * sy.y;
			next.x = curr.x - (next.y - curr.y) / k;
		}
		else
		{
			next_i = 2;
			next.x = sx.x >= 0 ? OX_DEVICE[sx.x] : C_HX * sx.x;
			next.y = curr.y - k * (next.x - curr.x);
		}
		if (next.x - uv.x < FLT_MIN)
		{
			// сюда попадаем и в случае когда треугольник полностью в одной ячейке лежит
			// и в случае когда прошлись по всем точкам...
			result += integrate_left_slant_chanel(prev_dens, curr, uv, (uv.x <= curr.x ? curr_i : 0) == 1, sx, sy, bv.x, ib);
			break;
		}
		result += integrate_left_slant_chanel(prev_dens, curr, next, (next.x <= curr.x ? curr_i : next_i) == 1, sx, sy, bv.x, ib);
		switch (next_i)
		{
		case 1:
			sy += 1;
			break;
		case 2:
			sx -= 1;
			break;
		}
		curr_i = next_i;
		curr = next;
	}
	return result;
}

__pure static double integrate_right_triangle_bottom_right(double* prev_dens, const c_dp_t& bv, const c_dp_t& uv)
{
	double k = 0;
	if (!try_get_slope_ratio(bv, uv, k)) return k;

	c_ip_t sx, sy;
	sx.x = static_cast<int>((bv.x + FLT_MIN) / C_HX);
	if (bv.x + FLT_MIN <= 0) sx.x -= 1;
	sx.y = sx.x + 1;
	sy.x = static_cast<int>((bv.y + FLT_MIN) / C_HY);
	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
	sy.y = sy.x + 1;

	c_ip_t ib(sx.x, sx.x + 1);
	double result = 0;
	int curr_i = 0, next_i;
	c_dp_t curr = bv, next;
	while (true)
	{
		double slope = sy.y >= 0 ? fabs(OY_DEVICE[sy.y] - curr.y) : fabs(C_HY * sy.y - curr.y);
		slope /= sx.y >= 0 ? fabs(OX_DEVICE[sx.y] - curr.x) : fabs(C_HX * sx.y - curr.x);
		if (slope <= k)
		{
			next_i = 1;
			next.y = sy.y >= 0 ? OY_DEVICE[sy.y] : C_HY * sy.y;
			next.x = bv.x + (next.y - bv.y) / k;
		}
		else
		{
			next_i = 2;
			next.x = sx.y >= 0 ? OX_DEVICE[sx.y] : C_HX * sx.y;
			next.y = bv.y + k * (next.x - bv.x);
		}
		if (next.x - uv.x > FLT_MIN)
		{
			result += integrate_right_slant_chanel(prev_dens, curr, uv, (uv.x <= curr.x ? 0 : curr_i) == 1, sx, bv.x, ib, sy);
			break;
		}
		result += integrate_right_slant_chanel(prev_dens, curr, next, (next.x <= curr.x ? next_i : curr_i) == 1, sx, bv.x, ib, sy);
		switch (next_i)
		{
		case 1:
			sy += 1;
			break;
		case 2:
			sx += 1;
			break;
		}
		curr_i = next_i;
		curr = next;
	}
	return result;
}

// __pure static double integrate_right_triangle_bottom_right_wall(const c_dp_t& bv, const c_dp_t& uv)
// {
// 	double k = 0;
// 	if (!try_get_slope_ratio(bv, uv, k)) return k;

// 	c_ip_t sx, sy;
// 	sx.x = static_cast<int>((bv.x + FLT_MIN) / C_HX);
// 	if (bv.x + FLT_MIN <= 0) sx.x -= 1;
// 	sx.y = sx.x + 1;
// 	sy.x = static_cast<int>((bv.y + FLT_MIN) / C_HY);
// 	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
// 	sy.y = sy.x + 1;

// 	c_ip_t ib(sx.x, sx.x + 1);
// 	double result = 0;
// 	int curr_i = 0, next_i;
// 	c_dp_t curr = bv, next;
// 	while (true)
// 	{
// 		double slope = sy.y >= 0 ? fabs(OY_DEVICE[sy.y] - curr.y) : fabs(C_HY * sy.y - curr.y);
// 		slope /= sx.y >= 0 ? fabs(OX_DEVICE[sx.y] - curr.x) : fabs(C_HX * sx.y - curr.x);
// 		if (slope <= k)
// 		{
// 			next_i = 1;
// 			next.y = sy.y >= 0 ? OY_DEVICE[sy.y] : C_HY * sy.y;
// 			next.x = bv.x + (next.y - bv.y) / k;
// 		}
// 		else
// 		{
// 			next_i = 2;
// 			next.x = sx.y >= 0 ? OX_DEVICE[sx.y] : C_HX * sx.y;
// 			next.y = bv.y + k * (next.x - bv.x);
// 		}
// 		if (next.x - uv.x > FLT_MIN)
// 		{
// 			result += integrate_right_slant_chanel(curr, uv, (uv.x <= curr.x ? 0 : curr_i) == 1, sx, bv.x, ib, sy);
// 			break;
// 		}
// 		result += integrate_right_slant_chanel(curr, next, (next.x <= curr.x ? next_i : curr_i) == 1, sx, bv.x, ib, sy);
// 		switch (next_i)
// 		{
// 		case 1:
// 			sy += 1;
// 			break;
// 		case 2:
// 			sx += 1;
// 			break;
// 		}
// 		curr_i = next_i;
// 		curr = next;
// 	}
// 	return result;
// }

__pure static double integrate_right_triangle_upper_left(double* prev_dens, const c_dp_t& bv, const c_dp_t& uv)
{
	double k = 0;
	if (!try_get_slope_ratio(bv, uv, k)) return k;

	c_ip_t sx, sy, ib;
	sx.x = static_cast<int>((bv.x + FLT_MIN) / C_HX); //   -  If bv.x is in grid edge I want it will be in the right side.
	if (bv.x + FLT_MIN <= 0) sx.x -= 1;
	sx.y = sx.x + 1;
	sy.x = static_cast<int>((bv.y + FLT_MIN) / C_HY); //   -  If bv.y is in grid edge I want it will be in the upper square.
	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
	sy.y = sy.x + 1;
	ib.x = static_cast<int>((uv.x - FLT_MIN) / C_HY); //   -  If uv.x is in grid edge I want it will be in the left side.
	if (uv.x - FLT_MIN <= 0) ib.x -= 1;
	ib.y = ib.x + 1;

	double result = 0;
	int curr_i = 0, next_i;
	c_dp_t curr = bv, next;
	while (true)
	{
		double slope = sy.y >= 0 ? OY_DEVICE[sy.y] - curr.y : fabs(C_HY * sy.y - curr.y);
		slope /= sx.y >= 0 ? OX_DEVICE[sx.y] - curr.x : fabs(C_HX * sx.y - curr.x);
		if (slope <= k)
		{
			next_i = 1;
			next.y = sy.y >= 0 ? OY_DEVICE[sy.y] : C_HY * sy.y;
			next.x = bv.x + (next.y - bv.y) / k;
		}
		else
		{
			next_i = 2;
			next.x = sx.y >= 0 ? OX_DEVICE[sx.y] : C_HX * sx.y;
			next.y = bv.y + k * (next.x - bv.x);
		}
		if (next.x - uv.x > FLT_MIN) // если следующая точка уже правее, чем наша граничная точка, то мы обработали канал
		{
			result += integrate_left_slant_chanel(prev_dens, curr, uv, (uv.x <= curr.x ? curr_i : 0) == 1, sx, sy, uv.x, ib);
			break;
		}
		result += integrate_left_slant_chanel(prev_dens, curr, next, (next.x <= curr.x ? curr_i : next_i) == 1, sx, sy, uv.x, ib);

		switch (next_i)
		{
		case 1:
			sy += 1;
			break;
		case 2:
			sx += 1;
			break;
		}
		curr_i = next_i;
		curr = next;
	}
	return result;
}

__pure static double integrate_right_triangle_upper_right(double* prev_dens, const c_dp_t& bv, const c_dp_t& uv)
{
	double k = 0;
	if (!try_get_slope_ratio(bv, uv, k)) return k;

	c_ip_t sx, sy, ib;
	sx.x = static_cast<int>((bv.x - FLT_MIN) / C_HX); //   -  If bv.x is in grid edge I want it will be between in the left side.
	if (bv.x - FLT_MIN <= 0) sx.x -= 1;
	sx.y = sx.x + 1;
	sy.x = static_cast<int>((bv.y + FLT_MIN) / C_HY); //   -  If bv.y is in grid edge I want it will be in the upper side.
	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
	sy.y = sy.x + 1;
	ib.x = static_cast<int>((uv.x + FLT_MIN) / C_HX);
	if (uv.x + FLT_MIN <= 0) ib.x -= 1;
	ib.y = ib.x + 1;

	double result = 0;
	int curr_i = 0, next_i;
	c_dp_t curr = bv, next;
	while (true)
	{
		double slope = sy.y >= 0 ? fabs(OY_DEVICE[sy.y] - curr.y) : fabs(C_HY * sy.y - curr.y);
		slope /= sx.x >= 0 ? fabs(curr.x - OX_DEVICE[sx.x]) : fabs(curr.x - C_HX * sx.x);
		if (slope <= k)
		{
			next_i = 1;
			next.y = sy.y >= 0 ? OY_DEVICE[sy.y] : C_HY * sy.y;
			next.x = bv.x - (next.y - bv.y) / k;
		}
		else
		{
			next_i = 2;
			next.x = sx.x >= 0 ? OX_DEVICE[sx.x] : C_HX * sx.x;
			next.y = bv.y - k * (next.x - bv.x);
		}
		if (next.x - uv.x < FLT_MIN)
		{
			result += integrate_right_slant_chanel(prev_dens, curr, uv, (uv.x <= curr.x ? 0 : curr_i) == 1, sx, uv.x, ib, sy);
			break;
		}
		result += integrate_right_slant_chanel(prev_dens, curr, next, (next.x <= curr.x ? next_i : curr_i) == 1, sx, uv.x, ib, sy);
		switch (next_i)
		{
		case 1:
			sy += 1;
			break;
		case 2:
			sx -= 1;
			break;
		}
		curr_i = next_i;
		curr = next;
	}
	return result;
}

// __pure static double integrate_right_triangle_upper_left_wall(const c_dp_t& bv, const c_dp_t& uv)
// {
// 	double k = 0;
// 	if (!try_get_slope_ratio(bv, uv, k)) return k;

// 	c_ip_t sx, sy, ib;
// 	sx.x = static_cast<int>((bv.x + FLT_MIN) / C_HX); //   -  If bv.x is in grid edge I want it will be in the right side.
// 	if (bv.x + FLT_MIN <= 0) sx.x -= 1;
// 	sx.y = sx.x + 1;
// 	sy.x = static_cast<int>((bv.y + FLT_MIN) / C_HY); //   -  If bv.y is in grid edge I want it will be in the upper square.
// 	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
// 	sy.y = sy.x + 1;
// 	ib.x = static_cast<int>((uv.x - FLT_MIN) / C_HY); //   -  If uv.x is in grid edge I want it will be in the left side.
// 	if (uv.x - FLT_MIN <= 0) ib.x -= 1;
// 	ib.y = ib.x + 1;

// 	double result = 0;
// 	int curr_i = 0, next_i;
// 	c_dp_t curr = bv, next;
// 	while (true)
// 	{
// 		double slope = sy.y >= 0 ? OY_DEVICE[sy.y] - curr.y : fabs(C_HY * sy.y - curr.y);
// 		slope /= sx.y >= 0 ? OX_DEVICE[sx.y] - curr.x : fabs(C_HX * sx.y - curr.x);
// 		if (slope <= k)
// 		{
// 			next_i = 1;
// 			next.y = sy.y >= 0 ? OY_DEVICE[sy.y] : C_HY * sy.y;
// 			next.x = bv.x + (next.y - bv.y) / k;
// 		}
// 		else
// 		{
// 			next_i = 2;
// 			next.x = sx.y >= 0 ? OX_DEVICE[sx.y] : C_HX * sx.y;
// 			next.y = bv.y + k * (next.x - bv.x);
// 		}
// 		if (next.x - uv.x > FLT_MIN) // если следующая точка уже правее, чем наша граничная точка, то мы обработали канал
// 		{
// 			result += integrate_left_slant_chanel(curr, uv, (uv.x <= curr.x ? curr_i : 0) == 1, sx, sy, uv.x, ib);
// 			break;
// 		}
// 		result += integrate_left_slant_chanel(curr, next, (next.x <= curr.x ? curr_i : next_i) == 1, sx, sy, uv.x, ib);

// 		switch (next_i)
// 		{
// 		case 1:
// 			sy += 1;
// 			break;
// 		case 2:
// 			sx += 1;
// 			break;
// 		}
// 		curr_i = next_i;
// 		curr = next;
// 	}
// 	return result;
// }

// __pure static double integrate_right_triangle_upper_right_wall(const c_dp_t& bv, const c_dp_t& uv)
// {
// 	double k = 0;
// 	if (!try_get_slope_ratio(bv, uv, k)) return k;

// 	c_ip_t sx, sy, ib;
// 	sx.x = static_cast<int>((bv.x - FLT_MIN) / C_HX); //   -  If bv.x is in grid edge I want it will be between in the left side.
// 	if (bv.x - FLT_MIN <= 0) sx.x -= 1;
// 	sx.y = sx.x + 1;
// 	sy.x = static_cast<int>((bv.y + FLT_MIN) / C_HY); //   -  If bv.y is in grid edge I want it will be in the upper side.
// 	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
// 	sy.y = sy.x + 1;
// 	ib.x = static_cast<int>((uv.x + FLT_MIN) / C_HX);
// 	if (uv.x + FLT_MIN <= 0) ib.x -= 1;
// 	ib.y = ib.x + 1;

// 	double result = 0;
// 	int curr_i = 0, next_i;
// 	c_dp_t curr = bv, next;
// 	while (true)
// 	{
// 		double slope = sy.y >= 0 ? fabs(OY_DEVICE[sy.y] - curr.y) : fabs(C_HY * sy.y - curr.y);
// 		slope /= sx.x >= 0 ? fabs(curr.x - OX_DEVICE[sx.x]) : fabs(curr.x - C_HX * sx.x);
// 		if (slope <= k)
// 		{
// 			next_i = 1;
// 			next.y = sy.y >= 0 ? OY_DEVICE[sy.y] : C_HY * sy.y;
// 			next.x = bv.x - (next.y - bv.y) / k;
// 		}
// 		else
// 		{
// 			next_i = 2;
// 			next.x = sx.x >= 0 ? OX_DEVICE[sx.x] : C_HX * sx.x;
// 			next.y = bv.y - k * (next.x - bv.x);
// 		}
// 		if (next.x - uv.x < FLT_MIN)
// 		{
// 			result += integrate_right_slant_chanel(curr, uv, (uv.x <= curr.x ? 0 : curr_i) == 1, sx, uv.x, ib, sy);
// 			break;
// 		}
// 		result += integrate_right_slant_chanel(curr, next, (next.x <= curr.x ? next_i : curr_i) == 1, sx, uv.x, ib, sy);
// 		switch (next_i)
// 		{
// 		case 1:
// 			sy += 1;
// 			break;
// 		case 2:
// 			sx -= 1;
// 			break;
// 		}
// 		curr_i = next_i;
// 		curr = next;
// 	}
// 	return result;
// }

__pure static double integrate_bottom_triangle(double* prev_dens, const c_dp_t& l, const c_dp_t& m, const c_dp_t& r)
{
	double result = 0;
	if (m.x == l.x)
	{
		result = integrate_right_triangle_bottom_right(prev_dens, m, r);
	}
	else if (m.x == r.x)
	{
		result = integrate_right_triangle_bottom_left(prev_dens, m, l);
	}
	else if (m.x < l.x)
	{
		result = integrate_right_triangle_bottom_right(prev_dens, m, r) - integrate_right_triangle_bottom_right(prev_dens, m, l);
	}
	else if (m.x > l.x && m.x < r.x)
	{
		result = integrate_right_triangle_bottom_left(prev_dens, m, l) + integrate_right_triangle_bottom_right(prev_dens, m, r);
	}
	else if (m.x > r.x)
	{
		result = integrate_right_triangle_bottom_left(prev_dens, m, l) - integrate_right_triangle_bottom_left(prev_dens, m, r);
	}
	return result;
}

__pure static double integrate_upper_triangle(double* prev_dens, const c_dp_t& l, const c_dp_t& m, const c_dp_t& r)
{
	double result = 0;
	if (m.x == l.x)
	{
		result = integrate_right_triangle_upper_right(prev_dens, r, m);
	}
	else if (m.x == r.x)
	{
		result = integrate_right_triangle_upper_left(prev_dens, l, m);
	}
	else if (m.x < l.x)
	{
		result = integrate_right_triangle_upper_right(prev_dens, r, m) - integrate_right_triangle_upper_right(prev_dens, l, m);
	}
	else if (m.x > l.x && m.x < r.x)
	{
		result = integrate_right_triangle_upper_left(prev_dens, l, m) + integrate_right_triangle_upper_right(prev_dens, r, m);
	}
	else if (m.x > r.x)
	{
		result = integrate_right_triangle_upper_left(prev_dens, l, m) - integrate_right_triangle_upper_left(prev_dens, r, m);
	}
	return result;
}

// x,y,z
__pure static double integrate_uniform_triangle(double* prev_dens, const c_dp_t& x, const c_dp_t& y, const c_dp_t& z)
{
	// точки должны идти в порядке возрастания y координаты, чтобы правильно отработала процедура интегрирования		

	//   a * x  +  b * y  = c.
	double a = z.y - x.y;
	if (fabs(a) < FLT_MIN) return FLT_MIN;
	double b = x.x - z.x;
	double c = b * x.y + a * x.x;
	c_dp_t ip((c - b * y.y) / a, y.y);

	//   Возможны 2 случая расположения точки пересечения относительно средней
	//   слева или справа.
	//   есди средняя точка справа от точки пересечения
	//   обменяем местами  X координаты, чтобы использовать один код для расчета
	c_dp_t t = y;
	if (t.x >= ip.x)
	{
		double tx = t.x;
		t.x = ip.x;
		ip.x = tx;
	}

	return integrate_upper_triangle(prev_dens, t, z, ip) + integrate_bottom_triangle(prev_dens, t, x, ip);
}

// __pure static double integrate_uniform_triangle_wall(const c_dp_t& x, const c_dp_t& y,
//                                               const c_dp_t& z, quad_type type)
// {
// 	// для точек YOt оси координат обозначены по другому
// 	// по OX будет откладываться значение y
// 	// по OY будет откладываться значение t
// 	// т.е. будет плоскость YOt
// 	// здесь y координата точки должна быть временем	
// 	switch (type)
// 	{
// 	case wall_1_middle_at:
// 	case wall_1_middle_in:
// 	case wall_1_middle_out:		
// 			// !phd\2014\fem\ggb\wa1\4.ggb
// 			// случай А
// 			if (x.x >= y.x)
// 			{
// 				double res = 0;
// 				double t = integrate_right_triangle_upper_left_wall(z, x);
// 				res += t;
// 				t = integrate_right_triangle_upper_left_wall(y, x);
// 				res += t;
// 				return res;
// 			}
// 			if (x.x < y.x && x.x > z.x) // случай B
// 			{
// 				double res = 0;
// 				double t = integrate_right_triangle_upper_left_wall(z, x);
// 				res += t;
// 				t = integrate_right_triangle_upper_right_wall(y, x);
// 				res += t;
// 				return res;
// 			}
// 			if (x.x <= z.x) // случай C
// 			{
// 				double res = 0;
// 				double t = integrate_right_triangle_upper_right_wall(y, x);
// 				res += t;
// 				t = integrate_right_triangle_upper_left_wall(z, x);
// 				res += t;
// 				return res;
// 			}
// 	case wall_2:
// 	{
// 		double t = 0;
// 		double res = 0;
// 		t = integrate_right_triangle_upper_right_wall(x, y);
// 		res += t;
// 		t = integrate_right_triangle_bottom_right_wall(y, z);
// 		res += t;
// 		return res;
// 	}
// 		break;

// 	default:
// 		return 0;
// 	}
// }

// __pure inline int get_wall_intersection_type_as_int(c_dp4_t* a)
// {
// 	int type = -1;
// 	bool is_four_point_on_the_wall = a[0].x <= 0 && a[1].x <= 0 && a[2].x <= 0 && a[3].x <= 0;
// 	bool is_three_point_on_the_wall = a[0].x <= 0 && a[1].x <= 0 && a[2].x <= 0 && a[3].x > 0;
// 	bool is_two_point_on_the_wall = a[0].x <= 0 && a[1].x <= 0 && a[2].x > 0 && a[3].x > 0;
// 	bool is_one_point_on_the_wall = a[0].x <= 0 && a[1].x > 0 && a[2].x > 0 && a[3].x > 0;
// 	if (is_four_point_on_the_wall)
// 	{
// 		type = 4;
// 	}
// 	else if (is_three_point_on_the_wall)
// 	{
// 		type = 3;
// 	}
// 	else if (is_two_point_on_the_wall)
// 	{
// 		type = 2;
// 	}
// 	else if (is_one_point_on_the_wall)
// 	{
// 		type = 1;
// 	}
// 	else
// 	{
// 		type = 0;
// 	}
// 	return type;
// }

// __pure inline static quad_type get_wall_intersection_type(c_dp4_t* a)
// {
// 	/*
// 	 формулы расчетов здесь http://www.pm298.ru/reshenie/fha0327.php
// 	 a[0] - alpha
// 	 a[1] - beta
// 	 a[2] - gamma
// 	 a[3] - theta
// 	 a[4] - mu
// 	 a[5] - nu

// 	 */

// 	int type = get_wall_intersection_type_as_int(a);
// 	switch (type)
// 	{
// 	case 4:
// 		return wall_4;
// 	case 3:
// 		{
// 			sort_by_x_asc(a);
// 			sort_by_y_desc_3(a);
// 			// рассчитаем точку пересечения OY и прямой a[0]:a[3]
// 			// тут не надо fabs, потому что a[3].x > a[0].x
// 			double y = a[3].x - a[0].x < FLT_MIN ? 0.5 * (a[0].y + a[3].y) : a[0].y - a[0].x * ((a[3].y - a[0].y) / (a[3].x - a[0].x));
// 			a[4] = c_dp4_t(0, y); // mu

// 			// рассчитаем точку пересечения OY и прямой a[2]:a[3]
// 			// тут не надо fabs, потому что a[3].x > a[2].x
// 			y = a[3].x - a[2].x < FLT_MIN ? 0.5 * (a[3].y + a[2].y) : a[2].y - a[2].x * ((a[3].y - a[2].y) / (a[3].x - a[2].x));
// 			a[5] = c_dp4_t(0, y); // nu

// 			if ((a[0].x - a[2].x) * (a[1].y - a[2].y) - (a[1].x - a[2].x) * (a[0].y - a[2].y) < FLT_MIN)
// 				return wall_3_middle_at;
// 			if (a[0].x < a[1].x && a[1].x > a[2].x)
// 				return wall_3_middle_out;
// 			if (a[1].x < a[0].x && a[1].x < a[2].x)
// 				return wall_3_middle_in;
// 		}
// 	case 2:
// 		{
// 			sort_by_xy_wall_2(a);
// 			double y = 0;
// 			if (a[2].x - a[1].x < FLT_MIN)
// 			{
// 				y = (a[1].y + a[2].y) * 0.5;
// 			}
// 			else
// 			{
// 				y = a[1].y - a[1].x * ((a[2].y - a[1].y) / (a[2].x - a[1].x));
// 			}
// 			a[4] = c_dp4_t(0, y); // mu

// 			if (a[3].x - a[0].x < FLT_MIN)
// 			{
// 				y = (a[0].y + a[3].y) * 0.5;
// 			}
// 			else
// 			{
// 				y = a[0].y - a[0].x * ((a[3].y - a[0].y) / (a[3].x - a[0].x));
// 			}
// 			a[5] = c_dp4_t(0, y); // nu
// 			return wall_2;
// 		}

// 	case 1: 
// 		{
// 			sort_by_x_asc(a);
// 			sort_by_y_desc_3(a);

// 			// для точек на стенке оси координат обозначены по другому
// 			// по OX будет откладываться значение y
// 			// по OY будет откладываться значение t
// 			// точка a[0] - точка, которая упала на стенку
// 			// считаем y компоненту
// 			a[0].x = a[0].y_initial - a[0].x_initial * func_v(C_UB, C_BB, C_LB, C_RB, C_TIME, a[0].x_initial, a[0].y_initial) / func_u(C_B, a[0].x_initial, a[0].y_initial);
// 			a[0].y = C_TIME - a[0].x_initial*(1 / func_u(C_B, a[0].x_initial, a[0].y_initial)); // здесь будет Mt = tk - Ax*tg(alpha); tg(alpha) = 1 / U(A)
			

// 			// рассчитаем точку пересечения OY и прямой a[0]:a[1]
// 			// тут не надо fabs, потому что a[1].x > a[0].x
// 			double y = a[1].x - a[0].x < FLT_MIN ? 0.5 * (a[0].y + a[1].y) : (a[0].y - a[0].x * ((a[1].y - a[0].y) / (a[1].x - a[0].x)));
// 			a[4] = c_dp4_t(0, y); // mu

// 			// рассчитаем точку пересечения OY и прямой a[0]:a[3]
// 			// тут не надо fabs, потому что a[3].x > a[0].x
// 			y = a[3].x - a[0].x < FLT_MIN ? 0.5 * (a[0].y + a[3].y) : (a[0].y - a[0].x * ((a[3].y - a[0].y) / (a[3].x - a[0].x)));
// 			a[5] = c_dp4_t(0, y); // nu

// 			if (is_points_belong_to_one_line(a[1], a[2], a[3]))
// 				return wall_1_middle_at;
// 			if (a[1].x < a[2].x && a[2].x > a[3].x)
// 				return wall_1_middle_out;
// 			if (a[2].x < a[1].x && a[2].x < a[3].x)
// 				return wall_1_middle_in;
// 			break;
// 		}
// 	default:
// 		return normal;
// 	}
// 	return normal;
// }

__pure static quad_type get_quadrangle_type(int i, int j,
                                     c_dp_t& a, c_dp_t& b, c_dp_t& c, c_dp_t& k, c_dp_t& m, c_dp_t& n, c_dp4_t* p)
{
	// TODO какой порядок тут все таки предполагется? против часовой начиная с верхней левой?	
	c_dp_t alpha((OX_DEVICE[i - 1] + OX_DEVICE[i]) * 0.5, (OY_DEVICE[j - 1] + OY_DEVICE[j]) * 0.5),
		beta((OX_DEVICE[i + 1] + OX_DEVICE[i]) * 0.5, (OY_DEVICE[j - 1] + OY_DEVICE[j]) * 0.5),
		gamma((OX_DEVICE[i + 1] + OX_DEVICE[i]) * 0.5, (OY_DEVICE[j + 1] + OY_DEVICE[j]) * 0.5),
		theta((OX_DEVICE[i - 1] + OX_DEVICE[i]) * 0.5, (OY_DEVICE[j + 1] + OY_DEVICE[j]) * 0.5);

	// get prev coordnates
	double u = func_u(C_B, alpha);
	double v = func_v(C_UB, C_BB, C_LB, C_RB, C_TIME, alpha);

	p[0].x = alpha.x - C_TAU * u;
	p[0].y = alpha.y - C_TAU * v;
	p[0].x_initial = alpha.x;
	p[0].y_initial = alpha.y;

	v = func_v(C_UB, C_BB, C_LB, C_RB, C_TIME, beta);
	u = func_u(C_B, beta);
	p[1].x = beta.x - C_TAU * u;
	p[1].y = beta.y - C_TAU * v;
	p[1].x_initial = beta.x;
	p[1].y_initial = beta.y;
	v = func_v(C_UB, C_BB, C_LB, C_RB, C_TIME, gamma);
	u = func_u(C_B, gamma);
	p[2].x = gamma.x - C_TAU * u;
	p[2].y = gamma.y - C_TAU * v;
	p[2].x_initial = gamma.x;
	p[2].y_initial = gamma.y;
	v = func_v(C_UB, C_BB, C_LB, C_RB, C_TIME, theta);
	u = func_u(C_B, theta);
	p[3].x = theta.x - C_TAU * u;
	p[3].y = theta.y - C_TAU * v;
	p[3].x_initial = theta.x;
	p[3].y_initial = theta.y;

	c_dp_t intersection = get_intersection_point(p[0], p[1], p[2], p[3]);
	if ((p[1].y - intersection.y) * (p[3].y - intersection.y) > 0) return pseudo; // ??
	if ((p[0].x - intersection.x) * (p[2].x - intersection.x) > 0) return pseudo; // ??	
	if (is_points_belong_to_one_line(p[0], p[1], p[3])) return pseudo;
	
	a = p[0];
	b = p[1];
	c = p[2];
	k = p[0];
	m = p[3];
	n = p[2];
return normal;
	//return get_wall_intersection_type(p);
}

// __pure static double integrate_wall_triangle(const c_dp_t wp, // wall point
//                                       double ly, // left y coordinate
//                                       double ry) // right y coordinate
// {
// 	return 0;
// }

// __pure static double integrate_wall_rectangle(const c_dp_t wp1, const c_dp_t wp2, const c_dp_t wp3, const c_dp_t wp4, double wp1y, double wp2y)
// {
// 	return 0;
// }

// __pure static double integrate_wall_rectangle(const c_dp_t wp1, const c_dp_t wp2, double wp1y, double wp2y)
// {
// 	return 0;
// }

// __pure static double integrate_wall_pentagon(const c_dp_t wp1, const c_dp_t wp2, const c_dp_t wp3, double y1, double y2)
// {
// 	return 0;
// }

// __pure static double integrate_pentagon(const c_dp_t x, const c_dp_t y, const c_dp_t z, double ly, double ry)
// {
// 	return 0;
// }

__pure static double integrate(double* prev_dens, int i, int j)
{
	c_dp_t a1, b1, c1, a2, b2, c2;
	c_dp4_t* p = new c_dp4_t[6];
	quad_type type = get_quadrangle_type(i, j, a1, b1, c1, a2, b2, c2, p);


	switch (type)
	{
	case wall_1_middle_in: // вообщем это один и тот же способ
	case wall_1_middle_out:
	case wall_1_middle_at:
		{
//			//тут получается всегда 3 треугольника
//			double result = 0;
//			double t = 0;			
//			c_dp_t v1 = c_dp_t(p[4].x, p[4].y);
//			c_dp_t v2 = c_dp_t(p[2].x, p[2].y);
//			c_dp_t v3 = c_dp_t(p[1].x, p[1].y);
//			sort_by_y_asc(v1, v2, v3);
//			t = integrate_uniform_triangle(v1, v2, v3);
//			result += t;
//										
//			v1 = c_dp_t(p[4].x, p[4].y);
//			v2 = c_dp_t(p[2].x, p[2].y);
//			v3 = c_dp_t(p[5].x, p[5].y);
//			sort_by_y_asc(v1, v2, v3);
//			t = integrate_uniform_triangle(v1, v2, v3); 
//			result += t;
//										
//			v1 = c_dp_t(p[3].x, p[3].y);
//			v2 = c_dp_t(p[2].x, p[2].y);
//			v3 = c_dp_t(p[5].x, p[5].y);
//			sort_by_y_asc(v1, v2, v3);
//			t = integrate_uniform_triangle(v1, v2, v3);
//			result += t;
//					
//			v1 = c_dp_t(p[0].x, p[0].y);
//			v2 = c_dp_t(p[4].x, p[4].y);
//			v3 = c_dp_t(p[5].x, p[5].y);
//			t = integrate_uniform_triangle_wall(v1, v2, v3, type);
//			result += t;
			//return result;
			//break;
		}
	case wall_2:
	{
		//double t = 0;
		//double result = 0;
			//// надо рассмотреть три случая
			//// 1. p2 внутри треугольника (p4,p5,p3) = итегрирование по 3 треугольникам
			//if (is_point_in_triangle(p[2], p[4], p[5], p[3]))
			//{
			//	c_dp_t v1 = p[5];
			//	c_dp_t v2 = p[3];
			//	c_dp_t v3 = p[4];
			//	//сразу их располагаем в порядке возростания y координаты
			//	t = integrate_uniform_triangle(v1, v2, v3);
			//	result += t;
			//				
			//	v1 = p[5];
			//	v2 = p[3];
			//	v3 = p[2];
			//	sort_by_y_asc(v1, v2, v3);
			//	t = integrate_uniform_triangle(v1, v2, v3);
			//	result += t;
			//				
			//	v1 = p[3];
			//	v2 = p[2];
			//	v3 = p[4];
			//	sort_by_y_asc(v1, v2, v3);
			//	t = integrate_uniform_triangle(v1, v2, v3);
			//	result += t;
			//}
			//// 2. p3 внутри треугольника (p4,p5,p2) =  итегрирование по 3 треугольникам
			//else if (is_point_in_triangle(p[3], p[4], p[5], p[2]))
			//{
			//	//сразу их располагаем в порядке возростания y координаты
			//	c_dp_t v1 = p[5];
			//	c_dp_t v2 = p[2];
			//	c_dp_t v3 = p[4];			
			//	t = integrate_uniform_triangle(v1, v2, v3);
			//	result += t;

			//	v1 = p[5];
			//	v2 = p[2];
			//	v3 = p[3];
			//	sort_by_y_asc(v1, v2, v3);
			//	t = integrate_uniform_triangle(v1, v2, v3);
			//	result += t;

			//	v1 = p[4];
			//	v2 = p[2];
			//	v3 = p[3];
			//	sort_by_y_asc(v1, v2, v3);
			//	t = integrate_uniform_triangle(v1, v2, v3);
			//	result += t;
			//}
			//// 3. ни 1 ни 2 условие - нормальный случай 4 угольник, интегрируем как обычно
			//else
			//{
			//	c_dp_t v1 = p[4];
			//	c_dp_t v2 = p[5];
			//	c_dp_t v3 = p[3];
			//	sort_by_y_asc(v1, v2, v3);
			//	t = integrate_uniform_triangle(v1, v2, v3);
			//	result += t;
			//	v1 = p[4];
			//	v2 = p[2];
			//	v3 = p[3];
			//	sort_by_y_asc(v1, v2, v3);
			//	t = integrate_uniform_triangle(v1, v2, v3);
			//	result += t;
			//}
//		if (type == wall_2){
//			c_dp_t v1 = c_dp_t(p[0].x, p[0].y);
//			c_dp_t v2 = c_dp_t(p[1].x, p[1].y);
//			c_dp_t v3 = c_dp_t(p[5].x, p[5].y);
//			sort_by_y_asc(v1, v2, v3);
//			t = integrate_uniform_triangle_wall(v1, v2, v3, type);
//			result += t;
//			v1 = c_dp_t(p[4].x, p[4].y);
//			v2 = c_dp_t(p[1].x, p[1].y);
//			v3 = c_dp_t(p[5].x, p[5].y);
//			sort_by_y_asc(v1, v2, v3);
//			t = integrate_uniform_triangle_wall(v1, v2, v3, type);
//			result += t;
//		}
			//return 0;
			//return result;
		}
	case wall_3_middle_in:
	case wall_3_middle_out:
	case wall_3_middle_at:
		{
			/*double result = 0;
		double t = 0;
		c_dp_t v1 = p[4];
		c_dp_t v2 = p[5];
		c_dp_t v3 = p[3];
		sort_by_y_asc(v1, v2, v3);
		t = integrate_uniform_triangle(a1, b1, c1);
		result += t;
		return result;*/
		}
	case wall_4:
		//return 0;
	case normal:
		{
			double result = 0;
			double t = 0;
			sort_by_y_asc(a1, b1, c1);
			if (i==1 && j == 1)
			{
				flag = 1;
			}
			else
			{
				flag = 0;
			}
			t = integrate_uniform_triangle(prev_dens, a1, b1, c1);
			// if (i==1 && j == 1)
			// {
			// 	printf("%s\n", "integrate");
			// 	printf("%lf\n", t);
			// }
			result += t;
			sort_by_y_asc(a2, b2, c2);
			t = integrate_uniform_triangle(prev_dens, a2, b2, c2);
			result += t;
			return result;
		}
	case concave:
	case convex:
	case pseudo:
		return -1;
	}
	return 0;
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
	delete [] OX;
	delete [] OY;
}

__global__ void kernel(double* prev_result, double* result)
{
	if (blockIdx.x * blockDim.x + threadIdx.x == 0)
	{
		printf("\nKERNEL PARAMS\n");
		printf("b = %f\n", C_B);
		printf("lbDom = %f\n", C_LB);
	 	printf("rbDom = %f\n", C_RB);
	 	printf("bbDom = %f\n", C_BB);
		printf("ubDom = %f\n", C_UB);
	 	printf("tau = %f\n", C_TAU);
	 	printf("ox length = %d\n", C_OX_LEN + 1);
	 	printf("oy length = %d\n", C_OX_LEN + 1);
	 	printf("PREV_TIME = %lf\n", C_PREV_TIME);
	 	printf("TIME = %lf\n", C_TIME);
	 	printf("C_XY_LEN = %d\n", C_XY_LEN);
	 	printf("%s\n", "OX DEVICE");
	 	for(int i = 0; i < C_OX_LEN + 1; i++)
	 	{
	 		printf("%f ", OX_DEVICE[i]);
	 	}
	 	printf("%s\n", "");
	 	printf("%s\n", "OY DEVICE");
	 	for(int i = 0; i < C_OY_LEN + 1; i++)
	 	{
	 		printf("%f ", OY_DEVICE[i]);
	 	}
	 	printf("%s\n", "");
	 	printf("%s\n", "PREV DENSITY");
	 	for(int i = 0; i < C_OX_LEN + 1; i++)
	 	{
	 		for(int j = 0; j < C_OY_LEN + 1; j++)
		 	{
		 		printf("%f ", prev_result[C_OX_LEN_1 * j + i]);
		 	}
		 	printf("%s\n", "");
	 	}
	 	
	}
	
	 for (int opt = blockIdx.x * blockDim.x + threadIdx.x; opt < C_XY_LEN; opt += blockDim.x * gridDim.x)
	 {		
	 	int i = opt % (C_OX_LEN + 1);
	 	int j = opt / (C_OY_LEN + 1);

	 	// расчет границы
	 	if (j == 0)  // bottom bound
	 	{
	 		result[ opt ]  = 1.1  +  sin( C_TIME * C_HX * j * C_BB );
	 	}
		else if (i == 0) // left bound
		{
			result[ opt ] = 1.1  +  sin( C_TIME * C_HX* i * C_LB );
		}
		else if (j == C_OY_LEN) // upper bound
		{ 
			result[ opt ] = 1.1  +  sin( C_TIME * C_HX * i * C_UB );
		}
		else if (i == C_OX_LEN) // right bound
		{ 
			result[ opt ] = 1.1  +  sin(  C_TIME * C_HX * j * C_RB );
		}
		else if (i > 0 && j > 0 && j != C_OY_LEN && i != C_OX_LEN)
		{                   
			double t = integrate(prev_result, i, j);	
			result[ opt ] =  t * C_INVERTED_HX_HY;
			if (opt == 13) 
			{
				printf("%s\n", "result = 13");
				printf("%lf\n", t);
				printf("%lf\n", result[ opt ]);
			}
			result[ opt ] += C_TAU * func_f(C_B, C_TIME, C_UB, C_BB, C_LB, C_RB, OX_DEVICE[i], OY_DEVICE[j]);
			if (opt == 13) 
			{
				printf("%s\n", " F result = 13");
				printf("%lf\n", C_TAU);
				printf("%lf\n", result[ opt ] );
				printf("%lf\n", func_f(C_B, C_TIME, C_UB, C_BB, C_LB, C_RB, OX_DEVICE[i], OY_DEVICE[j]) );
			}
		}
	 }
}

float solve_cuda(double* density)
{
//	const int gridSize = 256;
//	const int blockSize =  512; 
	const int gridSize = 1;
	const int blockSize =  1;
	double *result = NULL, *prev_result = NULL, *ox = NULL, *oy=NULL;
	int size = sizeof(double)*XY_LEN;
	double *prev_result_h = new double[XY_LEN];
	for (int j = 0; j < OY_LEN + 1; j++)
	{
		for (int i = 0; i < OX_LEN_1; i++)
		{
			prev_result_h[OX_LEN_1 * j + i] = analytical_solution(0, OX[i], OY[j]);
		}
	}

	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
    checkCuda(cudaMemcpyToSymbol(C_TAU, &TAU, sizeof(double)));	
	checkCuda(cudaMemcpyToSymbol(C_B, &B, sizeof(double)));
	checkCuda(cudaMemcpyToSymbol(C_LB, &LB, sizeof(double)));
	checkCuda(cudaMemcpyToSymbol(C_RB, &RB, sizeof(double)));
	checkCuda(cudaMemcpyToSymbol(C_BB, &BB, sizeof(double)));
	checkCuda(cudaMemcpyToSymbol(C_UB, &UB, sizeof(double)));
	checkCuda(cudaMemcpyToSymbol(C_INVERTED_HX_HY, &INVERTED_HX_HY, sizeof(double)));	
	checkCuda(cudaMemcpyToSymbol(C_HX, &HX, sizeof(int)));	
	checkCuda(cudaMemcpyToSymbol(C_HY, &HY, sizeof(int)));	
	checkCuda(cudaMemcpyToSymbol(C_OX_LEN_1, &OX_LEN_1, sizeof(int)));
	checkCuda(cudaMemcpyToSymbol(C_XY_LEN, &XY_LEN, sizeof(int)));
	checkCuda(cudaMemcpyToSymbol(C_OX_LEN, &OX_LEN, sizeof(int)));
	checkCuda(cudaMemcpyToSymbol(C_OY_LEN, &OY_LEN, sizeof(int)));
	
	checkCuda(cudaMalloc((void**)&(result), size) );
	checkCuda(cudaMemset(result, 0, size) );
	checkCuda(cudaMalloc((void**)&(prev_result), size) );
	checkCuda(cudaMalloc((void**)&(ox), sizeof(ox)*(OX_LEN+1)));
	checkCuda(cudaMalloc((void**)&(oy), sizeof(oy)*(OY_LEN+1)));
	checkCuda(cudaMemcpy(ox, OX, sizeof(ox)*(OX_LEN + 1), cudaMemcpyHostToDevice));	
	checkCuda(cudaMemcpy(oy, OY, sizeof(oy)*(OY_LEN + 1), cudaMemcpyHostToDevice));	
	checkCuda(cudaMemcpy(prev_result, prev_result_h, size, cudaMemcpyHostToDevice));	
	checkCuda(cudaMemcpyToSymbol(OX_DEVICE, &ox, sizeof(ox)));
	checkCuda(cudaMemcpyToSymbol(OY_DEVICE, &oy, sizeof(oy)));	

	cudaEventRecord(start, 0);   

	TIME = 0;
	int tl = 0;
	TIME_STEP_CNT = 2;
	int tempTl  = TIME_STEP_CNT - 1;
	while(tl < tempTl)
	{
		checkCuda(cudaMemcpyToSymbol(C_PREV_TIME, &TIME, sizeof(double)));
     	TIME = TAU * (tl+1);
	    checkCuda(cudaMemcpyToSymbol(C_TIME, &TIME, sizeof(double)));	
	    kernel<<<gridSize, blockSize>>>(prev_result, result);
		// checkCuda(cudaMemcpyToSymbol(C_PREV_TIME, &TIME, sizeof(double)));
  		// TIME = TAU * (tl+2);
		// checkCuda(cudaMemcpyToSymbol(C_TIME, &TIME, sizeof(double)));	
		// kernel<<<gridSize, blockSize>>>(result, prev_result); 
	    tl += 2;            
	}

	printf("%s\n", "");
	//checkCuda(cudaMemcpy(density, prev_result, size, cudaMemcpyDeviceToHost));	
	checkCuda(cudaMemcpy(density, result, size, cudaMemcpyDeviceToHost));	
	for (int j = 0; j < OX_LEN_1; j++)
	{
		for (int i = 0; i < OY_LEN + 1; i++)
		{
			printf("%lf ", density[OX_LEN_1 * i + j]);			
		}
		printf("%s\n", "");
	}
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	printf("Computation Time %f\n", time);
	cudaFree(result);
	cudaFree(prev_result);
	cudaFree(ox);
	cudaFree(oy);
	cudaDeviceReset();
	delete[] prev_result_h;
	return time;
}

inline void print_matrix11(double* a, int n, int m, int precision = 8) {
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

double* compute_density_cuda_internal(double b, double lb, double rb, double bb, double ub,
                        double tau, int time_step_count, int ox_length, int oy_length, double& norm, float& time)
{
#ifdef __NVCC__	
    init(b, lb, rb, bb, ub, tau, time_step_count, ox_length, oy_length);
	double* density = new double[XY_LEN];
	print_params(B, LB, RB, BB, UB, TAU, TIME_STEP_CNT, OX_LEN, OY_LEN);
	time = solve_cuda(density);
	//print_matrix11(density, 11, 11);
	norm = get_norm_of_error(density, TIME_STEP_CNT * TAU);
	clean();
	return density;
#else
        return NULL;
#endif
}
