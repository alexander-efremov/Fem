#include "common.h"
#include "point.h"
#include "utils.h"
#include <algorithm>

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
static int TL; //-V707
static double TIME;
static double PREV_TIME; // tau * (tl - 1)
static double INVERTED_HX_HY;

__pure inline static void sort_by_y_asc(dp_t& x, dp_t& y, dp_t& z)
{
	using std::swap; // TODO: переписать без swap (использовать std::move)
	if (x.y < y.y)
	{
		if (z.y < x.y) swap(x, z);
	}
	else
	{
		if (y.y < z.y) swap(x, y);
		else swap(x, z);
	}
	if (z.y < y.y) swap(y, z);
}

__pure inline static void sort_by_x(dp_t& x, dp_t& y, dp_t& z)
{
	double t;
	if (x.x < y.x)
	{
		if (z.x < x.x)
		{
			t = x.x;
			x.x = z.x;
			z.x = t;
			t = y.y;
			y.y = z.y;
			z.y = t;
			//swap(x, z);
		}
	}
	else
	{
		if (y.x < z.x)
		{
			t = x.x;
			x.x = y.x;
			y.x = t;
			t = x.y;
			x.y = y.y;
			y.y = t;
			//swap(x, y);
		}
		else
		{
			t = x.x;
			x.x = z.x;
			z.x = t;
			t = x.y;
			x.y = z.y;
			z.y = t;
			//swap(x, z);
		}
	}
	if (z.x < y.x)
	{
		t = y.x;
		y.x = z.x;
		z.x = t;
		t = y.y;
		y.y = z.y;
		z.y = t;
		//swap(y, z);
	}
}

//__pure inline void sort_by_y(dp_t& w, dp_t& x, dp_t& y, dp_t& z)
//{
//	if (w.y > x.y)  { double t = w.y; w.y = x.y; x.y = t; t = w.x; w.x = x.x; x.x = t; }
//	if (w.y > y.y)  { double t = w.y; w.y = y.y; y.y = t; t = w.x; w.x = y.x; y.x = t; }
//	if (w.y > z.y)  { double t = w.y; w.y = z.y; z.y = t; t = w.x; w.x = z.x; z.x = t; }
//	sort_by_y(x, y, z);
//}
//
//__pure inline void sort_by_x_asc(dp_t& w, dp_t& x, dp_t& y, dp_t& z)
//{
//	if (w.x > x.x)  { double t = w.y; w.y = x.y; x.y = t; t = w.x; w.x = x.x; x.x = t; }
//	if (w.x > y.x)  { double t = w.y; w.y = y.y; y.y = t; t = w.x; w.x = y.x; y.x = t; }
//	if (w.x > z.x)  { double t = w.y; w.y = z.y; z.y = t; t = w.x; w.x = z.x; z.x = t; }
//	sort_by_x_asc(x, y, z);
//}

inline void sort_by_y(dp_t* a)
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

inline void sort_by_y_desc_3(dp_t* a)
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

__pure inline void sort_by_x_asc(dp_t* a)
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
}

// получается порядок
/*

a[1]    a[2]
a[0]   a[3]
*/
__pure inline void sort_by_xy_wall_2(dp_t* a)
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

__pure inline static bool try_get_slope_ratio(const dp_t& bv, const dp_t& uv, double& value)
{
	if (fabs(bv.x - uv.x) < MIN_VALUE)
	{
		return false;
	}
	value = fabs((uv.y - bv.y) / (uv.x - bv.x)); // угловой коэффициент прямой
	if (value < MIN_VALUE)
	{
		return false;
	}
	return true;
}

__pure inline static dp_t get_intersection_point(const dp_t& alpha, const dp_t& beta, const dp_t& gamma, const dp_t& theta)
{
	double a1 = gamma.y - alpha.y;
	double b1 = alpha.x - gamma.x; //double b1 = -(gamma.x - alpha.x);
	double c1 = a1 * alpha.x + b1 * alpha.y;
	double a2 = theta.y - beta.y;
	double b2 = beta.x - theta.x; //double b2 = -(theta.x - beta.x);
	double c2 = a2 * beta.x + b2 * beta.y;
	return dp_t((b1 * c2 - b2 * c1) / (b1 * a2 - b2 * a1), (a1 * c2 - a2 * c1) / (-b1 * a2 + b2 * a1));
}

__pure inline static double sign(const dp_t& p1, const dp_t p2, const dp_t p3)
{
	return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

__pure inline static bool is_points_belong_to_one_line(const dp_t& p1, const dp_t p2, const dp_t p3)
{
	return sign(p1, p2, p3) == FLT_MIN ;
}

__pure inline static bool is_point_in_triangle(dp_t pt, dp_t v1, dp_t v2, dp_t v3)
{
	bool b1, b2, b3;
	b1 = sign(pt, v1, v2) < 0.0;
	b2 = sign(pt, v2, v3) < 0.0;
	b3 = sign(pt, v3, v1) < 0.0;
	return b1 == b2 && b2 == b3;
}

__pure inline static double analytical_solution(double t, double x, double y)
{
	return 1.1 + sin(t * x * y);
}

__pure inline static double func_u(double b, double x, double y)
{
	return b * y * (1 - y) * (M_PI_2 + atan(-x));
}

__pure inline static double func_u(double b, const dp_t& p)
{
	return func_u(b, p.x, p.y);
}

__pure inline static double func_v(double ub, double bb, double lb, double rb, double t, double x, double y)
{
	return atan(0.1 * (x - lb) * (x - rb) * (1 + t) * (y - ub) * (y - bb));
}

__pure inline static double func_v(double ub, double bb, double lb, double rb, double t, const dp_t& p)
{
	return func_v(ub, bb, lb, rb, t, p.x, p.y);
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

static double integrate_rectangle_one_cell(double py, double qy, double gx, double hx, const ip_t& sx, const ip_t& sy)
{
	double result, a, b;
	a = sx.y >= 0 && sy.y >= 0 ? OX[sx.y] : HX * sx.y;
	b = sx.y >= 0 && sy.y >= 0 ? OY[sy.y] : HY * sy.y; // ЭТО ПЛОТНОСТЬ С ПРЕДЫДУЩЕГО СЛОЯ ДЛЯ ДАННОЙ ЯЧЕЙКИ
	result = integrate_rectangle(py, qy, gx, hx, a, b) * (sx.x >= 0 && sy.x >= 0 ? PREV_DENSITY[OX_LEN_1 * sy.x + sx.x] : analytical_solution(PREV_TIME, sx.x * HX, sy.x * HY));
	a = sx.x >= 0 && sy.y >= 0 ? OX[sx.x] : HX * sx.x;
	b = sx.x >= 0 && sy.y >= 0 ? OY[sy.y] : HY * sy.y;
	result -= integrate_rectangle(py, qy, gx, hx, a, b) * (sx.x >= 0 && sy.x >= 0 ? PREV_DENSITY[OX_LEN_1 * sy.x + sx.y] : analytical_solution(PREV_TIME, sx.y * HX, sy.x * HY));
	a = sx.y >= 0 && sy.x >= 0 ? OX[sx.y] : HX * sx.y;
	b = sx.y >= 0 && sy.x >= 0 ? OY[sy.x] : HY * sy.x;
	result -= integrate_rectangle(py, qy, gx, hx, a, b) * (sx.x >= 0 && sy.x >= 0 ? PREV_DENSITY[OX_LEN_1 * sy.y + sx.x] : analytical_solution(PREV_TIME, sx.x * HX, sy.y * HY));
	a = sx.x >= 0 && sy.x >= 0 ? OX[sx.x] : HX * sx.x;
	b = sx.x >= 0 && sy.x >= 0 ? OY[sy.x] : HY * sy.x;
	result += integrate_rectangle(py, qy, gx, hx, a, b) * (sx.x >= 0 && sy.x >= 0 ? PREV_DENSITY[OX_LEN_1 * sy.y + sx.y] : analytical_solution(PREV_TIME, sx.y * HX, sy.y * HY));
	return result * INVERTED_HX_HY;
}

static double integrate_triangle_left_one_cell(const dp_t& bv, const dp_t& uv, double hx,
                                               const ip_t& sx, const ip_t& sy)
{
	double a_sl = (bv.x - uv.x) / (bv.y - uv.y); //   Coefficients of slant line: x = a_SL *y  +  b_SL.
	if (fabs(a_sl) <= FLT_MIN) return 0;
	double b_sl = uv.x - a_sl * uv.y;
	double result = 0, tmp, alpha, beta;
	alpha = sx.y >= 0 && sy.y >= 0 ? OY[sy.y] : HY * sy.y;
	beta = sx.y >= 0 && sy.y >= 0 ? OX[sx.y] : HX * sx.y;
	tmp = 0.25 * (sqr(uv.y - OY[sy.y]) - sqr(bv.y - OY[sy.y])) * sqr(hx - beta) - integrate_triangle(bv.y, uv.y, alpha, beta, a_sl, b_sl);
	result += tmp * (sx.x >= 0 && sx.y <= OX_LEN && sy.x >= 0 && sy.y <= OY_LEN ? PREV_DENSITY[OX_LEN_1 * sy.x + sx.x] : analytical_solution(PREV_TIME, sx.x * HX, sy.x * HY));
	beta = sx.x >= 0 && sy.y >= 0 ? OX[sx.x] : HX * sx.x;
	tmp = sqr(uv.y - OY[sy.y]) - sqr(bv.y - OY[sy.y]);
	tmp = -0.25 * tmp * sqr(hx - beta) + integrate_triangle(bv.y, uv.y, alpha, beta, a_sl, b_sl);
	result += tmp * (sx.x >= 0 && sx.y <= OX_LEN && sy.x >= 0 && sy.y <= OY_LEN ? PREV_DENSITY[OX_LEN_1 * sy.x + sx.y] : analytical_solution(PREV_TIME, sx.y * HX, sy.x * HY));
	alpha = sx.y >= 0 && sy.x >= 0 ? OY[sy.x] : HY * sy.x;
	beta = sx.y >= 0 && sy.x >= 0 ? OX[sx.y] : HX * sx.y;
	tmp = sqr(uv.y - OY[sy.x]) - sqr(bv.y - OY[sy.x]);
	tmp = -0.25 * tmp * sqr(hx - beta) + integrate_triangle(bv.y, uv.y, alpha, beta, a_sl, b_sl);
	result += tmp * (sx.x >= 0 && sx.y <= OX_LEN && sy.x >= 0 && sy.y <= OY_LEN ? PREV_DENSITY[OX_LEN_1 * sy.y + sx.x] : analytical_solution(PREV_TIME, sx.x * HX, sy.y * HY));
	alpha = sx.x >= 0 && sy.x >= 0 ? OY[sy.x] : HY * sy.x;
	beta = sx.x >= 0 && sy.x >= 0 ? OX[sx.x] : HX * sx.x;
	tmp = sqr(uv.y - OY[sy.x]) - sqr(bv.y - OY[sy.x]);
	tmp = 0.25 * tmp * sqr(hx - beta) - integrate_triangle(bv.y, uv.y, alpha, beta, a_sl, b_sl);
	result += tmp * (sx.x >= 0 && sx.y <= OX_LEN && sy.x >= 0 && sy.y <= OY_LEN ? PREV_DENSITY[OX_LEN_1 * sy.y + sx.y] : analytical_solution(PREV_TIME, sx.y * HX, sy.y * HY));
	return result * INVERTED_HX_HY;
}

static double integrate_right_slant_chanel(const dp_t& bv, const dp_t& uv, bool is_rect_truncated, const ip_t& sx, double b, const ip_t& sb, const ip_t& sy)
{
	if (fabs(uv.y - bv.y) <= FLT_MIN) return FLT_MIN ;
	double result = 0, gx = 0;
	double x = uv.x <= bv.x ? uv.x : bv.x;

	//   A. Under rectangle.
	result += -1 * integrate_triangle_left_one_cell(bv, uv, x, sx, sy);

	// case B: неполный прямоугольник    
	if (is_rect_truncated)
	{
		if (sx.x == sb.x) gx = b;
		if (sx.x > sb.x)
		{
			gx = sx.x >= 0 ? OX[sx.x] : HX * sx.x;
		}
		result += integrate_rectangle_one_cell(bv.y, uv.y, gx, x, sx, sy);
	}

	//   А теперь прибавим все прямоугольные куски, которые помещаются в ячейку
	ip_t ch_pos(sb.x, sb.x + 1);
	for (int j = sb.x; j < sx.x; j++)
	{
		if (j == sb.x) gx = b;
		else gx = ch_pos.x >= 0 ? OX[ch_pos.x] : HX * ch_pos.x;
		result += integrate_rectangle_one_cell(bv.y, uv.y, gx, ch_pos.x >= 0 ? OX[ch_pos.y] : HX * ch_pos.y, ch_pos, sy);
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

static double integrate_left_slant_chanel(const dp_t& bv, const dp_t& uv,
                                          bool is_rect_trunc, const ip_t& sx, const ip_t& sy,
                                          double b, const ip_t& sb)
{
	if (fabs(uv.y - bv.y) <= FLT_MIN) return FLT_MIN ;
	double result = 0, hx = 0; //   -  Left and right boundary for each integration.   
	double x = uv.x <= bv.x ? bv.x : uv.x;

	// case A: triangle
	result += integrate_triangle_left_one_cell(bv, uv, x, sx, sy);

	// case B: не полный прямоугольник
	if (is_rect_trunc)
	{ // это значит, что прямоугольник занимает не всю ячейку  
		hx = sx.x == sb.x ? b : (sx.y >= 0 ? OX[sx.y] : HX * sx.y);
		result += integrate_rectangle_one_cell(bv.y, uv.y, x, hx, sx, sy);
	}

	//   А теперь прибавим все прямоугольные куски, которые помещаются в ячейку
	ip_t ch_pos(sx.x + 1, sx.x + 2); //   - координаты канала
	for (int j = sx.x + 1; j < sb.x + 1; j++)
	{
		hx = ch_pos.y <= 0 ? HX * ch_pos.y : hx = OX[ch_pos.y];
		if (j == sb.x) hx = b;
		result += integrate_rectangle_one_cell(bv.y, uv.y, ch_pos.y <= 0 ? HX * ch_pos.x : OX[ch_pos.x], hx, ch_pos, sy);
		ch_pos.x += 1;
		ch_pos.y = ch_pos.x + 1;
	}
	return result;
}

// определим целочисленные индексы квадратов в которых лежат верхняя и нижняя точки треугольника
// sx = (x,y) координаты квадрата в которой лежит нижняя точка
// sy = (x,y) координаты квадрата в которой лежит верхняя точка
// в случае успешной проверки, k = будет  угловой коэфициент прямой

static double integrate_right_triangle_bottom_left(const dp_t& bv, const dp_t& uv)
{
	double k = 0;
	if (!try_get_slope_ratio(bv, uv, k)) return k;

	//   -  Index of current square by Ox and Oy axes. 
	ip_t sx, sy;
	sx.x = static_cast<int>((bv.x - FLT_MIN) / HX);
	if (bv.x - FLT_MIN <= 0) sx.x -= 1;
	sx.y = sx.x + 1;
	sy.x = static_cast<int>((bv.y + FLT_MIN) / HY);
	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
	sy.y = sy.x + 1;

	ip_t ib(sx.x, sx.x + 1);
	double result = 0;
	int curr_i = 0, next_i;
	dp_t curr = bv, next;
	while (true)
	{
		//TODO: sx.x и sx.y должны быть положительными всегда? Кажется для sx.x это всегда верно...
		double slope = sx.y >= 0 ? OY[sy.y] - curr.y : fabs(HY * sy.y - curr.y);
		slope /= sx.x >= 0 ? curr.x - OX[sx.x] : fabs(curr.x - HX * sx.x);
		if (slope <= k)
		{
			next_i = 1;
			next.y = sy.y >= 0 ? OY[sy.y] : HY * sy.y;
			next.x = curr.x - (next.y - curr.y) / k;
		}
		else
		{
			next_i = 2;
			next.x = sx.x >= 0 ? OX[sx.x] : HX * sx.x;
			next.y = curr.y - k * (next.x - curr.x);
		}
		if (next.x - uv.x < FLT_MIN)
		{
			// сюда попадаем и в случае когда треугольник полностью в одной ячейке лежит
			// и в случае когда прошлись по всем точкам...
			result += integrate_left_slant_chanel(curr, uv, (uv.x <= curr.x ? curr_i : 0) == 1, sx, sy, bv.x, ib);
			break;
		}
		result += integrate_left_slant_chanel(curr, next, (next.x <= curr.x ? curr_i : next_i) == 1, sx, sy, bv.x, ib);
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

static double integrate_right_triangle_bottom_right(const dp_t& bv, const dp_t& uv)
{
	double k = 0;
	if (!try_get_slope_ratio(bv, uv, k)) return k;

	ip_t sx, sy;
	sx.x = static_cast<int>((bv.x + FLT_MIN) / HX);
	if (bv.x + FLT_MIN <= 0) sx.x -= 1;
	sx.y = sx.x + 1;
	sy.x = static_cast<int>((bv.y + FLT_MIN) / HY);
	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
	sy.y = sy.x + 1;

	ip_t ib(sx.x, sx.x + 1);
	double result = 0;
	int curr_i = 0, next_i;
	dp_t curr = bv, next;
	while (true)
	{
		double slope = sy.y >= 0 ? fabs(OY[sy.y] - curr.y) : fabs(HY * sy.y - curr.y);
		slope /= sx.y >= 0 ? fabs(OX[sx.y] - curr.x) : fabs(HX * sx.y - curr.x);
		if (slope <= k)
		{
			next_i = 1;
			next.y = sy.y >= 0 ? OY[sy.y] : HY * sy.y;
			next.x = bv.x + (next.y - bv.y) / k;
		}
		else
		{
			next_i = 2;
			next.x = sx.y >= 0 ? OX[sx.y] : HX * sx.y;
			next.y = bv.y + k * (next.x - bv.x);
		}
		if (next.x - uv.x > FLT_MIN)
		{
			result += integrate_right_slant_chanel(curr, uv, (uv.x <= curr.x ? 0 : curr_i) == 1, sx, bv.x, ib, sy);
			break;
		}
		result += integrate_right_slant_chanel(curr, next, (next.x <= curr.x ? next_i : curr_i) == 1, sx, bv.x, ib, sy);
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

static double integrate_right_triangle_upper_left(const dp_t& bv, const dp_t& uv)
{
	double k = 0;
	if (!try_get_slope_ratio(bv, uv, k)) return k;

	ip_t sx, sy, ib;
	sx.x = static_cast<int>((bv.x + FLT_MIN) / HX); //   -  If bv.x is in grid edge I want it will be in the right side.
	if (bv.x + FLT_MIN <= 0) sx.x -= 1;
	sx.y = sx.x + 1;
	sy.x = static_cast<int>((bv.y + FLT_MIN) / HY); //   -  If bv.y is in grid edge I want it will be in the upper square.
	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
	sy.y = sy.x + 1;
	ib.x = static_cast<int>((uv.x - FLT_MIN) / HY); //   -  If uv.x is in grid edge I want it will be in the left side.
	if (uv.x - FLT_MIN <= 0) ib.x -= 1;
	ib.y = ib.x + 1;

	double result = 0;
	int curr_i = 0, next_i;
	dp_t curr = bv, next;
	while (true)
	{
		double slope = sy.y >= 0 ? OY[sy.y] - curr.y : fabs(HY * sy.y - curr.y);
		slope /= sx.y >= 0 ? OX[sx.y] - curr.x : fabs(HX * sx.y - curr.x);
		if (slope <= k)
		{
			next_i = 1;
			next.y = sy.y >= 0 ? OY[sy.y] : HY * sy.y;
			next.x = bv.x + (next.y - bv.y) / k;
		}
		else
		{
			next_i = 2;
			next.x = sx.y >= 0 ? OX[sx.y] : HX * sx.y;
			next.y = bv.y + k * (next.x - bv.x);
		}
		if (next.x - uv.x > FLT_MIN) // если следующая точка уже правее, чем наша граничная точка, то мы обработали канал
		{
			result += integrate_left_slant_chanel(curr, uv, (uv.x <= curr.x ? curr_i : 0) == 1, sx, sy, uv.x, ib);
			break;
		}
		result += integrate_left_slant_chanel(curr, next, (next.x <= curr.x ? curr_i : next_i) == 1, sx, sy, uv.x, ib);

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

static double integrate_right_triangle_upper_right(const dp_t& bv, const dp_t& uv)
{
	double k = 0;
	if (!try_get_slope_ratio(bv, uv, k)) return k;

	ip_t sx, sy, ib;
	sx.x = static_cast<int>((bv.x - FLT_MIN) / HX); //   -  If bv.x is in grid edge I want it will be between in the left side.
	if (bv.x - FLT_MIN <= 0) sx.x -= 1;
	sx.y = sx.x + 1;
	sy.x = static_cast<int>((bv.y + FLT_MIN) / HY); //   -  If bv.y is in grid edge I want it will be in the upper side.
	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
	sy.y = sy.x + 1;
	ib.x = static_cast<int>((uv.x + FLT_MIN) / HX);
	if (uv.x + FLT_MIN <= 0) ib.x -= 1;
	ib.y = ib.x + 1;

	double result = 0;
	int curr_i = 0, next_i;
	dp_t curr = bv, next;
	while (true)
	{
		double slope = sy.y >= 0 ? fabs(OY[sy.y] - curr.y) : fabs(HY * sy.y - curr.y);
		slope /= sx.x >= 0 ? fabs(curr.x - OX[sx.x]) : fabs(curr.x - HX * sx.x);
		if (slope <= k)
		{
			next_i = 1;
			next.y = sy.y >= 0 ? OY[sy.y] : HY * sy.y;
			next.x = bv.x - (next.y - bv.y) / k;
		}
		else
		{
			next_i = 2;
			next.x = sx.x >= 0 ? OX[sx.x] : HX * sx.x;
			next.y = bv.y - k * (next.x - bv.x);
		}
		if (next.x - uv.x < FLT_MIN)
		{
			result += integrate_right_slant_chanel(curr, uv, (uv.x <= curr.x ? 0 : curr_i) == 1, sx, uv.x, ib, sy);
			break;
		}
		result += integrate_right_slant_chanel(curr, next, (next.x <= curr.x ? next_i : curr_i) == 1, sx, uv.x, ib, sy);
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

static double integrate_bottom_triangle(const dp_t& l, const dp_t& m, const dp_t& r)
{
	double result = 0;
	if (m.x == l.x)
	{
		result = integrate_right_triangle_bottom_right(m, r);
	}
	else if (m.x == r.x)
	{
		result = integrate_right_triangle_bottom_left(m, l);
	}
	else if (m.x < l.x)
	{
		result = integrate_right_triangle_bottom_right(m, r) - integrate_right_triangle_bottom_right(m, l);
	}
	else if (m.x > l.x && m.x < r.x)
	{
		result = integrate_right_triangle_bottom_left(m, l) + integrate_right_triangle_bottom_right(m, r);
	}
	else if (m.x > r.x)
	{
		result = integrate_right_triangle_bottom_left(m, l) - integrate_right_triangle_bottom_left(m, r);
	}
	return result;
}

static double integrate_upper_triangle(const dp_t& l, const dp_t& m, const dp_t& r)
{
	double result = 0;
	if (m.x == l.x)
	{
		result = integrate_right_triangle_upper_right(r, m);
	}
	else if (m.x == r.x)
	{
		result = integrate_right_triangle_upper_left(l, m);
	}
	else if (m.x < l.x)
	{
		result = integrate_right_triangle_upper_right(r, m) - integrate_right_triangle_upper_right(l, m);
	}
	else if (m.x > l.x && m.x < r.x)
	{
		result = integrate_right_triangle_upper_left(l, m) + integrate_right_triangle_upper_right(r, m);
	}
	else if (m.x > r.x)
	{
		result = integrate_right_triangle_upper_left(l, m) - integrate_right_triangle_upper_left(r, m);
	}
	return result;
}

// x,y,z
static double integrate_uniform_triangle(const dp_t& x, const dp_t& y, const dp_t& z)
{
	// точки должны идти в порядке возрастания y координаты, чтобы правильно отработала процедура интегрирования		

	//   a * x  +  b * y  = c.
	double a = z.y - x.y;
	if (fabs(a) < FLT_MIN) return FLT_MIN ;
	double b = x.x - z.x;
	double c = b * x.y + a * x.x;
	dp_t ip((c - b * y.y) / a, y.y);

	//   Возможны 2 случая расположения точки пересечения относительно средней
	//   слева или справа.
	//   есди средняя точка справа от точки пересечения
	//   обменяем местами  X координаты, чтобы использовать один код для расчета
	dp_t t = y;
	if (t.x >= ip.x)
	{
		double tx = t.x;
		t.x = ip.x;
		ip.x = tx;
	}

	return integrate_upper_triangle(t, z, ip) + integrate_bottom_triangle(t, x, ip);
}

static double integrate_uniform_triangle_wall(const dp_t& x, const dp_t& y, const dp_t& z)
{
	return 0;
	//	dp_t cp(0, ); // из y до OY
	//	return integrate_upper_triangle(t, z, ip) + integrate_bottom_triangle(t, x, ip);
}

__pure inline int get_wall_intersection_type_as_int(dp_t* a)
{
	int type = -1;
	bool is_four_point_on_the_wall = a[0].x <= 0 && a[1].x <= 0 && a[2].x <= 0 && a[3].x <= 0;
	bool is_three_point_on_the_wall = a[0].x <= 0 && a[1].x <= 0 && a[2].x <= 0 && a[3].x > 0;
	bool is_two_point_on_the_wall = a[0].x <= 0 && a[1].x <= 0 && a[2].x > 0 && a[3].x > 0;
	bool is_one_point_on_the_wall = a[0].x <= 0 && a[1].x > 0 && a[2].x > 0 && a[3].x > 0;
	if (is_four_point_on_the_wall)
	{
		type = 4;
	}
	else if (is_three_point_on_the_wall)
	{
		type = 3;
	}
	else if (is_two_point_on_the_wall)
	{
		type = 2;
	}
	else if (is_one_point_on_the_wall)
	{
		type = 1;
	}
	else
	{
		type = 0;
	}
	return type;
}

__pure inline static quad_type get_wall_intersection_type(dp_t* a)
{
	/*
	 формулы расчетов здесь http://www.pm298.ru/reshenie/fha0327.php
	 a[0] - alpha
	 a[1] - beta
	 a[2] - gamma
	 a[3] - theta
	 a[4] - mu
	 a[5] - nu

	 */

	int type = get_wall_intersection_type_as_int(a);
	switch (type)
	{
	case 4:
		return wall_4;
	case 3:
		{
			sort_by_x_asc(a);
			sort_by_y_desc_3(a);
			// рассчитаем точку пересечения OY и прямой a[0]:a[3]
			// тут не надо fabs, потому что a[3].x > a[0].x
			double y = a[3].x - a[0].x < FLT_MIN ? 0.5 * (a[0].y + a[3].y) : a[0].y - a[0].x * ((a[3].y - a[0].y) / (a[3].x - a[0].x));
			a[4] = dp_t(0, y); // mu

			// рассчитаем точку пересечения OY и прямой a[2]:a[3]
			// тут не надо fabs, потому что a[3].x > a[2].x
			y = a[3].x - a[2].x < FLT_MIN ? 0.5 * (a[3].y + a[2].y) : a[2].y - a[2].x * ((a[3].y - a[2].y) / (a[3].x - a[2].x));
			a[5] = dp_t(0, y); // nu

			if ((a[0].x - a[2].x) * (a[1].y - a[2].y) - (a[1].x - a[2].x) * (a[0].y - a[2].y) < FLT_MIN)
				return wall_3_middle_at;
			if (a[0].x < a[1].x && a[1].x > a[2].x)
				return wall_3_middle_out;
			if (a[1].x < a[0].x && a[1].x < a[2].x)
				return wall_3_middle_in;
		}
	case 2:
		{
			sort_by_xy_wall_2(a);
			double y = 0;
			if (a[2].x - a[1].x < FLT_MIN)
			{
				y = (a[1].y + a[2].y) * 0.5;
			}
			else
			{
				y = a[1].y - a[1].x * ((a[2].y - a[1].y) / (a[2].x - a[1].x));
			}
			a[4] = dp_t(0, y); // mu

			if (a[3].x - a[0].x < FLT_MIN)
			{
				y = (a[0].y + a[3].y) * 0.5;
			}
			else
			{
				y = a[0].y - a[0].x * ((a[3].y - a[0].y) / (a[3].x - a[0].x));
			}
			a[5] = dp_t(0, y); // nu
			return wall_2;
		}

	case 1: //	wall_1.pdf 		
		{
			sort_by_x_asc(a);
			sort_by_y_desc_3(a);

			// точка a[0] - точка, которая упала на стенку, значит для нее x=t
			// считаем y компоненту

			// рассчитаем точку пересечения OY и прямой a[0]:a[1]
			// тут не надо fabs, потому что a[1].x > a[0].x
			double y = a[1].x - a[0].x < FLT_MIN ? 0.5 * (a[0].y + a[1].y) : (a[0].y - a[0].x * ((a[1].y - a[0].y) / (a[1].x - a[0].x)));
			a[4] = dp_t(0, y); // mu

			// рассчитаем точку пересечения OY и прямой a[0]:a[3]
			// тут не надо fabs, потому что a[3].x > a[0].x
			y = a[3].x - a[0].x < FLT_MIN ? 0.5 * (a[0].y + a[3].y) : (a[0].y - a[0].x * ((a[3].y - a[0].y) / (a[3].x - a[0].x)));
			a[5] = dp_t(0, y); // nu

			if (is_points_belong_to_one_line(a[1], a[2], a[3]))
				return wall_1_middle_at;
			if (a[1].x < a[2].x && a[2].x > a[3].x)
				return wall_1_middle_out;
			if (a[2].x < a[1].x && a[2].x < a[3].x)
				return wall_1_middle_in;
			break;
		}
	default:
		return normal;
	}
	return normal;
}

static quad_type get_quadrangle_type(int i, int j,
                                     dp_t& a, dp_t& b, dp_t& c, dp_t& k, dp_t& m, dp_t& n, dp_t* p)
{
	// TODO какой порядок тут все таки предполагется? против часовой начиная с верхней левой?	
	dp_t alpha((OX[i - 1] + OX[i]) * 0.5, (OY[j - 1] + OY[j]) * 0.5),
		beta((OX[i + 1] + OX[i]) * 0.5, (OY[j - 1] + OY[j]) * 0.5),
		gamma((OX[i + 1] + OX[i]) * 0.5, (OY[j + 1] + OY[j]) * 0.5),
		theta((OX[i - 1] + OX[i]) * 0.5, (OY[j + 1] + OY[j]) * 0.5);
	
	// get prev coordnates
	double u = func_u(B, alpha);
	double v = func_v(UB, BB, LB, RB, TIME, alpha);
	p[0].x = alpha.x - TAU * u;
	p[0].y = alpha.y - TAU * v;
	v = func_v(UB, BB, LB, RB, TIME, beta);
	u = func_u(B, beta);
	p[1].x = beta.x - TAU * u;
	p[1].y = beta.y - TAU * v;

	v = func_v(UB, BB, LB, RB, TIME, gamma);
	u = func_u(B, gamma);
	p[2].x= gamma.x - TAU * u;
	p[2].y= gamma.y - TAU * v;	
	v = func_v(UB, BB, LB, RB, TIME, theta);
	u = func_u(B, theta);
	p[3].x = theta.x - TAU * u;
	p[3].y = theta.y - TAU * v;	

	dp_t intersection = get_intersection_point(p[0], p[1], p[2], p[3]);
	if ((p[1].y - intersection.y) * (p[3].y - intersection.y) > 0) return pseudo; // ??
	if ((p[0].x - intersection.x) * (p[2].x - intersection.x) > 0) return pseudo; // ??	
	if (is_points_belong_to_one_line(p[0], p[1], p[3])) return pseudo;

	a = p[0];
	b = p[1];
	c = p[2];
	k = p[0];
	m = p[3];
	n = p[2];

	return get_wall_intersection_type(p);
}

static double integrate_wall_triangle(const dp_t wp, // wall point
                                      double ly, // left y coordinate
                                      double ry) // right y coordinate
{
	return 0;
}

static double integrate_wall_rectangle(const dp_t wp1, const dp_t wp2, const dp_t wp3, const dp_t wp4, double wp1y, double wp2y)
{
	return 0;
}

static double integrate_wall_rectangle(const dp_t wp1, const dp_t wp2, double wp1y, double wp2y)
{
	return 0;
}

static double integrate_wall_pentagon(const dp_t wp1, const dp_t wp2, const dp_t wp3, double y1, double y2)
{
	return 0;
}

static double integrate_pentagon(const dp_t x, const dp_t y, const dp_t z, double ly, double ry)
{
	return 0;
}

static double integrate(int i, int j)
{
	dp_t a1, b1, c1, a2, b2, c2;
	dp_t* p = new dp_t[6];
	quad_type type = get_quadrangle_type(i, j, a1, b1, c1, a2, b2, c2, p);


	switch (type)
	{
	case wall_1_middle_in: // вообщем это один и тот же способ
	case wall_1_middle_out:
	case wall_1_middle_at:
		//	{
		//	//тут получается всегда 3 треугольника
		//		double result = 0;
		//		double t = 0;			
		//		dp_t v1 = p[4];
		//		dp_t v2 = p[2];
		//		dp_t v3 = p[1];
		//		sort_by_y_asc(v1, v2, v3);
		//		t = integrate_uniform_triangle(v1, v2, v3);
		//		result += t;
		//					
		//		v1 = p[4];
		//		v2 = p[2];
		//		v3 = p[5];
		//		sort_by_y_asc(v1, v2, v3);
		//		t = integrate_uniform_triangle(v1, v2, v3); // почему то тут приходит отрицательный результат
		//		result += t;
		//					
		//		v1 = p[3];
		//		v2 = p[2];
		//		v3 = p[5];
		//		sort_by_y_asc(v1, v2, v3);
		//		t = integrate_uniform_triangle(v1, v2, v3);
		//		result += t;
		//
		//		v1 = p[0];
		//		v2 = p[4];
		//		v3 = p[5];
		//		sort_by_y_asc(v1, v2, v3);
		//		t = integrate_uniform_triangle_wall(v1, v2, v3);
		//		result += t;
		//		return result;
		//	}
	case wall_2:
		{
			//// надо рассмотреть три случая
			//// 1. p2 внутри треугольника (p4,p5,p3) = итегрирование по 3 треугольникам
			//if (is_point_in_triangle(p[2], p[4], p[5], p[3]))
			//{
			//	double t = 0;
			//	double result = 0;
			//	dp_t v1 = p[5];
			//	dp_t v2 = p[3];
			//	dp_t v3 = p[4];
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
			//	return result;
			//}
			//// 2. p3 внутри треугольника (p4,p5,p2) =  итегрирование по 3 треугольникам
			//else if (is_point_in_triangle(p[3], p[4], p[5], p[2]))
			//{
			//	double t = 0;
			//	double result = 0;
			//	//сразу их располагаем в порядке возростания y координаты
			//	dp_t v1 = p[5];
			//	dp_t v2 = p[2];
			//	dp_t v3 = p[4];			
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
			//	return result;
			//}
			//// 3. ни 1 ни 2 условие - нормальный случай 4 угольник, интегрируем как обычно
			//else
			//{
			//	double result = 0;
			//	double t = 0;
			//	dp_t v1 = p[4];
			//	dp_t v2 = p[5];
			//	dp_t v3 = p[3];
			//	sort_by_y_asc(v1, v2, v3);
			//	t = integrate_uniform_triangle(v1, v2, v3);
			//	result += t;
			//	v1 = p[4];
			//	v2 = p[2];
			//	v3 = p[3];
			//	sort_by_y_asc(v1, v2, v3);
			//	t = integrate_uniform_triangle(v1, v2, v3);
			//	result += t;
			//	return result;
			//}
		}
	case wall_3_middle_in:
	case wall_3_middle_out:
	case wall_3_middle_at:
		{
			/*double result = 0;
		double t = 0;
		dp_t v1 = p[4];
		dp_t v2 = p[5];
		dp_t v3 = p[3];
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
			t = integrate_uniform_triangle(a1, b1, c1);
			result += t;
			sort_by_y_asc(a2, b2, c2);
			t = integrate_uniform_triangle(a2, b2, c2);
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

static void solve(double* density)
{
	PREV_DENSITY = new double[XY_LEN];
	for (int j = 0; j < OY_LEN + 1; j++)
	{
		for (int i = 0; i < OX_LEN_1; i++)
		{
			PREV_DENSITY[OX_LEN_1 * j + i] = analytical_solution(0, OX[i], OY[j]);
		}
	}

	for (TL = 1; TL <= TIME_STEP_CNT; TL++)
	{
		PREV_TIME = TIME;
		TIME = TAU * TL;
		for (int i = 0; i <= OX_LEN; i++)
		{
			density[i] = analytical_solution(OX[i], BB, TIME);
			density[OX_LEN_1 * OY_LEN + i] = analytical_solution(OX[i], UB, TIME);
		}

		for (int i = 0; i <= OY_LEN; i++)
		{
			density[OX_LEN_1 * i] = analytical_solution(LB, OY[i], TIME);
			density[OX_LEN_1 * i + OX_LEN] = analytical_solution(RB, OY[i], TIME);
		}

		for (int j = 1; j < OY_LEN; j++)
		{
			for (int i = 1; i < OX_LEN; i++)
			{
				density[OX_LEN_1 * j + i] = integrate(i, j) * INVERTED_HX_HY;
				density[OX_LEN_1 * j + i] += TAU * func_f(B, TIME, UB, BB, LB, RB, OX[i], OY[j]);
			}
		}
		memcpy(PREV_DENSITY, density, XY_LEN * sizeof(double));// заменить на быструю версию из agnerasmlib
	}
	delete [] PREV_DENSITY;
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
	TL = 0;
	delete [] OX;
	delete [] OY;
}

double* compute_density(double b, double lb, double rb, double bb, double ub,
                        double tau, int time_step_count, int ox_length, int oy_length, double& norm)
{
	init(b, lb, rb, bb, ub, tau, time_step_count, ox_length, oy_length);
	double* density = new double[XY_LEN];
	print_params(B, LB, RB, BB, UB, TAU, TIME_STEP_CNT, OX_LEN, OY_LEN);
	solve(density);
	norm = get_norm_of_error(density, TIME_STEP_CNT * TAU);
	clean();
	return density;
}