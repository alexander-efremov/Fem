#include "common.h"
#include "point.h"
#include "utils.h"

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

static double B;
static double UB;
static double BB;
static double LB;
static double RB;
static double TAU;
static int OX_LEN;
static int OX_LEN_1; // OX_LEN_1
static int OY_LEN;
static int XY_LEN;
static int TIME_STEP_CNT;
static double HX;
static double HY;
static double* OX;
static double* OY;
static double* PREV_DENSITY;
static int TL;
static double TAU_TL;
static double TAU_TL_1; // tau * (tl - 1)
static double INVERTED_HX_HY;

__pure inline static void sort_by_y(dp_t& x, dp_t& y, dp_t& z)
{
	double t;
	if (x.y < y.y)
	{
		if (z.y < x.y)
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
		if (y.y < z.y)
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
			t = y.y;
			y.y = z.y;
			z.y = t;
			//swap(x, z);
		}
	}
	if (z.y < y.y)
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
	double b1 = -(gamma.x - alpha.x);
	double c1 = a1 * alpha.x + b1 * alpha.y;
	double a2 = theta.y - beta.y;
	double b2 = -(theta.x - beta.x);
	double c2 = a2 * beta.x + b2 * beta.y;
	return dp_t((b1 * c2 - b2 * c1) / (b1 * a2 - b2 * a1), (a1 * c2 - a2 * c1) / (-b1 * a2 + b2 * a1));
}

__pure inline static double get_vector_product(const dp_t& alpha, const dp_t beta, const dp_t theta)
{
	return (beta.x - alpha.x) * (theta.y - alpha.y) - (beta.y - alpha.y) * (theta.x - alpha.x);
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
	return atan((x - lb) * (x - rb) * (1 + t) * 0.1 * (y - ub) * (y - bb));
}

__pure inline static double func_v(double ub, double bb, double lb, double rb, double t, const dp_t& p)
{
	return func_v(ub, bb, lb, rb, t, p.x, p.y);
}

__pure inline static double func_f(double b, double tau_tl, double ub, double bb, double lb, double rb, double x, double y)
{
	double arg_v = (x - lb) * (x - rb) * (1 + tau_tl) * 0.1 * (y - ub) * (y - bb);
	double rho = analytical_solution(tau_tl, x, y);
	double drho_dt = x * y * cos(tau_tl * x * y);
	double drho_dx = tau_tl * y * cos(tau_tl * x * y);
	double dtho_dy = tau_tl * x * cos(tau_tl * x * y);
	double u = func_u(b, x, y);
	double v = func_v(ub, bb, lb, rb, tau_tl, x, y);
	double du_dx = -B * y * (1 - y) / (1 + x * x);
	double dv_dx = (x - lb) * (x - rb) * (1 + tau_tl) / 10 * (y - bb + y - ub);
	dv_dx /= (1 + arg_v * arg_v);
	double res = drho_dt + rho * du_dx + u * drho_dx + rho * dv_dx + v * dtho_dy;
	// print_f_params()...
	return res;
}

__pure inline static double integrate_rectangle(double py, double qy, double gx, double hx, double a, double b)
{
	return ((hx - a) * (hx - a) - (gx - a) * (gx - a)) * ((qy - b) * (qy - b) - (py - b) * (py - b)) * 0.25;
}

__pure inline static double integrate_triangle(double py, double qy, double alpha, double a, double b, double beta)
{
	return 0.5 * ((((qy - alpha) * (a * qy + b - beta) * (a * qy + b - beta) * (a * qy + b - beta)
		- (py - alpha) * (a * py + b - beta) * (a * py + b - beta) * (a * py + b - beta)) / (3 * a)) - ((a * qy + b - beta) * (a * qy + b - beta) * (a * qy + b - beta) * (a * qy + b - beta)
		- (a * py + b - beta) * (a * py + b - beta) * (a * py + b - beta) * (a * py + b - beta)) / (12 * a * a));
}

static double integrate_rectangle_one_cell(double py, double qy, double gx, double hx, const ip_t& sx, const ip_t& sy)
{
	double result, a, b;
	a = sx.y >= 0 && sy.y >= 0 ? OX[sx.y] : HX * sx.y;
	b = sx.y >= 0 && sy.y >= 0 ? OY[sy.y] : HY * sy.y; // ЭТО ПЛОТНОСТЬ С ПРЕДЫДУЩЕГО СЛОЯ ДЛЯ ДАННОЙ ЯЧЕЙКИ
	result = integrate_rectangle(py, qy, gx, hx, a, b) * (sx.x >= 0 && sy.x >= 0 ? PREV_DENSITY[OX_LEN_1 * sy.x + sx.x] : analytical_solution(TAU_TL_1, sx.x * HX, sy.x * HY));
	a = sx.x >= 0 && sy.y >= 0 ? OX[sx.x] : HX * sx.x;
	b = sx.x >= 0 && sy.y >= 0 ? OY[sy.y] : HY * sy.y;
	result -= integrate_rectangle(py, qy, gx, hx, a, b) * (sx.x >= 0 && sy.x >= 0 ? PREV_DENSITY[OX_LEN_1 * sy.x + sx.y] : analytical_solution(TAU_TL_1, sx.y * HX, sy.x * HY));
	a = sx.y >= 0 && sy.x >= 0 ? OX[sx.y] : HX * sx.y;
	b = sx.y >= 0 && sy.x >= 0 ? OY[sy.x] : HY * sy.x;
	result -= integrate_rectangle(py, qy, gx, hx, a, b) * (sx.x >= 0 && sy.x >= 0 ? PREV_DENSITY[OX_LEN_1 * sy.y + sx.x] : analytical_solution(TAU_TL_1, sx.x * HX, sy.y * HY));
	a = sx.x >= 0 && sy.x >= 0 ? OX[sx.x] : HX * sx.x;
	b = sx.x >= 0 && sy.x >= 0 ? OY[sy.x] : HY * sy.x;
	result += integrate_rectangle(py, qy, gx, hx, a, b) * (sx.x >= 0 && sy.x >= 0 ? PREV_DENSITY[OX_LEN_1 * sy.y + sx.y] : analytical_solution(TAU_TL_1, sx.y * HX, sy.y * HY));
	return result * INVERTED_HX_HY;
}

static double integrate_triangle_left_one_cell(const dp_t& bv, const dp_t& uv, double hx, // TODO: вычислять внутри!!!!!!
                                               const ip_t& sx, const ip_t& sy)
{
	double a_sl = (bv.x - uv.x) / (bv.y - uv.y); //   Coefficients of slant line: x = a_SL *y  +  b_SL.
	if (fabs(a_sl) <= FLT_MIN) return 0;
	double b_sl = uv.x - a_sl * uv.y;
	double result = 0, tmp, a, b;
	a = sx.y >= 0 && sy.y >= 0 ? OX[sx.y] : HX * sx.y;
	b = sx.y >= 0 && sy.y >= 0 ? OY[sy.y] : HY * sy.y;
	tmp = ((uv.y - OY[sy.y]) * (uv.y - OY[sy.y]) - (bv.y - OY[sy.y]) * (bv.y - OY[sy.y])) * (hx - a) * (hx - a) * 0.25 - integrate_triangle(bv.y, uv.y, b, a_sl, b_sl, a);
	result += tmp * (sx.x >= 0 && sx.y <= OX_LEN && sy.x >= 0 && sy.y <= OY_LEN ? PREV_DENSITY[OX_LEN_1 * sy.x + sx.x] : analytical_solution(TAU_TL_1, sx.x * HX, sy.x * HY));
	a = sx.x >= 0 && sy.y >= 0 ? OX[sx.x] : HX * sx.x;
	tmp = (uv.y - OY[sy.y]) * (uv.y - OY[sy.y]) - (bv.y - OY[sy.y]) * (bv.y - OY[sy.y]);
	tmp = tmp * (hx - a) * (hx - a) * -0.25 + integrate_triangle(bv.y, uv.y, b, a_sl, b_sl, a);
	result += tmp * (sx.x >= 0 && sx.y <= OX_LEN && sy.x >= 0 && sy.y <= OY_LEN ? PREV_DENSITY[OX_LEN_1 * sy.x + sx.y] : analytical_solution(TAU_TL_1, sx.y * HX, sy.x * HY));
	a = sx.y >= 0 && sy.x >= 0 ? OX[sx.y] : HX * sx.y;
	b = sx.y >= 0 && sy.x >= 0 ? OY[sy.x] : HY * sy.x;
	tmp = (uv.y - OY[sy.x]) * (uv.y - OY[sy.x]) - (bv.y - OY[sy.x]) * (bv.y - OY[sy.x]);
	tmp = tmp * (hx - a) * (hx - a) * -0.25 + integrate_triangle(bv.y, uv.y, b, a_sl, b_sl, a);
	result += tmp * (sx.x >= 0 && sx.y <= OX_LEN && sy.x >= 0 && sy.y <= OY_LEN ? PREV_DENSITY[OX_LEN_1 * sy.y + sx.x] : analytical_solution(TAU_TL_1, sx.x * HX, sy.y * HY));
	a = sx.x >= 0 && sy.x >= 0 ? OX[sx.x] : HX * sx.x;
	b = sx.x >= 0 && sy.x >= 0 ? OY[sy.x] : HY * sy.x;
	tmp = (uv.y - OY[sy.x]) * (uv.y - OY[sy.x]) - (bv.y - OY[sy.x]) * (bv.y - OY[sy.x]);
	tmp = tmp * (hx - a) * (hx - a) * 0.25 - integrate_triangle(bv.y, uv.y, b, a_sl, b_sl, a);
	result += tmp * (sx.x >= 0 && sx.y <= OX_LEN && sy.x >= 0 && sy.y <= OY_LEN ? PREV_DENSITY[OX_LEN_1 * sy.y + sx.y] : analytical_solution(TAU_TL_1, sx.y * HX, sy.y * HY));
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
	double x = uv.x <= bv.x ? bv.x : uv.x; // TODO: занести hx в integrate_triangle_left_one_cell

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
		{ //   Intersection with straight line parallel Ox axis.        
			next_i = 1;
			next.y = sy.y >= 0 ? OY[sy.y] : HY * sy.y;
			next.x = curr.x - (next.y - curr.y) / k;
		}
		else
		{ //   Intersection with straight line parallel Oy axis.            
			next_i = 2;
			next.x = sx.x >= 0 ? OX[sx.x] : HX * sx.x;
			next.y = curr.y - k * (next.x - curr.x);
		}
		if (next.x < (uv.x + FLT_MIN))
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
		{//   Intersection with straight line parallel Ox axis.            
			next_i = 1;
			next.y = sy.y >= 0 ? OY[sy.y] : HY * sy.y;
			next.x = bv.x + (next.y - bv.y) / k;
		}
		else
		{//   Intersection with straight line parallel OY axis.
			next_i = 2;
			next.x = sx.y >= 0 ? OX[sx.y] : HX * sx.y;
			next.y = bv.y + k * (next.x - bv.x);
		}
		if (next.x > (uv.x - FLT_MIN))
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
		{ //   intersection with straight line parallel Ox axis.
			next_i = 1;
			next.y = sy.y >= 0 ? OY[sy.y] : HY * sy.y;
			next.x = bv.x + (next.y - bv.y) / k;
		}
		else
		{//   intersection with straight line parallel OY axis.            
			next_i = 2;
			next.x = sx.y >= 0 ? OX[sx.y] : HX * sx.y;
			next.y = bv.y + k * (next.x - bv.x);
		}
		if (next.x > (uv.x - FLT_MIN)) // если следующая точка уже правее, чем наша граничная точка, то мы обработали канал
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
		{ //   Intersection with straight line parallel Ox axis.
			next_i = 1;
			next.y = sy.y >= 0 ? OY[sy.y] : HY * sy.y;
			next.x = bv.x - (next.y - bv.y) / k;
		}
		else
		{ //   Intersection with straight line parallel Oy axis.
			next_i = 2;
			next.x = sx.x >= 0 ? OX[sx.x] : HX * sx.x;
			next.y = bv.y - k * (next.x - bv.x);
		}
		if (next.x < uv.x + FLT_MIN)
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

static double integrate_uniform_triangle_wall(const dp_t& a, const dp_t& b, const dp_t& c)
{
	return 0;
}

static double integrate_uniform_triangle(const dp_t& x, dp_t& y, const dp_t& z)
{
	//   a * x  +  b * y  = c.
	double a = z.y - x.y;
	if (fabs(a) < FLT_MIN) return FLT_MIN ;
	double b = x.x - z.x;
	double c = b * x.y + a * x.x;
	dp_t ip((c - b * y.y) / a, y.y);

	//   Возможны 2 случая расположения точки перечеения относительно средней
	//   слева или справа.
	//   есди средняя точка справа от точки пересечения
	//   обменяем местами  X координаты, чтобы использовать один код для расчета
	if (y.x >= ip.x)
	{
		double tx = y.x;
		y.x = ip.x;
		ip.x = tx;
	}

	return integrate_upper_triangle(y, z, ip) + integrate_bottom_triangle(y, x, ip);
}

static quad_type get_coordinates_on_prev_layer(int ix, int iy,
                                               dp_t& alpha, dp_t& beta, dp_t& gamma, dp_t& theta)
{
	//   1 First of all let's compute coordinates of square vertexes.
	//  OX:
	if (ix == 0)
	{
		alpha.x = OX[ix];
		beta.x = (OX[ix] + OX[ix + 1]) * 0.5;
		gamma.x = (OX[ix] + OX[ix + 1]) * 0.5;
		theta.x = OX[ix];
	}
	else if (ix == OX_LEN)
	{
		alpha.x = (OX[ix - 1] + OX[ix]) * 0.5;
		beta.x = OX[ix];
		gamma.x = OX[ix];
		theta.x = (OX[ix - 1] + OX[ix]) * 0.5;
	}
	else
	{
		alpha.x = (OX[ix - 1] + OX[ix]) * 0.5;
		beta.x = (OX[ix + 1] + OX[ix]) * 0.5;
		gamma.x = (OX[ix + 1] + OX[ix]) * 0.5;
		theta.x = (OX[ix - 1] + OX[ix]) * 0.5;
	}

	//  OY:
	if (iy == 0)
	{
		alpha.y = OY[iy];
		beta.y = OY[iy];
		gamma.y = (OY[iy] + OY[iy + 1]) * 0.5;
		theta.y = (OY[iy] + OY[iy + 1]) * 0.5;
	}
	else if (iy == OY_LEN)
	{
		alpha.y = (OY[iy] + OY[iy - 1]) * 0.5;
		beta.y = (OY[iy] + OY[iy - 1]) * 0.5;
		gamma.y = OY[iy];
		theta.y = OY[iy];
	}
	else
	{
		alpha.y = (OY[iy] + OY[iy - 1]) * 0.5;
		beta.y = (OY[iy] + OY[iy - 1]) * 0.5;
		gamma.y = (OY[iy] + OY[iy + 1]) * 0.5;
		theta.y = (OY[iy] + OY[iy + 1]) * 0.5;
	}

	double u, v;

	// Now let's compute new coordinates on the previous time level of alpha, beta, gamma, theta points.
	u = func_u(B, alpha);
	v = func_v(UB, BB, LB, RB, TAU_TL, alpha);
	alpha.x -= TAU * u;
	alpha.y -= TAU * v;

	u = func_u(B, beta);
	v = func_v(UB, BB, LB, RB, TAU_TL, beta);
	beta.x -= TAU * u;
	beta.y -= TAU * v;

	u = func_u(B, gamma);
	v = func_v(UB, BB, LB, RB, TAU_TL, gamma);
	gamma.x -= TAU * u;
	gamma.y -= TAU * v;

	u = func_u(B, theta);
	v = func_v(UB, BB, LB, RB, TAU_TL, theta);
	theta.x -= TAU * u;
	theta.y -= TAU * v;

	dp_t intersection = get_intersection_point(alpha, beta, gamma, theta);
	if ((beta.y - intersection.y) * (theta.y - intersection.y) > 0) return pseudo; // ??
	if ((alpha.x - intersection.x) * (gamma.x - intersection.x) > 0) return pseudo; // ??
	double product = get_vector_product(alpha, beta, theta); // ?
	if (product < 0) return pseudo;

	// значит что точка улетела за левую границу
	if (theta.x < 0 || theta.y < 0 || beta.x < 0 || beta.y < 0 || gamma.x < 0 ||
		gamma.y < 0 || alpha.x < 0 || alpha.y < 0)
	{
		return normal;
		//return wall;
	}
	return normal;
}

// Type of quadrangle: 0 - pseudo; 1 - convex; 2 - concave;

static quad_type get_quadrangle_type(int ix, int iy,
                                     dp_t& a, //   -  First vertex of first triangle.
                                     dp_t& b, //   -  Second vertex of first triangle.
                                     dp_t& c, //   -  Third vertex of first triangle.
                                     dp_t& k, //   -  First vertex of second triangle.
                                     dp_t& m, //   -  Second vertex of second triangle.
                                     dp_t& n) //   -  Third vertex of second triangle.
{
	dp_t alpha, beta, gamma, theta; // coordinates on previous time layer
	quad_type type = get_coordinates_on_prev_layer(ix, iy, alpha, beta, gamma, theta);
	a = alpha;
	b = beta;
	c = gamma;
	k = alpha;
	m = theta;
	n = gamma;
	return type;
}

static double integrate(int ix, int iy)
{
	dp_t a1, b1, c1, a2, b2, c2;
	quad_type type = get_quadrangle_type(ix, iy, a1, b1, c1, a2, b2, c2);
	if (type != normal && type != wall)
	{
		return -1;
	}

	// чтобы правилно отработала процедура интегрирования
	// точки должны идти в порядке возрастания y координаты
	sort_by_y(a1, b1, c1);
	sort_by_y(a2, b2, c2);

	// check the type of triangle to select appropriate computation method
	double result = 0;
	switch (type)
	{
	case wall:
		result += integrate_uniform_triangle_wall(a1, b1, c1);
		result += integrate_uniform_triangle_wall(a2, b2, c2);
		return result;
	case normal:
		result += integrate_uniform_triangle(a1, b1, c1);
		result += integrate_uniform_triangle(a2, b2, c2);
		return result;
	case concave:
	case convex:
	case pseudo:
		return 0;
	}
	return 0;
}

static double get_norm_of_error(double* density, double ts_count_mul_steps)
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
		TAU_TL = TAU * TL;
		TAU_TL_1 = TAU * (TL - 1);
		for (int i = 0; i <= OX_LEN; i++)
		{
			density[i] = analytical_solution(OX[i], BB, TAU_TL);
			density[OX_LEN_1 * OY_LEN + i] = analytical_solution(OX[i], UB, TAU_TL);
		}

		for (int i = 0; i <= OY_LEN; i++)
		{
			density[OX_LEN_1 * i] = analytical_solution(LB, OY[i], TAU_TL);
			density[OX_LEN_1 * i + OX_LEN] = analytical_solution(RB, OY[i], TAU_TL);
		}

		for (int i = 1; i < OY_LEN; i++)
		{
			for (int j = 1; j < OX_LEN; j++)
			{
				density[OX_LEN_1 * i + j] = integrate(j, i) * INVERTED_HX_HY;
				density[OX_LEN_1 * i + j] += TAU * func_f(B, TAU_TL, UB, BB, LB, RB, OX[j], OY[i]);
			}
		}
		memcpy(PREV_DENSITY, density, XY_LEN * sizeof(double));// заменить на быструю версию из agnerasmlib
	}
	delete[] PREV_DENSITY;
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
	TAU_TL = 0;
	OX_LEN = 0;
	OY_LEN = 0;
	OX_LEN_1 = 0;
	TIME_STEP_CNT = 0;
	XY_LEN = 0;
	HX = 0;
	HY = 0;
	INVERTED_HX_HY = 0;
	TL = 0;
	delete[] OX;
	delete[] OY;
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