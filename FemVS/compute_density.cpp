#include "common.h"
#include "point.h"
#include "utils.h"

static int TMP_WALL_CNT = 0;
static double B;
static double UB;
static double BB;
static double LB;
static double RB;
static double TAU;
static int OX_LEN;
static int OY_LEN;
static int XY_LEN;
static int TIME_STEP_CNT;
static double HX;
static double HY;
static double *OX;
static double *OY;
static double *PREV_DENSITY;
static int TL;
static double TAU_TL;
static double TAU_TL_1; // tau * (tl - 1)
static double INVERTED_HX_HY;

inline static void init(double b, double lb, double rb, double bb, double ub,
	double tau, int time_step_count, int ox_length, int oy_length) {
	B = b;
	UB = ub;
	BB = bb;
	LB = lb;
	RB = rb;
	TAU = tau;
	TIME_STEP_CNT = time_step_count;
	XY_LEN = (ox_length + 1) * (oy_length + 1);
	OX_LEN = ox_length;
	OY_LEN = oy_length;
	OX = new double[OX_LEN + 1];
	OY = new double[OY_LEN + 1];
	for (int i = 0; i <= OX_LEN; ++i) OX[i] = lb + i * (rb - lb) / OX_LEN;
	for (int i = 0; i <= OY_LEN; ++i) OY[i] = bb + i * (ub - bb) / OY_LEN;
	HX = OX[1] - OX[0];
	HY = OY[1] - OY[0];
	INVERTED_HX_HY = 1 / HX / HY;
}

inline static void clean() {
	B = 0;
	UB = 0;
	BB = 0;
	LB = 0;
	RB = 0;
	TAU = 0;
	TAU_TL = 0;
	OX_LEN = 0;
	OY_LEN = 0;
	TIME_STEP_CNT = 0;
	TMP_WALL_CNT = 0;
	XY_LEN = 0;
	HX = 0;
	HY = 0;
	INVERTED_HX_HY = 0;
	TL = 0;
	delete[] OX;
	delete[] OY;
}

inline static void sort_by_y(dp_t& x, dp_t& y, dp_t& z) {
	if (x.y < y.y) {
		if (z.y < x.y)
			std::swap(x, z);
	}
	else {
		if (y.y < z.y)
			std::swap(x, y);
		else
			std::swap(x, z);
	}
	if (z.y < y.y) std::swap(y, z);
}

inline static bool try_get_slope_ratio(const dp_t &bv, const dp_t &uv, double &value) {
	if (fabs(bv.x - uv.x) < MIN_VALUE) {
		return false;
	}
	value = fabs((uv.y - bv.y) / (uv.x - bv.x)); // угловой коэффициент прямой
	if (value < MIN_VALUE) {
		return false;
	}
	return true;
}

inline static dp_t get_intersection_point(const dp_t& alpha, const dp_t& beta, const dp_t& gamma, const dp_t& theta) {
	dp_t result;
	dp_t alpha_to_gamma;
	dp_t beta_to_theta;
	double a_1LC, b_1LC, c_1LC;
	double a_2LC, b_2LC, c_2LC;
	alpha_to_gamma.x = gamma.x - alpha.x;
	alpha_to_gamma.y = gamma.y - alpha.y;
	a_1LC = alpha_to_gamma.y;
	b_1LC = -alpha_to_gamma.x;
	c_1LC = alpha_to_gamma.y * alpha.x - alpha_to_gamma.x * alpha.y;
	beta_to_theta.x = theta.x - beta.x;
	beta_to_theta.y = theta.y - beta.y;
	a_2LC = beta_to_theta.y;
	b_2LC = -beta_to_theta.x;
	c_2LC = beta_to_theta.y * beta.x - beta_to_theta.x * beta.y;
	result.x = (b_1LC * c_2LC - b_2LC * c_1LC) / (b_1LC * a_2LC - b_2LC * a_1LC);
	result.y = (a_1LC * c_2LC - a_2LC * c_1LC) / (-b_1LC * a_2LC + b_2LC * a_1LC);
	return result;
}

inline static double get_vector_product(const dp_t& alpha, const dp_t beta, const dp_t theta) {
	dp_t alpha_to_beta;
	dp_t alpha_to_theta;
	alpha_to_beta.x = beta.x - alpha.x;
	alpha_to_beta.y = beta.y - alpha.y;
	alpha_to_theta.x = theta.x - alpha.x;
	alpha_to_theta.y = theta.y - alpha.y;
	return alpha_to_beta.x * alpha_to_theta.y - alpha_to_beta.y * alpha_to_theta.x;
}

inline static double analytical_solution(double t, double x, double y) {
	return 1.1 + sin(t * x * y);
}

inline static double init_side(double x, double y, double t) {
	return analytical_solution(t, x, y);
}

inline static double func_u(double x, double y) {
	return B * y * (1 - y) * (M_PI_2 + atan(-x));
}

inline static double func_u(const dp_t &p) {
	return func_u(p.x, p.y);
}

inline static double func_v(double t, double x, double y) {
	return atan((x - LB) * (x - RB) * (1 + t) / 10 * (y - UB) * (y - BB));
}

inline static double func_v(double t, const dp_t &p) {
	return func_v(t, p.x, p.y);
}

inline static double func_f(double x, double y) {
	double arg_v = (x - LB) * (x - RB) * (1 + TAU_TL) / 10 * (y - UB) * (y - BB);
	double rho = analytical_solution(TAU_TL, x, y);
	double drho_dt = x * y * cos(TAU_TL * x * y);
	double drho_dx = TAU_TL * y * cos(TAU_TL * x * y);
	double dtho_dy = TAU_TL * x * cos(TAU_TL * x * y);
	double u = func_u(x, y);
	double v = func_v(TAU_TL, x, y);
	double du_dx = -B * y * (1 - y) / (1 + x * x);
	double dv_dx = (x - LB) * (x - RB) * (1 + TAU_TL) / 10 * (y - BB + y - UB);
	dv_dx /= (1 + arg_v * arg_v);
	double res = drho_dt + rho * du_dx + u * drho_dx + rho * dv_dx + v * dtho_dy;
	// print_f_params()...
	return res;
}

inline static double integrate_rectangle(double py, double qy, double gx, double hx, double a, double b) {
	return ((hx - a) * (hx - a) - (gx - a) * (gx - a)) * ((qy - b) * (qy - b) - (py - b) * (py - b)) * 0.25;
}

inline static double integrate_triangle(double py, double qy, double alpha, double a, double b,
	double beta) {
	return (((qy - alpha) * (a * qy + b - beta) * (a * qy + b - beta) * (a * qy + b - beta)
		- (py - alpha) * (a * py + b - beta) * (a * py + b - beta) * (a * py + b - beta)) / (3 * a)) - ((a * qy + b - beta) * (a * qy + b - beta) * (a * qy + b - beta) * (a * qy + b - beta)
		- (a * py + b - beta) * (a * py + b - beta) * (a * py + b - beta) * (a * py + b - beta)) / (12 * a * a);
}

double integrate_rectangle_one_cell(double py, double qy, double gx, double hx,
	const ip_t &sx, const ip_t &sy) {
	double result, tmp = 0;
	double rho[2][2];
	if (sx.x >= 0 && sy.x >= 0) {
		rho[0][0] = PREV_DENSITY[(OX_LEN + 1) * sy.x + sx.x];
		rho[0][1] = PREV_DENSITY[(OX_LEN + 1) * sy.y + sx.x];
		rho[1][0] = PREV_DENSITY[(OX_LEN + 1) * sy.x + sx.y];
		rho[1][1] = PREV_DENSITY[(OX_LEN + 1) * sy.y + sx.y];
	}
	else {
		// TODO: убрать потому что это неверно (надо расчитывать граничные условия)
		rho[0][0] = analytical_solution(TAU_TL_1, sx.x * HX, sy.x * HY);
		rho[0][1] = analytical_solution(TAU_TL_1, sx.x * HX, sy.y * HY);
		rho[1][0] = analytical_solution(TAU_TL_1, sx.y * HX, sy.x * HY);
		rho[1][1] = analytical_solution(TAU_TL_1, sx.y * HX, sy.y * HY);
		TMP_WALL_CNT++;
	}

	if (sx.y >= 0 && sy.y >= 0) {
		tmp = integrate_rectangle(py, qy, gx, hx, OX[sx.y], OY[sy.y]);
	}
	else {
		tmp = integrate_rectangle(py, qy, gx, hx, HX * sx.y, HY * sy.y);
	}
	tmp = tmp * INVERTED_HX_HY;
	result = tmp * rho[0][0];
	if (sx.x >= 0 && sy.y >= 0) {
		tmp = integrate_rectangle(py, qy, gx, hx, OX[sx.x], OY[sy.y]);
	}
	else {
		tmp = integrate_rectangle(py, qy, gx, hx, HX * sx.x, HY * sy.y);
	}
	tmp = tmp * INVERTED_HX_HY;
	result -= tmp * rho[1][0];
	if (sx.y >= 0 && sy.x >= 0) {
		tmp = integrate_rectangle(py, qy, gx, hx, OX[sx.y], OY[sy.x]);
	}
	else {
		tmp = integrate_rectangle(py, qy, gx, hx, HX * sx.y, HY * sy.x);
	}
	tmp = tmp * INVERTED_HX_HY;
	result -= tmp * rho[0][1];
	if (sx.x >= 0 && sy.x >= 0) {
		tmp = integrate_rectangle(py, qy, gx, hx, OX[sx.x], OY[sy.x]);
	}
	else {
		tmp = integrate_rectangle(py, qy, gx, hx, HX * sx.x, HY * sy.x);
	}
	tmp = tmp * INVERTED_HX_HY;
	return result + tmp * rho[1][1];
}

double integrate_triangle_left_one_cell(const dp_t &bv, const dp_t &uv, double hx,
	const ip_t &sx, const ip_t &sy) {
	if (fabs(bv.y - uv.y) <= FLT_MIN) return 0;
	double a_sl = (bv.x - uv.x) / (bv.y - uv.y); //   Coefficients of slant line: x = a_SL *y  +  b_SL.
	if (fabs(a_sl) <= FLT_MIN) return 0;
	double b_sl = uv.x - a_sl * uv.y;

	double result, tmp, tmp_integral;
	double rho[2][2];

	if (sx.x >= 0 && sx.y <= OX_LEN && sy.x >= 0 && sy.y <= OY_LEN) {
		rho[0][0] = PREV_DENSITY[(OX_LEN + 1) * sy.x + sx.x];
		rho[0][1] = PREV_DENSITY[(OX_LEN + 1) * sy.y + sx.x];
		rho[1][0] = PREV_DENSITY[(OX_LEN + 1) * sy.x + sx.y];
		rho[1][1] = PREV_DENSITY[(OX_LEN + 1) * sy.y + sx.y];
	}
	else {
		// TODO: убрать потому что это неверно (надо расчитывать граничные условия)
		// норма должна уменьшиться
		rho[0][0] = analytical_solution(TAU_TL_1, sx.x * HX, sy.x * HY);
		rho[0][1] = analytical_solution(TAU_TL_1, sx.x * HX, sy.y * HY);
		rho[1][0] = analytical_solution(TAU_TL_1, sx.y * HX, sy.x * HY);
		rho[1][1] = analytical_solution(TAU_TL_1, sx.y * HX, sy.y * HY);
		TMP_WALL_CNT++;
	}

	//   1
	tmp = (uv.y - OY[sy.y]) * (uv.y - OY[sy.y]) - (bv.y - OY[sy.y]) * (bv.y - OY[sy.y]);
	if (sx.y >= 0 && sy.y >= 0) {
		tmp *= (hx - OX[sx.y]) * (hx - OX[sx.y]) * 0.25;
		tmp_integral = integrate_triangle(bv.y, uv.y, OY[sy.y], a_sl, b_sl, OX[sx.y]);
	}
	else {
		tmp *= (hx - HX * sx.y) * (hx - HX * sx.y) * 0.25;
		tmp_integral = integrate_triangle(bv.y, uv.y, HY * sy.y, a_sl, b_sl, HX * sx.y);
	}
	tmp -= tmp_integral * 0.5;
	result = tmp * rho[0][0] * INVERTED_HX_HY;

	//   2
	tmp = (uv.y - OY[sy.y]) * (uv.y - OY[sy.y]) - (bv.y - OY[sy.y]) * (bv.y - OY[sy.y]);
	if (sx.x >= 0 && sy.y >= 0) {
		tmp *= (hx - OX[sx.x]) * (hx - OX[sx.x]) * -0.25;
		tmp_integral = integrate_triangle(bv.y, uv.y, OY[sy.y], a_sl, b_sl, OX[sx.x]);
	}
	else {
		tmp *= (hx - HX * sx.x) * (hx - HX * sx.x) * -0.25;
		tmp_integral = integrate_triangle(bv.y, uv.y, HY * sy.y, a_sl, b_sl, HX * sx.x);
	}
	tmp += tmp_integral * 0.5;
	result += tmp * rho[1][0] * INVERTED_HX_HY;

	//   3
	tmp = (uv.y - OY[sy.x]) * (uv.y - OY[sy.x]) - (bv.y - OY[sy.x]) * (bv.y - OY[sy.x]);
	if (sx.y >= 0 && sy.x >= 0) {
		tmp *= (hx - OX[sx.y]) * (hx - OX[sx.y]) * -0.25;
		tmp_integral = integrate_triangle(bv.y, uv.y, OY[sy.x], a_sl, b_sl, OX[sx.y]);
	}
	else {
		tmp *= (hx - HX * sx.y) * (hx - HX * sx.y) * -0.25;
		tmp_integral = integrate_triangle(bv.y, uv.y, HY * sy.x, a_sl, b_sl, HX * sx.y);
	}
	tmp += tmp_integral * 0.5;
	result += tmp * rho[0][1] * INVERTED_HX_HY;

	//   4
	tmp = (uv.y - OY[sy.x]) * (uv.y - OY[sy.x]) - (bv.y - OY[sy.x]) * (bv.y - OY[sy.x]);
	if (sx.x >= 0 && sy.x >= 0) {
		tmp *= (hx - OX[sx.x]) * (hx - OX[sx.x]) * 0.25;
		tmp_integral = integrate_triangle(bv.y, uv.y, OY[sy.x], a_sl, b_sl, OX[sx.x]);
	}
	else {
		tmp *= (hx - HX * sx.x) * (hx - HX * sx.x) * 0.25;
		tmp_integral = integrate_triangle(bv.y, uv.y, HY * sy.x, a_sl, b_sl, HX * sx.x);
	}
	tmp -= tmp_integral * 0.5;
	result += tmp * rho[1][1] * INVERTED_HX_HY;
	return result;
}

double integrate_chanel_slant_right(const dp_t& bv, const dp_t& uv,
	int curr_i, int next_i, const ip_t &sx, double b, const ip_t &sb,
	const ip_t &sy) {
	if (fabs(uv.y - bv.y) <= FLT_MIN) return fabs(uv.y - bv.y);

	double result = 0, gx = 0, hx = 0;
	dp_t mv, rv;
	int m_i;
	if (uv.x <= bv.x) {
		mv = uv;
		m_i = next_i;
	}
	else {
		mv = bv;
		m_i = curr_i;
	}

	//   A. Under rectangle.
	result += -1 * integrate_triangle_left_one_cell(bv, uv, mv.x, sx, sy);

	// case B: не полный прямоугольник    
	if (m_i == 1) {
		if (sx.x == sb.x) gx = b;
		if (sx.x > sb.x) {
			gx = sx.x >= 0 ? OX[sx.x] : HX * sx.x;
		}
		result += integrate_rectangle_one_cell(bv.y, uv.y, gx, mv.x, sx, sy);
	}

	//   А теперь прибавим все прямоугольные куски, которые помещаются в ячейку
	ip_t ch_pos(sb.x, sb.x + 1);
	for (int j = sb.x; j < sx.x; j++) {
		if (j == sb.x) gx = b;
		else gx = ch_pos.x >= 0 ? OX[ch_pos.x] : HX * ch_pos.x;
		hx = ch_pos.x >= 0 ? OX[ch_pos.y] : HX * ch_pos.y;
		result += integrate_rectangle_one_cell(bv.y, uv.y, gx, hx, ch_pos, sy);
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

double integrate_chanel_slant_left(const dp_t& bv, const dp_t& uv,
	int curr_i, int next_i, const ip_t &sx, const ip_t &sy,
	double b, const ip_t &sb) {
	if (fabs(uv.y - bv.y) <= FLT_MIN) return fabs(uv.y - bv.y);

	dp_t lv, mv; //   -  Left and middle vertices.
	int m_i = 0; //   -  Where middle vertex is.        
	double result = 0, gx = 0, hx = 0; //   -  Left and right boundary for each integration.   

	// зачем то определили среднюю точку
	if (uv.x <= bv.x) {
		mv = bv;
		m_i = curr_i;
	}
	else {
		mv = uv;
		m_i = next_i;
	}

	// case A: triangle
	result += integrate_triangle_left_one_cell(bv, uv, mv.x, sx, sy);

	// case B: не полный прямоугольник
	if (m_i == 1) { // это значит, что прямоугольник занимает не всю ячейку  
		hx = sx.x == sb.x ? b : (sx.y >= 0 ? OX[sx.y] : HX * sx.y);
		gx = mv.x;
		result += integrate_rectangle_one_cell(bv.y, uv.y, gx, hx, sx, sy);
	}

	//   А теперь прибавим все прямоугольные куски, которые помещаются в ячейку
	ip_t ch_pos(sx.x + 1, sx.x + 2); //   - координаты канала
	for (int j = sx.x + 1; j < sb.x + 1; j++) {
		hx = ch_pos.y <= 0 ? HX * ch_pos.y : hx = OX[ch_pos.y];
		if (j == sb.x) hx = b;
		gx = ch_pos.y <= 0 ? HX * ch_pos.x : OX[ch_pos.x];
		result += integrate_rectangle_one_cell(bv.y, uv.y, gx, hx, ch_pos, sy);
		ch_pos.x += 1;
		ch_pos.y = ch_pos.x + 1;
	}
	return result;
}

// определим целочисленные индексы квадратов в которых лежат верхняя и нижняя точки треугольника
// sx = (x,y) координаты квадрата в которой лежит нижняя точка
// sy = (x,y) координаты квадрата в которой лежит верхняя точка
// в случае успешной проверки, k = будет  угловой коэфициент прямой

double integrate_right_triangle_bottom_left(const dp_t& bv, const dp_t& uv) {
	double k = 0;
	if (!try_get_slope_ratio(bv, uv, k)) return k;

	//   -  Index of current square by Ox and Oy axes. 
	ip_t sx, sy;
	sx.x = static_cast<int> ((bv.x - FLT_MIN) / HX);
	if (bv.x - FLT_MIN <= 0) sx.x -= 1;
	sx.y = sx.x + 1;
	sy.x = static_cast<int> ((bv.y + FLT_MIN) / HY);
	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
	sy.y = sy.x + 1;

	ip_t ib(sx.x, sx.x + 1); //   -  Index of right boundary.   

	double result = 0;
	int curr_i = 0, next_i = 0;
	dp_t curr = bv, next;
	while (true) {
		//TODO: sx.x и sx.y должны быть положительными всегда? Кажется для sx.x это всегда верно...
		double slope = sx.y >= 0 ? OY[sy.y] - curr.y : fabs(HY * sy.y - curr.y);
		slope /= sx.x >= 0 ? curr.x - OX[sx.x] : fabs(curr.x - HX * sx.x);
		if (slope <= k) { //   Intersection with straight line parallel Ox axis.        
			next_i = 1;
			next.y = sy.y >= 0 ? OY[sy.y] : HY * sy.y;
			next.x = curr.x - (next.y - curr.y) / k;
		}
		else { //   Intersection with straight line parallel Oy axis.            
			next_i = 2;
			next.x = sx.x >= 0 ? OX[sx.x] : HX * sx.x;
			next.y = curr.y - k * (next.x - curr.x);
		}
		if (next.x < (uv.x + FLT_MIN)) {
			// сюда попадаем и в случае когда треугольник полностью в одной ячейке лежит
			// и в случае когда прошлись по всем точкам...
			next_i = 0;
			next = uv;
			result += integrate_chanel_slant_left(curr, next, curr_i, next_i, sx, sy, bv.x, ib);
			break;
		}
		result += integrate_chanel_slant_left(curr, next, curr_i, next_i, sx, sy, bv.x, ib);
		switch (next_i) {
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

double integrate_right_triangle_bottom_right(const dp_t& bv, const dp_t& uv) {
	double k = 0;
	if (!try_get_slope_ratio(bv, uv, k)) return k;

	ip_t sx, sy;
	sx.x = static_cast<int> ((bv.x + FLT_MIN) / HX);
	if (bv.x + FLT_MIN <= 0) sx.x -= 1;
	sx.y = sx.x + 1;
	sy.x = static_cast<int> ((bv.y + FLT_MIN) / HY);
	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
	sy.y = sy.x + 1;

	ip_t ib(sx.x, sx.x + 1);
	double result = 0;
	int curr_i = 0, next_i = 0;
	dp_t curr = bv, next;
	while (true) {
		double slope = sy.y >= 0 ? fabs(OY[sy.y] - curr.y) : fabs(HY * sy.y - curr.y);
		slope /= sx.y >= 0 ? fabs(OX[sx.y] - curr.x) : fabs(HX * sx.y - curr.x);
		if (slope <= k) {//   Intersection with straight line parallel Ox axis.            
			next_i = 1;
			next.y = sy.y >= 0 ? OY[sy.y] : HY * sy.y;
			next.x = bv.x + (next.y - bv.y) / k;
		}
		else {//   Intersection with straight line parallel Oy axis.
			next_i = 2;
			next.x = sx.y >= 0 ? OX[sx.y] : HX * sx.y;
			next.y = bv.y + k * (next.x - bv.x);
		}
		if (next.x > (uv.x - FLT_MIN)) {
			next = uv;
			next_i = 0;
			result += integrate_chanel_slant_right(curr, next, curr_i, next_i, sx, bv.x, ib, sy);
			break;
		}
		result += integrate_chanel_slant_right(curr, next, curr_i, next_i, sx, bv.x, ib, sy);
		switch (next_i) {
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

double integrate_right_triangle_upper_left(const dp_t& bv, const dp_t& uv) {
	double k = 0;
	if (!try_get_slope_ratio(bv, uv, k)) return k;

	ip_t sx, sy, ib;
	sx.x = static_cast<int> ((bv.x + FLT_MIN) / HX); //   -  If bv.x is in grid edge I want it will be in the right side.
	if (bv.x + FLT_MIN <= 0) sx.x -= 1;
	sx.y = sx.x + 1;
	sy.x = static_cast<int> ((bv.y + FLT_MIN) / HY); //   -  If bv.y is in grid edge I want it will be in the upper square.
	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
	sy.y = sy.x + 1;
	ib.x = static_cast<int> ((uv.x - FLT_MIN) / HY); //   -  If uv.x is in grid edge I want it will be in the left side.
	if (uv.x - FLT_MIN <= 0) ib.x -= 1;
	ib.y = ib.x + 1;

	double result = 0;
	int curr_i = 0, next_i = 0;
	dp_t curr = bv, next;
	while (true) {
		double slope = sy.y >= 0 ? OY[sy.y] - curr.y : fabs(HY * sy.y - curr.y);
		slope /= sx.y >= 0 ? OX[sx.y] - curr.x : fabs(HX * sx.y - curr.x);
		if (slope <= k) { //   intersection with straight line parallel Ox axis.
			next_i = 1;
			next.y = sy.y >= 0 ? OY[sy.y] : HY * sy.y;
			next.x = bv.x + (next.y - bv.y) / k;
		}
		else {//   intersection with straight line parallel Oy axis.            
			next_i = 2;
			next.x = sx.y >= 0 ? OX[sx.y] : HX * sx.y;
			next.y = bv.y + k * (next.x - bv.x);
		}
		if (next.x > (uv.x - FLT_MIN)) {
			next_i = 0;
			next = uv;
			result += integrate_chanel_slant_left(curr, next, curr_i, next_i, sx, sy, uv.x, ib);
			break;
		}
		result += integrate_chanel_slant_left(curr, next, curr_i, next_i, sx, sy, uv.x, ib);

		switch (next_i) {
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

double integrate_right_triangle_upper_right(const dp_t& bv, const dp_t& uv) {
	double k = 0;
	if (!try_get_slope_ratio(bv, uv, k)) return k;

	ip_t sx, sy, ib;
	sx.x = static_cast<int> ((bv.x - FLT_MIN) / HX); //   -  If bv.x is in grid edge I want it will be between in the left side.
	if (bv.x - FLT_MIN <= 0) sx.x -= 1;
	sx.y = sx.x + 1;
	sy.x = static_cast<int> ((bv.y + FLT_MIN) / HY); //   -  If bv.y is in grid edge I want it will be in the upper side.
	if (bv.y + FLT_MIN <= 0) sy.x -= 1;
	sy.y = sy.x + 1;
	ib.x = static_cast<int> ((uv.x + FLT_MIN) / HX);
	if (uv.x + FLT_MIN <= 0) ib.x -= 1;
	ib.y = ib.x + 1;

	double result = 0;
	int curr_i = 0, next_i = 0;
	dp_t curr = bv, next;
	while (true) {
		double slope = sy.y >= 0 ? fabs(OY[sy.y] - curr.y) : fabs(HY * sy.y - curr.y);
		slope /= sx.x >= 0 ? fabs(curr.x - OX[sx.x]) : fabs(curr.x - HX * sx.x);
		if (slope <= k) { //   Intersection with straight line parallel Ox axis.
			next_i = 1;
			next.y = sy.y >= 0 ? OY[sy.y] : HY * sy.y;
			next.x = bv.x - (next.y - bv.y) / k;
		}
		else { //   Intersection with straight line parallel Oy axis.
			next_i = 2;
			next.x = sx.x >= 0 ? OX[sx.x] : HX * sx.x;
			next.y = bv.y - k * (next.x - bv.x);
		}
		if (next.x < uv.x + FLT_MIN) {
			next_i = 0;
			next = uv;
			result += integrate_chanel_slant_right(curr, next, curr_i, next_i, sx, uv.x, ib, sy);
			break;
		}
		result += integrate_chanel_slant_right(curr, next, curr_i, next_i, sx, uv.x, ib, sy);
		switch (next_i) {
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

double integrate_bottom_triangle(const dp_t& l, const dp_t& r, const dp_t& m) {
	double result = 0;
	if (m.x == l.x) {
		result = integrate_right_triangle_bottom_right(m, r);
	}
	else if (m.x == r.x) {
		result = integrate_right_triangle_bottom_left(m, l);
	}
	else if (m.x < l.x) {
		result = integrate_right_triangle_bottom_right(m, r);
		result -= integrate_right_triangle_bottom_right(m, l);
	}
	else if (m.x > l.x && m.x < r.x) {
		result = integrate_right_triangle_bottom_left(m, l);
		result += integrate_right_triangle_bottom_right(m, r);
	}
	else if (m.x > r.x) {
		result = integrate_right_triangle_bottom_left(m, l);
		result -= integrate_right_triangle_bottom_left(m, r);
	}
	return result;
}

double integrate_upper_triangle(const dp_t& l, const dp_t& r, const dp_t& m) {
	double result = 0;
	if (m.x == l.x) {
		result = integrate_right_triangle_upper_right(r, m);
	}
	else if (m.x == r.x) {
		result = integrate_right_triangle_upper_left(l, m);
	}
	else if (m.x < l.x) {
		result = integrate_right_triangle_upper_right(r, m);
		result -= integrate_right_triangle_upper_right(l, m);
	}
	else if (m.x > l.x && m.x < r.x) {
		result = integrate_right_triangle_upper_left(l, m);
		result += integrate_right_triangle_upper_right(r, m);
	}
	else if (m.x > r.x) {
		result = integrate_right_triangle_upper_left(l, m);
		result -= integrate_right_triangle_upper_left(r, m);
	}
	return result;
}

double integrate_uniform_triangle_wall(const dp_t& a, const dp_t& b, const dp_t& c) {
	return 0;
}

double integrate_uniform_triangle(const dp_t& x, dp_t& y, const dp_t& z) {
	//   a * x  +  b * y  = c.
	double a = z.y - x.y;
	if (fabs(a) < FLT_MIN) return FLT_MIN;
	double b = x.x - z.x;
	double c = b * x.y + a * x.x;
	dp_t ip((c - b * y.y) / a, y.y);

	//   Возможны 2 случая расположения точки перечеения относительно средней
	//   слева или справа.
	//   есди средняя точка справа от точки пересечения
	//   обменяем местами  X координаты, чтобы использовать один код для расчета
	if (y.x >= ip.x) {
		double tx = y.x;
		y.x = ip.x;
		ip.x = tx;
	}

	return integrate_upper_triangle(y, ip, z)
		+ integrate_bottom_triangle(y, ip, x);
}

quad_type get_coordinates_on_prev_layer(int ix, int iy,
	dp_t& alpha, dp_t& beta, dp_t& gamma, dp_t& theta) {
	//   1 First of all let's compute coordinates of square vertexes.
	//  OX:
	if (ix == 0) {
		alpha.x = OX[ix];
		beta.x = (OX[ix] + OX[ix + 1]) * 0.5;
		gamma.x = (OX[ix] + OX[ix + 1]) * 0.5;
		theta.x = OX[ix];
	}
	else if (ix == OX_LEN) {
		alpha.x = (OX[ix - 1] + OX[ix]) * 0.5;
		beta.x = OX[ix];
		gamma.x = OX[ix];
		theta.x = (OX[ix - 1] + OX[ix]) * 0.5;
	}
	else {
		alpha.x = (OX[ix - 1] + OX[ix]) * 0.5;
		beta.x = (OX[ix + 1] + OX[ix]) * 0.5;
		gamma.x = (OX[ix + 1] + OX[ix]) * 0.5;
		theta.x = (OX[ix - 1] + OX[ix]) * 0.5;
	}

	//  OY:
	if (iy == 0) {
		alpha.y = OY[iy];
		beta.y = OY[iy];
		gamma.y = (OY[iy] + OY[iy + 1]) * 0.5;
		theta.y = (OY[iy] + OY[iy + 1]) * 0.5;
	}
	else if (iy == OY_LEN) {
		alpha.y = (OY[iy] + OY[iy - 1]) * 0.5;
		beta.y = (OY[iy] + OY[iy - 1]) * 0.5;
		gamma.y = OY[iy];
		theta.y = OY[iy];
	}
	else {
		alpha.y = (OY[iy] + OY[iy - 1]) * 0.5;
		beta.y = (OY[iy] + OY[iy - 1]) * 0.5;
		gamma.y = (OY[iy] + OY[iy + 1]) * 0.5;
		theta.y = (OY[iy] + OY[iy + 1]) * 0.5;
	}

	double u, v;

	// Now let's compute new coordinates on the previous time level of alpha, beta, gamma, theta points.
	u = func_u(alpha);
	v = func_v(TAU_TL, alpha);
	alpha.x -= TAU * u;
	alpha.y -= TAU * v;

	u = func_u(beta);
	v = func_v(TAU_TL, beta);
	beta.x -= TAU * u;
	beta.y -= TAU * v;

	u = func_u(gamma);
	v = func_v(TAU_TL, gamma);
	gamma.x -= TAU * u;
	gamma.y -= TAU * v;

	u = func_u(theta);
	v = func_v(TAU_TL, theta);
	theta.x -= TAU * u;
	theta.y -= TAU * v;

	dp_t intersection = get_intersection_point(alpha, beta, gamma, theta);
	if ((beta.y - intersection.y) * (theta.y - intersection.y) > 0) return pseudo; // ??
	if ((alpha.x - intersection.x) * (gamma.x - intersection.x) > 0) return pseudo; // ??
	double product = get_vector_product(alpha, beta, theta); // ?
	if (product < 0) return pseudo;

	// значит что точка улетела за левую границу
	if (theta.x < 0 || theta.y < 0 || beta.x < 0 || beta.y < 0 || gamma.x < 0 ||
		gamma.y < 0 || alpha.x < 0 || alpha.y < 0) {
		return normal;
		//return wall;
	}
	return normal;
}

// Type of quadrangle: 0 - pseudo; 1 - convex; 2 - concave;

quad_type get_quadrangle_type(int ix, int iy,
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

double integrate(int ix, int iy) {
	dp_t a1, b1, c1, a2, b2, c2;
	quad_type type = get_quadrangle_type(ix, iy, a1, b1, c1, a2, b2, c2);
	if (type != normal && type != wall) {
		return -1;
	}

	// чтобы правилно отработала процедура интегрирования
	// точки должны идти в порядке возрастания y координаты
	sort_by_y(a1, b1, c1);
	sort_by_y(a2, b2, c2);

	// check the type of triangle to select appropriate computation method
	double result = 0;
	switch (type) {
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

double get_norm_of_error(double* density, double ts_count_mul_steps) {
	double r = 0;
	for (int k = 1; k < OY_LEN; ++k) {
		for (int j = 1; j < OX_LEN; ++j) {
			r += fabs(analytical_solution(ts_count_mul_steps, OX[j], OY[k])
				- density[(OY_LEN + 1) * k + j]);
		}
	}
	return HX * HY * r;
}

void solve(double* density) {
	PREV_DENSITY = new double[XY_LEN];
	for (int j = 0; j < OY_LEN + 1; j++) {
		for (int i = 0; i < OX_LEN + 1; i++) {
			PREV_DENSITY[(OX_LEN + 1) * j + i] = analytical_solution(0, OX[i], OY[j]);
		}
	}

	for (TL = 1; TL <= TIME_STEP_CNT; TL++) {
		TAU_TL = TAU * TL;
		TAU_TL_1 = TAU * (TL - 1);
		for (int i = 0; i <= OX_LEN; i++) {
			density[i] = init_side(OX[i], BB, TAU_TL);
			density[(OX_LEN + 1) * OY_LEN + i] = init_side(OX[i], UB, TAU_TL);
		}

		for (int i = 0; i <= OY_LEN; i++) {
			density[(OX_LEN + 1) * i] = init_side(LB, OY[i], TAU_TL);
			density[(OX_LEN + 1) * i + OX_LEN] = init_side(RB, OY[i], TAU_TL);
		}

		for (int i = 1; i < OY_LEN; i++) {
			for (int j = 1; j < OX_LEN; j++) {
				density[(OX_LEN + 1) * i + j] = integrate(j, i) * INVERTED_HX_HY;
				density[(OX_LEN + 1) * i + j] += TAU * func_f(OX[j], OY[i]);
			}
		}
		memcpy(PREV_DENSITY, density, XY_LEN * sizeof(double));
	}
	delete[] PREV_DENSITY;
}


double* compute_density(double b, double lb, double rb, double bb, double ub,
	double tau, int time_step_count, int ox_length, int oy_length, double &norm) {
	init(b, lb, rb, bb, ub, tau, time_step_count, ox_length, oy_length);
	double* density = new double[XY_LEN];
	print_params(B, LB, RB, BB, UB, TAU, TIME_STEP_CNT, OX_LEN, OY_LEN);
	solve(density);
	norm = get_norm_of_error(density, TIME_STEP_CNT * TAU);
	printf("%d x %d wall count = %d\n", OX_LEN + 1, OY_LEN + 1, TMP_WALL_CNT);
	clean();
	return density;
}