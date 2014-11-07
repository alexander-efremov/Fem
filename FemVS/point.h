#ifndef MISC_H
#define	MISC_H

static const double MIN_VALUE = 1.e-12;

struct dp_t {
	double x;
	double y;

	double operator[](int i) {
		return (i == 0) ? x : y;
	}

	dp_t(double x = 0, double y = 0)
		: x(x), y(y) {
	}

	dp_t& operator=(const dp_t& p) {
		if (this != &p) {
			x = p.x;
			y = p.y;
		}

		return *this;
	}

	dp_t operator+(const dp_t& p) const {
		return dp_t(p.x + x, p.y + y);
	}

	bool operator==(const dp_t& p) const {
		return (x == p.x && y == p.y);
	}

	int operator<(const dp_t& p) {
		return ((x < p.x) || ((x == p.x) && (y < p.y)));
	}

	int operator>(const dp_t& p) {
		return ((x > p.x) || ((x == p.x) && (y > p.y)));
	}
};

struct ip_t {
	short x;
	short y;

	double operator[](short i) {
		return (i == 0) ? x : y;
	}

	ip_t(short x = 0, short y = 0)
		: x(x), y(y) {
	}

	ip_t& operator=(const ip_t& p) {
		if (this != &p) {
			x = p.x;
			y = p.y;
		}

		return *this;
	}

	ip_t operator+(const ip_t& p) const {
		return ip_t(p.x + x, p.y + y);
	}

	ip_t& operator+=(short a) {
		x += a;
		y += a;
		return *this;
	}

	ip_t& operator-=(short a) {
		x -= a;
		y -= a;
		return *this;
	}

	ip_t& operator-(short a) {
		x -= a;
		y -= a;
		return *this;
	}

	ip_t& operator+(short a) {
		x += a;
		y += a;
		return *this;
	}

	bool operator==(const ip_t& p) const {
		return (x == p.x && y == p.y);
	}

	short operator<(const ip_t& p) {
		return ((x < p.x) || ((x == p.x) && (y < p.y)));
	}

	short operator>(const ip_t& p) {
		return ((x > p.x) || ((x == p.x) && (y > p.y)));
	}
};

inline void sort_by_y(dp_t& x, dp_t& y, dp_t& z) {
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

inline bool try_get_slope_ratio(const dp_t &bv, const dp_t &uv, double &value) {
	if (fabs(bv.x - uv.x) < MIN_VALUE) {
		return false;
	}
	value = fabs((uv.y - bv.y) / (uv.x - bv.x)); // угловой коэффициент прямой
	if (value < MIN_VALUE) {
		return false;
	}
	return true;
}

inline dp_t get_intersection_point(const dp_t& alpha, const dp_t& beta, const dp_t& gamma, const dp_t& theta) {
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

inline double get_vector_product(const dp_t& alpha, const dp_t beta, const dp_t theta) {
	dp_t alpha_to_beta;
	dp_t alpha_to_theta;
	alpha_to_beta.x = beta.x - alpha.x;
	alpha_to_beta.y = beta.y - alpha.y;
	alpha_to_theta.x = theta.x - alpha.x;
	alpha_to_theta.y = theta.y - alpha.y;
	return alpha_to_beta.x * alpha_to_theta.y - alpha_to_beta.y * alpha_to_theta.x;
}

#endif /* MISC_H */