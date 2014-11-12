#ifndef POINT_H
#define	POINT_H

static const double MIN_VALUE = 1.e-12;

struct dp_t {
	double x;
	double y;	

	dp_t(double x, double y)
		: x(x), y(y) {
	}

	dp_t()
		: x(0), y(0) {
	}

	inline double operator[](int i) {
		return (i == 0) ? x : y;
	}

	inline dp_t& operator=(const dp_t& p) {
		if (this != &p) {
			x = p.x;
			y = p.y;
		}

		return *this;
	}

	inline dp_t operator+(const dp_t& p) const {
		return dp_t(p.x + x, p.y + y);
	}

	inline bool operator==(const dp_t& p) const {
		return (x == p.x && y == p.y);
	}

	inline int operator<(const dp_t& p) {
		return ((x < p.x) || ((x == p.x) && (y < p.y)));
	}

	inline int operator>(const dp_t& p) {
		return ((x > p.x) || ((x == p.x) && (y > p.y)));
	}
};

struct ip_t {
	int x;
	int y;


	ip_t(int x, int y)
		: x(x), y(y) {
	}

	ip_t()
		: x(0), y(0) {
	}


	inline double operator[](int i) {
		return (i == 0) ? x : y;
	}

	ip_t& operator=(const ip_t& p) {
		if (this != &p) {
			x = p.x;
			y = p.y;
		}

		return *this;
	}

	inline ip_t operator+(const ip_t& p) const {
		return ip_t(p.x + x, p.y + y);
	}

	inline ip_t& operator+=(int a) {
		x += a;
		y += a;
		return *this;
	}

	inline ip_t& operator-=(int a) {
		x -= a;
		y -= a;
		return *this;
	}

	inline ip_t& operator-(int a) {
		x -= a;
		y -= a;
		return *this;
	}

	inline ip_t& operator+(int a) {
		x += a;
		y += a;
		return *this;
	}

	inline bool operator==(const ip_t& p) const {
		return (x == p.x && y == p.y);
	}

	inline int operator<(const ip_t& p) {
		return ((x < p.x) || ((x == p.x) && (y < p.y)));
	}

	inline int operator>(const ip_t& p) {
		return ((x > p.x) || ((x == p.x) && (y > p.y)));
	}
};

#endif /* POINT_H */