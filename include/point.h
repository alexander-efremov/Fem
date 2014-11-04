
#ifndef MISC_H
#define	MISC_H


struct point_t {
    double x;
    double y;

    double operator[](int i) {
        return (i == 0) ? x : y;
    }

    point_t(point_t&& other)
    : x(std::move(other.x)), y(std::move(other.y)) {
    }

    point_t(double x = 0, double y = 0)
    : x(x), y(y) {
    }

    point_t& operator=(const point_t& p) {
        if (this != &p) {
            x = p.x;
            y = p.y;
        }

        return *this;
    }

    point_t operator+(const point_t& p) const {
        return point_t(p.x + x, p.y + y);
    }

    bool operator==(const point_t& p) const {
        return (x == p.x && y == p.y);
    }

    int operator<(const point_t &p) {
        return ((x < p.x) || ((x == p.x) && (y < p.y)));
    }

    int operator>(const point_t &p) {
        return ((x > p.x) || ((x == p.x) && (y > p.y)));
    }
};

void sort_by_y(point_t& x, point_t& y, point_t& z) {
    if (x.y < y.y) {
        if (z.y < x.y) 
            std::swap(x, z);
    } else {
        if (y.y < z.y) 
            std::swap(x, y);
        else 
            std::swap(x, z);
    }
    if (z.y < y.y) std::swap(y, z);
}

point_t get_intersection_point(point_t *alpha, point_t *beta, point_t *gamma, point_t *theta) {
    point_t result;
    point_t alpha_to_gamma;
    point_t beta_to_theta;
    double a_1LC, b_1LC, c_1LC;
    double a_2LC, b_2LC, c_2LC;
    alpha_to_gamma.x = gamma->x - alpha->x;
    alpha_to_gamma.y = gamma->y - alpha->y;
    a_1LC = alpha_to_gamma.y;
    b_1LC = -alpha_to_gamma.x;
    c_1LC = alpha_to_gamma.y * alpha->x - alpha_to_gamma.x * alpha->y;
    beta_to_theta.x = theta->x - beta->x;
    beta_to_theta.y = theta->y - beta->y;
    a_2LC = beta_to_theta.y;
    b_2LC = -beta_to_theta.x;
    c_2LC = beta_to_theta.y * beta->x - beta_to_theta.x * beta->y;
    result.x = (b_1LC * c_2LC - b_2LC * c_1LC) / (b_1LC * a_2LC - b_2LC * a_1LC);
    result.y = (a_1LC * c_2LC - a_2LC * c_1LC) / (-b_1LC * a_2LC + b_2LC * a_1LC);
    return result;
}

inline void ptcpy(point_t *dst, const point_t *source) {
    dst->x = source->x;
    dst->y = source->y;
}

double get_vector_product(point_t *alpha, point_t *beta, point_t *theta) {
    point_t alpha_to_beta;
    point_t alpha_to_theta;
    alpha_to_beta.x = beta->x - alpha->x;
    alpha_to_beta.y = beta->y - alpha->y;
    alpha_to_theta.x = theta->x - alpha->x;
    alpha_to_theta.y = theta->y - alpha->y;
    return alpha_to_beta.x * alpha_to_theta.y - alpha_to_beta.y * alpha_to_theta.x;
}

#ifdef	__cplusplus
extern "C" {
#endif




#ifdef	__cplusplus
}
#endif

#endif	/* MISC_H */

