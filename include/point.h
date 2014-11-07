#ifndef MISC_H
#define	MISC_H

static const double MIN_VALUE = 1.e-12;

struct dp_t {
    double x;
    double y;

    double operator[](int i) {
        return (i == 0) ? x : y;
    }

    dp_t(double x, double y)
    : x(x), y(y) {
    }
    
     dp_t()
    : x(0), y(0) {
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
    int x;
    int y;

    double operator[](int i) {
        return (i == 0) ? x : y;
    }
    
    ip_t(int x, int y)
    : x(x), y(y) {
    }

    ip_t()
    : x(0), y(0) {
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
    
    ip_t& operator+=(int a) {
        x += a;
        y += a;
        return *this;
    }
    
    ip_t& operator-=(int a) {
        x -= a;
        y -= a;
        return *this;
    }
    
    ip_t& operator-(int a) {
        x -= a;
        y -= a;
        return *this;
    }
    
    ip_t& operator+(int a) {
        x += a;
        y += a;
        return *this;
    }

    bool operator==(const ip_t& p) const {
        return (x == p.x && y == p.y);
    }

    int operator<(const ip_t& p) {
        return ((x < p.x) || ((x == p.x) && (y < p.y)));
    }

    int operator>(const ip_t& p) {
        return ((x > p.x) || ((x == p.x) && (y > p.y)));
    }
};

#endif /* MISC_H */