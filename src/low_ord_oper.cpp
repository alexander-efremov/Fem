#include "common.h"
#include "point.h"
#include "utils.h"

static int TMP_WALL_CNT = 0;
static const double _PI = 3.14159265358979323846264338327;
static const double _MIN_VALUE = 1.e-12;

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

double analytical_solution(double t, double x, double y) {
    return 1.1 + sin(t * x * y);
}

double init_side(double x, double y, double t, bound bound) {
    switch (bound) {
        case up:
            return analytical_solution(t, x, y);
        case bottom:
            return analytical_solution(t, x, y);
        case left:
            return analytical_solution(t, x, y);
        case right:
            return analytical_solution(t, x, y);
    }
    return 0;
}

double integrate_first_type(double py,
        double qy,
        double gx,
        double hx,
        double a,
        double b) {
    double integ = (hx - a) * (hx - a) - (gx - a) * (gx - a);
    integ *= (qy - b) * (qy - b) - (py - b) * (py - b);
    return integ / 4;
}

double integrate_second_type(double py,
        double qy,
        double alpha,
        double a,
        double b,
        double beta) {
    double tmp, integ;
    tmp = (qy - alpha) * (a * qy + b - beta) * (a * qy + b - beta) * (a * qy + b - beta);
    tmp -= (py - alpha) * (a * py + b - beta) * (a * py + b - beta) * (a * py + b - beta);
    integ = tmp / (3 * a);
    tmp = (a * qy + b - beta) * (a * qy + b - beta) * (a * qy + b - beta) * (a * qy + b - beta);
    tmp -= (a * py + b - beta) * (a * py + b - beta) * (a * py + b - beta) * (a * py + b - beta);
    return integ - tmp / (12 * a * a);
}

double integrate_rectangle_one_cell(double Py,
        double Qy,
        double Gx,
        double Hx,
        int tl,
        const ip_t &indCurSqOx, //   -  Index of current square by Ox axis.
        const ip_t &indCurSqOy,
        const double* ox,
        const double* oy,
        double* density) {    
    double result;
    double tmp;
    double rho[2][2];
    double t = TAU * (tl - 1.);
    double x, y;
    if (indCurSqOx.x >= 0 && indCurSqOy.x >= 0) {
        rho[0][0] = density[(OX_LEN + 1) * indCurSqOy.x + indCurSqOx.x];
        rho[0][1] = density[(OX_LEN + 1) * indCurSqOy.y + indCurSqOx.x];
        rho[1][0] = density[(OX_LEN + 1) * indCurSqOy.x + indCurSqOx.y];
        rho[1][1] = density[(OX_LEN + 1) * indCurSqOy.y + indCurSqOx.y];
    } else {
        // TODO: убрать потому что это неверно (надо расчитывать граничные условия)
        x = indCurSqOx.x * HX;
        y = indCurSqOy.x * HY;
        rho[0][0] = analytical_solution(t, x, y);
        x = indCurSqOx.x * HX;
        y = indCurSqOy.y * HY;
        rho[0][1] = analytical_solution(t, x, y);
        x = indCurSqOx.y * HX;
        y = indCurSqOy.x * HY;
        rho[1][0] = analytical_solution(t, x, y);
        x = indCurSqOx.y * HX;
        y = indCurSqOy.y * HY;
        rho[1][1] = analytical_solution(t, x, y);

        TMP_WALL_CNT++;
    }

    if (indCurSqOx.y >= 0 && indCurSqOy.y >= 0) {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, ox[indCurSqOx.y], oy[indCurSqOy.y]);
    } else {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, HX * indCurSqOx.y, HY * indCurSqOy.y);
    }
    tmp = tmp / HX / HY;
    result = tmp * rho[0][0]; //   rhoInPrevTL[ indCurSqOx.x ][ indCurSqOy.x ];
    if (indCurSqOx.x >= 0 && indCurSqOy.y >= 0) {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, ox[indCurSqOx.x], oy[indCurSqOy.y]);
    } else {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, HX * indCurSqOx.x, HY * indCurSqOy.y);
    }
    tmp = tmp / HX / HY;
    result = result - tmp * rho[1][0]; //   rhoInPrevTL[ sx.y ][ indCurSqOy.x ];
    if (indCurSqOx.y >= 0 && indCurSqOy.x >= 0) {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, ox[indCurSqOx.y], oy[indCurSqOy.x]);
    } else {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, HX * indCurSqOx.y, HY * indCurSqOy.x);
    }
    tmp = tmp / HX / HY;
    result -= tmp * rho[0][1]; //   rhoInPrevTL[ indCurSqOx.x ][ indCurSqOy.y ];
    if (indCurSqOx.x >= 0 && indCurSqOy.x >= 0) {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, ox[indCurSqOx.x], oy[indCurSqOy.x]);
    } else {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, HX * indCurSqOx.x, HY * indCurSqOy.x);
    }
    tmp = tmp / HX / HY;
    return result + tmp * rho[1][1]; //   rhoInPrevTL[ sx.y ][ indCurSqOy.y ];
}

double integrate_triangle_left_one_cell(
        double Py,
        double Qy,
        double a_SL,
        double b_SL,
        double Hx,
        int tl,
        const ip_t &sox, //   -  Index of current square by Ox axis.
        const ip_t &soy, //   -  Index of current square by Oy axis.        
        const double* ox,
        const double* oy,
        double* density) {
    double result;
    double tmp, bufInteg_D;
    double rho[2][2];
    
    if (sox.x >= 0 && sox.y <= OX_LEN && soy.x >= 0 && soy.y <= OY_LEN) {        
        rho[0][0] = density[(OX_LEN + 1) * soy.x + sox.x];
        rho[0][1] = density[(OX_LEN + 1) * soy.y + sox.x];
        rho[1][0] = density[(OX_LEN + 1) * soy.x + sox.y];
        rho[1][1] = density[(OX_LEN + 1) * soy.y + sox.y];        
    }    
    else {
        // TODO: убрать потому что это неверно (надо расчитывать граничные условия)
        // норма должна уменьшиться
        rho[0][0] = analytical_solution(TAU * (tl - 1), sox.x * HX, soy.x * HY);
        rho[0][1] = analytical_solution(TAU * (tl - 1), sox.x * HX, soy.y * HY);
        rho[1][0] = analytical_solution(TAU * (tl - 1), sox.y * HX, soy.x * HY);
        rho[1][1] = analytical_solution(TAU * (tl - 1), sox.y * HX, soy.y * HY);
        TMP_WALL_CNT++;
    }

    //   1.
    tmp = (Qy - oy[soy.y]) * (Qy - oy[soy.y]) - (Py - oy[soy.y]) * (Py - oy[soy.y]);
    if ((sox.y >= 0) && (soy.y >= 0)) {
        tmp = tmp * (Hx - ox[sox.y]) * (Hx - ox[sox.y]) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, oy[soy.y], a_SL, b_SL, ox[sox.y]);
    } else {
        tmp = tmp * (Hx - HX * sox.y) * (Hx - HX * sox.y) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, HY * soy.y, a_SL, b_SL, HX * sox.y);
    }
    tmp -= bufInteg_D / 2.;
    result = tmp * rho[0][0] / HX / HY;

    //   2.
    tmp = (Qy - oy[soy.y]) * (Qy - oy[soy.y]) - (Py - oy[soy.y]) * (Py - oy[soy.y]);
    if ((sox.x >= 0) && (soy.y >= 0)) {
        tmp = -1. * tmp * (Hx - ox[sox.x]) * (Hx - ox[sox.x]) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, oy[soy.y], a_SL, b_SL, ox[sox.x]);
    } else {
        tmp = -1. * tmp * (Hx - HX * sox.x) * (Hx - HX * sox.x) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, HY * soy.y, a_SL, b_SL, HX * sox.x);
    }
    tmp += bufInteg_D / 2.;
    result += tmp * rho[1][0] / HX / HY;

    //   3.
    tmp = (Qy - oy[soy.x]) * (Qy - oy[soy.x]) - (Py - oy[soy.x]) * (Py - oy[soy.x]);
    if ((sox.y >= 0) && (soy.x >= 0)) {
        tmp = -1. * tmp * (Hx - ox[sox.y]) * (Hx - ox[sox.y]) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, oy[soy.x], a_SL, b_SL, ox[sox.y]);
    } else {
        tmp = -1. * tmp * (Hx - HX * sox.y) * (Hx - HX * sox.y) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, HY * soy.x, a_SL, b_SL, HX * sox.y);
    }
    tmp += bufInteg_D / 2.;
    result += tmp * rho[0][1] / HX / HY;

    //   4.
    tmp = (Qy - oy[soy.x]) * (Qy - oy[soy.x]) - (Py - oy[soy.x]) * (Py - oy[soy.x]);
    if ((sox.x >= 0) && (soy.x >= 0)) {
        tmp = tmp * (Hx - ox[sox.x]) * (Hx - ox[sox.x]) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, oy[soy.x], a_SL, b_SL, ox[sox.x]);
    } else {
        tmp = tmp * (Hx - HX * sox.x) * (Hx - HX * sox.x) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, HY * soy.x, a_SL, b_SL, HX * sox.x);
    }
    tmp -= bufInteg_D / 2.;
    result += tmp * rho[1][1] / HX / HY;

    return result;
}

double integrate_triangle_right_one_cell(double py,
        double qy,
        double a_SL,
        double b_SL,
        double gx,
        int tl,
        const ip_t &indCurSqOx, //   -  Index of current square by Ox axis.
        const ip_t &indCurSqOy, //   -  Index of current square by Oy axis.
        const double* ox,
        const double* oy,
        double* density) {
    return -1. * integrate_triangle_left_one_cell(
            py, qy,
            a_SL, b_SL, gx,
            tl,
            indCurSqOx, //   -  Index of current square by Ox axis.
            indCurSqOy, //   -  Index of current square by Oy axis.            
            ox, oy,
            density);
}

double integrate_chanel_slant_right(int tl,
        const dp_t& bv, int wTrPCI, //   -  Where travel point current (botton vertex) is.
        const dp_t& uv, int wTrPNI, //   -  Where travel point next (upper vertex) is.
        //
        const ip_t &indCurSqOx, //   -  Index by OX axis where bv and uv are.
        //
        double lb, const ip_t &indLB, //   -  Left boundary by Ox. Index by OX axis where lb is.
        //
        const ip_t &indCurSqOy, //   -  Index of current square by Oy axis.
        //
        const double* ox,
        const double* oy,
        double* density) {
    dp_t mv, rv; //   -  Middle and right vertices.
    int wMvI = 0; //   -  Where middle vertex is.
    ip_t indCurSqOxToCh; //   -  Indices of current square by Ox axis to be changed. Under which we want to integrate.
    double h = ox[1] - ox[0];
    double a_SL, b_SL; //   -  Coefficients of slant line: x = a_SL *y  +  b_SL.
    double Gx = 0, Hx = 0; //   -  Left boundary for each integration.
    double result = 0.;
    double tmp;

    //   Let's compute helpful values.
    if (uv.x <= bv.x) {
        mv = uv;
        wMvI = wTrPNI;
        rv = bv;
    } else {
        mv = bv;
        wMvI = wTrPCI;
        rv = uv;
    }

    if ((fabs(uv.y - bv.y)) <= _MIN_VALUE) {
        //   Computation is impossible. Too smale values. Let's return some approximate value.
        //   buf_D  =  (uv.y - bv.y)  *  ((uv.x + bv.x) /2.  -  lb) * rhoInPrevTL[ indCurSqOx.x ][ indCurSqOy.x ];
        return fabs(uv.y - bv.y); //   fabs(uv.y - bv.y);
    }


    //   First step: from "lb" to "masOX[ indCurSqOx.x ]" by iteration.
    //   integ  += fabs( mv[0] - lb) * fabs(uv.y - bv.y);

    indCurSqOxToCh.x = indLB.x;
    indCurSqOxToCh.y = indCurSqOxToCh.x + 1;

    for (int j = indLB.x; j < indCurSqOx.x; j++) {
        //   If this is first cell we should integrate under rectangle only.
        if (indCurSqOxToCh.x >= 0) {
            Gx = ox[indCurSqOxToCh.x];
            Hx = ox[indCurSqOxToCh.y];
        }


        if (indCurSqOxToCh.x < 0) {
            Gx = h * indCurSqOxToCh.x;
            Hx = h * indCurSqOxToCh.y;
        }

        if (j == indLB.x) {
            Gx = lb;
        }

        tmp = integrate_rectangle_one_cell(
                bv.y, //   -  double Py,
                uv.y, //   -  double Qy,
                Gx, //   -  double Gx,
                Hx, //   -  double Hx,
                //
                tl,
                indCurSqOxToCh, //   -  Index of current square by Ox axis.
                indCurSqOy, //   -  Index of current square by Oy axis.
                //
                ox,
                oy,
                density);

        result += tmp;

        indCurSqOxToCh.x += 1;
        indCurSqOxToCh.y = indCurSqOxToCh.x + 1;
    }

    //   Integration. Second step: under [ indCurSqOx.x; sx.y ] square.

    //   A. Under rectangle.
    if (wMvI == 1) {
        if (indCurSqOx.x == indLB.x) {
            Gx = lb;
        }

        if (indCurSqOx.x > indLB.x) {
            if (indCurSqOx.x >= 0) {
                Gx = ox[indCurSqOx.x];
            }

            if (indCurSqOx.x < 0) {
                Gx = h * indCurSqOx.x;
            }
        }

        tmp = integrate_rectangle_one_cell(bv.y, //   -  double Py,
                uv.y, //   -  double Qy,
                //
                Gx, //   -  double Gx,
                mv[0], //   -  double Hx,
                //
                tl,
                indCurSqOx, //   -  Index of current square by Ox axis.
                indCurSqOy, //   -  Index of current square by Oy axis.
                //
                ox,
                oy,
                density);

        result += tmp;
    }

    //   B. Under triangle.

    if (fabs(uv.y - bv.y) > _MIN_VALUE) {
        //   integ += fabs(uv.y - bv.y) * (rv[0] - mv[0]) /2.;
        //   Coefficients of slant line: x = a_SL *y  +  b_SL.
        a_SL = (uv.x - bv.x) / (uv.y - bv.y);
        b_SL = bv.x - a_SL * bv.y;


        //   Integration under one cell triangle.

        if (fabs(a_SL) > _MIN_VALUE) {
            tmp = integrate_triangle_right_one_cell(
                    bv.y, //   -  double Py,
                    uv.y, //   -  double Qy,
                    //
                    a_SL,
                    b_SL,
                    mv[0], //   -  double Gx,
                    //
                    tl,
                    indCurSqOx, //   -  Index of current square by Ox axis.
                    indCurSqOy, //   -  Index of current square by Oy axis.
                    //
                    ox,
                    oy,
                    density);

            result += tmp;
        }
    }

    return result;
}

double integrate_chanel_slant_left(int tl,
        const dp_t& bv, int curr_i, //   -  Where travel point current (bottom vertex) is.
        const dp_t& uv, int next_i, //   -  Where travel point next (upper vertex) is.
        //
        const ip_t &sox, //   -  Index by OX axis where bv and uv are.
        const ip_t &soy, //   -  Index of current square by Oy axis.
        //
        double rb, const ip_t &irb, //   -  Right boundary by Ox. Index by OX axis where rb is.        
        const double* ox,
        const double* oy,
        double* density) {

    if (fabs(uv.y - bv.y) <= _MIN_VALUE) return fabs(uv.y - bv.y);


    dp_t lv, mv; //   -  Left and middle vertices.
    int wMvI = 0; //   -  Where middle vertex is.    
    double a_SL, b_SL; //   -  Coefficients of slant line: x = a_SL *y  +  b_SL.
    double gx = 0, hx = 0; //   -  Left and right boundary for each integration.
    double result = 0.;

    //   Let's compute helpful values.
    if (uv.x <= bv.x) {
        lv = uv;
        mv = bv;
        wMvI = curr_i;
    } else {
        lv = bv;
        mv = uv;
        wMvI = next_i;
    }

    //   Integration. First step: under [ indCurSqOx.x; sx.y ] square.
    //   A. Under triangle.
    if (fabs(uv.y - bv.y) > _MIN_VALUE) {
        //   Coefficients of slant line: x = a_SL *y  +  b_SL.
        a_SL = (uv.x - bv.x) / (uv.y - bv.y);
        b_SL = bv.x - a_SL * bv.y;
        //   Integration under one cell triangle.
        if (fabs(a_SL) > _MIN_VALUE) {
            result += integrate_triangle_left_one_cell(bv.y, uv.y, a_SL, b_SL, mv[0],
                    tl, sox, soy, ox, oy, density);
        }
    }

    //   B. Under rectangle. Need to be checking.
    if (wMvI == 1) {
        if (sox.x == irb.x) {
            hx = rb;
        } else if (sox.x < irb.x) {
            hx = sox.y >= 0 ? ox[sox.y] : HX * sox.y;
        }
        result += integrate_rectangle_one_cell(bv.y, uv.y, mv[0], hx, tl,
                sox, soy,
                ox, oy, density);
    }

    //   Second step: from "masOX[ sx.y ]" to "rb" by iteration.
    ip_t sox_ch(sox.x + 1, sox.x + 2); //   -  Indices of current square by Ox axis to be changed. Under which we want to integrate.    

    for (int j = sox.x + 1; j < irb.x + 1; j++) {
        //   If this is first cell we should integrate under triangle only.
        gx = sox_ch.y <= 0 ? HX * sox_ch.x : ox[sox_ch.x];
        hx = sox_ch.y <= 0 ? HX * sox_ch.y : hx = ox[sox_ch.y];
        if (j == irb.x) hx = rb;
        result += integrate_rectangle_one_cell(bv.y, uv.y, gx, hx, tl,
                sox_ch, soy,
                ox, oy,
                density);
        sox_ch.x += 1;
        sox_ch.y = sox_ch.x + 1;
    }

    return result;
}

// определим целочисленные индексы квадратов в которых лежат верхняя и нижняя точки треугольника
// sx = (x,y) координаты квадрата в которой лежит нижняя точка
// sy = (x,y) координаты квадрата в которой лежит верхняя точка
// в случае успешной проверки, k = будет  угловой коэфициент прямой

double integrate_right_triangle_bottom_left(const dp_t& bv, const dp_t& uv, int tl,
        const double* ox, const double* oy, double* density) {
    double k = 0.;
    if (!try_get_slope_ratio(bv, uv, k)) return k;

    //   -  Index of current square by Ox and Oy axes. 
    ip_t sx, sy;
    sx.x = static_cast<short> ((bv.x - _MIN_VALUE) / HX);
    if (bv.x - _MIN_VALUE <= 0) sx.x -= 1;
    sx.y = sx.x + 1;
    sy.x = static_cast<short> ((bv.y + _MIN_VALUE) / HY);
    if (bv.y + _MIN_VALUE <= 0) sy.x -= 1;
    sy.y = sy.x + 1;

    ip_t ib(sx.x, sx.x + 1); //   -  Index of right boundary.   

    double result = 0.;
    short curr_i = 0, next_i = 0;
    dp_t curr = bv, next;
    while (true) {
        //TODO: sx.x и sx.y должны быть положительными всегда? Кажется для sx.x это всегда верно...
        double dx = sx.x >= 0 ? curr.x - ox[sx.x] : fabs(curr.x - HX * sx.x); // Distance to nearest Ox and Oy straight lines.
        double dy = sx.y >= 0 ? oy[sy.y] - curr.y : fabs(HY * sy.y - curr.y);
        if (dy / dx <= k) { //   Intersection with straight line parallel Ox axis.        
            next_i = 1;
            next.y = sy.y >= 0 ? oy[sy.y] : HY * sy.y;
            next.x = curr.x - (next.y - curr.y) / k;
        } else { //   Intersection with straight line parallel Oy axis.            
            next_i = 2;
            next.x = sx.x >= 0 ? ox[sx.x] : HX * sx.x;
            next.y = curr.y - k * (next.x - curr.x);
        }
        if (next.x < (uv.x + _MIN_VALUE)) {
            next_i = 0;
            next = uv;
            result += integrate_chanel_slant_left(tl, curr, curr_i, next, next_i,
                    sx, sy, bv.x, ib, ox, oy, density);
            break;
        }

        result += integrate_chanel_slant_left(tl, curr, curr_i, next, next_i,
                sx, sy, bv.x, ib, ox, oy, density);

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

double integrate_right_triangle_bottom_right(const dp_t& bv, const dp_t& uv, int tl,
        const double* ox, const double* oy, double* density) {
    double k = 0.;
    if (!try_get_slope_ratio(bv, uv, k)) return k;
    ip_t sx, sy;

    sx.x = static_cast<short> ((bv.x + _MIN_VALUE) / HX);
    if (bv.x + _MIN_VALUE <= 0) sx.x -= 1;
    sx.y = sx.x + 1;
    sy.x = static_cast<short> ((bv.y + _MIN_VALUE) / HY);
    if (bv.y + _MIN_VALUE <= 0) sy.x -= 1;
    sy.y = sy.x + 1;

    ip_t ib(sx.x, ib.x + 1);

    double result = 0.;
    short curr_i = 0, next_i = 0;
    dp_t curr = bv, next;
    while (true) {
        double dox = sx.y >= 0 ? fabs(ox[sx.y] - curr.x) : fabs(HX * sx.y - curr.x);
        double doy = sy.y >= 0 ? fabs(oy[sy.y] - curr.y) : fabs(HY * sy.y - curr.y);
        if (doy / dox <= k) {//   Across with straight line parallel Ox axis.            
            next_i = 1;
            next.y = sy.y >= 0 ? oy[sy.y] : HY * sy.y;
            next.x = bv.x + (next.y - bv.y) / k;
        } else {//   Across with straight line parallel Oy axis.
            next_i = 2;
            next.x = sx.y >= 0 ? ox[sx.y] : HX * sx.y;
            next.y = bv.y + k * (next.x - bv.x);
        }
        if (next.x > (uv.x - _MIN_VALUE)) {
            next = uv;
            next_i = 0;
            result += integrate_chanel_slant_right(tl, curr, curr_i,
                    next, next_i,
                    sx, bv.x, ib,
                    sy, ox, oy, density);
            break;
        }
        result += integrate_chanel_slant_right(tl, curr, curr_i,
                next, next_i,
                sx, bv.x, ib,
                sy, ox, oy, density);

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

double integrate_right_triangle_upper_left(const dp_t& bv, const dp_t& uv, int tl,
        const double* ox, const double* oy, double* density) {
    double k = 0.;
    if (!try_get_slope_ratio(bv, uv, k)) return k;

    ip_t sx, sy, ib;
    sx.x = static_cast<short> ((bv.x + _MIN_VALUE) / HX); //   -  If bv.x is in grid edge I want it will be in the right side.
    if (bv.x + _MIN_VALUE <= 0) sx.x -= 1;
    sx.y = sx.x + 1;
    sy.x = static_cast<short> ((bv.y + _MIN_VALUE) / HY); //   -  If bv.y is in grid edge I want it will be in the upper square.
    if (bv.y + _MIN_VALUE <= 0) sy.x -= 1;
    sy.y = sy.x + 1;
    ib.x = static_cast<short> ((uv.x - _MIN_VALUE) / HY); //   -  If uv.x is in grid edge I want it will be in the left side.
    if (uv.x - _MIN_VALUE <= 0) ib.x -= 1;
    ib.y = ib.x + 1;

    short curr_i = 0, next_i = 0;
    dp_t curr = bv, next;
    double result = 0.;
    while (true) {
        double dox = sx.y >= 0 ? ox[sx.y] - curr.x : fabs(HX * sx.y - curr.x);
        double doy = sy.y >= 0 ? oy[sy.y] - curr.y : fabs(HY * sy.y - curr.y);
        if (doy / dox <= k) { //   intersection with straight line parallel Ox axis.
            next_i = 1;
            next.y = sy.y >= 0 ? oy[sy.y] : HY * sy.y;
            next.x = bv.x + (next.y - bv.y) / k;
        } else {//   intersection with straight line parallel Oy axis.            
            next_i = 2;
            next.x = sx.y >= 0 ? ox[sx.y] : HX * sx.y;
            next.y = bv.y + k * (next.x - bv.x);
        }
        if (next.x > (uv.x - _MIN_VALUE)) {
            next_i = 0;
            next = uv;
            result += integrate_chanel_slant_left(tl, curr, curr_i, next, next_i,
                    sx, sy, uv.x, ib, ox, oy, density);
            break;
        }
        result += integrate_chanel_slant_left(tl, curr, curr_i, next, next_i,
                sx, sy, uv.x, ib, ox, oy, density);

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

double integrate_right_triangle_upper_right(const dp_t& bv, const dp_t& uv, int tl,
        const double* ox, const double* oy, double* density) {
    double k = 0.;
    if (!try_get_slope_ratio(bv, uv, k)) return k;

    ip_t sx, sy, ib;

    sx.x = static_cast<short> ((bv.x - _MIN_VALUE) / HX); //   -  If bv.x is in grid edge I want it will be between in the left side.
    if (bv.x - _MIN_VALUE <= 0) sx.x -= 1;
    sx.y = sx.x + 1;

    ib.x = static_cast<short> ((uv.x + _MIN_VALUE) / HX);
    if (uv.x + _MIN_VALUE <= 0) ib.x -= 1;
    ib.y = ib.x + 1;
    sy.x = static_cast<short> ((bv.y + _MIN_VALUE) / HY); //   -  If bv.y is in grid edge I want it will be in the upper side.
    if (bv.y + _MIN_VALUE <= 0) sy.x -= 1;
    sy.y = sy.x + 1;

    double result = 0.;
    short curr_i = 0, next_i = 0;
    dp_t curr = bv, next;
    while (true) {
        double distOx = sx.x >= 0 ? fabs(curr.x - ox[sx.x]) : fabs(curr.x - HX * sx.x);
        double distOy = sy.y >= 0 ? fabs(oy[sy.y] - curr.y) : fabs(HY * sy.y - curr.y);
        if (distOy / distOx <= k) { //   Intersection with straight line parallel Ox axis.
            next_i = 1;
            next.y = sy.y >= 0 ? oy[sy.y] : HY * sy.y;
            next.x = bv.x - (next.y - bv.y) / k;
        } else { //   Intersection with straight line parallel Oy axis.
            next_i = 2;
            next.x = sx.x >= 0 ? ox[sx.x] : HX * sx.x;
            next.y = bv.y - k * (next.x - bv.x);
        }
        if (next.x < uv.x + _MIN_VALUE) {
            next_i = 0;
            next = uv;
            result += integrate_chanel_slant_right(tl, curr, curr_i, next, next_i,
                    sx, uv.x, ib, sy, ox, oy, density);
            break;
        }
        result += integrate_chanel_slant_right(tl, curr, curr_i, next, next_i,
                sx, uv.x, ib, sy, ox, oy, density);

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

double integrate_bottom_triangle(int tl, const dp_t& l, const dp_t& r, const dp_t& m,
        const double* ox, const double* oy, double* density) {
    double result = 0.;
    if (m.x == l.x) {
        result = integrate_right_triangle_bottom_right(m, r, tl, ox, oy, density);
    } else if (m.x == r.x) {
        result = integrate_right_triangle_bottom_left(m, l, tl, ox, oy, density);
    } else if (m.x < l.x) {
        result = integrate_right_triangle_bottom_right(m, r, tl, ox, oy, density);
        result -= integrate_right_triangle_bottom_right(m, l, tl, ox, oy, density);
    } else if (m.x > l.x && m.x < r.x) {
        result = integrate_right_triangle_bottom_left(m, l, tl, ox, oy, density);
        result += integrate_right_triangle_bottom_right(m, r, tl, ox, oy, density);
    } else if (m.x > r.x) {
        result = integrate_right_triangle_bottom_left(m, l, tl, ox, oy, density);
        result -= integrate_right_triangle_bottom_left(m, r, tl, ox, oy, density);
    }
    return result;
}

double integrate_upper_triangle(int tl, const dp_t& l, const dp_t& r, const dp_t& m,
        const double* ox, const double* oy, double* density) {
    double result = 0.;
    if (m.x == l.x) {
        result = integrate_right_triangle_upper_right(r, m, tl, ox, oy, density);
    } else if (m.x == r.x) {
        result = integrate_right_triangle_upper_left(l, m, tl, ox, oy, density);
    } else if (m.x < l.x) {
        result = integrate_right_triangle_upper_right(r, m, tl, ox, oy, density);
        result -= integrate_right_triangle_upper_right(l, m, tl, ox, oy, density);
    } else if (m.x > l.x && m.x < r.x) {
        result = integrate_right_triangle_upper_left(l, m, tl, ox, oy, density);
        result += integrate_right_triangle_upper_right(r, m, tl, ox, oy, density);
    } else if (m.x > r.x) {
        result = integrate_right_triangle_upper_left(l, m, tl, ox, oy, density);
        result -= integrate_right_triangle_upper_left(r, m, tl, ox, oy, density);
    }
    return result;
}

double integrate_uniform_triangle_wall(int tl, const dp_t& a, const dp_t& b, const dp_t& c,
        const double* ox, const double* oy, double* density) {
    return 0;
}

double integrate_uniform_triangle(int tl,
        dp_t& x,
        dp_t& y,
        dp_t& z,
        const double* ox,
        const double* oy,
        double* density) {
    //   a * x  +  b * y  = c.
    double a = z.y - x.y;
    if (fabs(a) < _MIN_VALUE) return _MIN_VALUE;
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

    return integrate_upper_triangle(tl, y, ip, z, ox, oy, density)
            + integrate_bottom_triangle(tl, y, ip, x, ox, oy, density);
}

inline double func_u(double x, double y) {
    return B * y * (1. - y) * (_PI / 2. + atan(-x));
}

inline double func_v(double t, double x, double y) {
    return atan((x - LB) * (x - RB) * (1. + t) / 10. * (y - UB) * (y - BB));
}

double func_f(
        double tl_on_tau,
        double x,
        double y) {
    double arg_v = (x - LB) * (x - RB) * (1. + tl_on_tau) / 10. * (y - UB) * (y - BB);
    double rho = analytical_solution(tl_on_tau, x, y);
    double dRhoDT = x * y * cos(tl_on_tau * x * y);
    double dRhoDX = tl_on_tau * y * cos(tl_on_tau * x * y);
    double dRhoDY = tl_on_tau * x * cos(tl_on_tau * x * y);
    double u = func_u(x, y);
    double v = func_v(tl_on_tau, x, y);
    double duDX = -B * y * (1. - y) / (1. + x * x);
    double dvDY = (x - LB) * (x - RB) * (1. + tl_on_tau) / 10. * (y - BB + y - UB);
    dvDY /= (1. + arg_v * arg_v);
    double res = dRhoDT + rho * duDX + u * dRhoDX + rho * dvDY + v * dRhoDY;

    //  printf("x = %f\n", x);
    //  printf("y = %f\n", y);
    //  printf("arg_v = %f\n", arg_v);
    //  printf("rho = %f\n", rho);
    //  printf("dRhoDT = %f\n", dRhoDT);
    //  printf("dRhoDX = %f\n", dRhoDX);
    //  printf("dRhoDY = %f\n", dRhoDY);
    //  printf("u = %f\n", u);
    //  printf("duDX = %f\n", duDX);
    //  printf("v = %f\n", v);
    //  printf("dvDY = %f\n", dvDY);
    //  printf("res = %f\n", res);

    return res;
}

quad_type get_coordinates_on_prev_layer(int cur_tl, int ix, int iy,
        const double* ox,
        const double* oy,
        dp_t& alpha, dp_t& beta, dp_t& gamma, dp_t& theta) {
    //   1. First of all let's compute coordinates of square vertexes.
    //  OX:
    if (ix == 0) {
        alpha.x = ox[ix];
        beta.x = (ox[ix] + ox[ix + 1]) / 2.;
        gamma.x = (ox[ix] + ox[ix + 1]) / 2.;
        theta.x = ox[ix];
    } else if (ix == OX_LEN) {
        alpha.x = (ox[ix - 1] + ox[ix]) / 2.;
        beta.x = ox[ix];
        gamma.x = ox[ix];
        theta.x = (ox[ix - 1] + ox[ix]) / 2.;
    } else {
        alpha.x = (ox[ix - 1] + ox[ix]) / 2.;
        beta.x = (ox[ix + 1] + ox[ix]) / 2.;
        gamma.x = (ox[ix + 1] + ox[ix]) / 2.;
        theta.x = (ox[ix - 1] + ox[ix]) / 2.;
    }

    //  OY:
    if (iy == 0) {
        alpha.y = oy[iy];
        beta.y = oy[iy];
        gamma.y = (oy[iy] + oy[iy + 1]) / 2.;
        theta.y = (oy[iy] + oy[iy + 1]) / 2.;
    } else if (iy == OY_LEN) {
        alpha.y = (oy[iy] + oy[iy - 1]) / 2.;
        beta.y = (oy[iy] + oy[iy - 1]) / 2.;
        gamma.y = oy[iy];
        theta.y = oy[iy];
    } else {
        alpha.y = (oy[iy] + oy[iy - 1]) / 2.;
        beta.y = (oy[iy] + oy[iy - 1]) / 2.;
        gamma.y = (oy[iy] + oy[iy + 1]) / 2.;
        theta.y = (oy[iy] + oy[iy + 1]) / 2.;
    }

    double u, v;

    // Now let's compute new coordinates on the previous time level of alpha, beta, gamma, theta points.
    u = func_u(alpha.x, alpha.y);
    v = func_v(TAU * cur_tl, alpha.x, alpha.y);
    alpha.x -= TAU * u;
    alpha.y -= TAU * v;

    u = func_u(beta.x, beta.y);
    v = func_v(TAU * cur_tl, beta.x, beta.y);
    beta.x -= TAU * u;
    beta.y -= TAU * v;

    u = func_u(gamma.x, gamma.y);
    v = func_v(TAU * cur_tl, gamma.x, gamma.y);
    gamma.x -= TAU * u;
    gamma.y -= TAU * v;

    u = func_u(theta.x, theta.y);
    v = func_v(TAU * cur_tl, theta.x, theta.y);
    theta.x -= TAU * u;
    theta.y -= TAU * v;

    dp_t intersection = get_intersection_point(alpha, beta, gamma, theta);
    if ((beta.y - intersection.y) * (theta.y - intersection.y) > 0.) return pseudo; // ??
    if ((alpha.x - intersection.x) * (gamma.x - intersection.x) > 0.) return pseudo; // ??
    double product = get_vector_product(alpha, beta, theta); // ?
    if (product < 0.) return pseudo;

    // значит что точка улетела за левую границу
    if (theta.x < 0 ||
            theta.y < 0 ||
            beta.x < 0 ||
            beta.y < 0 ||
            gamma.x < 0 ||
            gamma.y < 0 ||
            alpha.x < 0 ||
            alpha.y < 0) {
        return normal;
        //return wall;
    }
    return normal;
}


// Type of quadrangle: 0 - pseudo; 1 - convex; 2 - concave;

quad_type get_quadrangle_type(int tl, int ix, int iy,
        const double* ox,
        const double* oy,
        dp_t& a, //   -  First vertex of first triangle.
        dp_t& b, //   -  Second vertex of first triangle.
        dp_t& c, //   -  Third vertex of first triangle.
        dp_t& k, //   -  First vertex of second triangle.
        dp_t& m, //   -  Second vertex of second triangle.
        dp_t& n) //   -  Third vertex of second triangle.
{
    dp_t alpha, beta, gamma, theta; // coordinates on previous time layer
    quad_type type = get_coordinates_on_prev_layer(tl, ix, iy, ox, oy, alpha, beta, gamma, theta);
    a = alpha;
    b = beta;
    c = gamma;
    k = alpha;
    m = theta;
    n = gamma;
    return type;
}

double integrate(double tl, int ix, int iy,
        const double* ox,
        const double* oy,
        double* density) {
    dp_t a1, b1, c1, a2, b2, c2;
    quad_type type = get_quadrangle_type(tl, ix, iy, ox, oy, a1, b1, c1, a2, b2, c2);
    if (type != normal && type != wall) {
        return -1.;
    }

    // чтобы правилно отработала процедура интегрирования
    // точки должны идти в порядке возрастания y координаты
    sort_by_y(a1, b1, c1);
    sort_by_y(a2, b2, c2);

    // check the type of triangle to select appropriate computation method
    double result = 0.;
    switch (type) {
        case wall:
            result += integrate_uniform_triangle_wall(tl, a1, b1, c1, ox, oy,
                    density);
            result += integrate_uniform_triangle_wall(tl, a2, b2, c2, ox, oy,
                    density);
            return result;
        case normal:
            result += integrate_uniform_triangle(tl, a1, b1, c1, ox, oy, density);
            result += integrate_uniform_triangle(tl, a2, b2, c2, ox, oy, density);
            return result;
        case concave:
        case convex:
        case pseudo:
            return 0.;
    }
    return 0;
}

double get_norm_of_error(double* density, int x_length, int y_length, double* ox,
        double* oy,
        double ts_count_mul_steps) {
    double result = 0.;
    for (int k = 1; k < y_length; ++k) {
        for (int j = 1; j < x_length; ++j) {
            result += fabs(analytical_solution(ts_count_mul_steps, ox[j], oy[k])
                    - density[(x_length + 1) * k + j]);
        }
    }
    double hx = ox[1] - ox[0];
    double hy = oy[1] - oy[0];
    return hx * hy * result;
}

double solve(const double* ox, const double* oy, double* density) {
    double* prev_density = new double [ XY_LEN ];
    for (int iy = 0; iy < OY_LEN + 1; iy++) {
        for (int ix = 0; ix < OX_LEN + 1; ix++) {
            prev_density[(OX_LEN + 1) * iy + ix] = analytical_solution(0., ox[ix], oy[iy]);
        }
    }

    for (int it = 1; it < TIME_STEP_CNT + 1; it++) {
        for (int i = 0; i <= OX_LEN; i++) {
            density[i] = init_side(ox[i], BB, TAU * it, bottom);
            density[(OX_LEN + 1) * OY_LEN + i] = init_side(ox[i], UB, TAU * it, up);
        }

        for (int i = 0; i <= OY_LEN; i++) {
            density[(OX_LEN + 1) * i] = init_side(LB, oy[i], TAU * it, left);
            density[(OX_LEN + 1) * i + OX_LEN] = init_side(RB, oy[i], TAU * it, right);
        }

        for (int iy = 1; iy < OY_LEN; iy++) {
            for (int ix = 1; ix < OX_LEN; ix++) {
                int index = (OX_LEN + 1) * iy + ix;

                double value = integrate(it, ix, iy, ox, oy, prev_density);

                double h = (ox[ix + 1] - ox[ix - 1]) / 2.;
                value /= h;
                h = (oy[iy + 1] - oy[iy - 1]) / 2.;
                value /= h;

                double rp = func_f(TAU * it, ox[ix], oy[iy]);
                density[index] = value + TAU * rp;
            }
        }
        memcpy(prev_density, density, (OX_LEN + 1) * (OY_LEN + 1) * sizeof (double));
    }

    delete[] prev_density;
    return 0;
}

double* compute_density(double b,
        double lb,
        double rb,
        double bb,
        double ub,
        double tau,
        int time_step_count,
        int ox_length,
        int oy_length,
        double* norm) {
    B = b;
    UB = ub;
    BB = bb;
    LB = lb;
    RB = rb;
    TAU = tau;
    OX_LEN = ox_length;
    OY_LEN = oy_length;
    TIME_STEP_CNT = time_step_count;
    XY_LEN = (ox_length + 1) * (oy_length + 1);

    double* density = new double [ XY_LEN ];
    double* ox = new double [ OX_LEN + 1 ];
    double* oy = new double [ OY_LEN + 1 ];

    for (int i = 0; i <= OX_LEN; i++) {
        ox[i] = lb + i * (rb - lb) / OX_LEN;
    }

    for (int i = 0; i <= OY_LEN; i++) {
        oy[i] = bb + i * (ub - bb) / OY_LEN;
    }

    HX = oy[1] - oy[0];
    HY = oy[1] - oy[0];

    print_params(B, LB, RB, BB, UB, TAU, time_step_count, OX_LEN, OY_LEN);

    solve(ox, oy, density);

    *norm = get_norm_of_error(density, OX_LEN, OY_LEN, ox, oy,
            time_step_count * TAU);
    //  printf("Norm L1 = %f\n", *norm);
    printf("%d x %d wall count = %d\n", ox_length + 1, oy_length + 1, TMP_WALL_CNT);

    delete[] ox;
    delete[] oy;
    return density;
}