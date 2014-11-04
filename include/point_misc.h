/* 
 * File:   misc.h
 * Author: alexanderefremov
 *
 * Created on November 4, 2014, 12:06 AM
 */

#ifndef MISC_H
#define	MISC_H

#include "common.h"


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

inline void ptcpy(point_t *dst, const point_t *source)
{
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

