#ifndef QUADRATIC_ORDERS_H
#define QUADRATIC_ORDERS_H
#include <pari/pari.h>

GEN rqoinit(GEN d);
/* 
 * Real quadratic order initialization.
 * Input: Nonsquare t_INT d.
 * Output: [d, [d_0, f], s, D], where
 *          d = d_0*f^2 and 
 * D is the discriminant of the order with conductor f in Q(sqrt(d_0));
 * s = 2 if d_0 = 1 (mod 4) and s = 1 otherwise.
 */

GEN rqiinit(GEN S, GEN Q, GEN P);
/* 
 * Real quadratic ideal initalization.
 * Input: t_INT S, Q, P.
 * Output: [S, [Q, P]].
 * (Does not perform a check that this representation is valid.)
 */

void checkrqi(GEN O, GEN a, const char *f);
/*
 * Check real quadratic ideal.
 * Input: Real quadratic order O (as output by rqoinit), ideal (as output by rqiinit), name of function calling.
 * Throws a type error if a is not a valid O-ideal.
 */

GEN inucomp(GEN O, GEN a, GEN b, long flag);
/* 
 * Ideal NUCOMP.
 * Input: Order O (as output by rqoinit), reduced ideals a and b (as output by rqiinit).
 * Output: [[1, [Q, P]], [A, B, C]], where [1, [Q, P]] = a*b (computed using a version of Shanks' NUCOMP), coefficients A, B, C (as t_INT) of relative generator (A+B*sqrt(d))/C.
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
 */

#endif