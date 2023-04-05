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

GEN imultiply(GEN O, GEN a, GEN b);
/*
 * Ideal multiply.
 * Input: Real quadratic order (as output by rqoinit), reduced real quadratic ideals (as output by rqiinit).
 * Output: Real quadratic ideal a*b (not necessarily primitive or reduced).
*/

GEN qiimultiply(GEN O, GEN qi, GEN a);
/*
 * Quadratic integer ideal multiply.
 * Input: Quadratic integer qi in the form [A, B, C] = (A + B*sqrt(d))/C (as output for example by inucomp), primitive ideal a = [1, [Q, P]] (as for example output by imultiply or rqiinit).
 * Output: The primitive ideal qi*a.
*/

GEN inucomp(GEN O, GEN a, GEN b, long flag);
/* 
 * Ideal NUCOMP.
 * Input: Order O (as output by rqoinit), reduced ideals a and b (as output by rqiinit).
 * Output: [[1, [Q, P]], [A, B, C]], where [1, [Q, P]] = a*b (computed using a version of Shanks' NUCOMP), coefficients A, B, C (as t_INT) of relative generator (A+B*sqrt(d))/C.
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
 */

GEN regulatorcf(GEN O, long prec, long flag);
/*
 * Regulator continued fraction
 * Input: Order O (as output by rqoinit), precision prec, flag = 0, 1, 2.
 * Output: Regulator of O as t_REAL to precision prec;
 *         flag controls whether the infrastructure is printed:
 *         0 to print nothing,
 *         1 to print only the ideals,
 *         2 to print ideals and distances.
*/

GEN regulatorshanks(GEN O, long prec, long flag);
/*
 * Input: Order O (as output by rqoinit), precision prec, flag = 0, 1.
 * Output: Regulator of O as t_REAL to precision prec;
 *         flag controls whether the baby/giant steps are printed:
 *         0 to print nothing,
 *         1 to print the steps.
 * Known Bugs:
 * - In some cases the giant step is too small and computation might return ~0 erroneously, because it hits a baby step immediately.
*/

#endif