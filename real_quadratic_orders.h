#ifndef QUADRATIC_ORDERS_H
#define QUADRATIC_ORDERS_H
#include <pari/pari.h>

GEN rqoinit(GEN d);
/* Real quadratic order initialization.
 * Input:   Nonsquare t_INT d.
 * Output:  [d, [d_0, f], s, D], where
 *          d = d_0*f^2 and 
 *          D is the discriminant of the order with conductor f in Q(sqrt(d_0));
 *          s = 2 if d_0 = 1 (mod 4) and s = 1 otherwise.
 */

GEN rqiinit(GEN S, GEN Q, GEN P);
/* Real quadratic ideal initialization.
 * Input:   Integers S, Q, P.
 * Output:  [S, [Q, P]].
 */

void checkrqi(GEN O, GEN a, const char *f);
/* Check real quadratic ideal.
 * Input:   Real quadratic order O (as output by rqoinit);
            real quadratic ideal (as output by rqiinit); 
            name of function calling.
 * Throws an appropriate error if a is not a valid O-ideal, given in the form (S)(Q, P).
 */

GEN imultiply(GEN O, GEN a, GEN b); // DEPRECATED
/* Ideal multiply.
 * Input:   Real quadratic order O (as output by rqoinit);
            reduced real quadratic ideals a and b (as output by rqiinit).
 * Output:  [S, [Q, P]], where (S)(Q, P) = a*b is not necessarily reduced and not necessarily primitive.
*/

GEN qiimultiply(GEN O, GEN qi, GEN a); // DEPRECATED
/* Quadratic integer ideal multiply.
 * Input:   Quadratic integer qi in the form [A, B, C] = (A + B*sqrt(d))/C (as output for example by inucomp);
            primitive ideal a (as for example output by imultiply or rqiinit).
 * Output:  [1, [Q, P]] = qi*a.
*/

GEN inucomp(GEN O, GEN a, GEN b, long flag);
/* Ideal NUCOMP.
 * Input:   Real quadratic order O (as output by rqoinit);
            reduced ideals a and b (as output by rqiinit);
            flag = {0, 1}.
 * Output:  [[1, [Q, P]], [A, B, C]], 
            where [1, [Q, P]] = a*b (computed using a version of Shanks' NUCOMP), 
            A, B, C (are integer coefficients of a relative generator (A+B*sqrt(d))/C.
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
 */

GEN regulatorcf(GEN O, long prec, long flag);
/* Regulator (via) continued fraction.
 * Input:   Order O (as output by rqoinit); 
            precision prec;
            flag = {0, 1, 2}.
 * Output:  [R, [x, y, o]], where
            R is the regulator of O to precision prec and
            [x, y] are the coefficients of the fundamental unit, so (x + y*sqrt(d))/2 is of norm o = {+-1} and generates O;
            this means x^2 - d*y^2 = 4o is the least solution to X^2 - d*Y^2 = +-4.
 * The flag controls whether the infrastructure is printed:
 *      0 to print nothing,
 *      1 to print only the ideals (do marginally less computations compared to 2),
 *      2 to print ideals and distances.
*/

GEN regulatorshanks(GEN O, long prec, long flag);
/* Regulator (via) Shanks' baby-step giant-step method.
 * Input:   Real quadratic order O (as output by rqoinit);
            precision prec;
            flag = {0, 1}.
 * Output:  R, the regulator of O to precision prec.
 * The flag controls whether the baby/giant steps are printed:
 *      0 to print nothing,
 *      1 to print the steps.
 * Known Bugs:
 *      - In some cases the giant step is too small and computation might return ~0 erroneously, because it hits a baby step immediately.
*/

#endif