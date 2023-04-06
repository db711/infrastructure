#ifndef FP_REPRESENTATIONS_H
#define FP_REPRESENTATIONS_H
#include "real_quadratic_orders.h"

GEN fprepinit (GEN f, GEN p, GEN b, GEN d, GEN k);
/*
 * (f, p) representation initialization.
 * Input:   t_REAL f >= 1 (or NULL for unknown), t_INT p, d and k, real quadratic ideal b (as output by rqiinit).
 * Output:  [[f, p], [b, d, k]] = [[f, p], [[S, [Q, P]], d, k]].
 * (Representation is checked for validity, ideal b is not.)
 */

GEN tofprep(GEN O, GEN b, GEN theta, GEN prec);
/*
 * To (f, p) representation.
 * Input:   Real quadratic order O (as output by rqoinit), ideal b (as outpute by rqiinit), quadratic integer [A, B, C] = (A + B*sqrt(d))/C.
 * Output:  (1, p) representation (b, d, k) of (1/theta)*b.
*/

GEN fpremove(GEN fprep, GEN T, GEN C, GEN s);
/*
 * Input:   An (f, p) representation fprep (as output by fprepinit) = (mb, d, k) of a; t_INT T, C and s,
 *          such that m = |(A + B*sqrt(d))/C| and T = 2^s*A + B*floor(2^s*sqrt(d)) and s such that 2^s*|C| > 2^(p+4)*|B|
 * Output:  [[f + 9/8, p], [b, d, k]], an (f + 9/8, p) representation (b, d', k') of a.
 */

GEN numult(GEN O, GEN fprep1, GEN fprep2, long flag); 
/*
 * Input:   Real quadratic order O (as output by rqoinit);
 *          Reduced (f', p) representation fprep1 = (b', d', k') of invertible ideal a';
 *          Reduced (f'', p) representation fprep2 = (b'', d'', k'') of invertible ideal a''.
 * Output:  Reduced (f, p) representation (b, d, k) of the product a = a'*a'',
 *          where b = [1, [Q, P]] with (P + sqrt(d))/Q > 1, -1 < (P - sqrt(d))/Q < 0,
 *          k <= k' + k'' + 1,
 *          f = f* + 17/8, where f* = f' + f'' + 2^(-p)*f'*f''.
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
 */

GEN wnear(GEN O, GEN fprep, GEN w);
/*
 * Input:   Real quadratic order O (as output by rqoinit);
 *          A reduced (f, p) representation fprep of ideal a, t_INT w.
 * Output:  A w-near (f + 9/8, p) representation of the ideal a.
 */

GEN wmult(GEN O, GEN fprep1, GEN fprep2, GEN w, long flag);
/*
 * Input:   Real quadratic order O (as output by rqoinit);
 *          Reduced w-near (f', p) representation fprep1 = (b', d', k') of invertible ideal a';
 *          Reduced w-near (f'', p) representation fprep2 = (b'', d'', k'') of invertible ideal a''.
 *          (w-nearness is not, and in fact cannot, be checked.)
 * Output:  Reduced w-near (f* + 13/4, p) representation (b, d, k) of the product a = a'*a'',
 *          where f* = f' + f'' + 2^(-p)*f'*f''.
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
*/

#endif