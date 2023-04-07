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
 * Output:  [[[f, p], [b_, d, k]], [a, b]], where
 *          (b_, d, k) is a reduced (f, p) representation of the product a = a'*a'',
 *          where b_ = [1, [Q, P]] with (P + sqrt(d))/Q > 1, -1 < (P - sqrt(d))/Q < 0,
 *          k <= k' + k'' + 1,
 *          f = f* + 17/8, where f* = f' + f'' + 2^(-p)*f'*f'' and
 *          a, b are integers such that 
 *          b_ = ((a + b*sqrt(d))/s)/(N(b')*N(b''))*b'*b''. 
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
 */

GEN wnear(GEN O, GEN fprep, GEN w);
/*
 * Input:   Real quadratic order O (as output by rqoinit);
 *          A reduced (f, p) representation fprep = (b, d, k) of ideal a,
 *          where b = (Q, P) with P + floor(sqrt(d)) >= Q and 0 <= floor(sqrt(d)) - P <= Q,
 *          positive integer w.
 * Output:  A w-near (f + 9/8, p) representation of the ideal a.
 */

GEN ewnear(GEN O, GEN fprep, GEN w);
/*
 * Extended wnear.
 * Input:   Real quadratic order O (as output by rqoinit);
 *          A reduced (f, p) representation fprep = (b, d, k) of ideal a, 
 *          where b = (Q, P) with  with P + floor(sqrt(d)) >= Q and 0 <= floor(sqrt(d)) - P <= Q,
 *          positive integer w with k < w.
 * Output:  [[[f, p], [c, g, h]], [a, b]] with integers a, b, where
 *          (c, g, h) is a w-near (f + 9/8, p) representation of the ideal a and
 *          c = ((a + b*sqrt(d))/Q)b.
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

GEN iexp(GEN O, GEN fprep, GEN w, GEN n, long flag);
/*
 * Ideal exponentiation.
 * Input:   Real quadratic order O (as output by rqoinit),
 *          w-near (f', p representation) (b', d', k') of invertible real quadratic ideal a (this isn't checked),
 *          positive integer n.
 * Output:  A w-near (f, p)  representation (b, d, k) of a^n for some suitable f.
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
*/

GEN addxy(GEN O, GEN fprep1, GEN fprep2, GEN x, GEN y, long flag);
/*
 * Input:   Real quadratic order O (as output by rqoinit),
 *          x-near (f', p) representation (a[x], d', k') of the ideal (1),
 *          y-near (f'', p) rerepsentation (a[y], d'', k'') of the ideal (1).
 * Output:  An (x+y)-near (f, p) representation (a[x+y], d, k) of (1), where
 *          f = f' + f'' + (f'*f'')/2^p + 13/4.
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
*/

GEN ax(GEN O, GEN x, GEN p);
/*
 * Input:   Real quadratic order O (as output by rqoinit), positive integers x and p.
 * Output:  An x-near (f, p) representation (a[x], d, k) of the ideal (1) in O for some f \in [1, 2^p).
*/

GEN eaddxy(GEN O, GEN fprep1, GEN fprep2, GEN x, GEN y, long flag);
/*
 * Extended addxy.
 * Input:   Real quadratic order O (as output by rqoinit),
 *          x-near (f', p) representation (a[x], d', k') of the ideal (1),
 *          y-near (f'', p) rerepsentation (a[y], d'', k'') of the ideal (1).
 * Output:  [[[f, p], [a[x+y], d, k]], [a,b]], where
 *          (a[x+y], d, k) is an (x+y)-near (f, p) representation of (1) with f = f' + f'' + (f'*f'')/2^p + 13/4 and
 *          (a, b) are integers such that 
 *          a[x+y] = ((lambda*theta'*theta'')/(N(a[x])*N(a[y])))(1), where
 *          lambda = (a + b*sqrt(d))/s and a[x] = theta'*a, a[y] = theta''*a.
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
*/


#endif