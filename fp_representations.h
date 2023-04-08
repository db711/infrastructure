#ifndef FP_REPRESENTATIONS_H
#define FP_REPRESENTATIONS_H
#include "real_quadratic_orders.h"

GEN fprepinit (GEN f, GEN p, GEN b, GEN d, GEN k);
/* (f, p) representation initialization.
 * Input:   Real number f >= 1 (or NULL for unknown/unimportnat); 
            integers p, d and k;
            real quadratic ideal b = [S, [Q, P]] (as output by rqiinit).
 * Output:  [[f, p], [b, d, k]] = [[f, p], [[S, [Q, P]], d, k]].
 * (Representation is checked for validity, ideal b is not.)
 */

GEN fpremove(GEN fprep, GEN T, GEN C, GEN s);
/* Remove factor from (f, p) representations.
 * Input:   (f, p) representation fprep = [[f, p], [mb, d, k]] of some real quadratic ideal a;
            integers T, C and s, such that 
            if m = |(A + B*sqrt(d))/C|, then T = 2^s*A + B*floor(2^s*sqrt(d)) and s such that 2^s*|C| > 2^(p+4)*|B|.
 * Output:  [[f + 9/8, p], [b, d, k]], an (f + 9/8, p) representation (b, d', k') of a.
 */

GEN numult(GEN O, GEN fprep1, GEN fprep2, long flag); 
/* NUMULT: NUCOMP for (f, p) representations.
 * Input:   Real quadratic order O (as output by rqoinit);
            Reduced (f', p) representation fprep1 = [[f', p], [b', d', k']] of invertible ideal a';
            Reduced (f'', p) representation fprep2 = [[f'', p], [b'', d'', k'']] of invertible ideal a''.
 * Output:  [[[f, p], [b_, d, k]], [a, b]], where
            (b_, d, k) is a reduced (f, p) representation of the product a = a'*a'',
            where b_ = [1, [Q, P]] with (P + sqrt(d))/Q > 1, -1 < (P - sqrt(d))/Q < 0,
            k <= k' + k'' + 1,
            f = f* + 17/8, where f* = f' + f'' + 2^(-p)*f'*f'' and
            a, b are integers such that 
            b_ = ((a + b*sqrt(d))/s)/(N(b')*N(b''))*b'*b''. 
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
 */

GEN wnear(GEN O, GEN fprep, GEN w);
/* (Find) w-near (f, p) representation.
 * Input:   Real quadratic order O (as output by rqoinit);
            reduced (f, p) representation fprep = [[f, p], [b', d', k']] of ideal a,
            where b' = (Q, P) with P + floor(sqrt(d)) >= Q and 0 <= floor(sqrt(d)) - P <= Q;
            positive integer w.
 * Output:  [[f + 9/8, p], [b, d, k]], a w-near (f + 9/8, p) representation of the ideal a.
 */

GEN ewnear(GEN O, GEN fprep, GEN w);
/* Extended wnear.
 * Input:   Real quadratic order O (as output by rqoinit);
            reduced (f, p) representation fprep = [[f, p], [b', d', k']] of ideal a, 
            where b = (Q, P) with  with P + floor(sqrt(d)) >= Q and 0 <= floor(sqrt(d)) - P <= Q;
            positive integer w with k < w.
 * Output:  [[[f + 9/8, p], [c, g, h]], [a, b]] with integers a, b, where
            (c, g, h) is a w-near (f + 9/8, p) representation of the ideal a and
            c = ((a + b*sqrt(d))/Q)b.
 */

GEN wmult(GEN O, GEN fprep1, GEN fprep2, GEN w, long flag);
/* w-near multiplication.
 * Input:   Real quadratic order O (as output by rqoinit);
            Reduced w-near (f', p) representation fprep1 = [[f', p], [b', d', k']] of invertible ideal a';
            Reduced w-near (f'', p) representation fprep2 = [[f'', p], [b'', d'', k'']] of invertible ideal a''.
            (w-nearness is not, and in fact cannot, be checked.)
 * Output:  Reduced w-near (f* + 13/4, p) representation (b, d, k) of the product a = a'*a'',
            where f* = f' + f'' + 2^(-p)*f'*f''.
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
*/

GEN iexp(GEN O, GEN fprep, GEN w, GEN n, long flag);
/* Ideal exponentiation.
 * Input:   Real quadratic order O (as output by rqoinit);
            w-near (f', p representation) [[f', p], [b', d', k']] = fprep of invertible real quadratic ideal a;
            positive integer n.
 * Output:  [[f, p], [b, d, k]], a w-near (f, p)  representation (b, d, k) of a^n for some suitable f.
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
*/

GEN addxy(GEN O, GEN fprep1, GEN fprep2, GEN x, GEN y, long flag);
/* Add x(-near) y(-near).
 * Input:   Real quadratic order O (as output by rqoinit);
            [[f', p], [a[x], d', k'], an x-near (f', p) representation (a[x], d', k') of the ideal (1),
            [[f'', p], [a[y], d'', k''], an y-near (f'', p) representation (a[y], d'', k'') of the ideal (1).
 * Output:  [[f, p], [a[x+y], d, k], an (x+y)-near (f, p) representation (a[x+y], d, k) of (1), where
            f = f' + f'' + (f'*f'')/2^p + 13/4.
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
*/

GEN ax(GEN O, GEN x, GEN p);
/* (Find) a[x].
 * Input:   Real quadratic order O (as output by rqoinit);
            positive integers x and p.
 * Output:  [[f, p], [a[x], d, k]], where (a[x], d, k) is an (f, p) representation of the ideal (1) in O for some f \in [1, 2^p).
*/

GEN eaddxy(GEN O, GEN fprep1, GEN fprep2, GEN x, GEN y, long flag);
/* Extended addxy.
 * Input:   Real quadratic order O (as output by rqoinit);
            [[f', p], [a[x], d', k'], an x-near (f', p) representation (a[x], d', k') of the ideal (1),
            [[f'', p], [a[y], d'', k''], an y-near (f'', p) rerepsentation (a[y], d'', k'') of the ideal (1).
 * Output:  [[[f, p], [a[x+y], d, k]], [a,b]], where
            (a[x+y], d, k) is an (x+y)-near (f, p) representation of (1) with f = f' + f'' + (f'*f'')/2^p + 13/4 and
            (a, b) are integers such that 
            a[x+y] = ((lambda*theta'*theta'')/(N(a[x])*N(a[y])))(1), where
            lambda = (a + b*sqrt(d))/s and a[x] = theta'*a, a[y] = theta''*a.
 * Set flag = 0 to skip some tests if you are sure that your input is correct.
*/

GEN crax(GEN O, GEN x, GEN p);
/* Compact representation ax.
 * Input:   Real quadratic order O (as output by rqoinit);
            positive integers x and p,
            that satisfy 2^p > 11.2x*max(16,log_2(x)).      
 * Output:  [[f, p], [a[x], d, k]], [[[m_0, n_0], L_0], [[m_1, n_1], L_1], ..., [[m_l, n_l], L_l]].
            an x-near (f, p) representation (a[x], d, k) of the ideal (1) in O, where f < 2^(p-4)
            as well as integer pairs (m_i, n_i) and positive integers L_i for i = 0, 1, ..., l = floor(log_2(x)).
 * When this algorithm terminates we have a[x] = (theta) for
        theta = lambda_l * \prod_{i=0}^{l-1}(lambda_i/L_{i+1})^{2^(l-i)}, where lambda_i = (m_i + n_i*sqrtd(d))/s.
*/

#endif