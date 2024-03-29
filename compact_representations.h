#ifndef COMPACT_REPRESENTATIONS_H
#define COMPACT_REPRESENTATIONS_H
#include "fp_representations.h"

GEN expandcr(GEN O, GEN cr);
/* Expand compact representation.
 * Input:   Real quadratic order O (as output by rqoinit);
            list [[m_1, n_1], L_1], ..., [[m_l, n_l], L_l], [[m_{l+1}, n_{l+1}], L_{l+1}]] (as output for example by cr), such that
                theta = lambda_{l+1} * \prod_{i=0}^l (lambda_i/L_{i+1})^{2^(l-i)}, where lambda_i = (m_i + n_i*sqrtd(d))/s
            is a compact representation of some theta in O.
 * Output:  [a, b], such that theta = (a + b*sqrt(d))/s.
*/

GEN crmodm(GEN O, GEN cr, GEN m);
/* Compact representation p-adic valuation.
 * Input:   Real quadratic order O (as output by rqoinit);
            a compact representation [[m_1, n_1], L_1], ..., [[m_l, n_l], L_l], [[m_{l+1}, n_{l+1}], L_{l+1}]] (as output for example by cr),
            where m divides at most L_{l+1};
            integer m.
 * Output:  [x, y], with 0 <= x, y < m, where
                (x + y*sqrt(d))/s = theta 
                = lambda_{l+1} * \prod_{i=0}^l (lambda_i/L_{i+1})^{2^(l-i)}, where lambda_i = (m_i + n_i*sqrtd(d))/s.
*/

GEN crmodm_alt(GEN O, GEN cr, GEN m);
/* Compact representation p-adic valuation.
 * Input:   Real quadratic order O (as output by rqoinit);
            a compact representation [[m_1, n_1], L_1], ..., [[m_l, n_l], L_l], [[m_{l+1}, n_{l+1}], L_{l+1}]] (as output for example by cr),
            integer m.
 * Output:  [x, y], with 0 <= x, y < m, where
                (x + y*sqrt(d))/s = theta 
                = lambda_{l+1} * \prod_{i=0}^l (lambda_i/L_{i+1})^{2^(l-i)}, where lambda_i = (m_i + n_i*sqrtd(d))/s.
*/

ulong crpval(GEN O, GEN cr, ulong p, ulong c);
/* Compact representation p-adic valuation.
 * Input:   Real quadratic order O (as output by rqoinit);
            a compact representation [[m_1, n_1], L_1], ..., [[m_l, n_l], L_l], [[m_{l+1}, n_{l+1}], L_{l+1}]] (as output for example by cr),
            where p divides at most L_{l+1};
            integer p;
            c = {1, 2}.
 * Output:  Integer n, the p-adic valuation of the c'th coordinate of [x, y] = (x + y*sqrt(d))/s = theta
                = lambda_{l+1} * \prod_{i=0}^l (lambda_i/L_{i+1})^{2^(l-i)}, where lambda_i = (m_i + n_i*sqrtd(d))/s.
*/

ulong crpval_alt(GEN O, GEN cr, ulong p, ulong c);
/* Compact representation p-adic valuation alternative.
 * Input:   Real quadratic order O (as output by rqoinit);
            a compact representation [[m_1, n_1], L_1], ..., [[m_l, n_l], L_l], [[m_{l+1}, n_{l+1}], L_{l+1}]] (as output for example by cr),
            integer p;
            c = {1, 2}.
 * Output:  Integer n, the p-adic valuation of the c'th coordinate of [x, y] = (x + y*sqrt(d))/s = theta
                = lambda_{l+1} * \prod_{i=0}^l (lambda_i/L_{i+1})^{2^(l-i)}, where lambda_i = (m_i + n_i*sqrtd(d))/s.
*/

GEN crsmoothpart(GEN O, GEN cr, GEN S, ulong c); // O(logR*B*max{B, logD})
/* Compact representation smooth part.
 * Input:   Real quadratic order O (as output by rqoinit);
            a compact representation [[m_1, n_1], L_1], ..., [[m_l, n_l], L_l], [[m_{l+1}, n_{l+1}], L_{l+1}]] (as output for example by cr),
            where for every prime p <= B, p divides at most L_{i+1};
            vector S of primes below a bound B;
            c = {1, 2}.
 * Output:  Integer n, the B-smooth part of the c'th coordinate of [x, y] = (x + y*sqrt(d))/s = theta
                = lambda_{l+1} * \prod_{i=0}^l (lambda_i/L_{i+1})^{2^(l-i)}, where lambda_i = (m_i + n_i*sqrtd(d))/s;
            that is the c'th coordinate can be written as n*r, where r is only divisible by primes p > B.
*/

GEN crsmoothpart2(GEN O, GEN b, GEN y, GEN q, GEN S, ulong c); // O(logR*B*max{logB, logD}) (crude bound; probably lower in practice)
/* Compact representation smooth part, alternative.
 * Input:   Real quadratic order O (as output by rqoinit);
            reduced principal ideal b = [1, [Q, P]] with (P + sqrt(d))/Q > 1 and -1 < (P - sqrt(d))/Q < 0;
            y, q positive rational number, such that
                |log_2(theta) - y| < q.
            where (theta) = b;
            vector S of primes below a bound B;
            c = {1, 2}.
 * Output:  Integer n, the B-smooth part of the c'th coordinate of [x, y] = (x + y*sqrt(d))/s = theta,
                (theta) = b.
*/

GEN crsmoothpart_alt(GEN O, GEN cr, GEN S, ulong c);
/* Compact rerpresentation smooth part.
 * Input:   Real quadratic order O (as output by rqoinit);
            a compact representation [[m_1, n_1], L_1], ..., [[m_l, n_l], L_l], [[m_{l+1}, n_{l+1}], L_{l+1}]] (as output for example by cr);
            vector S of primes below a bound B;
            c = {1, 2}.
 * Output:  Integer n, the B-smooth part of the c'th coordinate of [x, y] = (x + y*sqrt(d))/s = theta
                = lambda_{l+1} * \prod_{i=0}^l (lambda_i/L_{i+1})^{2^(l-i)}, where lambda_i = (m_i + n_i*sqrtd(d))/s;
            that is the c'th coordinate can be written as n*r, where r is only divisible by primes p > B.
*/

long crnorm_sign(GEN O, GEN cr);
/* Compact representation norm sign.
 * Input:   Real quadratic order O (as output by rqoinit);
            compact representation cr;
* Output:   Sign of the norm of cr: either -1 or +1.
*/

#endif