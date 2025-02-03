#ifndef UTILITY_H
#define UTILITY_H
#include <pari/pari.h>

ulong sigbits(GEN x);
/* Significant bits.
 * Input:   Integer x (assumed to be normalized).
 * Output:  The number of significant bits in x, i.e floor(log_2(x)) + 1.
*/

GEN rqfumodm(GEN bnf, GEN m); // This function is very bad.
/* Real quadratic fundamental unit mod m.
 * Input:   bnf (as output by bnfinit0);
            integer m.
 * Output:  [x, y], the coefficients of the fundamental unit modulo some m' | m.
*/

GEN rqfumodm_fixed(GEN bnf, GEN m); // This function is very bad.
/* Real quadratic fundamental unit mod m fixed.
 * Input:   bnf (as output by bnfinit0);
            integer m.
 * Output:  [x, y], the coefficients of the fundamental unit modulo m.
*/

ulong rqfupval(GEN bnf, ulong p, ulong c); // This function is very bad.
/* Real quadratic fundamental unit p-adic valuation.
 * Input:  bnf (as output by bnfinit0);
           prime number p;
           c = {1, 2}.
 * Output: Integer n, the p-adic valuation of the c'th coordinate of the fundamental unit.
*/

GEN rqfusmoothpart(GEN bnf, ulong B, ulong c); // This function is very bad.
/* Real quadratic fundamental unit smooth part.
 * Input:   bnf (as output by bnfinit0);
            integer B;
            c = {1, 2}.
 * Output:  Integer n, the B-smooth part of the c'th coordinate of the fundamental unit.
*/


GEN mulqidivn (GEN O, GEN qi1, GEN qi2, GEN N); // Returns [x, y] = (x + y*sqrt(d))/(s*N) = ((a + b*sqrt(d))/s * (a_ + b_*sqrt(d))/s)/N, where qi1 = [a, b] and qi2 = [a_, b_].
GEN mulqig (GEN O, GEN qi1, GEN qi2); // Returns [x, y] = (x + y*sqrt(d))/s = (a + b*sqrt(d))/s * (a_ + b_*sqrt(d))/s, where qi1 = [a, b] and qi2 = [a_, b_].
GEN invqi(GEN qi); // Returns [a, -b], where qi = [a, b].

#endif