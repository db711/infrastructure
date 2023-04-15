#ifndef UTILITY_H
#define UTILITY_H
#include <pari/pari.h>

ulong sigbits(GEN x);
/* Significant bits.
 * Input:   Integer x (assumed to be normalized).
 * Output:  The number of significant bits in x, i.e floor(log_2(x)) + 1.
*/

GEN split(GEN h, GEN m);
/* Split a number.
 * Input:   Integer h; positive integer m.
 * Output:  [r, s], such that r and s are positive integers with
            m = r*s, gcd(r, h) = 1 and any prime divisor of s divides h.
*/

#endif