#ifndef UTILITY_H
#define UTILITY_H
#include <pari/pari.h>

ulong sigbits(GEN x);
/* Significant bits.
 * Input:   Integer x (assumed to be normalized).
 * Output:  The number of significant bits in x, i.e floor(log_2(x)) + 1.
*/

GEN mulqidivn (GEN O, GEN qi1, GEN qi2, GEN N); // Returns [x, y] = (x + y*sqrt(d))/(s*N) = ((a + b*sqrt(d))/s * (a_ + b_*sqrt(d))/s)/N, where qi1 = [a, b] and qi2 = [a_, b_].
GEN mulqig (GEN O, GEN qi1, GEN qi2); // Returns [x, y] = (x + y*sqrt(d))/s = (a + b*sqrt(d))/s * (a_ + b_*sqrt(d))/s, where qi1 = [a, b] and qi2 = [a_, b_].
GEN invqi(GEN qi); // Returns [a, -b], where qi = [a, b].

#endif