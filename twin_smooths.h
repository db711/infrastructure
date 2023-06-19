#ifndef TWIN_SMOOTHS_H
#define TWIN_SMOOTHS_H
#include "compact_representations.h"

GEN twin_smooth_d(GEN S, GEN d, ulong m, GEN P);
/* Twin smooths (for value d)
 * Input:   Vector S of primes below a bound B;;
            positive, even, squarefree integer d;
            upper bound m.
 * Output:  Vector containing all x, such that x(x+1) is B-smooth, corresponding to d or NULL if no such exist.
*/

GEN twin_smooth(ulong B);
/* Twin smooths.
 * Input:   Positive integer B;
 * Output:  Vector containing all x, such that x(x+1) is B-smooth.
*/

GEN regulator_range(GEN O, ulong A, ulong B);
/* Regulator in range.
 * Input:   Real quadratic order O (as output by rqoinit);
            positive integers A < B.
 * Output:  If the regulator R is in [A, B]:
                [[f, p], [b, d, k]], where (b, d, k) is an (f, p) representation of (1) in O with b = (1) and k \approx R/log2. 
            Otherwise:
                NULL.
*/

#endif