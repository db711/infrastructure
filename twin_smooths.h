#ifndef TWIN_SMOOTHS_H
#define TWIN_SMOOTHS_H
#include "compact_representations.h"
#define LOWER_BOUND_I 240 // lower bound (in bits) for twin smooths NIST-I
#define UPPER_BOUND_I 260 // upper bound (in bits) for twin smooths NIST-I
#define LOWER_BOUND_III 370 // lower bound (in bits) for twin smooths NIST-III
#define UPPER_BOUND_III 390 // upper bound (in bits) for twin smooths NIST-III
#define LOWER_BOUND_V 500 // lower bound (in bits) for twin smooths NIST-V
#define UPPER_BOUND_V 520 // upper bound (in bits) for twin smooths NIST-V

GEN twin_smooth_d(GEN S, GEN d, ulong m, GEN P);
/* Twin smooths (for value d)
 * Input:   Vector S of primes below a bound B;
            positive, even, squarefree integer d;
            upper bound m.
 * Output:  Vector containing all x, such that x(x+1) is B-smooth, corresponding to d.
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
 * Output:  If the regulator R (or an integer multiple) is in [A, B]:
                [[f, p], [b, d, k]], where (b, d, k) is an (f, p) representation of (1) in O with b = (1) and k \approx R/log2. 
            Otherwise:
                NULL.
*/

GEN twin_smooth_range_d_small(ulong B, GEN d, ulong bot, ulong top, ulong m);
/* Twin smooths (for value d) in range, small.
 * Input:   Smoothness bound B;
            positive, even, squarefree integer d;
            range [bot, top];
            upper bound m.
 * Output:  Vector containing all x, such that x(x+1) is B-smooth corresponding to d.
            Here we only perform a check if the regulator (or a multiple of it) corresponding to d is in [bot, top] and otherwise skip the discriminant, this ensures that x does not have much more than top bits.
            The value m decides how many solutions (including the fundamental) are checked, the algorithm is supposed to be run with m < 5.
*/

GEN regulator_cryptographic(GEN O);
/* Regulator of cryptographic size.
 * Input:   Real quadratic order O (as output by rqoinit);
 * Output:  Vector containing the x-coordinate of units in O having the correct number of bits (per NIST security level)
            and the chance to be prime (may be empty).
*/

#endif