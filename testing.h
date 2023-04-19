#ifndef TESTING_H
#define TESTING_H
#include <pari/pari.h>
#include "real_quadratic_orders.h"
#include "fp_representations.h"
#include "compact_representations.h"
#include "utility.h"

ulong timedtest(GEN (*f)(GEN O, ulong prec, long flag), ulong i, ulong n);
/* Timed tests.
 * Input:   Function to time f (regulatorcf/regulatorshanks);
            integers i and n.
 * Output:  The time it took in ms.
 * The function tests n random squarefree d in the range [10^i, 10^(i+1)).
 * Assumes that the RNG of PARI is already set.
*/

ulong testcr(ulong i, ulong n, GEN m);
/* Test compact representations.
 * Input:   Integers i and n as well as m.
 * Output:  Number of times that the compact representation was smaller than the standard decimal representation.
 * The function tests n random squarefree d in the range [10^i, 10^(i+1))
   and outputs an error message if a compact representation is computed wrong.
   The integer m is used to compute a representation where at most the factor that appears with to the power 1 shares any factors with m.
 * Assumes that the RNG of PARI is already set.
*/

ulong testcrsmoothpart(ulong i, ulong n, ulong B, ulong c, ulong sw);
/* Test compact representation smooth part.
 * Input:   Integers i and n as well as B;
            c = {1, 2};
            sw = {1, 2, 3}.
 * Output:  The time it took in ms.
 * The function tests either crsmoothpart or crsmoothpart2 or crsmoothpart_alt (depending on wheter sw = 1 or sw = 2 or sw = 3)
   by computing the B-smooth part of the c'th coordinate of of a compact representation of the fundamental unit (x + y*sqrt(d))/s,
   which is randomly chose in [10^i, 10^(i+1)); this is repeated n times.
 * Assumes that the RNG of PARI is already set.
*/

#endif