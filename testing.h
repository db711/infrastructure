#ifndef TESTING_H
#define TESTING_H
#include <pari/pari.h>
#include "real_quadratic_orders.h"
#include "fp_representations.h"
#include "compact_representations.h"

long timedtest(GEN (*f)(GEN O, long prec, long flag), int i, int n);
/* Timed tests.
 * Input:   Function to time f (regulatorcf/regulatorshanks);
            integers i and n.
 * Output:  The time it took in ms.
 * The function tests n random squarefree d in the range [10^i, 10^(i+1)).
 * Assumes that the RNG of PARI is already set.
*/

void testcr(int i, int n);
/* Test compact representations.
 * Input:   Integers i and n.
 * Output:  Nothing.
 * The function tests n random squarefree d in the range [10^i, 10^(i+1))
   and outputs an error message if a compact representation is computed wrong.
 * Assumes that the RNG of PARI is already set.
*/

#endif