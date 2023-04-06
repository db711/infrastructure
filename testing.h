#ifndef TESTING_H
#define TESTING_H
#include <pari/pari.h>

long timedtest(GEN (*f)(GEN O, long prec, long flag), int i, int n);
/*
 * Timed tests.
 * Input:   Function to time f (regulatorcf/regulatorshanks), int i, int n.
 * Output:  The time it took in ms.
 * The function tests n random squarefree d in the range [10^i, 10^(i+1));
 * assumes that the RNG of pari is already set.
*/

#endif