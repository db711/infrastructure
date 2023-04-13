#ifndef UTILITY_H
#define UTILITY_H
#include <pari/pari.h>

long sigbits(GEN x);
/* Significant bits.
 * Input:   Integer x (assumed to be normalized).
 * Output:  The number of significant bits in x, i.e floor(log_2(x)) + 1.
*/

#endif