#ifndef TWIN_SMOOTHS_H
#define TWIN_SMOOTHS_H
#include "compact_representations.h"

GEN twin_smooth_d(ulong B, GEN d, ulong m);
/* Twin smooths (for value d)
 * Input:   Positive integer B;
            positive, even, squarefree integer d;
            upper bound m.
 * Output:  Vector containing all x, such that x(x+1) is B-smooth, corresponding to d or NULL if no such exist.
*/

GEN twin_smooth(ulong B);
/* Twin smooths.
 * Input:   Positive integer B;
            filename filename.
 * Output:  Vector containing all x, such that x(x+1) is B-smooth.
*/

#endif