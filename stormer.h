#ifndef STORMER_H
#define STORMER_H
#include <pari/pari.h>

GEN stormer_gen(long length, GEN d, GEN lb, GEN ub, GEN bv, long* h, long l);
/* Stormer generator. 
 * Input:   length (of the considered list of primes);
            d as a starting value in the nodes;
            lb (lower bound) for d in the nodes;
            ub (upper bound) for d in the nodes;
            bv (bit vector) as a starting value, having the same length as lop (or NULL for the first one);
            h, bv's effective height or 0;
            l (least effective height).
 * Output:  A singly linked list starting at [bv, d, prev] (if bv != NULL, this d is different than the input)
            with backlinks up to the root.
            It is ensured that each d stored in a node is not smaller than lb and not greater than ub.
            It is ensured that the first bv has effective height at least l.
            (if bv is supplied, it is not checked whether it is malformed with respect to lb, ub, h, l).
            Returns NULL if d > ub.
*/

GEN stormer_next(GEN node, long length, GEN ub, long* h, long l, long m);
/* Stormer next.
 * Input:   node (as returned by stormer_gen or stormer_next);
            length (of the considered list of primes);
            ub (upperbound) for sol in the nodes;
            h, storing the effective height of current leaf;
            l, m bounding the effective height: l <= h <= m is maintained (i.e. needs to be true when first calling this).
 * Output:  A singly linked list [bv, d, prev], starting at the next node, going back to root;
            This function overwrites the part of the PARI stack taken up by the singly linked list returned by stormer_gen
            and restores avma before returning.
            For some reason if you want to create (and keep) data between calls of stormer_next,
            you have to leave space between the singly linked list and the data, i.e. use it something like this:
                stormer = stormer_gen(...);
                set_avma (avma - 64); // possibly adjust
                // create data
                stormer_next(stormer, ...);
                // create more data 
                stormer_next(stormer, ...);
                // ...
            Returns NULL if there are no more leaves to return.
*/

#endif