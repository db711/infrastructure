#ifndef STORMER_H
#define STORMER_H
#include <pari/pari.h>

GEN logsofprimes(ulong B, long prec);
/* Logarithms of primes.
 * Input:   Smoothness bound B;
            precision prec.
 * Output:  log(primes(primepi(B))).
*/

GEN stormer_gen(GEN lop, GEN sol, GEN ub, GEN bv);
/* Stormer generator. 
 * Input:   lop (as returned by logsofprimes);
            sol (sum of logs) as a starting value in the nodes;
            ub (upperbound) for sol in the nodes (may be t_INFINITY);
            bv (bit vector) as a starting value, having the same length as lop (NULL is treated as [0, ..., 0]). 
 * Output:  A singly linked list starting at [bv, sol, prev]
            with backlinks up to the root.
            It is ensured that each sol stored in a node is not greater than ub.
            (bv has to have this property, this is not checked)
            Returns NULL if sol > ub.
*/

GEN stormer_next(GEN node, GEN lop, GEN ub);
/* Stormer next.
 * Input:   node with (as returned by stormer_gen or stormer_next);
            lop (as returned by logsofprimes);
            ub (upperbound) for sol in the nodes (may be t_INFINITY);
 * Output:  A singly linked list [bv, sol, prev], starting at the next node, going back to root;
            This function overwrites the part of the PARI stack taken up by the singly linked list returned by stormer_gen
            and restores avma before returning.
            Returns NULL if there are no more leaves to return.
*/

#endif