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
            (bv has to have this property)
            Returns NULL if sol > ub.
*/

GEN stormer_next(GEN node, GEN lop, GEN ub, GEN *old);
/* Stormer next.
 * Input:   node with (as returned by stormer_gen or stormer_next);
            lop (as returned by logsofprimes);
            ub (upperbound) for sol in the nodes (may be t_INFINITY);
            *old (just some GEN *, that will be overwritten).
 * Output:  A singly linked list [node, sol, prev], starting at the next node, going back to root;
            old now points to newest node in that list, that was already on the PARI stack.
            If the return value and *old differ after this function terminates, 
            new nodes have been added to the singly linked list.
            To ensure that the list remains continuous in memory, 
            all data that was added to the PARI stack before the call of this function has to be copied (and pointers updated) and memory has to be cleaned.
            To ensure this, use the function like this:
                GEN node, old;
                // call of node = stormer_gen or previous call of node = stormer_next is here
                // data is created
                pari_sp lbot = avma;
                node = stormer_next(node, lop, &old);
                if (stormer != old)
                {
                    // gcopy all the data you want to keep, updating the pointers
                    // alternatively if you know what you are doing, 
                    // you can use gerepileall instead of the following gerepile
                }
                node = gerepile((pari_sp)gel(old,2), lbot, node);
            It is ensured that each sol stored in a node is not greater than ub.
            If all valid leaves have already been returned, this returns NULL.
*/

#endif