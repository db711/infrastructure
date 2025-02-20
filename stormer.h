#ifndef STORMER_H
#define STORMER_H
#include <pari/pari.h>

typedef struct {
    pari_sp bot; // current node
    pari_sp top; // root node
} sstack; // stormer stack
/* We think of this as a (separate) stack storing nodes.
 * This stack is initialized with stormer_gen, 
 * further elements are generated (i.e. bot is changed) by running stormer_next.
 * Currently it is required that the PARI stack is not modified 
 * (or at least cleaned, meaning avma = gel(bot,2))
 * between running stormer_gen and calling stormer_next 
 * (respectively between any subsequent calls of stormer_next).
 * Further calls to stormer_next after all leaves were generated are undefined.
*/ 

GEN logsofprimes(ulong B, long prec);
/* Logarithms of primes.
 * Input:   Smoothness bound B;
            precision prec.
 * Output:  log(primes(primepi(B))).
*/

GEN createnode(GEN bv, GEN sol);
/* Create node.
 * Input:   bit vector bv \in {0, 1}^*;
            real number sol.
 * Output:  [bv, sol].
*/

int isleaf (GEN node, GEN lop);
/* is leaf?
 * Input:   node [bv, sol] as created by createnode or leftchild/rightchild,
            lop (as returned by logsofprimes).
 * Output:  1 if the node is a leaf, 0 otherwise.
*/

GEN leftchild(GEN node);
/* Left child.
 * Input:   Node [bv, sol] (as returned by createnode or *child).
 * Output:  Node [vec_append(bv,0), sol]
            undefined behavior if node is already a leaf.
*/

GEN rightchild(GEN node, GEN lop);
/* Right child.
 * Input:   Node [bv, sol] (as returned by createnode or *child);
            lop (as returned by logsofprimes).
 * Output:  Node [vec_append(bv,1), sol+lop(length(bv))];
            undefined behavior if node is already a leaf.
*/

sstack stormer_gen(GEN lop, GEN sol);
/* Stormer generator. 
 * Input:   lop (as returned by logsofprimes);
            sol (sum of logs) as a starting value in the nodes. 
 * Output:  A sstack structure with bot set to the node [[0,...,0], sol].
*/

sstack stormer_next(sstack stormer, GEN lop);
/* Stormer next.
 * Input:   sstack structure stormer (as returned by stormer_gen or stormer_next);
            lop (as returned by logsofprimes). 
 * Output:  A sstack structure with bot set to the next leaf node.
            (Undefined behavior if the sstack structure already returned all leaves.)
*/

#endif