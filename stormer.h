#ifndef STORMER_H
#define STORMER_H
#include <pari/pari.h>

typedef struct {
    pari_sp bot; // points to root
    pari_sp top; // points to current element
    pari_ulong *sizes;
    int length;
} sstack; // stormer stack
/* We think of this struct as its own stack of nodes, i.e. [bv, sol].
 * bot points to the bottom (so always the root node),
 * top points to a (the next) leaf node.
 * sizes is a list, which stores the size of the nodes on the stack (in order),
 * length is the length of this list (equal to 2*(length of logs of primes)-1).
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
 * Output:  A sstack structure with top set to the node [[0,...,0], sol].
*/

sstack stormer_next(sstack stormer, GEN lop);
/* Stormer next.
 * Input:   sstack structure stormer (as returned by stormer_gen or stormer_next);
            lop (as returned by logsofprimes). 
 * Output:  A sstack structure with top set to the next leaf node.
*/

//void stormeri_write_txt(GEN sol, ulong B, ulong top, long prec, FILE* fptr);
/* Stormer indices write txt.
 * Input:   sum of logs sol (as a starting value);
            smoothness bound B;
            upper bound (in bits) top;
            precision prec;
            path path.
 * Output:  Writes the corresponding indices of the Stormer discriminants associated to (B, top) to the txt file path.
*/

#endif