#ifndef STORMER_H
#define STORMER_H
#include <pari/pari.h>

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

GEN leftchild(GEN node, GEN lop);
/* Left child.
 * Input:   Node [bv, sol] (as returned by createnode or *child);
            lop (as returned by logsofprimes).
 * Output:  Node [vec_append(bv,0), sol]
            or NULL, if node is already a leaf.
*/

GEN rightchild(GEN node, GEN lop);
/* Right child.
 * Input:   Node [bv, sol] (as returned by createnode or *child);
            lop (as returned by logsofprimes).
 * Output:  Node [vec_append(bv,1), sol+lop(length(bv))]
            or NULL, if node is already a leaf.
*/

void stormeri_write_txt(ulong B, ulong top, long prec, char* path);
/* Stormer indices write txt.
 * Input:   Smoothness bound B;
            upper bound (in bits) top;
            precision prec;
            path path.
 * Output:  Writes the corresponding indices of the Stormer discriminants associated to (B, top) to the txt file path.
*/

#endif