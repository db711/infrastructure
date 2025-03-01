#include "stormer.h"

GEN
logsofprimes(ulong B, long prec)
{
    return glog(primes((long)uprimepi(B)),prec);
}

static inline GEN
createnode(GEN bv, GEN sol, GEN prev)
/* Create node.
 * Input:   bit vector bv \in {0, 1}^*;
            real number sol;
            other node prev or NULL.
 * Output:  [bv, sol, prev].
            (bv and sol are copied, prev is stored as a pointer only)
*/
{
    if (typ(bv) != t_VECSMALL) pari_err_TYPE("createnode",bv);
    if (typ(sol) != t_REAL) pari_err_TYPE("createnode",sol);
    if (prev != NULL && typ(prev) != t_VEC) pari_err_TYPE("createnode",prev);
    GEN res = cgetg(4,t_VEC);
    gel(res,1) = gcopy(bv);
    gel(res,2) = gcopy(sol);
    gel(res,3) = prev;
    return res;
}

static inline int
isleaf (GEN node, GEN lop)
/* is leaf?
 * Input:   node [bv, sol, prev] as created by createnode or leftchild/rightchild,
            lop (as returned by logsofprimes).
 * Output:  1 if the node is a leaf, 0 otherwise.
*/
{
    if (lg(gel(node,1)) >= lg(lop)) return 1;
    return 0;
}

static inline int
isleftchild (GEN node)
/* is left child?
 * Input:   node [bv, sol, prev] as created by createnode or leftchild/rightchild,
 * Output:  1 if the node is a left child, 0 otherwise.
*/
{
    if (lg(gel(node,1)) > 1 && (long)gel(gel(node,1),lg(gel(node,1))-1) == 0) return 1;
    return 0;
}

static inline GEN
leftchild(GEN node)
/* Left child.
 * Input:   node [bv, sol, prev] (as returned by createnode or *child);
 * Output:  Node [vec_append(bv,0), sol, node].
*/
{
    GEN res = cgetg(4,t_VEC);
    gel(res,1) = vecsmall_append(gel(node,1),0);
    gel(res,2) = gcopy(gel(node,2));
    gel(res,3) = node;
    return res;
}

static inline GEN
rightchild(GEN node, GEN lop)
/* Right child.
 * Input:   node [bv, sol, prev] (as returned by createnode or *child);
            lop (as returned by logsofprimes);
 * Output:  Node [vec_append(bv,1), sol+lop(length(lop)-length(bv)), node];
*/
{
    GEN res = cgetg(4,t_VEC);
    pari_sp ltop;
    if (lg(gel(node,1)) >= lg(lop)) pari_err_COMPONENT("rightchild",">=",gel(node,1),lop);
    gel(res,1) = vecsmall_append(gel(node,1),1);
    ltop = avma; gel(res,2) = gerepileupto(ltop,addrr(gel(node,2),gel(lop,lg(lop)+1-lg(gel(res,1)))));
    gel(res,3) = node;
    return res;
}

static inline int
checkrightchild(GEN node, GEN lop, GEN ub)
/* Check right child.
 * Input:   node [bv, sol, prev] (as returned by createnode or *child);
            lop (as returned by logsofprimes);
            ub (upperbound) for sol in the nodes (may be t_INFINITY);
 * Output:  1 if the rightchild of node is within the ub, 0 otherwise.
*/
{
    pari_sp ltop = avma;
    int ret;
    if (gcmp(addrr(gel(node,2),gel(lop,lg(lop)-lg(gel(node,1)))),ub) > 0) ret = 0;
    else ret = 1;
    set_avma(ltop); return ret;
}

GEN
stormer_gen(GEN lop, GEN sol, GEN ub, GEN bv)
{ 
    if (typ(lop) != t_VEC) pari_err_TYPE("stormer_gen", lop);
    if (typ(sol) != t_REAL) pari_err_TYPE("stormer_gen",sol);
    if (bv != NULL && typ(bv) != t_VECSMALL) pari_err_TYPE("stormer_gen",bv);
    GEN node;
    pari_sp av = avma;
    ulong i;
    if (gcmp(sol,ub) > 0) return NULL;
    node = gerepileupto(av,createnode(cgetg(1,t_VECSMALL),sol,NULL));
    if (bv == NULL) while (!isleaf(node,lop)) node = leftchild(node);
    else
    {
        for (i = 1; i < lg(bv); i++)
        {
            if ((long)gel(bv,i) == 1) node = rightchild(node,lop);
            else node = leftchild(node);
        }
    }
    return node;
}

GEN 
stormer_next(GEN node, GEN lop, GEN ub)
{ 
    GEN prev;
    pari_sp ltop = avma;
    do
    {
        prev = node;
        node = gel(node,3);
        if (node == NULL) return NULL;
        if (isleftchild(prev) && checkrightchild(node,lop,ub))
        {
            set_avma((pari_sp)gel(node,2)); 
            node = rightchild(node,lop);
            while (!isleaf(node,lop)) node = leftchild(node);
            set_avma(ltop);
            return node;
        }
    } while (1);
}