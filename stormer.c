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

static inline GEN
leftchild(GEN node, GEN prev)
/* Left child.
 * Input:   node [bv, sol, prev] (as returned by createnode or *child);
            node prev.
 * Output:  Node [vec_append(bv,0), sol, prev].
*/
{
    GEN res = cgetg(4,t_VEC);
    gel(res,1) = vecsmall_append(gel(node,1),0);
    gel(res,2) = gcopy(gel(node,2));
    gel(res,3) = prev;
    return res;
}

static inline GEN
rightchild(GEN node, GEN lop, GEN prev)
/* Right child.
 * Input:   node [bv, sol, prev] (as returned by createnode or *child);
            lop (as returned by logsofprimes);
            node prev.
 * Output:  Node [vec_append(bv,1), sol+lop(length(lop)-length(bv)), prev];
*/
{
    GEN res = cgetg(4,t_VEC);
    pari_sp ltop;
    if (lg(gel(node,1)) >= lg(lop)) pari_err_COMPONENT("rightchild",">=",gel(node,1),lop);
    gel(res,1) = vecsmall_append(gel(node,1),1);
    ltop = avma; gel(res,2) = gerepileupto(ltop,addrr(gel(node,2),gel(lop,lg(lop)+1-lg(gel(res,1)))));
    gel(res,3) = prev;
    return res;
}

GEN
stormer_gen(GEN lop, GEN sol, GEN ub, GEN bv)
{ 
    if (typ(lop) != t_VEC) pari_err_TYPE("stormer_gen", lop);
    if (typ(sol) != t_REAL) pari_err_TYPE("stormer_gen",sol);
    if (bv != NULL && typ(bv) != t_VECSMALL) pari_err_TYPE("stormer_gen",bv);
    GEN rc, node;
    pari_sp av = avma;
    ulong i;
    if (gcmp(sol,ub) > 0) return NULL;
    node = gerepileupto(avma,createnode(cgetg(1,t_VECSMALL),sol,NULL));
    if (bv == NULL)
    {
        while (!isleaf(node,lop)) 
        {
            av = avma; rc = rightchild(node,lop,node);
            if (gcmp(gel(rc,2),ub) > 0) node = gerepileupto(av,leftchild(node,node));
            else node = leftchild(node,rc);
        }
    }
    else
    {
        pari_printf("\ttrying to walk to %Ps\n",bv);
        for (i = 1; i < lg(bv); i++)
        {
            if ((long)gel(bv,i) == 1) node = rightchild(node,lop,node);
            else // code duplication...
            {
                av = avma; rc = rightchild(node,lop,node);
                if (gcmp(gel(rc,2),ub) > 0) node = gerepileupto(av,leftchild(node,node));
                else node = leftchild(node,rc);
            }
        }
    }
    return node;
}

GEN 
stormer_next(GEN node, GEN lop, GEN ub, GEN *old)
{ 
    GEN rc;
    pari_sp av;
    *old = gel(node,3);
    while (lg(gel(node,1)) > lg(gel(*old,1)))
    {
        node = *old;
        *old = gel(node,3);
        if (*old == NULL) return NULL;
    }
    node = *old;
    while (!isleaf(node,lop)) // code duplication...
    {
        av = avma; rc = rightchild(node,lop,node);
        if (gcmp(gel(rc,2),ub) > 0) node = gerepileupto(av,leftchild(node,node));
        else node = leftchild(node,rc);
    }
    return node;
}