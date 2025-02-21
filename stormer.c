#include "stormer.h"

GEN
logsofprimes(ulong B, long prec)
{
    return glog(primes((long)uprimepi(B)),prec);
}

GEN
createnode(GEN bv, GEN sol, GEN prev)
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

int
isleaf (GEN node, GEN lop)
{
    if (lg(gel(node,1)) >= lg(lop)) return 1;
    return 0;
}

GEN
leftchild(GEN node, GEN prev)
{
    GEN res = cgetg(4,t_VEC);
    gel(res,1) = vecsmall_append(gel(node,1),0);
    gel(res,2) = gcopy(gel(node,2));
    gel(res,3) = prev;
    return res;
}

GEN
rightchild(GEN node, GEN lop, GEN prev)
{
    GEN res = cgetg(4,t_VEC);
    pari_sp ltop;
    gel(res,1) = vecsmall_append(gel(node,1),1);
    ltop = avma; gel(res,2) = gerepileupto(ltop,addrr(gel(node,2),gel(lop,lg(lop)+1-lg(gel(res,1)))));
    gel(res,3) = prev;
    return res;
}

GEN
stormer_gen(GEN lop, GEN sol, GEN ub, GEN bv)
{ // TODO: incorporate bv
    if (typ(lop) != t_VEC) pari_err_TYPE("stormer_gen", lop);
    if (typ(sol) != t_REAL) pari_err_TYPE("stormer_gen",sol);
    if (bv != NULL && typ(bv) != t_VECSMALL) pari_err_TYPE("stormer_gen",bv);
    GEN rc;
    pari_sp av;
    if (gcmp(sol,ub) > 0) return NULL;
    pari_sp ltop = avma;
    GEN node = gerepileupto(ltop,createnode(cgetg(1,t_VECSMALL),sol,NULL));
    //need cases depending on bv = NULL or not
    while (!isleaf(node,lop)) 
    {
        av = avma;
        rc = rightchild(node,lop,node);
        if (gcmp(gel(rc,2),ub) > 0)
        {
            set_avma(av);
            node = leftchild(node,node);
        }
        else node = leftchild(node,rc);
    }
    return node;
}

GEN 
stormer_next(GEN node, GEN lop, GEN ub, GEN *old)
{ // TODO: how to recognize we reached the end?
    GEN rc;
    pari_sp av;
    *old = gel(node,3);
    while (lg(gel(node,1)) > lg(gel(*old,1)))
    {
        node = *old;
        *old = gel(node,3);
    }
    node = *old;
    while (!isleaf(node,lop)) // code duplication...
    {
        av = avma;
        rc = rightchild(node,lop,node);
        if (gcmp(gel(rc,2),ub) > 0)
        {
            set_avma(av);
            node = leftchild(node,node);
        }
        else node = leftchild(node,rc);
    }
    return node;
}