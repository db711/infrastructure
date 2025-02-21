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
stormer_gen(GEN lop, GEN sol)
{ // TODO: supply a starting bv and incorporate sol
    pari_sp ltop = avma;
    GEN curr = gerepileupto(ltop,createnode(cgetg(1,t_VECSMALL),sol,NULL));
    while (!isleaf(curr,lop)) curr = leftchild(curr,rightchild(curr,lop,curr));
    return curr;
}

GEN 
stormer_next(GEN node, GEN lop, GEN *old)
{ // TODO: incorporate sol
    if ((long)gel(gel(node,1),lg(gel(node,1))-1) == 1) 
    {
        do { node = gel(node,3); } while ((long)gel(gel(node,1),lg(gel(node,1))-1) == 1);
        *old = node = gel(node,3);
        while (!isleaf(node,lop)) node = leftchild(node,rightchild(node,lop,node));
    }
    else *old = node = gerepile((pari_sp)gel(gel(node,3),2),(pari_sp)gel(node,2),gel(node,3));
    return node;
}