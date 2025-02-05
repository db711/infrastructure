#include "stormer.h"

GEN
logsofprimes(ulong B, ulong prec)
{
    return glog(primes((long)uprimepi(B)),prec);
}

GEN
createnode(GEN bv, GEN sol)
{
    //TODO: error handling
    GEN res = cgetg(3,t_VEC);
    gel(res,1) = shallowcopy(bv);
    gel(res,2) = gcopy(sol);
    return res;
}

GEN
leftchild(GEN node)
{
    //TODO: error handling
    GEN res = cgetg(3,t_VEC);
    gel(res,1) = vec_append(gel(node,1),gen_0);
    gel(res,2) = gcopy(gel(node,2));
    return res;
}

GEN
rightchild(GEN node, GEN lop)
{
    //TODO: error handling
    GEN res = cgetg(3,t_VEC);
    pari_sp ltop;
    gel(res,1) = vec_append(gel(node,1),gen_1);
    ltop = avma; gel(res,2) = gerepileupto(ltop,addrr(gel(node,2),gel(lop,lg(gel(res,1))-1)));
    return res;
}