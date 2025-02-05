#include "stormer.h"

GEN
logsofprimes(ulong B, long prec)
{
    return glog(primes((long)uprimepi(B)),prec);
}

GEN
createnode(GEN bv, GEN sol)
{
    if (typ(bv) != t_VEC) pari_err_TYPE("createnode",bv);
    if (typ(sol) != t_REAL) pari_err_type("createnode",sol);
    GEN res = cgetg(3,t_VEC);
    gel(res,1) = shallowcopy(bv);
    gel(res,2) = gcopy(sol);
    return res;
}

GEN
leftchild(GEN node, GEN lop)
{
    if (lg(gel(node,1)) > lg(lop)) return NULL; // node is a leaf
    GEN res = cgetg(3,t_VEC);
    gel(res,1) = vec_append(gel(node,1),gen_0);
    gel(res,2) = gcopy(gel(node,2));
    return res;
}

GEN
rightchild(GEN node, GEN lop)
{
    if (lg(gel(node,1)) > lg(lop)) return NULL; // node is a leaf
    GEN res = cgetg(3,t_VEC);
    pari_sp ltop;
    gel(res,1) = vec_append(gel(node,1),gen_1);
    ltop = avma; gel(res,2) = gerepileupto(ltop,addrr(gel(node,2),gel(lop,lg(gel(res,1))-1)));
    return res;
}

void
printfailures(ulong B, ulong top, long prec)
{
    GEN ub,lop;
    pari_sp ltop = avma, av;
    av = avma; ub = gerepileupto(av,logr_abs(mulir(gen_2,addir(gen_1,gcosh(mulis(gen_2,top),prec))))); //can this be simplified?
    lop = logsofprimes(B,prec);
    set_avma(ltop);
    return;
}