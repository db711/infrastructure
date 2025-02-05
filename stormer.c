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
    if (typ(sol) != t_REAL) pari_err_TYPE("createnode",sol);
    GEN res = cgetg(3,t_VEC);
    gel(res,1) = shallowcopy(bv);
    gel(res,2) = gcopy(sol);
    return res;
}

GEN
leftchild(GEN node, GEN lop)
{
    if (lg(gel(node,1)) >= lg(lop)) return NULL; // node is a leaf
    GEN res = cgetg(3,t_VEC);
    gel(res,1) = vec_append(gel(node,1),gen_0);
    gel(res,2) = gcopy(gel(node,2));
    return res;
}

GEN
rightchild(GEN node, GEN lop)
{
    if (lg(gel(node,1)) >= lg(lop)) return NULL; // node is a leaf
    GEN res = cgetg(3,t_VEC);
    pari_sp ltop;
    gel(res,1) = vec_append(gel(node,1),gen_1);
    ltop = avma; gel(res,2) = gerepileupto(ltop,addrr(gel(node,2),gel(lop,lg(gel(res,1))-1)));
    return res;
}

static void
branch_print(GEN node, GEN lop, GEN ub)
{
    if (node == NULL) return;
    if(cmprr(gel(node,2),ub) > 0)  pari_printf("%Ps\n",gel(node,1)); //need to get rid of nodes after this
    else
    {
        branch_print(rightchild(node,lop),lop,ub);
        branch_print(leftchild(node,lop),lop,ub);
    }
    return;
}

void
printfailures(ulong B, ulong top, long prec)
{
    GEN node, ub, lop, vec;
    pari_sp ltop = avma, av;
    av = avma; ub = gerepileupto(av,logr_abs(mulir(gen_2,addir(gen_1,gcosh(mulis(gen_2,top),prec))))); //can this be simplified?
    lop = logsofprimes(B,prec);
    av = avma;
    vec = cgetg(1,t_VEC);
    node = gerepileupto(av,createnode(vec,itor(gen_0,prec)));
    pari_printf("%Ps\n",lop);
    pari_printf("upper bound: %Ps\n\n",ub);
    branch_print(node,lop,ub);
    set_avma(ltop);
    return;
}