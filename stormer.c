#include "stormer.h"

GEN
logsofprimes(ulong B, long prec)
{
    return glog(primes((long)uprimepi(B)),prec);
}

GEN
createnode(GEN bv, GEN sol)
{
    if (typ(bv) != t_VECSMALL) pari_err_TYPE("createnode",bv);
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
    gel(res,1) = vecsmall_append(gel(node,1),0);
    gel(res,2) = gcopy(gel(node,2));
    return res;
}

GEN
rightchild(GEN node, GEN lop)
{
    if (lg(gel(node,1)) >= lg(lop)) return NULL; // node is a leaf
    GEN res = cgetg(3,t_VEC);
    pari_sp ltop;
    gel(res,1) = vecsmall_append(gel(node,1),1);
    ltop = avma; gel(res,2) = gerepileupto(ltop,addrr(gel(node,2),gel(lop,lg(gel(res,1))-1)));
    return res;
}

static int
stormer_branch(GEN node, GEN lop, GEN ub, int ft)
{
    pari_sp av = avma;
    if (cmprr(gel(node,2),ub) > 0) return 1; // stops computing this branch
    else 
    {
        if (lg(gel(node,1)) >= lg(lop)) pari_printf("%Ps\n",node); // leaf node
        else
        {
            if (ft == 0) av = avma; ft = stormer_branch(rightchild(node,lop),lop,ub,ft); set_avma(av);
            av = avma; stormer_branch(leftchild(node,lop),lop,ub,ft); set_avma(av);
        }
    }
    return 0;
}

void
stormer_print(ulong B, ulong top, long prec)
{
    GEN node, ub, lop, vec;
    pari_sp ltop = avma, av;
    av = avma; ub = gerepileupto(av,logr_abs(mulir(gen_2,addir(gen_1,gcosh(mulis(gen_2,top),prec))))); //can this be simplified?
    lop = logsofprimes(B,prec);
    av = avma;
    vec = cgetg(1,t_VECSMALL);
    node = gerepileupto(av,createnode(vec,itor(gen_0,prec)));
    pari_printf("%Ps\n",lop);
    pari_printf("upper bound: %Ps\n\n",ub);
    stormer_branch(node,lop,ub,0);
    set_avma(ltop);
    return;
}