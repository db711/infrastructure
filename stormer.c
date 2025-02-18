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

static inline int
stormeri_write_txt_branch(GEN node, GEN lop, GEN ub, int ft, FILE* fptr)
{
    pari_sp ltop = avma;
    if (cmprr(gel(node,2),ub) > 0) return 1; // stops computing this branch
    else 
    {
        if (lg(gel(node,1)) >= lg(lop)) pari_fprintf(fptr,"%Ps\n",bits_to_int(gel(node,1),lg(gel(node,1))-1));
        else
        {
            stormeri_write_txt_branch(leftchild(node,lop),lop,ub,ft,fptr); set_avma(ltop);
            if (ft == 0) ft = stormeri_write_txt_branch(rightchild(node,lop),lop,ub,ft,fptr); set_avma(ltop);
        }
    }
    set_avma(ltop);
    return 0;
}

void
stormeri_write_txt(GEN sol, ulong B, ulong top, long prec, FILE* fptr)
{
    GEN node, ub, lop, vec;
    pari_sp ltop = avma, av;
    av = avma; ub = gerepileupto(av,logr_abs(mulir(gen_2,addir(gen_1,gcosh(mulis(gen_2,top),prec))))); //can this be simplified?
    lop = logsofprimes(B,prec);
    av = avma;
    vec = cgetg(1,t_VECSMALL);
    node = gerepileupto(av,createnode(vec,sol));
    stormeri_write_txt_branch(node,lop,ub,0,fptr);
    set_avma(ltop);
    return;
}