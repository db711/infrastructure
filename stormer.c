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
    gel(res,1) = gcopy(bv);
    gel(res,2) = gcopy(sol);
    return res;
}

int
isleaf (GEN node, GEN lop)
{
    if (lg(gel(node,1)) >= lg(lop)) return 1;
    return 0;
}

GEN
leftchild(GEN node)
{
    GEN res = cgetg(3,t_VEC);
    gel(res,1) = vecsmall_append(gel(node,1),0);
    gel(res,2) = gcopy(gel(node,2));
    return res;
}

GEN
rightchild(GEN node, GEN lop)
{
    GEN res = cgetg(3,t_VEC);
    pari_sp ltop;
    gel(res,1) = vecsmall_append(gel(node,1),1);
    ltop = avma; gel(res,2) = gerepileupto(ltop,addrr(gel(node,2),gel(lop,lg(lop)+1-lg(gel(res,1)))));
    return res;
}

sstack
stormer_gen(GEN lop, GEN sol)
{ // TODO: supply a starting bv and incorporate sol
    sstack stormer = {0};
    stormer.top = avma;
    stormer.bot= (pari_sp)gerepileupto(stormer.top,createnode(cgetg(1,t_VECSMALL),sol));
    while (!isleaf((GEN)stormer.bot,lop))
    {
        rightchild((GEN)stormer.bot,lop);
        stormer.bot = (pari_sp)leftchild((GEN)stormer.bot);
    }
    return stormer;
}

sstack
stormer_next(sstack stormer, GEN lop)
{ // TODO: incorporate sol
    GEN prev, curr = (GEN)stormer.bot;
    if ((long)gel(gel(curr,1),lg(gel(curr,1))-1) == 1) 
    {
        prev = curr;
        do { prev = (GEN)((pari_sp)prev+gsizebyte(prev)-sizeof(long)); } while ((long)gel(gel(prev,1),lg(gel(prev,1))-1) == 1);
        prev = (GEN)((pari_sp)prev+gsizebyte(prev));
        stormer.bot = (pari_sp)gerepile((pari_sp)(gel(prev,2)),(pari_sp)gel(curr,2),prev);
        while (!isleaf((GEN)stormer.bot,lop))
        {
            rightchild((GEN)stormer.bot,lop);
            stormer.bot = (pari_sp)leftchild((GEN)stormer.bot);
        }
    }
    else
    {
        prev = (GEN)(stormer.bot+gsizebyte((GEN)stormer.bot));
        stormer.bot = (pari_sp)gerepile((pari_sp)(gel(prev,2)),(pari_sp)gel(curr,2),prev);
    }
    return stormer;
}