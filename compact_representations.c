#include "compact_representations.h"
#include "utility.h"

GEN 
expandcr(GEN O, GEN cr)
{
    pari_sp ltop = avma, av, av2;
    GEN theta, fac;
    long i, j, l;
    l = lg(cr)-1;
    theta = gmael(cr,l,1);
    for (i = l-1; i > 0; i--)
    {
        av2 = avma;
        fac = cgetg(3,t_VEC);
        av = avma; gel(fac,1) = gerepileupto(av,gdiv(gmael3(cr,i,1,1),gmael(cr,i+1,2)));
        av = avma; gel(fac,2) = gerepileupto(av,gdiv(gmael3(cr,i,1,2),gmael(cr,i+1,2)));
        for (j = 1; j < l-i; j++) fac = gerepileupto(av2,mulqig(O,fac,fac));
        theta = gerepileupto(ltop,mulqig(O,theta,fac));
    }
    return gerepileupto(ltop,theta);
}

GEN 
crmodm (GEN O, GEN cr, GEN m)
{
    if (typ(m) != t_INT) pari_err_TYPE("crmodm",m);
    pari_sp ltop = avma, av;
    GEN D, G, res, tmp, m_;
    m_ = mulii(gen_2,m);
    ulong l = lg(cr), i, sw = 0; 
    av = avma;
    D = mkintmod(modii(gmael(cr,2,2),m),m);
    G = cgetg(3,t_VEC);
    if (!cmpii(gel(O,3),gen_1) || Mod2(m)) sw = 1;
    if (sw)
    {
        gel(G,1) = mkintmod(modii(gmael3(cr,1,1,1),m),m);
        gel(G,2) = mkintmod(modii(gmael3(cr,1,1,2),m),m);
    }
    else
    {
        gel(G,1) = modii(gmael3(cr,1,1,1),m_);
        gel(G,2) = modii(gmael3(cr,1,1,2),m_);
    }
    for (i = 2; i < l-2; i++)
    {
        D = gmul(gmael(cr,i+1,2),gsqr(D));
        if (sw) G = mulqig(O,gmael(cr,i,1),mulqig(O,G,G));
        else
        {
            tmp = mulqig(O,G,G);
            gel(tmp,1) = modii(gel(tmp,1),m_);
            gel(tmp,2) = modii(gel(tmp,2),m_);
            tmp = mulqig(O,tmp,gmael(cr,i,1));
            G = cgetg(3,t_VEC);
            gel(G,1) = modii(gel(tmp,1),m_);
            gel(G,2) = modii(gel(tmp,2),m_);
        }
        gerepileall(av,2,&D,&G);
    }
    D = gsqr(D);
    if (sw) G = mulqig(O,G,G);
    else
    {
        tmp = mulqig(O,G,G);
        G = cgetg(3,t_VEC);
        gel(G,1) = modii(gel(tmp,1),m_);
        gel(G,2) = modii(gel(tmp,2),m_);
    }
    tmp = mulqig(O,G,mulqidivn(O,gmael(cr,l-1,1),gmael(cr,l-2,1),gmael(cr,l-1,2)));
    if (sw) 
    {
        res = cgetg(3,t_VEC);
        av = avma; gel(res,1) = gerepileupto(av,lift(gdiv(gel(tmp,1),D)));
        av = avma; gel(res,2) = gerepileupto(av,lift(gdiv(gel(tmp,2),D)));
    }
    else
    {
        res = cgetg(3,t_VEC);
        av = avma; gel(res,1) = gerepileupto(av,modii(mulii(gel(tmp,1),ginvmod(lift(D),m)),m));
        av = avma; gel(res,2) = gerepileupto(av,modii(mulii(gel(tmp,2),ginvmod(lift(D),m)),m));
    }
    return gerepileupto(ltop,res);
}

GEN 
crmodm_alt(GEN O, GEN cr, GEN m)
{
    pari_sp ltop = avma, av, av2, av3;
    GEN D, G, G_, r, r_, d, n, tmp, res;
    ulong l = lg(cr), i;
    av2 = avma;
    D = mulii(powii(gel(O,3),mulis(gen_2,l-2)),mulii(m,gmael(cr,l-1,2)));
    for (i = 2; i < l-1; i++) D = gerepileupto(av2,mulii(D,sqri(gmael(cr,i,2))));
    G = cgetg(3,t_VEC);
    gel(G,1) = modii(gmael3(cr,1,1,1),D);
    gel(G,2) = modii(gmael3(cr,1,1,2),D);
    r = gen_1;
    for (i = 2; i < l-1; i++) 
    {
        tmp = cgetg(3,t_VEC);
        av = avma; gel(tmp,1) = gerepileupto(av,addii(sqri(gel(G,1)),mulii(sqri(gel(G,2)),gel(O,1))));
        av = avma; gel(tmp,2) = gerepileupto(av,mulii(gen_2,mulii(gel(G,1),gel(G,2))));
        G_ = cgetg(3,t_VEC);
        av = avma; gel(G_,1) = gerepileupto(av,addii(mulii(gel(tmp,1),gmael3(cr,i,1,1)),mulii(mulii(gel(tmp,2),gmael3(cr,i,1,2)),gel(O,1))));
        av = avma; gel(G_,2) = gerepileupto(av,addii(mulii(gel(tmp,1),gmael3(cr,i,1,2)),mulii(gel(tmp,2),gmael3(cr,i,1,1))));
        r = sqri(r);
        av = avma; tmp = gerepileupto(av,mulii(sqri(gel(O,3)),mulii(sqri(gmael(cr,i,2)),r)));
        av = avma; d = gerepileupto(av,gcdii(tmp,gcdii(gel(G_,1),gel(G_,2))));
        r = diviiexact(tmp,d);
        av = avma; D = gerepileupto(av,diviiexact(D,gcdii(D,d)));
        G = cgetg(3,t_VEC);
        gel(G,1) = modii(diviiexact(gel(G_,1),d),D);
        gel(G,2) = modii(diviiexact(gel(G_,2),d),D);
        av3 = avma;
        av = avma; r_ = gerepileupto(av,diviiexact(r,gcdii(r,D)));
        while (cmpii(gcdii(r_,D),gen_1) > 1) r_ = gerepileupto(av3,diviiexact(r_,gcdii(r_,D)));
        if (cmpii(r_,gen_1) > 1)
        {
            n = ginvmod(r_,D);
            r = diviiexact(r,r_);
            G_ = cgetg(3,t_VEC);
            av = avma; gel(G_,1) = gerepileupto(av,modii(mulii(gel(G,1),n),D));
            av = avma; gel(G_,1) = gerepileupto(av,modii(mulii(gel(G,1),n),D));
            G = G_;
        }
        gerepileall(av2,3,&G,&r,&D);
    }
    av = avma; D = gerepileupto(av,diviiexact(D,gcdii(D,mulii(sqri(gel(O,3)),gmael(cr,l-1,2)))));
    if (cmpii(D,m)) pari_err_BUG("crmodm_alt"); //crmodm_alt(O,cr,mulii(m,diviiexact(m,D)));
    G = mulqig(O,G,gmael(cr,l-1,1));
    res = cgetg(3,t_VEC);
    av = avma; gel(res,1) = gerepileupto(av,modii(diviiexact(gel(G,1),gmael(cr,l-1,2)),D));
    av = avma; gel(res,2) = gerepileupto(av,modii(diviiexact(gel(G,2),gmael(cr,l-1,2)),D));
    return gerepileupto(ltop,res);
}

ulong
crpval(GEN O, GEN cr, ulong p, ulong c)
{
    pari_sp ltop = avma;
    ulong k = 16, k_ = 0, m; //starting value of 16 good?
    GEN tmp;
    tmp = gel(crmodm(O,cr,powuu(p,k)),c);
    while (!cmpii(tmp,gen_0))
    {
        tmp = gel(crmodm(O,cr,powuu(p,k)),c);
        k_ = k;
        k *= 2;
        set_avma(ltop);
    }
    while (1) // v_p is now in [k_,k); do a binary search
    {
        m = k_ + (k-k_)/2;
        tmp = gel(crmodm(O,cr,powuu(p,m)),c);
        if (!cmpii(tmp,gen_0)) k_ = m;
        else k = m;
        set_avma(ltop);
        if (k_+1 == k) break;
    }
    return gc_long(ltop,k_);
}

ulong
crpval_alt(GEN O, GEN cr, ulong p, ulong c)
{ //code duplication...
    pari_sp ltop = avma;
    ulong k = 16, k_ = 0, m; //starting value of 16 good?
    GEN tmp;
    tmp = gel(crmodm_alt(O,cr,powuu(p,k)),c);
    while (!cmpii(tmp,gen_0))
    {
        tmp = gel(crmodm(O,cr,powuu(p,k)),c);
        k_ = k;
        k *= 2;
        set_avma(ltop);
    }
    while (1) // v_p is now in [k_,k); do a binary search
    {
        m = k_ + (k-k_)/2;
        tmp = gel(crmodm_alt(O,cr,powuu(p,m)),c);
        if (!cmpii(tmp,gen_0)) k_ = m;
        else k = m;
        set_avma(ltop);
        if (k_+1 == k) break;
    }
    return gc_long(ltop,k_);
}

GEN 
crsmoothpart(GEN O, GEN cr, GEN S, ulong c)
{
    pari_sp ltop = avma;
    GEN res = gen_1;
    ulong i;
    for (i = 1; i < lg(S); i++) res = gerepileupto(ltop,mulii(res,powuu((ulong)gel(S,i),crpval(O,cr,(ulong)gel(S,i),c))));
    return gerepileupto(ltop,res);
}

GEN
crsmoothpart2(GEN O, GEN b, GEN y, GEN q, GEN S, ulong c)
{
    pari_sp ltop = avma;
    GEN res = gen_1;
    ulong i;
    for (i = 1; i < lg(S); i++) res = gerepileupto(ltop,mulii(res,powuu((ulong)gel(S,i),crpval(O,gel(cr(O,b,y,q,stoi((ulong)gel(S,i))),2),(ulong)gel(S,i),c))));
    return gerepileupto(ltop,res);
}

GEN 
crsmoothpart_alt(GEN O, GEN cr, GEN S, ulong c)
{ //code duplication...
    pari_sp ltop = avma;
    GEN res = gen_1;
    ulong i;
    for (i = 1; i < lg(S); i++) res = gerepileupto(ltop,mulii(res,powuu((ulong)gel(S,i),crpval_alt(O,cr,(ulong)gel(S,i),c))));
    return gerepileupto(ltop,res);
}

long 
crnorm_sign(GEN O, GEN cr)
{
    pari_sp ltop = avma;
    return gc_long(ltop,signe(subii(sqri(gmael3(cr,lg(cr)-2,1,1)),mulii(gel(O,1),sqri(gmael3(cr,lg(cr)-2,1,2)))))*signe(gmael(cr,lg(cr)-1,2))*signe(subii(sqri(gmael3(cr,lg(cr)-1,1,1)),mulii(gel(O,1),sqri(gmael3(cr,lg(cr)-1,1,2))))));
}