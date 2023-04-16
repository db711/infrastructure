#include "compact_representations.h"
#include "utility.h"

GEN expandcr(GEN O, GEN cr)
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
    pari_sp ltop = avma, av;
    GEN D, G, res, tmp;
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
        gel(G,1) = modii(gmael3(cr,1,1,1),m);
        gel(G,2) = modii(gmael3(cr,1,1,2),m);
    }
    for (i = 2; i < l-2; i++)
    {
        D = gmul(gmael(cr,i+1,2),gsqr(D));
        if (sw) G = mulqig(O,gmael(cr,i,1),mulqig(O,G,G));
        else
        {
            tmp = mulqig(O,G,G);
            gel(tmp,1) = modii(gel(tmp,1),m);
            gel(tmp,2) = modii(gel(tmp,2),m);
            tmp = mulqig(O,tmp,gmael(cr,i,1));
            G = cgetg(3,t_VEC);
            gel(G,1) = modii(gel(tmp,1),m);
            gel(G,2) = modii(gel(tmp,2),m);
        }
        gerepileall(av,2,&D,&G);
    }
    D = gsqr(D);
    if (sw) G = mulqig(O,G,G);
    else
    {
        tmp = mulqig(O,G,G);
        G = cgetg(3,t_VEC);
        gel(G,1) = modii(gel(tmp,1),m);
        gel(G,2) = modii(gel(tmp,2),m);
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
