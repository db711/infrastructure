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
    ulong l = lg(cr), i; 
    if (!cmpii(gel(O,3),gen_2) && !Mod2(m)) m = mulii(gen_2,m); //this case requires extra care
    av = avma;
    D = mkintmod(modii(gmael(cr,2,2),m),m);
    G = cgetg(3,t_VEC);
    if (Mod2(m))
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
        if (Mod2(m)) G = mulqig(O,gmael(cr,i,1),mulqig(O,G,G));
        else
        {
            tmp = mulqig(O,G,G);
            gel(tmp,1) = modii(gel(tmp,1),m);
            gel(tmp,2) = modii(gel(tmp,2),m);
            tmp = mulqig(O,tmp,G);
            G = cgetg(3,t_VEC);
            gel(G,1) = modii(gel(tmp,1),m);
            gel(G,2) = modii(gel(tmp,2),m);
        }
        gerepileall(av,2,&D,&G);
    }
    D = gsqr(D);
    if (Mod2(m)) G = mulqig(O,G,G);
    else
    {
        tmp = mulqig(O,G,G);
        G = cgetg(3,t_VEC);
        gel(G,1) = modii(gel(tmp,1),m);
        gel(G,2) = modii(gel(tmp,2),m);
    }
    tmp = mulqig(O,G,mulqidivn(O,gmael(cr,l-1,1),gmael(cr,l-2,1),gmael(cr,l-1,2)));
    if (Mod2(m)) 
    {
        res = cgetg(3,t_VEC);
        gel(res,1) = gdiv(gel(tmp,1),D);
        gel(res,2) = gdiv(gel(tmp,2),D);
    }
    else
    {
        m = diviiexact(m,gen_2);
        res = cgetg(3,t_VEC);
        av = avma; gel(res,1) = gerepileupto(av,modii(mulii(gel(tmp,1),ginvmod(lift(D),m)),m));
        av = avma; gel(res,2) = gerepileupto(av,modii(mulii(gel(tmp,2),ginvmod(lift(D),m)),m));
    }
    return gerepileupto(ltop,res);
}
