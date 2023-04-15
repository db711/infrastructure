#include "compact_representations.h"

static inline GEN
mulqifrac (GEN O, GEN qi1, GEN qi2)
{ // Returns [x, y] = (x + y*sqrt(d))/s = (a + b*sqrt(d))/s * (a_ + b_*sqrt(d))/s, where qi1 = [a, b] and qi2 = [a_, b_].
    pari_sp ltop = avma, av;
    GEN res = cgetg(3,t_VEC);
    av = avma; gel(res,1) = gerepileupto(av,gdiv(gadd(gmul(gel(qi1,1),gel(qi2,1)),gmul(gmul(gel(qi1,2),gel(qi2,2)),gel(O,1))),gel(O,3)));
    av = avma; gel(res,2) = gerepileupto(av,gdiv(gadd(gmul(gel(qi1,1),gel(qi2,2)),gmul(gel(qi1,2),gel(qi2,1))),gel(O,3)));
    return gerepileupto(ltop,res);
}

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
        for (j = 1; j < l-i; j++) fac = gerepileupto(av2,mulqifrac(O,fac,fac));
        theta = gerepileupto(ltop,mulqifrac(O,theta,fac));
    }
    return gerepileupto(ltop,theta);
}
