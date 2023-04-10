#include "compact_representations.h"

static inline GEN
mulqi (GEN O, GEN a, GEN b, GEN a_, GEN b_)
{ // Returns [x, y] = (x + y*sqrt(d))/s = (a + b*sqrt(d))/s * (a_ + b_*sqrt(d))/2.
    pari_sp ltop = avma, av;
    GEN res = cgetg(3,t_VEC);
    av = avma; gel(res,1) = gerepileupto(av,gdiv(gadd(gmul(a,a_),gmul(gmul(b,b_),gel(O,1))),gel(O,3)));
    av = avma; gel(res,2) = gerepileupto(av,gdiv(gadd(gmul(a,b_),gmul(b,a_)),gel(O,3)));
    return gerepileupto(ltop, res);
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
        for (j = 1; j < l-i; j++)
        {
            fac = gerepileupto(av2,mulqi(O,gel(fac,1),gel(fac,2),gel(fac,1),gel(fac,2)));
        }
        theta = gerepileupto(ltop,mulqi(O,gel(theta,1),gel(theta,2),gel(fac,1),gel(fac,2)));
    }
    return gerepileupto(ltop,theta);
}
