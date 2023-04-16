#include "utility.h"

ulong
sigbits (GEN x)
{
    if (typ(x) != t_INT) pari_err_TYPE("sigbits",x);
    ulong i = 0, x_;
    x_ = *int_MSW(x);
    while (x_ >>= 1) i++;
    return (lgefint(x)-3)*BITS_IN_LONG+i+1;
}

GEN
mulqidivn (GEN O, GEN qi1, GEN qi2, GEN N)
{ // Returns [x, y] = (x + y*sqrt(d))/(s*N) = ((a + b*sqrt(d))/s * (a_ + b_*sqrt(d))/s)/N, where qi1 = [a, b] and qi2 = [a_, b_].
    pari_sp ltop = avma, av;
    GEN res = cgetg(3,t_VEC);
    av = avma; gel(res,1) = gerepileupto(av,diviiexact(addii(mulii(gel(qi1,1),gel(qi2,1)),mulii(mulii(gel(qi1,2),gel(qi2,2)),gel(O,1))),mulii(gel(O,3),N)));
    av = avma; gel(res,2) = gerepileupto(av,diviiexact(addii(mulii(gel(qi1,1),gel(qi2,2)),mulii(gel(qi1,2),gel(qi2,1))),mulii(gel(O,3),N)));
    return gerepileupto(ltop,res);
}

GEN
mulqig (GEN O, GEN qi1, GEN qi2)
{ // Returns [x, y] = (x + y*sqrt(d))/s = (a + b*sqrt(d))/s * (a_ + b_*sqrt(d))/s, where qi1 = [a, b] and qi2 = [a_, b_].
    pari_sp ltop = avma, av;
    GEN res = cgetg(3,t_VEC);
    av = avma; gel(res,1) = gerepileupto(av,gdiv(gadd(gmul(gel(qi1,1),gel(qi2,1)),gmul(gmul(gel(qi1,2),gel(qi2,2)),gel(O,1))),gel(O,3)));
    av = avma; gel(res,2) = gerepileupto(av,gdiv(gadd(gmul(gel(qi1,1),gel(qi2,2)),gmul(gel(qi1,2),gel(qi2,1))),gel(O,3)));
    return gerepileupto(ltop,res);
}

GEN 
invqi(GEN qi)
{ // Returns [a, -b], where qi = [a, b].
    pari_sp ltop = avma;
    GEN res;
    res = cgetg(3,t_VEC);
    gel(res,1) = gcopy(gel(qi,1));
    gel(res,2) = negi(gel(qi,2));
    return gerepileupto(ltop,res);
}