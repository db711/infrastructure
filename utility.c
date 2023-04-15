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
split(GEN h, GEN m)
{
    if (typ(h) != t_INT) pari_err_TYPE("split",h);
    if (typ(m) != t_INT) pari_err_TYPE("split",m);
    if (cmpii(m,gen_1) < 0) pari_err_DOMAIN("split",itostr(m),"<",gen_1,m);
    GEN g, r, s, res;
    pari_sp ltop = avma;
    g = gcdii(m,h);
    r = diviiexact(m,g);
    s = g;
    while (cmpii(g,gen_1) > 0)
    {
        g = gcdii(r,g);
        r = diviiexact(r,g);
        s = mulii(g,s);
        gerepileall(ltop,3,&g,&r,&s);
    }
    res = cgetg(3,t_VEC);
    gel(res,1) = gcopy(r);
    gel(res,2) = gcopy(s);
    return gerepileupto(ltop,res);
}