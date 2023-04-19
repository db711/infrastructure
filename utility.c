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
rqfumodm(GEN bnf, GEN m)
{
    pari_sp ltop = avma, av, av2, av3;
    GEN fu, zk, c, c_, fac, faci = gen_0, D, R, R_, P, N, d, res, res_, tmp, tmp_, tmp__;
    ulong i, j;
    fu = gel(bnf_compactfu(bnf),1);
    fac = nf_get_zkden(bnf_get_nf(bnf));
    P = nf_get_pol(bnf_get_nf(bnf));
    zk = nf_get_zk(bnf_get_nf(bnf));
    c = constant_coeff(gmul(fac,gel(zk,2)));
    c_ = leading_coeff(gmul(fac,gel(zk,2)));
    av3 = avma;
    D = m;
    R = mkintmod(gen_1,m);
    if (Mod4(nf_get_disc(bnf_get_nf(bnf)))) R_ = gen_2;
    else R_ = gen_1;
    for (i = 1; i < lg(gel(fu,1)); i++)
    {
        if (typ(gmael(fu,1,i)) == t_INT) //this is always prime
        {
            if (cmpii(gmael(fu,2,i),gen_0) > 0) //exponent positive
            {
                if (cmpii(gcdii(gmael(fu,1,i),m),gen_1) > 0) //not invertible
                {
                    av = avma; R_ = gerepileupto(av,mulii(R_,powii(gmael(fu,1,i),(gmael(fu,2,i)))));
                }
                else //invertible
                {
                    av = avma; R = gerepileupto(av,gmul(R,gpow(mkintmod(modii(gmael(fu,1,i),m),m),gmael(fu,2,i),0)));
                }
            }
            else //exponent negative
            {
                if (cmpii(gcdii(gmael(fu,1,i),m),gen_1) > 0) //not invertible
                {
                    av = avma; D = gerepileupto(av,mulii(D,powii(gmael(fu,1,i),negi(gmael(fu,2,i)))));
                }
                else //invertible
                {
                    av = avma; R = gerepileupto(av,gmul(R,gpow(mkintmod(ginvmod(gmael(fu,1,i),m),m),negi(gmael(fu,2,i)),0)));
                }
            }
        }
        else
        {
            faci = addii(faci,gmael(fu,2,i));
            if (cmpii(gmael(fu,2,i),gen_0) < 0)
            {
                av = avma; N = gerepileupto(av,quadnorm(mkquad(P,addii(mulii(fac,gel(gmael(fu,1,i),1)),mulii(c,gel(gmael(fu,1,i),2))),mulii(c_,gel(gmael(fu,1,i),2)))));
                av2 = avma;
                d = gcdii(N,m);
                while (cmpii(d,gen_1) > 0)
                {
                    av = avma; D = gerepileupto(av,mulii(D,powii(d,negi(gmael(fu,2,i)))));
                    N = diviiexact(N,d);
                    gerepileall(av2,2,&D,&N);
                    d = gcdii(N,m);
                }
                av = avma; R = gerepileupto(av,gmul(R,gpow(mkintmod(ginvmod(N,m),m),negi(gmael(fu,2,i)),0)));
            }
        }
    }
    if(cmpii(faci,gen_0) > 0) D = mulii(D,powii(fac,faci));
    else R_ = mulii(R_,powii(fac,negi(faci)));
    d = gcdii(D,R_);
    R_ = diviiexact(R_,d);
    D = diviiexact(D,d);
    gerepileall(av3,3,&D,&R,&R_);
    res_ = mkquad(P,gen_1,gen_0);
    for (i = 1; i < lg(gel(fu,1)); i++)
    {
        if (typ(gmael(fu,1,i)) == t_COL) //integers already taken care of
        {
            if (cmpii(gmael(fu,2,i),gen_0) > 0) tmp = mkquad(P,addii(mulii(fac,gel(gmael(fu,1,i),1)),mulii(c,gel(gmael(fu,1,i),2))),mulii(c_,gel(gmael(fu,1,i),2)));
            else tmp = mkquad(P,addii(mulii(fac,gel(gmael(fu,1,i),1)),mulii(c,gel(gmael(fu,1,i),2))),negi(mulii(c_,gel(gmael(fu,1,i),2))));
            av2 = avma; 
            tmp_ = shallowcopy(tmp);
            for (j = 1; abscmpui(j,gmael(fu,2,i)) < 0; j++)
            {
                tmp__ = gmul(tmp,tmp_);
                d = gcdii(D,gcdii(gel(tmp__,2),gel(tmp__,3)));
                D = diviiexact(D,d);
                tmp_ = mkquad(P,modii(diviiexact(gel(tmp__,2),d),D),modii(diviiexact(gel(tmp__,3),d),D));
                gerepileall(av2,3,&tmp_,&D,&res_);
            }
            tmp__ = gmul(res_,tmp_);
            d = gcdii(D,gcdii(gel(tmp__,2),gel(tmp__,3)));
            D = diviiexact(D,d);
            res_ = mkquad(P,modii(diviiexact(gel(tmp__,2),d),D),modii(diviiexact(gel(tmp__,3),d),D));
        }
    }
    R_ = mkintmod(modii(lift(R_),D),D);
    res = cgetg(3,t_COL);
    av = avma; gel(res,1) = gerepileupto(av,gmul(gmul(gel(res_,2),R),R_));
    av = avma; gel(res,2) = gerepileupto(av,gmul(gmul(gel(res_,3),R),R_));
    return gerepileupto(ltop,res);
}

GEN
rqfumodm_fixed(GEN bnf, GEN m)
{
    pari_sp ltop = avma;
    GEN res = rqfumodm(bnf,m);
    while (cmpii(gel(gel(res,1),1),m)) res = gerepileupto(ltop,rqfumodm(bnf,mulii(m,diviiexact(m,gel(gel(res,1),1)))));
    return res;
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