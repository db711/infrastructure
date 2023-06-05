#include "twin_smooths.h"

GEN 
twin_smooth_d(ulong B, GEN d, ulong m)
{
    GEN O, b, e, e_, I, R, ret = NULL, r, tmp;
    long i, j, ri = 0;
    pari_sp ltop = avma, av;
    O = rqoinit(d);
    b = pci(O);
    d = mulsi(4,d);
    I = const_vecsmall(m+1,1);
    if(cmprr(glog(d,DEFAULTPREC),dbltor(27.631021115928548208215897456212370491213217863545275712399934811)) < 0)
    {
        e = quadunit0(d,-1); 
        if (gsigne(gnorm(e)) == -1) e = gsqr(e);
        e_ = e;
        for (i = 1; i < lg(I); i++)
        {
            if (!gel(I,i)) continue;
            if (Z_issmooth(gel(e,3),B))
            {
                if (!ri)
                {
                    ret = mkvec(diviiexact(subii(gel(e,2),gen_1),gen_2));
                    ri = 1;
                }
                else ret = vec_append(ret,diviiexact(subii(gel(e,2),gen_1),gen_2));
            }
            else for (j = i; j < lg(I); j += j) gel(I,j) = 0;
            e = gmul(e,e_);
        }
    }
    else
    {
        R = gel(quadclassunit0(d,0,NULL,DEFAULTPREC),4);
        if (cmpri(R,stoi(25)) < 0) // in this case the continued fraction method is still faster
        { //code duplication...
            e = quadunit0(d,-1);
            if (gsigne(gnorm(e)) == -1) e = gsqr(e);
            e_ = e;
            for (i = 1; i < lg(I); i++)
            {
                if (!gel(I,i)) continue;
                if (Z_issmooth(gel(e,3),B))
                {
                    if (!ri)
                    {
                        ret = mkvec(diviiexact(subii(gel(e,2),gen_1),gen_2));
                        ri = 1;
                    }
                    else ret = vec_append(ret,diviiexact(subii(gel(e,2),gen_1),gen_2));
                }
                else for (j = i; j < lg(I); j += i) gel(I,j) = 0;
                e = gmul(e,e_);
            }
        }
        else
        {
            av = avma; tmp = gerepileupto(av,addrr(glog(gen_2,DEFAULTPREC),glog(gsqrt(gel(O,1),DEFAULTPREC),DEFAULTPREC)));
            av = avma; r = gerepileupto(av,roundr(gdiv(R,mplog2(DEFAULTPREC))));
            if (crnorm_sign(O,gel(cr(O,b,r,ghalf,gen_0),2)) == -1) 
            {
                R = mulir(gen_2,R);
                r = mulii(gen_2,r);
            }
            av = avma;
            if (gcmp(absr(subrr(glog(crsmoothpart2(O,b,r,ghalf,B,2),DEFAULTPREC),subrr(R,tmp))),ghalf) <= 0)
            {
                set_avma(av);
                e = quadunit0(d,-1);
                if (cmpii(gnorm(e),gen_1)) e = gsqr(e);
                ret = mkvec(diviiexact(subii(gel(e,2),gen_1),gen_2));
                for (i = 2; i < lg(I); i++)
                {
                    if (!gel(I,i)) continue;
                    av = avma;
                    if (gcmp(absr(subrr(glog(crsmoothpart2(O,b,mulsi(i,r),ghalf,B,2),DEFAULTPREC),subrr(mulsr(i,R),tmp))),ghalf) <= 0)
                    {
                        set_avma(av);
                        ret = vec_append(ret,diviiexact(subii(gel(gpow(e,stoi(i),DEFAULTPREC),2),gen_1),gen_2));
                    }
                    else for (j = i; j < lg(I); j += i) gel(I,j) = 0;
                }
            }
        }
    }
    if (ret != NULL) return gerepileupto(ltop,gcopy(ret));
    else 
    {
        set_avma(ltop);
        return NULL;
    }
}

GEN
twin_smooth(ulong B)
{
    GEN S, s, D, d, ret, res = mkvec(gen_1);
    pari_sp ltop = avma, av, av2;
    long i, m;
    forsubset_t T;
    S = primes_interval_zv(3,B);
    m = gel(S,lg(S)-1);
    forallsubset_init(&T,lg(S)-1);
    av2 = avma;
    s = forsubset_next(&T);
    while ((s = forsubset_next(&T)) != NULL)
    {
        D = shallowextract(S,s);
        d = gen_2;
        for (i = 1; i < lg(D); i++)
        {
            av = avma; d = gerepileupto(av,mulis(d,gel(D,i)));
        }
        ret = twin_smooth_d(B,d,m+1);
        if (ret != NULL) res = gerepileupto(av2,gcopy(concat(res,ret)));
    }
    return gerepileupto(ltop,gcopy(res));
}