#include "twin_smooths.h"

GEN 
twin_smooth_d(GEN S, GEN d, ulong m, GEN P)
{
    if (typ(d) != t_INT) pari_err_TYPE("twin_smooth_d",d);
    if (typ(P) != t_INT) pari_err_TYPE("twin_smooth_d",P);
    if (typ(S) != t_VECSMALL) pari_err_TYPE("twin_smooth_d",S);
    GEN O, b, e, e_, I, R, ret = NULL, r, tmp, smooth_part = gen_1, ret2;
    long i, j;
    pari_sp ltop = avma, av;
    O = rqoinit(d);
    b = pci(O);
    d = mulsi(4,d);
    I = const_vecsmall(m+1,1);
    if(cmprr(glog(d,DEFAULTPREC),dbltor(27.631021115928548208215897456212370491213217863545275712399934811)) < 0) // if d >= 10^12, using compact representations is worth it
    {
        e = quadunit0(d,-1); 
        if (gsigne(gnorm(e)) == -1) e = gsqr(e);
        e_ = e;
        for (i = 1; i < lg(I); i++)
        {
            if (!I[i]) continue;
            if (Z_issmooth(gel(e,3),S[lg(S)-1]))
            {
                if (ret == NULL) ret = mkvec(diviiexact(subii(gel(e,2),gen_1),gen_2));
                else ret = vec_append(ret,diviiexact(subii(gel(e,2),gen_1),gen_2));
            }
            else for (j = i; j < lg(I); j += i) I[j] = 0;
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
                if (!I[i]) continue;
                if (Z_issmooth(gel(e,3),S[lg(S)-1]))
                {
                    if (ret == NULL) ret = mkvec(diviiexact(subii(gel(e,2),gen_1),gen_2));
                    else ret = vec_append(ret,diviiexact(subii(gel(e,2),gen_1),gen_2));
                }
                else for (j = i; j < lg(I); j += i) I[j] = 0;
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
            if (S[lg(S)-1] < 2000) smooth_part = crsmoothpart(O,gel(cr(O,b,r,ghalf,P),2),S,2);
            else
            {
                if(cmprr(glog(d,DEFAULTPREC),dbltor(39.143946580898776628305854729634191529218725306689140592566574316)) < 0) smooth_part = crsmoothpart_alt(O,gel(cr(O,b,r,ghalf,gen_0),2),S,2); // d < 10^17
                else
                {
                    if (S[lg(S)-1] < 5200) smooth_part = crsmoothpart(O,gel(cr(O,b,r,ghalf,P),2),S,2);
                    else smooth_part = crsmoothpart2(O,b,r,ghalf,S,2);
                }
            }
            if (gcmp(absr(subrr(glog(smooth_part,DEFAULTPREC),subrr(R,tmp))),ghalf) <= 0)
            {
                set_avma(av);
                e = quadunit0(d,-1);
                if (cmpii(gnorm(e),gen_1)) e = gsqr(e);
                ret = mkvec(diviiexact(subii(gel(e,2),gen_1),gen_2));
                for (i = 2; i < lg(I); i++)
                {
                    if (!I[i]) continue;
                    av = avma;
                    if (S[lg(S)-1] < 2000) smooth_part = crsmoothpart(O,gel(cr(O,b,r,ghalf,P),2),S,2);
                    else
                    {
                        if(cmprr(glog(d,DEFAULTPREC),dbltor(39.143946580898776628305854729634191529218725306689140592566574316)) < 0) smooth_part = crsmoothpart_alt(O,gel(cr(O,b,r,ghalf,gen_0),2),S,2); // d < 10^17
                        else
                        {
                            if (S[lg(S)-1] < 5200) smooth_part = crsmoothpart(O,gel(cr(O,b,r,ghalf,P),2),S,2);
                            else smooth_part = crsmoothpart2(O,b,r,ghalf,S,2);
                        }
                    }
                    if (gcmp(absr(subrr(glog(smooth_part,DEFAULTPREC),subrr(mulsr(i,R),tmp))),ghalf) <= 0)
                    {
                        set_avma(av);
                        ret = vec_append(ret,diviiexact(subii(gel(gpow(e,stoi(i),DEFAULTPREC),2),gen_1),gen_2));
                    }
                    else for (j = i; j < lg(I); j += i) I[j] = 0;
                }
            }
        }
    }
    if (ret != NULL)
    {
        ret2 = cgetg(lg(ret),t_VEC);
        for (i = 1; i < lg(ret); i++) gel(ret2,i) = gcopy(gel(ret,i));
        return gerepileupto(ltop,ret2);
    }
    else 
    {
        set_avma(ltop);
        return NULL;
    }
}

GEN
twin_smooth(ulong B)
{
    GEN S, s, D, d, ret, res = NULL, res2, P = gen_1;
    pari_sp ltop = avma, av, av2;
    long i, m, k;
    forsubset_t T;
    S = primes_interval_zv(2,B);
    m = ((S[lg(S)-1])-1)/2;
    forallsubset_init(&T,lg(S)-1);
    forsubset_next(&T); forsubset_next(&T);
    for (k = 1; k < lg(S); k++) P = mulis(P,S[k]);
    res = twin_smooth_d(S,gen_2,m,P);
    av2 = avma;
    while ((s = forsubset_next(&T)) != NULL)
    {
        D = shallowextract(S,s);
        d = gen_2;
        for (i = 1; i < lg(D); i++)
        {
            av = avma; d = gerepileupto(av,mulis(d,(ulong)gel(D,i)));
        }
        ret = twin_smooth_d(S,d,m,P);
        if (ret != NULL) res = gerepileupto(av2,concat(res,ret));
    }
    res2 = cgetg(lg(res),t_VEC);
    for (i = 1; i < lg(res); i++) gel(res2,i) = gcopy(gel(res,i));
    return gerepileupto(ltop,res2);
}