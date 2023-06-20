#include "twin_smooths.h"
#include "utility.h"

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
    d = gel(O,4);
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
            if(equalii(gel(O,3),gen_2)) smooth_part = diviiexact(smooth_part,gen_2);
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
        ret2 = cgetg(1,t_VEC);
        return ret2;
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

GEN
regulator_range(GEN O, ulong A, ulong B)
{
    if (A > B) pari_err_DOMAIN("regulator_range",itostr(stoi(A)),">",stoi(B),stoi(A));
    GEN sqrtd_, fprep, b, w, w_, p, s, tmp, tmp2, twossqrtd, Q_1, P_1, M, P_m1, P_0, Q_m1, T_m1, Q_0, T_0, T_1, q, c, e, t, g, h;
    pari_sp ltop = avma, av, av2;
    b = pci(O);
    sqrtd_ = sqrti(gel(O,1));
    w = stoi(A);
    w_ = stoi(B);
    av = avma; p = gerepileupto(av,addsi(4.48543+sigbits(w),gmax_shallow(stoi(4),stoi(sigbits(stoi(sigbits(w)+1))+1)))); //generous upper bound
    fprep = ax(O,w,p);
    s = gmax_shallow(gen_0,addis(gmael(fprep,1,2),5-sigbits(gmael4(fprep,2,1,2,1))));
    av = avma; if (cmpii(gmael4(fprep,2,1,2,1),powii(gen_2,subii(addis(gmael(fprep,1,2),4),s))) <= 0) s = gerepileupto(av,addii(s,gen_1));
    else set_avma(av);
    tmp = powii(gen_2,s);
    twossqrtd = sqrti(mulii(sqri(tmp),gel(O,1)));
    Q_1 = gmael4(fprep,2,1,2,1);
    P_1 = gmael4(fprep,2,1,2,2);
    av = avma; M = gerepileupto(av,gceil(gmul(powii(gen_2,addii(subii(addii(gmael(fprep,1,2),s),gmael(fprep,2,3)),w_)),gdiv(gmael4(fprep,2,1,2,1),gmael(fprep,2,2)))));
    av2 = avma;
    P_m1 = P_0 = Q_m1 = T_m1 = gen_0; //dummy values
    av = avma; Q_0 = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(gmael4(fprep,2,1,2,2))),gmael4(fprep,2,1,2,1)));
    av = avma; T_0 = gerepileupto(av,addii(mulii(negi(tmp),P_1),twossqrtd)); 
    T_1 = mulii(tmp,gmael4(fprep,2,1,2,1));
    while (cmpii(T_1,M) <= 0)
    {
        if (!cmpii(Q_1,gmael2(b,2,1)) && !cmpii(P_1,gmael2(b,2,2))) break;
        av = avma; q = gerepileupto(av,divii(addii(P_1,sqrtd_),Q_1));
        P_m1 = P_0; P_0 = P_1;
        av = avma; P_1 = gerepileupto(av,subii(mulii(q,Q_1),P_1));
        Q_m1 = Q_0; Q_0 = Q_1;
        av = avma; Q_1 = gerepileupto(av,subii(Q_m1,mulii(q,subii(P_1,P_0))));
        T_m1 = T_0; T_0 = T_1;
        av = avma; T_1 = gerepileupto(av,addii(mulii(q,T_1),T_m1));
        gerepileall(av2,9,&P_m1,&P_0,&P_1,&Q_m1,&Q_0,&Q_1,&T_m1,&T_0,&T_1);
    }
    if (!cmpii(Q_1,gmael2(b,2,1)) && !cmpii(P_1,gmael2(b,2,2)))
    {
        return gen_1;
        c = rqiinit(gen_1,Q_1,P_1);
        av = avma; e = gerepileupto(av,gceil(gmul(powii(gen_2,addis(subii(gmael(fprep,1,2),s),3)),gdiv(T_0,gmael4(fprep,2,1,2,1)))));
        av = avma;
        if (cmpii(mulii(gmael(fprep,2,2),e),powii(gen_2,addis(addii(subii(mulii(gen_2,gmael(fprep,1,2)),gmael(fprep,2,3)),w),3))) <= 0)
        {
            set_avma(av);
            av = avma; e = gerepileupto(av,gceil(gmul(powii(gen_2,addis(subii(gmael(fprep,1,2),s),3)),gdiv(T_1,gmael4(fprep,2,1,2,1)))));
        }
        else set_avma(av);
        av = avma; t = gerepileupto(av,subsi(sigbits(e)+sigbits(gmael(fprep,2,2))-6,mulii(gen_2,gmael(fprep,1,2))));
        tmp = mulii(e,gmael(fprep,2,2)); tmp2 = powii(gen_2,addii(addis(t,4),mulii(gen_2,gmael(fprep,1,2))));
        if (gcmp(tmp2,tmp) < 0)
        {
            if (cmpii(mulii(tmp2,gen_2),tmp) < 0) t = gerepileupto(av,addii(t,gen_2));
            else t = gerepileupto(av,addii(t,gen_1));
        }
        else set_avma(av);
        av = avma; g = gerepileupto(av,gceil(gdiv(mulii(e,gmael(fprep,2,2)),powii(gen_2,addis(addii(gmael(fprep,1,2),t),3)))));
        h = addii(gmael(fprep,2,3),t);
        if (gmael(fprep,1,1) == NULL) return gerepileupto(ltop,fprepinit(NULL,gmael(fprep,1,2),c,g,h));
        else return gerepileupto(ltop,fprepinit(gadd(gmael(fprep,1,1),gdiv(powii(stoi(3),gen_2),powis(gen_2,3))),gmael(fprep,1,2),c,g,h));
    }
    else
    {
        set_avma(ltop);
        return NULL;
    }
}

GEN
twin_smooth_range_d_small(ulong B, GEN d, ulong bot, ulong top, ulong m)
{
    GEN O, e, e_, I, ret = NULL, ret2;
    ulong i;
    pari_sp ltop = avma;
    O = rqoinit(d);
    if (regulator_range(O,bot,top))
    {
        d = gel(O,4);
        I = const_vecsmall(m,1);
        e = quadunit0(d,-1); 
        if (gsigne(gnorm(e)) == -1) e = gsqr(e);
        e_ = e;
        for (i = 1; i < lg(I); i++)
        {
            if (Z_issmooth(gel(e,3),B))
            {
                if (ret == NULL) ret = mkvec(diviiexact(subii(gel(e,2),gen_1),gen_2));
                else ret = vec_append(ret,diviiexact(subii(gel(e,2),gen_1),gen_2));
            }
            e = gmul(e,e_);
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
        ret2 = cgetg(1,t_VEC);
        return ret2;
    }
}