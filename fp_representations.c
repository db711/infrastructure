#include "fp_representations.h"
#include "utility.h"

static inline GEN 
split(GEN h, GEN m)
{ // Split integer m into [r, s] with m = r*s such that gcd(r, s) = 1 and any prime dividing s also divides h.
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

GEN
fprepinit (GEN f, GEN p, GEN b, GEN d, GEN k)
{
    GEN res, res1, res2;
    pari_sp ltop = avma;
    if (f != NULL && typ(f) != t_REAL) pari_err_TYPE("fprepinit",f); //could be fixed
    if (typ(p) != t_INT) pari_err_TYPE("fprepinit",p);
    if (typ(d) != t_INT) pari_err_TYPE("fprepinit",d);
    if (typ(k) != t_INT) pari_err_TYPE("fprepinit",k);
    if (f != NULL && gcmp(f,gen_1) < 0) pari_err_DOMAIN("fprepinit",GENtostr(f),"<",gen_1,f);
    if (cmpii(p,gen_0) < 0) pari_err_DOMAIN("fprepinit",itostr(p),"<",gen_0,p);
    if (cmpii(powii(gen_2,p),d) >= 0 || cmpii(powii(gen_2,addii(p,gen_1)),d) < 0) pari_err_DOMAIN("fprepinit",itostr(d),"",NULL,d);
    set_avma(ltop);
    res = cgetg(3,t_VEC);
    res1 = cgetg(3,t_VEC);
    gel(res1,1) = (f == NULL) ? NULL : gcopy(f);
    gel(res1,2) = gcopy(p); //should really every (f, p) representation carry this?
    gel(res,1) = res1;
    res2 = cgetg(4,t_VEC);
    gel(res2,1) = gcopy(b); 
    gel(res2,2) = gcopy(d);
    gel(res2,3) = gcopy(k);
    gel(res,2) = res2;
    return res;
}

GEN
fpremove (GEN fprep, GEN T, GEN C, GEN s)
{
    if (typ(C) != t_INT) pari_err_TYPE("fpremove",C);
    if (typ(T) != t_INT) pari_err_TYPE("fpremove",T);
    if (typ(s) != t_INT) pari_err_TYPE("fpremove",s);
    if (equalii(C,gen_0)) pari_err_DOMAIN("fpremove","C","==",gen_0,C);
    pari_sp ltop = avma, av;
    GEN gen_3 = addii(gen_2,gen_1), e, d, k;
    ulong t;
    av = avma; e = gerepileupto(av,ground(gmul(powii(gen_2,subii(addii(gmael(fprep,1,2),gen_3),s)),gabs(gdiv(T,C),DEFAULTPREC)))); //precision should be ignored
    t = sigbits(e)-sigbits(gmael(fprep,2,2))-3;
    av = avma; if (cmpii(e,mulii(powis(gen_2,t+3),gmael(fprep,2,2))) >= 0) t++; set_avma(av);
    av = avma; d = gerepileupto(av,gceil(gmul(powii(gen_2,addis(addii(gmael(fprep,1,2),gen_3),t)),gdiv(gmael(fprep,2,2),e))));
    av = avma; k = gerepileupto(av,subis(gmael(fprep,2,3),t));
    if (gmael(fprep,1,1) == NULL) return gerepileupto(ltop,fprepinit(NULL,gmael(fprep,1,2),gmael(fprep,2,1),d,k));
    else return gerepileupto(ltop,fprepinit(gadd(gmael(fprep,1,1),gdiv(powii(gen_3,gen_2),powii(gen_2,gen_3))),gmael(fprep,1,2),gmael(fprep,2,1),d,k));
}

GEN
numult (GEN O, GEN fprep1, GEN fprep2, long flag)
{
    if (flag && cmpii(gmael(fprep1,1,2),gmael(fprep2,1,2))) pari_err_DOMAIN("numult",itostr(gmael(fprep2,1,2)),"!=",gmael(fprep1,1,2),gmael(fprep2,1,2));
    GEN b, e, h, s, T, res, res2;
    pari_sp ltop = avma, av;
    b = inucomp(O,gmael(fprep1,2,1),gmael(fprep2,2,1),flag);
    av = avma;
    if (cmpii(mulii(gmael(fprep1,2,2),gmael(fprep2,2,2)),powii(gen_2,addii(mulii(gen_2,gmael(fprep1,1,2)),gen_1))) <= 0)
    {
        e = gerepileupto(av,divii(mulii(gmael(fprep1,2,2),gmael(fprep2,2,2)),powii(gen_2,gmael(fprep1,1,2))));
        h = addii(gmael(fprep1,2,3),gmael(fprep2,2,3));
    }
    else
    {
        e = gerepileupto(av,gceil(gdiv(mulii(gmael(fprep1,2,2),gmael(fprep2,2,2)),powii(gen_2,addii(gmael(fprep1,1,2),gen_1)))));
        av = avma; h = gerepileupto(av,addii(addii(gmael(fprep1,2,3),gmael(fprep2,2,3)),gen_1));
    }
    if (cmpii(gmael(b,2,2),gen_0) > 0)
    {
        s = addis(gmael(fprep1,1,2),sigbits(gmael(b,2,2))-sigbits(gmael3(b,1,2,1))+4);
        av = avma; if (cmpii(mulii(powii(gen_2,subii(addis(gmael(fprep1,1,2),4),s)),gmael(b,2,2)),gmael3(b,1,2,1)) >= 0) s = gerepileupto(av,addii(s, gen_1));
        else set_avma(av);
    }
    else s = gen_0;
    av = avma; T = gerepileupto(av,addii(mulii(powii(gen_2,s),gmael(b,2,1)),mulii(gmael(b,2,2),gfloor(gmul(powii(gen_2,s),gsqrt(gel(O,1),DEFAULTPREC)))))); //should throw error if precision is too low
    res = cgetg(3,t_VEC);
    if (gmael(fprep1,1,1) == NULL || gmael(fprep2,1,1) == NULL) 
    {
        av = avma; gel(res,1) = gerepileupto(av,fpremove(fprepinit(NULL,gmael(fprep1,1,2),gel(b,1),e,h),T,gmael(b,2,3),s));
    }
    else
    {
        av = avma; gel(res,1) = gerepileupto(av,fpremove(fprepinit(addri(addrr(addrr(gmael(fprep1,1,1),gmael(fprep2,1,1)),gmul(powii(gen_2,negi(gmael(fprep1,1,2))),mulrr(gmael(fprep1,1,1),gmael(fprep2,1,1)))),gen_1),gmael(fprep1,1,2),gel(b,1),e,h),T,gmael(b,2,3),s));
    }
    res2 = cgetg(3,t_VEC);
    gel(res2,1) = gcopy(gmael(b,2,1));
    gel(res2,2) = negi(gmael(b,2,2));
    gel(res,2) = res2;
    return gerepileupto(ltop,res);
}

GEN 
wnear(GEN O, GEN fprep, GEN w)
{
    if (typ(w) != t_INT) pari_err_TYPE("wnear",w);
    if (cmpii(w,gen_1) < 0) pari_err_DOMAIN("wnear","w","<",gen_1,w);
    pari_sp ltop = avma, av, av2;
    if (cmpii(addii(gmael4(fprep,2,1,2,2),sqrti(gel(O,1))),gmael4(fprep,2,1,2,1)) < 0) pari_err_DOMAIN("wnear",itostr(addii(gmael4(fprep,2,1,2,2),sqrti(gel(O,1)))),"<",gmael4(fprep,2,1,2,1),gmael(fprep,2,1));
    if (cmpii(gen_0,subii(sqrti(gel(O,1)),gmael4(fprep,2,1,2,2))) > 0 || cmpii(subii(sqrti(gel(O,1)),gmael4(fprep,2,1,2,2)),gmael4(fprep,2,1,2,1)) > 0) pari_err_DOMAIN("wnear",itostr(subii(sqrti(gel(O,1)),gmael4(fprep,2,1,2,2))),"",NULL,gmael(fprep,2,1));
    set_avma(ltop);
    GEN s, Q_0, Q_1, Q_m1, Q_m2, Q_m3, P_0, P_1, P_m1, M, T_0, T_1, T_m1, T_m2, T_m3, q, twossqrtd, sqrtd_, tmp, tmp2, gen_3 = addii(gen_2,gen_1), e, e_, c, t, g, h;
    sqrtd_ = sqrti(gel(O,1));
    if (cmpii(gmael(fprep,2,3),w) < 0)
    {
        s = addis(gmael(fprep,1,2),5-sigbits(gmael4(fprep,2,1,2,1)));
        av = avma; if (cmpii(gmael4(fprep,2,1,2,1),powii(gen_2,subii(addis(gmael(fprep,1,2),4),s))) <= 0) s = gerepileupto(av,addii(s,gen_1));
        else set_avma(av);
        tmp = powii(gen_2,s);
        twossqrtd = sqrti(mulii(sqri(tmp),gel(O,1)));
        Q_1 = gmael4(fprep,2,1,2,1);
        P_1 = gmael4(fprep,2,1,2,2);
        av = avma; M = gerepileupto(av,gceil(gmul(powii(gen_2,addii(subii(addii(gmael(fprep,1,2),s),gmael(fprep,2,3)),w)),gdiv(gmael4(fprep,2,1,2,1),gmael(fprep,2,2)))));
        av2 = avma;
        P_m1 = P_0 = Q_m1 = T_m1 = gen_0; //dummy values
        av = avma; Q_0 = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(gmael4(fprep,2,1,2,2))),gmael4(fprep,2,1,2,1)));
        av = avma; T_0 = gerepileupto(av,addii(mulii(negi(tmp),P_1),twossqrtd)); 
        T_1 = mulii(tmp,gmael4(fprep,2,1,2,1));
        while (cmpii(T_1,M) <= 0)
        {
            av = avma; q = gerepileupto(av,divii(addii(P_1,sqrtd_),Q_1));
            P_m1 = P_0; P_0 = P_1;
            av = avma; P_1 = gerepileupto(av,subii(mulii(q,Q_1),P_1));
            Q_m1 = Q_0; Q_0 = Q_1;
            av = avma; Q_1 = gerepileupto(av,subii(Q_m1,mulii(q,subii(P_1,P_0))));
            T_m1 = T_0; T_0 = T_1;
            av = avma; T_1 = gerepileupto(av,addii(mulii(q,T_1),T_m1));
            gerepileall(av2,9,&P_m1,&P_0,&P_1,&Q_m1,&Q_0,&Q_1,&T_m1,&T_0,&T_1);
        }
        av = avma; e = gerepileupto(av,gceil(gmul(powii(gen_2,addii(subii(gmael(fprep,1,2),s),gen_3)),gdiv(T_0,gmael4(fprep,2,1,2,1)))));
        av = avma;
        if (cmpii(mulii(gmael(fprep,2,2),e),powii(gen_2,addii(addii(subii(mulii(gen_2,gmael(fprep,1,2)),gmael(fprep,2,3)),w),gen_3))) <= 0)
        {
            set_avma(av); c = rqiinit(gen_1,Q_0,P_0);
            //also set next element
        }
        else
        {
            set_avma(av); c = rqiinit(gen_1,Q_m1,P_m1);
            av = avma; e = gerepileupto(av,gceil(gmul(powii(gen_2,addii(subii(gmael(fprep,1,2),s),gen_3)),gdiv(T_m1,gmael4(fprep,2,1,2,1)))));
            //also set next element
        }
        av = avma; t = gerepileupto(av,subsi(sigbits(e)+sigbits(gmael(fprep,2,2))-6,mulii(gen_2,gmael(fprep,1,2))));
        tmp = mulii(e,gmael(fprep,2,2)); tmp2 = powii(gen_2,addii(addis(t,4),mulii(gen_2,gmael(fprep,1,2))));
        if (gcmp(tmp2,tmp) < 0)
        {
            if (cmpii(mulii(tmp2,gen_2),tmp) < 0) t = gerepileupto(av,addii(t,gen_2));
            else t = gerepileupto(av,addii(t,gen_1));
        }
        else set_avma(av);
        //also set for next element and initialize that   
        av = avma; g = gerepileupto(av,gceil(gdiv(mulii(e,gmael(fprep,2,2)),powii(gen_2,addii(addii(gmael(fprep,1,2),t),gen_3)))));
        h = addii(gmael(fprep,2,3),t);
    }
    else
    {
        av = avma; s = gerepileupto(av,addii(gmael(fprep,1,2),addii(gen_2,gen_2)));
        tmp = powii(gen_2,s);
        twossqrtd = sqrti(mulii(sqri(tmp),gel(O,1)));
        Q_0 = gmael4(fprep,2,1,2,1);
        P_1 = gmael4(fprep,2,1,2,2); 
        av = avma; M = gerepileupto(av,mulii(gmael(fprep,2,2),powii(gen_2,addii(subii(gmael(fprep,2,3),w),addii(gen_2,gen_2)))));
        av2 = avma;
        Q_m3 = Q_m2 = Q_m1 = T_m3 = T_m2 = T_m1 = gen_0; //dummy values
        av = avma; Q_1 = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(gmael4(fprep,2,1,2,2))),gmael4(fprep,2,1,2,1)));
        av = avma; T_0 = gerepileupto(av,mulii(tmp,Q_0));
        av = avma; T_1 = gerepileupto(av,addii(mulii(tmp,P_1),twossqrtd)); 
        while (cmpii(T_1,mulii(Q_1,M)) < 0)
        {
            av = avma; q = gerepileupto(av,divii(addii(P_1,sqrtd_),Q_1));
            P_0 = P_1;
            av = avma; P_1 = gerepileupto(av,subii(mulii(q,Q_1),P_1));
            Q_m3 = Q_m2; Q_m2 = Q_m1; Q_m1 = Q_0;
            Q_0 = Q_1;
            av = avma; Q_1 = gerepileupto(av,subii(Q_m1,mulii(q,subii(P_1,P_0))));
            T_m3 = T_m2; T_m2 = T_m1; T_m1 = T_0;
            T_0 = T_1;
            av = avma; T_1 = gerepileupto(av,addii(mulii(q,T_1),T_m1));
            gerepileall(av2,12,&P_0,&P_1,&Q_m3,&Q_m2,&Q_m1,&Q_0,&Q_1,&T_m3,&T_m2,&T_m1,&T_0,&T_1);
        }
        av = avma; q = gerepileupto(av,divii(addii(P_1,sqrtd_),Q_1));
        P_0 = P_1;
        av = avma; P_1 = gerepileupto(av,subii(mulii(q,Q_1),P_1));
        av = avma; e = gerepileupto(av,gceil(gdiv(T_1,mulii(gen_2,Q_1))));
        av = avma; e_ = gerepileupto(av,gceil(gdiv(T_0,mulii(gen_2,Q_0))));
        av2 = avma;
        av = avma; tmp = gerepileupto(av,mulii(gmael(fprep,2,2),powii(gen_2,addii(subii(gmael(fprep,2,3),w),gen_3))));
        if (cmpii(e_,tmp) >= 0)
        {
            e = e_;
            av = avma; e_ = gerepileupto(av,gceil(gdiv(T_m2,mulii(gen_2,Q_m2))));
            if (cmpii(e_,tmp) >= 0)
            {
                e = e_;
                av = avma; e_ = gerepileupto(av,gceil(gdiv(T_m3,mulii(gen_2,Q_m3))));
                if (cmpii(e_,tmp) >= 0) pari_err_IMPL("this");
            }
            c = rqiinit(gen_1,Q_0,P_0);
            //also set next element
        }
        else c = rqiinit(gen_1,Q_1,P_1); //also set next element
        av = avma; t = gerepileupto(av,subsi(sigbits(e)+sigbits(gmael(fprep,2,2))-6,mulii(gen_2,gmael(fprep,1,2))));
        tmp = mulii(e,gmael(fprep,2,2)); tmp2 = powii(gen_2,addii(addis(t,4),mulii(gen_2,gmael(fprep,1,2))));
        if (gcmp(tmp2,tmp) < 0)
        {
            if (cmpii(mulii(tmp2,gen_2),tmp) < 0) t = gerepileupto(av,addii(t,gen_2));
            else t = gerepileupto(av,addii(t,gen_1));
        }
        else set_avma(av);
        av = avma; t = gerepileupto(av,stoi(sigbits(e)-sigbits(gmael(fprep,2,2))-3));
        av = avma; if (cmpii(e,mulii(powii(gen_2,t+3),gmael(fprep,2,2))) >= 0) t = gerepileupto(av,addii(t,gen_1));
        else set_avma(av);
        //also set for next element and initialize that
        av = avma; g = gerepileupto(av,gceil(gmul(powii(gen_2,addii(addii(gmael(fprep,1,2),gen_3),t)),gdiv(gmael(fprep,2,2),e))));
        h = subii(gmael(fprep,2,3),t);
    }
    if (gmael(fprep,1,1) == NULL) return gerepileupto(ltop,fprepinit(NULL,gmael(fprep,1,2),c,g,h));
    else return gerepileupto(ltop,fprepinit(gadd(gmael(fprep,1,1),gdiv(powii(gen_3,gen_2),powii(gen_2,gen_3))),gmael(fprep,1,2),c,g,h));
}

GEN 
ewnear(GEN O, GEN fprep, GEN w)
{
    if (typ(w) != t_INT) pari_err_TYPE("ewnear",w);
    if (cmpii(w,gen_1) < 0) pari_err_DOMAIN("ewnear","w","<",gen_1,w);
    if (cmpii(gmael(fprep,2,3),w) >= 0) pari_err_DOMAIN("ewnear","w","<=",gmael(fprep,2,3),w); //TODO: Implement other case
    pari_sp ltop = avma, av, av2;
    if (cmpii(addii(gmael4(fprep,2,1,2,2),sqrti(gel(O,1))),gmael4(fprep,2,1,2,1)) < 0) pari_err_DOMAIN("ewnear",itostr(addii(gmael4(fprep,2,1,2,2),sqrti(gel(O,1)))),"<",gmael4(fprep,2,1,2,1),gmael(fprep,2,1));
    if (cmpii(gen_0,subii(sqrti(gel(O,1)),gmael4(fprep,2,1,2,2))) > 0 || cmpii(subii(sqrti(gel(O,1)),gmael4(fprep,2,1,2,2)),gmael4(fprep,2,1,2,1)) > 0) pari_err_DOMAIN("ewnear",itostr(subii(sqrti(gel(O,1)),gmael4(fprep,2,1,2,2))),"",NULL,gmael(fprep,2,1));
    set_avma(ltop);
    GEN B_0, B_1, B_m1, a, b, s, Q_0, Q_1, Q_m1, P_0, P_1, P_m1, M, T_0, T_1, T_m1, q, twossqrtd, sqrtd_, tmp, tmp2, gen_3 = addii(gen_2,gen_1), e, c, t, g, h, res;
    sqrtd_ = sqrti(gel(O,1));
    B_0 = gen_1;
    B_1 = gen_0;
    s = addis(gmael(fprep,1,2),5-sigbits(gmael4(fprep,2,1,2,1)));
    av = avma; if (cmpii(gmael4(fprep,2,1,2,1),powii(gen_2,subii(addis(gmael(fprep,1,2),4),s))) <= 0) s = gerepileupto(av,addii(s,gen_1));
    else set_avma(av);
    tmp = powii(gen_2,s);
    twossqrtd = sqrti(mulii(sqri(tmp),gel(O,1)));
    Q_1 = gmael4(fprep,2,1,2,1);
    P_1 = gmael4(fprep,2,1,2,2);
    av = avma; M = gerepileupto(av,gceil(gmul(powii(gen_2,addii(subii(addii(gmael(fprep,1,2),s),gmael(fprep,2,3)),w)),gdiv(gmael4(fprep,2,1,2,1),gmael(fprep,2,2)))));
    av2 = avma;
    P_m1 = P_0 = Q_m1 = T_m1 = gen_0; //dummy values
    av = avma; Q_0 = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(gmael4(fprep,2,1,2,2))),gmael4(fprep,2,1,2,1)));
    av = avma; T_0 = gerepileupto(av,addii(mulii(negi(tmp),P_1),twossqrtd)); 
    T_1 = mulii(tmp,gmael4(fprep,2,1,2,1));
    while (cmpii(T_1,M) <= 0)
    {
        av = avma; q = gerepileupto(av,divii(addii(P_1,sqrtd_),Q_1));
        P_m1 = P_0; P_0 = P_1;
        av = avma; P_1 = gerepileupto(av,subii(mulii(q,Q_1),P_1));
        Q_m1 = Q_0; Q_0 = Q_1;
        av = avma; Q_1 = gerepileupto(av,subii(Q_m1,mulii(q,subii(P_1,P_0))));
        T_m1 = T_0; T_0 = T_1;
        av = avma; T_1 = gerepileupto(av,addii(mulii(q,T_1),T_m1));
        B_m1 = B_0; B_0 = B_1;
        av = avma; B_1 = gerepileupto(av,addii(mulii(q,B_1),B_m1));
        gerepileall(av2,12,&P_m1,&P_0,&P_1,&Q_m1,&Q_0,&Q_1,&T_m1,&T_0,&T_1,&B_m1,&B_0,&B_1);
    }
    av = avma; e = gerepileupto(av,gceil(gmul(powii(gen_2,addii(subii(gmael(fprep,1,2),s),gen_3)),gdiv(T_0,gmael4(fprep,2,1,2,1)))));
    av = avma;
    if (cmpii(mulii(gmael(fprep,2,2),e),powii(gen_2,addii(addii(subii(mulii(gen_2,gmael(fprep,1,2)),gmael(fprep,2,3)),w),gen_3))) <= 0)
    {
        set_avma(av); c = rqiinit(gen_1,Q_0,P_0);
        av = avma; a = gerepileupto(av,diviiexact(subii(T_0,mulii(twossqrtd,B_0)),tmp)); 
        b = B_0;
    }
    else
    {
        set_avma(av); c = rqiinit(gen_1,Q_m1,P_m1);
        av = avma; e = gerepileupto(av,gceil(gmul(powii(gen_2,addii(subii(gmael(fprep,1,2),s),gen_3)),gdiv(T_m1,gmael4(fprep,2,1,2,1)))));
        av = avma; a = gerepileupto(av,diviiexact(subii(T_m1,mulii(twossqrtd,B_m1)),tmp));
        b = B_m1;
    }
    av = avma; t = gerepileupto(av,subsi(sigbits(e)+sigbits(gmael(fprep,2,2))-6,mulii(gen_2,gmael(fprep,1,2))));
    tmp = mulii(e,gmael(fprep,2,2)); tmp2 = powii(gen_2,addii(addis(t,4),mulii(gen_2,gmael(fprep,1,2))));
    if (gcmp(tmp2,tmp) < 0)
    {
        if (cmpii(mulii(tmp2,gen_2),tmp) < 0) t = gerepileupto(av,addii(t,gen_2));
        else t = gerepileupto(av,addii(t,gen_1));
    }
    else set_avma(av);
    av = avma; g = gerepileupto(av,gceil(gdiv(mulii(e,gmael(fprep,2,2)),powii(gen_2,addii(addii(gmael(fprep,1,2),t),gen_3)))));
    h = addii(gmael(fprep,2,3),t);
    res = cgetg(3,t_VEC);
    if (gmael(fprep,1,1) == NULL)
    {
        av = avma; gel(res,1) = gerepileupto(av,fprepinit(NULL,gmael(fprep,1,2),c,g,h));
    }
    else 
    {
        av = avma; gel(res,1) = gerepileupto(av,fprepinit(gadd(gmael(fprep,1,1),gdiv(powii(gen_3,gen_2),powii(gen_2,gen_3))),gmael(fprep,1,2),c,g,h));
    }
    gel(res,2) = mkvec2copy(a,b);
    return gerepileupto(ltop,res);
}

GEN 
wmult(GEN O, GEN fprep1, GEN fprep2, GEN w, long flag)
{
    pari_sp ltop = avma;
    return gerepileupto(ltop,wnear(O,gel(numult(O,fprep1, fprep2, flag),1),w));
}

GEN 
iexp(GEN O, GEN fprep, GEN w, GEN n, long flag)
{
    pari_sp ltop = avma, av;
    GEN ben, fprep_;
    ulong i, l;
    ben = binary_zv(n); l = lg(ben);
    fprep_ = shallowcopy(fprep);
    av = avma;
    for (i = 2; i < l; i++)
    {
        fprep_ = gerepileupto(av,wmult(O,fprep,fprep,w,flag));
        if (ben[i] == 1) fprep_ = gerepileupto(av,wmult(O,fprep_,fprep,w,flag));
    }
    return gerepileupto(ltop,fprep_); //value of f could be better approximated
}

GEN 
addxy(GEN O, GEN fprep1, GEN fprep2, GEN x, GEN y, long flag)
{
    if (typ(x) != t_INT) pari_err_TYPE("ax",x);
    if (typ(y) != t_INT) pari_err_TYPE("ax",y);
    if (cmpii(x,gen_0) <= 0) pari_err_DOMAIN("ax",itostr(x),"<=",gen_0,x);
    if (cmpii(y,gen_0) <= 0) pari_err_DOMAIN("ax",itostr(y),"<=",gen_0,y);
    pari_sp ltop = avma;
    return gerepileupto(ltop,wmult(O,fprep1,fprep2,addii(x,y),flag));
}

GEN 
ax(GEN O, GEN x, GEN p)
{
    if (typ(x) != t_INT) pari_err_TYPE("ax",x);
    if (typ(p) != t_INT) pari_err_TYPE("ax",p);
    if (cmpii(x,gen_0) <= 0) pari_err_DOMAIN("ax",itostr(x),"<=",gen_0,x);
    if (cmpii(p,gen_0) <= 0) pari_err_DOMAIN("ax",itostr(p),"<=",gen_0,p);
    pari_sp ltop = avma, av, av2;
    GEN bex, fprep, s;
    ulong i, l;
    bex = binary_zv(x); l = lg(bex);
    s = gen_1;
    av = avma; fprep = gerepileupto(av,wnear(O,fprepinit(itor(gen_1,DEFAULTPREC),p,pci(O),addii(powii(gen_2,p),gen_1),gen_0),gen_1)); //maybe change the gen_1,DEFAULTPREC)?
    av2 = avma;
    for (i = 2; i < l; i++)
    {
        av = avma; fprep = gerepileupto(av,addxy(O,fprep,fprep,s,s,0));
        s = mulii(gen_2,s);
        if (bex[i] == 1)
        {
            s = addii(s,gen_1);
            av = avma; fprep = gerepileupto(av,wnear(O,fprep,s));
        }
        gerepileall(av2,2,&fprep,&s);
    }
    return gerepileupto(ltop,fprep);
}

GEN 
eaddxy(GEN O, GEN fprep1, GEN fprep2, GEN x, GEN y, long flag)
{
    if (typ(x) != t_INT) pari_err_TYPE("eaddxy",x);
    if (typ(y) != t_INT) pari_err_TYPE("eaddxy",y);
    if (cmpii(x,gen_0) <= 0) pari_err_DOMAIN("eaddxy",itostr(x),"<=",gen_0,x);
    if (cmpii(y,gen_0) <= 0) pari_err_DOMAIN("eaddxy",itostr(y),"<=",gen_0,y);
    pari_sp ltop = avma, av;
    GEN fprep, fprep_, res, res2, r_, r_2;
    fprep = numult(O,fprep1,fprep2,flag);
    av = avma; fprep_ = gerepileupto(av,ewnear(O,gel(fprep,1),addii(x,y)));
    res = cgetg(3,t_VEC);
    if (gmael3(fprep_,1,1,1) == NULL) //in this case gcopy fails...
    {
        r_ = cgetg(3,t_VEC);
        r_2 = cgetg(3,t_VEC);
        gel(r_2,1) = NULL;
        gel(r_2,2) = gcopy(gmael3(fprep_,1,1,2));
        gel(r_,1) = r_2;
        gel(r_,2) = gcopy(gmael(fprep_,1,2));
        gel(res,1) = r_;
    }
    else gel(res,1) = gcopy(gel(fprep_,1));
    res2 = cgetg(3,t_VEC);
    av = avma; gel(res2,1) = gerepileupto(av,diviiexact(addii(mulii(gmael(fprep,2,1),gmael(fprep_,2,1)),mulii(gel(O,1),mulii(gmael(fprep,2,2),gmael(fprep_,2,2)))),gmael5(fprep,1,2,1,2,1)));
    av = avma; gel(res2,2) = gerepileupto(av,diviiexact(addii(mulii(gmael(fprep,2,1),gmael(fprep_,2,2)),mulii(gmael(fprep,2,2),gmael(fprep_,2,1))),gmael5(fprep,1,2,1,2,1)));
    gel(res,2) = res2;
    return gerepileupto(ltop,res);
}

GEN
crax(GEN O, GEN x, GEN p, GEN m)
{
    if (typ(x) != t_INT) pari_err_TYPE("crax",x);
    if (typ(p) != t_INT) pari_err_TYPE("crax",p);
    if (typ(m) != t_INT) pari_err_TYPE("crax",m);
    if (cmpii(x,gen_0) <= 0) pari_err_DOMAIN("crax",itostr(x),"<=",gen_0,x);
    if (cmpii(p,gen_0) <= 0) pari_err_DOMAIN("crax",itostr(p),"<=",gen_0,p);
    pari_sp ltop = avma, av, av2, av3;
    if (gcmp(powii(gen_2,p),gmax_shallow(powii(gen_2,addii(gen_2,gen_1)),gdiv(glog(x,DEFAULTPREC),mplog2(DEFAULTPREC)))) <= 0) pari_err_PREC("p"); //different error?
    set_avma(ltop);
    GEN bex, s, fprep, fprep_, list, el, el2, el2_, N, eli, res, beta, beta_;
    GEN h_1, h_2, h_3, b, r_1, r_2, r_3, s_1, s_2, spl, x_, y_, N_, Nc, tmp;
    res = cgetg(3,t_VEC);
    ulong i, l;
    av3 = avma;
    bex = binary_zv(x); l = lg(bex);
    av2 = avma;
    list = cgetg(l,t_VEC);
    for (i = 1; i < l; i++) gel(list,i) = gen_0; //initialize list so we can gerepile it
    av = avma; fprep = gerepileupto(av,ewnear(O,fprepinit(itor(gen_1,DEFAULTPREC),p,pci(O),addii(powii(gen_2,p),gen_1),gen_0),gen_1));
    s = gen_1;
    el = cgetg(3,t_VEC);
    gel(el,1) = mkvec2copy(gmael(fprep,2,1),gmael(fprep,2,2));
    gel(el,2) = gen_1;
    beta = mkvec2(gel(O,3),gen_0);
    Nc = gen_1;
    for (i = 2; i < l; i++)
    {
        N = diviiexact(gmael5(fprep,1,2,1,2,1),gel(O,3));
        fprep = eaddxy(O,gel(fprep,1),gel(fprep,1),s,s,0);
        s = mulii(gen_2,s);
        if (bex[i] == 1)
        {
            s = addii(s,gen_1);
            fprep_ = ewnear(O,gel(fprep,1),s);
            el2 = cgetg(3,t_VEC);
            eli = cgetg(3,t_VEC);
            av = avma; gel(eli,1) = gerepileupto(av,diviiexact(addii(mulii(gmael(fprep,2,1),gmael(fprep_,2,1)),mulii(gel(O,1),mulii(gmael(fprep,2,2),gmael(fprep_,2,2)))),gmael5(fprep,1,2,1,2,1)));
            av = avma; gel(eli,2) = gerepileupto(av,diviiexact(addii(mulii(gmael(fprep,2,1),gmael(fprep_,2,2)),mulii(gmael(fprep,2,2),gmael(fprep_,2,1))),gmael5(fprep,1,2,1,2,1)));
            gel(el2,1) = eli;
            fprep = fprep_;
        }
        else
        {
            el2 = cgetg(3,t_VEC);
            gel(el2,1) = mkvec2copy(gmael(fprep,2,1),gmael(fprep,2,2));
        }
        gel(el2,2) = gcopy(N);
        if (cmpii(m,gen_0))
        {
            if(i > 2)
            {
                beta_ = beta;
                N_ = diviiexact(gmael5(fprep,1,2,1,2,1),gel(O,3));
                if (i < l-1 && cmpii(gcdii(N_,m),gen_1) > 0)
                {
                    h_1 = N_;
                    av = avma; tmp = gerepileupto(av,mulii(gmael(O,2,2),subii(gel(O,3),gen_1)));
                    av = avma; b = gerepileupto(av,modii(diviiexact(subii(gmael5(fprep,1,2,1,2,2),tmp),gel(O,3)),h_1));
                    av = avma; h_3 = gerepileupto(av,diviiexact(subii(sqri(addii(mulii(gel(O,3),b),tmp)),gel(O,1)),mulii(h_1,sqri(gel(O,3)))));
                    av = avma; h_2 = gerepileupto(av,addii(diviiexact(mulii(gen_2,addii(mulii(gel(O,3),b),tmp)),gel(O,3)),addii(h_1,h_3)));
                    spl = split(h_1,m); r_1 = gel(spl,1); s_1 = gel(spl,2);
                    spl = split(h_2,s_1); r_2 = gel(spl,1); s_2 = gel(spl,2);
                    spl = split(h_3,s_2); r_3 = gel(spl,1); // s_3 not needed
                    x_ = lift(chinese1_coprime_Z(mkvec2(mkintmod(gen_1,mulii(r_1,r_2)),mkintmod(gen_0,r_3))));
                    y_ = lift(chinese1_coprime_Z(mkvec2(mkintmod(gen_1,mulii(r_2,r_3)),mkintmod(gen_0,r_1))));
                    beta = cgetg(3,t_VEC);
                    av = avma; gel(beta,1) = gerepileupto(av,addii(mulii(gel(O,3),addii(mulii(x_,h_1),mulii(b,y_))),mulii(y_,tmp)));
                    gel(beta,2) = gcopy(y_);
                }
                else
                {
                    beta = mkvec2(mulii(N_,gel(O,3)),gen_0);
                } 
                el2_ = cgetg(3,t_VEC);
                av = avma; gel(el2_,1) = gerepileupto(av,mulqidivn(O,mulqidivn(O,invqi(beta),gel(el2,1),gen_1),mulqidivn(O,beta_,beta_,gen_1),mulii(N_,sqri(N))));
                gel(el2_,2) = Nc;
                av = avma; Nc = gerepileupto(av,diviiexact(subii(sqri(gel(beta,1)),mulii(sqri(gel(beta,2)),gel(O,1))),mulii(sqri(gel(O,3)),N_)));
                el2 = el2_;
            }
            gerepileall(av2,7,&el2,&s,&fprep,&el,&list,&beta,&Nc);
        }
        else gerepileall(av2,5,&el2,&s,&fprep,&el,&list);
        gel(list,i-1) = el;
        el = el2;
    }
    fprep = gel(fprep,1);
    gel(list,i-1) = el;
    gerepileall(av3,2,&fprep,&list);
    gel(res,1) = fprep;
    gel(res,2) = list;
    return gerepileupto(ltop,res); 
}

GEN 
find(GEN O, GEN fprep, GEN c)
{
    pari_sp ltop = avma, av, av2;
    GEN s, Q_0, Q_1, P_0, P_1, G_0, G_1, B_0, B_1, swap, sqrtd_, q, T, e, tmp, tmp2, gen_3, t, res;
    tmp = gsqrt(gel(O,1),DEFAULTPREC);
    if (cmpri(addir(gmael4(fprep,2,1,2,2),tmp),gmael4(fprep,2,1,2,1)) <= 0) pari_err_DOMAIN("find",itostr(gmael4(fprep,2,1,2,1)),">=",addir(gmael4(fprep,2,1,2,2),tmp),gmael(fprep,2,1));
    tmp = subir(gmael4(fprep,2,1,2,2),tmp);
    if (cmpir(negi(gmael4(fprep,2,1,2,1)),tmp) >= 0 || cmpri(tmp,gen_0) >= 0) pari_err_DOMAIN("find",GENtostr(tmp),"",NULL,gmael(fprep,2,1));
    set_avma(ltop);
    sqrtd_ = sqrti(gel(O,1));
    s = addis(gmael(fprep,1,2),5-sigbits(gmael4(fprep,2,1,2,1)));
    av = avma; if (cmpii(gmael4(fprep,2,1,2,1),powii(gen_2,subii(addis(gmael(fprep,1,2),4),s))) <= 0) s = gerepileupto(av,addii(s,gen_1));
    else set_avma(av);
    Q_1 = gmael4(fprep,2,1,2,1);
    P_1 = gmael4(fprep,2,1,2,2);
    av2 = avma;
    av = avma; Q_0 = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(gmael4(fprep,2,1,2,2))),gmael4(fprep,2,1,2,1)));
    G_0 = negi(gmael4(fprep,2,1,2,2));
    G_1 = gmael4(fprep,2,1,2,1);
    B_0 = gen_1;
    B_1 = gen_0;
    while (cmpii(Q_1,gmael(c,2,1)) || cmpii(P_1,gmael(c,2,2)))
    {
        av = avma; q = gerepileupto(av,divii(addii(P_1,sqrtd_),Q_1));
        P_0 = P_1;
        av = avma; P_1 = gerepileupto(av,subii(mulii(q,Q_1),P_1));
        swap = Q_1;
        av = avma; Q_1 = gerepileupto(av,subii(Q_0,mulii(q,subii(P_1,P_0))));
        Q_0 = swap;
        swap = B_1;
        av = avma; B_1 = gerepileupto(av,addii(mulii(q,B_1),B_0));
        B_0 = swap;
        swap = G_1;
        av = avma; G_1 = gerepileupto(av,addii(mulii(q,G_1),G_0));
        G_0 = swap;
        gerepileall(av2,8,&P_0,&P_1,&Q_0,&Q_1,&B_0,&B_1,&G_0,&G_1);
    }
    tmp = powii(gen_2,s);
    gen_3 = addii(gen_2,gen_1);
    av = avma; T = gerepileupto(av,addii(mulii(tmp,G_1),mulii(B_1,gfloor(gmul(tmp,gsqrt(gel(O,1),DEFAULTPREC)))))); //should throw error if precision is too low
    av = avma; e = gerepileupto(av,gceil(gmul(powii(gen_2,subii(addii(gmael(fprep,1,2),gen_3),s)),gdiv(T,gmael4(fprep,2,1,2,1)))));
    av = avma; t = gerepileupto(av,subsi(sigbits(e)+sigbits(gmael(fprep,2,2))-6,mulii(gen_2,gmael(fprep,1,2))));
    tmp = mulii(e,gmael(fprep,2,2)); tmp2 = powii(gen_2,addii(addis(t,4),mulii(gen_2,gmael(fprep,1,2))));
    if (gcmp(tmp2,tmp) < 0)
    {
        if (cmpii(mulii(tmp2,gen_2),tmp) < 0) t = gerepileupto(av,addii(t,gen_2));
        else t = gerepileupto(av,addii(t,gen_1));
    }
    else set_avma(av);
    res = cgetg(3,t_VEC);
    if (gmael(fprep,1,1) == NULL)
    {
        av = avma; gel(res,1) = gerepileupto(av,fprepinit(NULL,gmael(fprep,1,2),rqiinit(gen_1,Q_1,P_1),gceil(gdiv(mulii(e,gmael(fprep,2,2)),powii(gen_2,addii(addii(gmael(fprep,1,2),t),gen_3)))),addii(gmael(fprep,2,3),t)));
    }
    else
    {
        av = avma; gel(res,1) = gerepileupto(av,fprepinit(gadd(gmael(fprep,1,1),ginv(addii(gen_2,gen_2))),gmael(fprep,1,2),rqiinit(gen_1,Q_1,P_1),gceil(gdiv(mulii(e,gmael(fprep,2,2)),powii(gen_2,addii(addii(gmael(fprep,1,2),t),gen_3)))),addii(gmael(fprep,2,3),t)));
    }
    gel(res,2) = mkvec2copy(G_1,B_1);
    return gerepileupto(ltop,res);
}

GEN 
cr(GEN O, GEN b, GEN y, GEN q, GEN m)
{
    if (typ(y) != t_FRAC && typ(y) != t_INT) pari_err_TYPE("cr",y);
    if (typ(q) != t_FRAC && typ(q) != t_INT) pari_err_TYPE("cr",q);
    if (gcmp(y,gen_0) <= 0) pari_err_DOMAIN("cr",GENtostr(y),"<=",gen_0,y);
    if (gcmp(q,gen_0) <= 0) pari_err_DOMAIN("cr",GENtostr(y),"<=",gen_0,q);
    pari_sp ltop = avma, av;
    GEN x, p, fprep, fprep_, res, res_, list;
    ulong i;
    //GEN c, gen_3 = addii(gen_2,gen_1);
    //av = avma; c = gerepileupto(av,addii(addii(gen_3,gen_3),gceil(gmul(gen_3,q)))); //this is an upper bound on the steps performed by find
    av = avma; x = gerepileupto(av,subii(gfloor(gsub(y,q)),gen_1));
    av = avma; p = gerepileupto(av,addsi(4.48543+sigbits(x),gmax_shallow(stoi(4),stoi(sigbits(stoi(sigbits(x)+1))+1)))); //generous upper bound
    fprep = crax(O,x,p,m);
    fprep_ = find(O,gel(fprep,1),b);
    res_ = cgetg(3,t_VEC);
    gel(res_,1) = gel(fprep_,2);
    gel(res_,2) = diviiexact(gmael5(fprep,1,2,1,2,1),gel(O,3));
    res = cgetg(3,t_VEC);
    gel(res,1) = gcopy(gel(fprep_,1));
    list = cgetg(lg(gel(fprep,2))+1,t_VEC);
    for (i = 1; i < lg(gel(fprep,2)); i++) //there must be a better way...
    {
        gel(list,i) = gcopy(gel(gel(fprep,2),i));
    }
    gel(list,i) = gcopy(res_);
    gel(res,2) = list;
    return gerepileupto(ltop,res);
}
