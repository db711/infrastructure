#include "fp_representations.h"

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
    long t;
    av = avma; e = gerepileupto(av,ground(gmul(powii(gen_2,subii(addii(gmael(fprep,1,2),gen_3),s)),gabs(gdiv(T,C),DEFAULTPREC)))); //precision should be ignored
    av = avma; t = dbllog2r(itor(e,DEFAULTPREC))-2-dbllog2r(itor(gmael(fprep,2,2),DEFAULTPREC)); set_avma(av); //is this enough?
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
        e = gerepileupto(av,divii(mulii(gmael(fprep1,2,2),gmael(fprep2,2,2)),powii(gen_2,addii(gmael(fprep1,1,2),gen_1))));
        av = avma; h = gerepileupto(av,addii(addii(gmael(fprep1,2,3),gmael(fprep2,2,3)),gen_1));
    }
    av = avma; s = gerepileupto(av,ceilr(addir(addii(gmael(fprep1,1,2),addii(gen_2,gen_2)),dbltor(dbllog2r(itor(gmael(b,2,2),DEFAULTPREC))-dbllog2r(itor(gmael3(b,1,2,1),DEFAULTPREC)))))); //is this enough?
    if ((cmpii(s,gen_0) < 0)) s = gen_0;
    av = avma; T = gerepileupto(av,addii(mulii(powii(gen_2,s),gmael(b,2,1)),mulii(gmael(b,2,2),floorr(mulir(powii(gen_2,s),gsqrt(gel(O,1),DEFAULTPREC)))))); //is this enough?
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
    GEN s, Q_0, Q_1, Q_m1, Q_m2, Q_m3, P_0, P_1, P_m1, M, T_0, T_1, T_m1, T_m2, T_m3, q, sqrtd, sqrtd_, tmp, gen_3 = addii(gen_2,gen_1), e, e_, c, t, g, h;
    sqrtd = gsqrt(gel(O,1),DEFAULTPREC); //is this enough?
    sqrtd_ = sqrti(gel(O,1));
    if (cmpii(gmael(fprep,2,3),w) < 0)
    {
        av = avma; s = gerepileupto(av,addis(gmael(fprep,1,2),5-dbllog2r(itor(gmael4(fprep,2,1,2,1),DEFAULTPREC)))); //is this enough?
        tmp = powii(gen_2,s);
        Q_1 = gmael4(fprep,2,1,2,1);
        P_1 = gmael4(fprep,2,1,2,2);
        av = avma; M = gerepileupto(av,gceil(gmul(powii(gen_2,addii(subii(addii(gmael(fprep,1,2),s),gmael(fprep,2,3)),w)),gdiv(gmael4(fprep,2,1,2,1),gmael(fprep,2,2)))));
        av2 = avma;
        P_m1 = P_0 = Q_m1 = T_m1 = gen_0; //dummy values
        av = avma; Q_0 = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(gmael4(fprep,2,1,2,2))),gmael4(fprep,2,1,2,1)));
        av = avma; T_0 = gerepileupto(av,addii(mulii(negi(tmp),P_1),floorr(mulir(tmp,sqrtd))));
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
        av = avma; t = gerepileupto(av,subsi(dbllog2r(itor(mulii(e,gmael(fprep,2,2)),DEFAULTPREC))-3,mulii(gen_2,gmael(fprep,1,2)))); //is this enough?
        //also set for next element and initialize that   
        av = avma; g = gerepileupto(av,gceil(gdiv(mulii(e,gmael(fprep,2,2)),powii(gen_2,addii(addii(gmael(fprep,1,2),t),gen_3)))));
        h = addii(gmael(fprep,2,3),t);
    }
    else
    {
        av = avma; s = gerepileupto(av,addii(gmael(fprep,1,2),addii(gen_2,gen_2)));
        tmp = powii(gen_2,s);
        Q_0 = gmael4(fprep,2,1,2,1);
        P_1 = gmael4(fprep,2,1,2,2); 
        av = avma; M = gerepileupto(av,mulii(gmael(fprep,2,2),powii(gen_2,addii(subii(gmael(fprep,2,3),w),addii(gen_2,gen_2)))));
        av2 = avma;
        Q_m3 = Q_m2 = Q_m1 = T_m3 = T_m2 = T_m1 = gen_0; //dummy values
        av = avma; Q_1 = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(gmael4(fprep,2,1,2,2))),gmael4(fprep,2,1,2,1)));
        av = avma; T_0 = gerepileupto(av,mulii(tmp,Q_0));
        av = avma; T_1 = gerepileupto(av,addii(mulii(tmp,P_1),floorr(mulir(tmp,sqrtd))));
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
        t = stoi(dbllog2r(itor(e,DEFAULTPREC))-2-dbllog2r(itor(gmael(fprep,2,2),DEFAULTPREC))); //is this enough?
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
    if (cmpii(gmael(fprep,2,3),w) >= 0) pari_err_DOMAIN("ewnear","w","<=",gmael(fprep,2,3),w); //TOD: Implement other case
    pari_sp ltop = avma, av, av2;
    if (cmpii(addii(gmael4(fprep,2,1,2,2),sqrti(gel(O,1))),gmael4(fprep,2,1,2,1)) < 0) pari_err_DOMAIN("ewnear",itostr(addii(gmael4(fprep,2,1,2,2),sqrti(gel(O,1)))),"<",gmael4(fprep,2,1,2,1),gmael(fprep,2,1));
    if (cmpii(gen_0,subii(sqrti(gel(O,1)),gmael4(fprep,2,1,2,2))) > 0 || cmpii(subii(sqrti(gel(O,1)),gmael4(fprep,2,1,2,2)),gmael4(fprep,2,1,2,1)) > 0) pari_err_DOMAIN("ewnear",itostr(subii(sqrti(gel(O,1)),gmael4(fprep,2,1,2,2))),"",NULL,gmael(fprep,2,1));
    set_avma(ltop);
    GEN B_0, B_1, B_m1, a, b, s, Q_0, Q_1, Q_m1, P_0, P_1, P_m1, M, T_0, T_1, T_m1, q, sqrtd, sqrtd_, tmp, gen_3 = addii(gen_2,gen_1), e, c, t, g, h, res;
    sqrtd = gsqrt(gel(O,1),DEFAULTPREC); //is this enough?
    sqrtd_ = floorr(sqrtd);
    B_0 = gen_1;
    B_1 = gen_0;
    av = avma; s = gerepileupto(av,addis(gmael(fprep,1,2),5-dbllog2r(itor(gmael4(fprep,2,1,2,1),DEFAULTPREC)))); //is this enough?
    tmp = powii(gen_2,s);
    Q_1 = gmael4(fprep,2,1,2,1);
    P_1 = gmael4(fprep,2,1,2,2);
    av = avma; M = gerepileupto(av,gceil(gmul(powii(gen_2,addii(subii(addii(gmael(fprep,1,2),s),gmael(fprep,2,3)),w)),gdiv(gmael4(fprep,2,1,2,1),gmael(fprep,2,2)))));
    av2 = avma;
    P_m1 = P_0 = Q_m1 = T_m1 = gen_0; //dummy values
    av = avma; Q_0 = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(gmael4(fprep,2,1,2,2))),gmael4(fprep,2,1,2,1)));
    av = avma; T_0 = gerepileupto(av,addii(mulii(negi(tmp),P_1),floorr(mulir(tmp,sqrtd))));
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
        av = avma; a = gerepileupto(av,diviiexact(subii(T_0,mulii(floorr(mulir(tmp,sqrtd)),B_0)),tmp));
        b = B_0;
    }
    else
    {
        set_avma(av); c = rqiinit(gen_1,Q_m1,P_m1);
        av = avma; e = gerepileupto(av,gceil(gmul(powii(gen_2,addii(subii(gmael(fprep,1,2),s),gen_3)),gdiv(T_m1,gmael4(fprep,2,1,2,1)))));
        av = avma; a = gerepileupto(av,diviiexact(subii(T_m1,mulii(floorr(mulir(tmp,sqrtd)),B_m1)),tmp));
        b = B_m1;
    }
    av = avma; t = gerepileupto(av,subsi(dbllog2r(itor(mulii(e,gmael(fprep,2,2)),DEFAULTPREC))-3,mulii(gen_2,gmael(fprep,1,2)))); //is this enough? no...
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
    long i, l;
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
    long i, l;
    bex = binary_zv(x); l = lg(bex);
    s = gen_1;
    av = avma; fprep = gerepileupto(av,wnear(O,fprepinit(itor(gen_1,DEFAULTPREC),p,rqiinit(gen_1,gel(O,3),subii(addii(mulii(gel(O,3),divii(addii(subii(sqrti(gel(O,1)),gel(O,3)),gen_1),gel(O,3))),gel(O,3)),gen_1)),addii(powii(gen_2,p),gen_1),gen_0),gen_1)); //maybe change the gen_1,DEFAULTPREC)?
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
    if (typ(x) != t_INT) pari_err_TYPE("ax",x);
    if (typ(y) != t_INT) pari_err_TYPE("ax",y);
    if (cmpii(x,gen_0) <= 0) pari_err_DOMAIN("ax",itostr(x),"<=",gen_0,x);
    if (cmpii(y,gen_0) <= 0) pari_err_DOMAIN("ax",itostr(y),"<=",gen_0,y);
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
    av = avma; gel(res2,1) = gerepileupto(av,diviiexact(addii(mulii(gmael(fprep,2,1),gmael(fprep_,2,1)),mulii(gel(O,1),mulii(gmael(fprep,2,2),gmael(fprep_,2,2)))),diviiexact(gmael5(fprep,1,2,1,2,1),gel(O,3))));
    av = avma; gel(res2,2) = gerepileupto(av,diviiexact(addii(mulii(gmael(fprep,2,1),gmael(fprep_,2,2)),mulii(gmael(fprep,2,2),gmael(fprep_,2,1))),diviiexact(gmael5(fprep,1,2,1,2,1),gel(O,3))));
    gel(res,2) = res2;
    return gerepileupto(ltop,res);
}

GEN
crax(GEN O, GEN x, GEN p)
{
    if (typ(x) != t_INT) pari_err_TYPE("crax",x);
    if (typ(p) != t_INT) pari_err_TYPE("crax",p);
    if (cmpii(x,gen_0) <= 0) pari_err_DOMAIN("crax",itostr(x),"<=",gen_0,x);
    if (cmpii(p,gen_0) <= 0) pari_err_DOMAIN("crax",itostr(p),"<=",gen_0,p);
    pari_sp ltop = avma, av, av2, av3;
    if (gcmp(powii(gen_2,p),gmax_shallow(powii(gen_2,addii(gen_2,gen_1)),gdiv(glog(x,DEFAULTPREC),mplog2(DEFAULTPREC)))) <= 0) pari_err_PREC("p"); //different error?
    set_avma(ltop);
    GEN bex, s, fprep, fprep_, list, el, el2, N, eli, res;
    res = cgetg(3,t_VEC);
    long i, l;
    av3 = avma;
    bex = binary_zv(x); l = lg(bex);
    av2 = avma;
    list = cgetg(l,t_VEC);
    for (i = 1; i < l; i++) gel(list,i) = gen_0; //initialize list so we can gerepile it
    av = avma; fprep = gerepileupto(av,ewnear(O,fprepinit(itor(gen_1,DEFAULTPREC),p,rqiinit(gen_1,gel(O,3),addii(mulii(gel(O,3),divii(addii(subii(sqrti(gel(O,1)),gel(O,3)),gen_1),gel(O,3))),subii(gel(O,3),gen_1))),addii(powii(gen_2,p),gen_1),gen_0),gen_1));
    s = gen_1;
    el = cgetg(3,t_VEC);
    gel(el,1) = mkvec2copy(gmael(fprep,2,1),gmael(fprep,2,2));
    gel(el,2) = gen_1;
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
        gerepileall(av2,5,&el2,&s,&fprep,&el,&list);
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