#include "fp_representations.h"

GEN
fprepinit (GEN f, GEN p, GEN b, GEN d, GEN k)
{
    GEN res = cgetg(3,t_VEC), res1 = cgetg(3,t_VEC), res2 = cgetg(4,t_VEC);
    pari_sp ltop = avma;
    if (f != NULL && typ(f) != t_REAL) pari_err_TYPE("fprepinit",f); //could be fixed
    if (typ(p) != t_INT) pari_err_TYPE("fprepinit",p);
    if (typ(d) != t_INT) pari_err_TYPE("fprepinit",d);
    if (typ(k) != t_INT) pari_err_TYPE("fprepinit",k);
    if (f != NULL && gcmp(f,gen_1) < 0) pari_err_DOMAIN("fprepinit",GENtostr(f),"<",gen_1,f);
    if (cmpii(p,gen_0) < 0) pari_err_DOMAIN("fprepinit",itostr(p),"<",gen_0,p);
    if (cmpii(powii(gen_2,p),d) >= 0 || cmpii(powii(gen_2,addii(p,gen_1)),d) < 0) pari_err_DOMAIN("fprepinit",itostr(d),"",NULL,d);
    set_avma(ltop);
    gel(res1,1) = (f == NULL) ? NULL : gcopy(f);
    gel(res1,2) = gcopy(p);
    gel(res2,1) = gcopy(b);
    gel(res2,2) = gcopy(d);
    gel(res2,3) = gcopy(k);
    gel(res,1) = res1;
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
    av = avma; e = gerepileupto(av,roundr(mulir(powii(gen_2,subii(addii(gmael(fprep,1,2),gen_3),s)),absr(rdivii(T,C,DEFAULTPREC)))));
    av = avma; t = dbllog2r(itor(e,DEFAULTPREC))-2-dbllog2r(itor(gmael(fprep,2,2),DEFAULTPREC)); set_avma(av);
    av = avma; d = gerepileupto(av,ceilr(mulir(powii(gen_2,addis(addii(gmael(fprep,1,2),gen_3),t)),rdivii(gmael(fprep,2,2),e,DEFAULTPREC))));
    av = avma; k = gerepileupto(av,subis(gmael(fprep,2,3),t));
    if (gmael(fprep,1,1) == NULL) return gerepileupto(ltop,fprepinit(NULL,gmael(fprep,1,2),gmael(fprep,2,1),d,k));
    else return gerepileupto(ltop,fprepinit(addrr(gmael(fprep,1,1),rdivii(powii(gen_3,gen_2),powii(gen_2,gen_3),DEFAULTPREC)),gmael(fprep,1,2),gmael(fprep,2,1),d,k));
}

GEN
numult (GEN O, GEN fprep1, GEN fprep2, long flag)
{
    if (flag && cmpii(gmael(fprep1,1,2),gmael(fprep2,1,2))) pari_err_DOMAIN("numult",itostr(gmael(fprep2,1,2)),"!=",gmael(fprep1,1,2),gmael(fprep2,1,2));
    GEN b, e, h, s, T;
    pari_sp ltop = avma, av;
    b = inucomp(O,gmael(fprep1,2,1),gmael(fprep2,2,1),flag);
    av = avma;
    if (cmpii(mulii(gmael(fprep1,2,2),gmael(fprep2,2,2)),powii(gen_2,addii(mulii(gen_2,gmael(fprep1,1,2)),gen_1))) <= 0)
    {
        e = gerepileupto(av,floorr(rdivii(mulii(gmael(fprep1,2,2),gmael(fprep2,2,2)),powii(gen_2,gmael(fprep1,1,2)),DEFAULTPREC)));
        h = addii(gmael(fprep1,2,3),gmael(fprep2,2,3));
    }
    else
    {
        e = gerepileupto(av,floorr(rdivii(mulii(gmael(fprep1,2,2),gmael(fprep2,2,2)),powii(gen_2,addii(gmael(fprep1,1,2),gen_1)),DEFAULTPREC)));
        av = avma; h = gerepileupto(av,addii(addii(gmael(fprep1,2,3),gmael(fprep2,2,3)),gen_1));
    }
    av = avma; s = gerepileupto(av,ceilr(addir(addii(gmael(fprep1,1,2),addii(gen_2,gen_2)),dbltor(dbllog2r(itor(gmael(b,2,2),DEFAULTPREC))-dbllog2r(itor(gmael3(b,1,2,1),DEFAULTPREC))))));
    if ((cmpii(s,gen_0) < 0)) s = gen_0;
    av = avma; T = gerepileupto(av,addii(mulii(powii(gen_2,s),gmael(b,2,1)),mulii(gmael(b,2,2),floorr(mulir(powii(gen_2,s),gsqrt(gel(O,1),DEFAULTPREC))))));
    if (gmael(fprep1,1,1) == NULL || gmael(fprep2,1,1) == NULL) return gerepileupto(ltop,fpremove(fprepinit(NULL,gmael(fprep1,1,2),gel(b,1),e,h),T,gmael(b,2,3),s));
    else return gerepileupto(ltop,fpremove(fprepinit(addri(addrr(addrr(gmael(fprep1,1,1),gmael(fprep2,1,1)),gmul(powii(gen_2,negi(gmael(fprep1,1,2))),mulrr(gmael(fprep1,1,1),gmael(fprep2,1,1)))),gen_1),gmael(fprep1,1,2),gel(b,1),e,h),T,gmael(b,2,3),s));
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
    sqrtd = gsqrt(gel(O,1),DEFAULTPREC);
    sqrtd_ = floorr(sqrtd);
    if (cmpii(gmael(fprep,2,3),w) < 0)
    {
        av = avma; s = gerepileupto(av,addis(gmael(fprep,1,2),5-dbllog2r(itor(gmael4(fprep,2,1,2,1),DEFAULTPREC))));
        tmp = powii(gen_2,s);
        Q_1 = gmael4(fprep,2,1,2,1);
        P_1 = gmael4(fprep,2,1,2,2);
        av = avma; M = gerepileupto(av,ceilr(mulir(powii(gen_2,addii(subii(addii(gmael(fprep,1,2),s),gmael(fprep,2,3)),w)),rdivii(gmael4(fprep,2,1,2,1),gmael(fprep,2,2),DEFAULTPREC))));
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
        av = avma; t = gerepileupto(av,subsi(dbllog2r(itor(mulii(e,gmael(fprep,2,2)),DEFAULTPREC))-3,mulii(gen_2,gmael(fprep,1,2))));
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
        av = avma; e = gerepileupto(av,ceilr(rdivii(T_1,mulii(gen_2,Q_1),DEFAULTPREC)));
        av = avma; e_ = gerepileupto(av,ceilr(rdivii(T_0,mulii(gen_2,Q_0),DEFAULTPREC)));
        av2 = avma;
        av = avma; tmp = gerepileupto(av,mulii(gmael(fprep,2,2),powii(gen_2,addii(subii(gmael(fprep,2,3),w),gen_3))));
        if (cmpii(e_,tmp) >= 0)
        {
            e = e_;
            av = avma; e_ = gerepileupto(av,ceilr(rdivii(T_m2,mulii(gen_2,Q_m2),DEFAULTPREC)));
            if (cmpii(e_,tmp) >= 0)
            {
                e = e_;
                av = avma; e_ = gerepileupto(av,ceilr(rdivii(T_m3,mulii(gen_2,Q_m3),DEFAULTPREC)));
                if (cmpii(e_,tmp) >= 0) pari_err_IMPL("this");
            }
            c = rqiinit(gen_1,Q_0,P_0);
            //also set next element
        }
        else c = rqiinit(gen_1,Q_1,P_1); //also set next element
        t = stoi(dbllog2r(itor(e,DEFAULTPREC))-2-dbllog2r(itor(gmael(fprep,2,2),DEFAULTPREC)));
        //also set for next element and initialize that
        av = avma; g = gerepileupto(av,gceil(gmul(powii(gen_2,addii(addii(gmael(fprep,1,2),gen_3),t)),gdiv(gmael(fprep,2,2),e))));
        h = subii(gmael(fprep,2,3),t);
    }
    if (gmael(fprep,1,1) == NULL) return gerepileupto(ltop,fprepinit(NULL,gmael(fprep,1,2),c,g,h));
    else return gerepileupto(ltop,fprepinit(addrr(gmael(fprep,1,1),rdivii(powii(gen_3,gen_2),powii(gen_2,gen_3),DEFAULTPREC)),gmael(fprep,1,2),c,g,h));
}

GEN 
wmult(GEN O, GEN fprep1, GEN fprep2, GEN w, long flag)
{
    pari_sp ltop = avma;
    return gerepileupto(ltop,wnear(O,numult(O,fprep1,fprep2,flag),w));
}