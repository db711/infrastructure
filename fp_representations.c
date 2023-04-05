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
    av = avma; t = dbllog2r(itor(e,DEFAULTPREC)) - 2.5 - dbllog2r(itor(gmael(fprep,2,2),DEFAULTPREC)); set_avma(av);
    av = avma; d = gerepileupto(av,ceilr(mulir(powii(gen_2,addis(addii(gmael(fprep,1,2),gen_3),t)),rdivii(gmael(fprep,2,2),e,DEFAULTPREC))));
    av = avma; k = gerepileupto(av,subis(gmael(fprep,2,3),t));
    if (gmael(fprep,1,1) == NULL) return gerepileupto(ltop,fprepinit(NULL,gmael(fprep,1,2),gmael(fprep,2,1),d,k));
    else return gerepileupto(ltop,fprepinit(gadd(gmael(fprep,1,1),rdivii(powii(gen_3,gen_2),powii(gen_2,gen_3),DEFAULTPREC)),gmael(fprep,1,2),gmael(fprep,2,1),d,k));
}

GEN
numult (GEN O, GEN fprep1, GEN fprep2, long flag)
{
    if (flag && cmpii(gmael(fprep1,1,2),gmael(fprep2,1,2))) pari_err_DOMAIN("numult",itostr(gmael(fprep2,1,2)),"!=",gmael(fprep1,1,2),gmael(fprep2,1,2));
    GEN b, e, h, s, T;
    pari_sp ltop = avma, av;
    b = inucomp(O,gmael(fprep1,2,1),gmael(fprep2,2,1),flag);
    pari_printf("%Ps\n",b);
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