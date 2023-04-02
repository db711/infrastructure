#include "fp_representations.h"

struct fp_rep 
fprepinit (GEN f, GEN p, GEN b, GEN d, GEN k)
{
    struct fp_rep a;
    pari_sp av;
    if (f != NULL && typ(f) != t_REAL) pari_err_TYPE("fprepinit",f); //could be fixed
    if (typ(p) != t_INT) pari_err_TYPE("fprepinit",p);
    if (typ(d) != t_INT) pari_err_TYPE("fprepinit",d);
    if (typ(k) != t_INT) pari_err_TYPE("fprepinit",k);
    if (f != NULL && gcmp(f,gen_1) < 0) pari_err_DOMAIN("fprepinit",GENtostr(f),"<",gen_1,f);
    if (cmpii(p,gen_0) < 0) pari_err_DOMAIN("fprepinit",itostr(p),"<",gen_0,p);
    av = avma;
    if (cmpii(powii(gen_2,p),d) >= 0 || cmpii(d,powii(gen_2,addii(p,gen_1))) < 0) pari_err_DOMAIN("fprepinit",itostr(d),"",NULL,d);
    set_avma(av);
    a.f = (f == NULL) ? NULL : gcopy(f);
    a.p = p;
    a.b = gcopy(b);
    a.d = gcopy(d);
    a.k = gcopy(k);
    return a;
}

struct fp_rep 
fpremove (struct fp_rep fprep, GEN T, GEN C, GEN s)
{
    if (typ(C) != t_INT) pari_err_TYPE("fpremove",C);
    if (typ(T) != t_INT) pari_err_TYPE("fpremove",T);
    if (typ(s) != t_INT) pari_err_TYPE("fpremove",s);
    if (equalii(C, gen_0)) pari_err_DOMAIN("fpremove","C","==",gen_0,C);
    pari_sp ltop = avma;
    GEN gen_3 = addii(gen_2,gen_1), e, t, d, k, f;
    struct fp_rep fprep_;
    e = roundr(mulir(powii(gen_2,subii(addii(fprep.p,gen_3),s)),absr(rdivii(T,C,DEFAULTPREC))));
    t = ceilr(divrr(mplog(rdivii(e,mulii(powii(gen_2,gen_3),fprep.d),DEFAULTPREC)),mplog2(DEFAULTPREC)));
    d = ceilr(mulir(powii(gen_2,addii(addii(fprep.p,gen_3),t)),rdivii(fprep.d,e,DEFAULTPREC)));
    k = divii(fprep.k,t);
    if (fprep.f == NULL)
    {
        f = NULL;
        gerepileall(ltop,2,&d,&k);
    }
    else
    {
        f = addrr(fprep.f,rdivii(powii(gen_3,gen_2),powii(gen_2,gen_3),DEFAULTPREC));
        gerepileall(ltop,3,&d,&k,&f);
    }
    fprep_.f = f;
    fprep_.p = fprep.p;
    fprep_.b = gcopy(fprep.b);
    fprep_.d = d;
    fprep_.k = k;
    return fprep_;
}

struct fp_rep 
numult (GEN O, struct fp_rep fprep1, struct fp_rep fprep2)
{
    if (cmpii(fprep1.p,fprep2.p)) pari_err_DOMAIN("numult",itostr(fprep2.p),"!=",fprep1.p,fprep2.p);
    struct fp_rep fprep;
    fprep.p = fprep1.p;
    pari_sp ltop = avma;
    GEN b, e, h, s, T, a, b;
    b = inucomp(O,fprep1.b,fprep2.b,0);
    if (cmpii(mulii(fprep1.d,fprep2.d),powii(gen_2,addii(mulii(gen_2,fprep.p),gen_1))) <= 0)
    {
        e = floorr(rdivii(mulii(fprep1.d,fprep2.d),powii(gen_2,fprep.p),DEFAULTPREC));
        h = addii(fprep1.k,fprep2.k);
    }
    else
    {
        e = floorr(rdivii(mulii(fprep1.d,fprep2.d),powii(gen_2,addii(fprep.p,gen_1)),DEFAULTPREC));
        h = addii(addii(fprep1.k,fprep2.k),gen_1);
    }
    s = gmax(gen_0,ceilr(subrr(addir(addii(fprep.p,addii(gen_2,gen_2)),divrr(mplog(itor(gmael(b,2,2),DEFAULTPREC)),mplog2(DEFAULTPREC))),divrr(mplog(itor(gmael3(b,1,1,1),DEFAULTPREC)),mplog2(DEFAULTPREC)))));
    T = addii(mulii(powii(gen_2,s),gmael(b,2,1)),mulii(gmael(b,2,2),floorr(mulir(powii(gen_2,s),gsqrt(gel(O,1),DEFAULTPREC)))));
    a = gcopy(gmael(b,2,1));
    b = negi(gmael(b,2,2));
    fprep.b = gel(b,1);
    fprep.d = e;
    fprep.k = h;
    fprep = fpremove(fprep,T,gmael(b,2,3),s);
    if (fprep1.f != NULL && fprep2.f != NULL)
    {
        fprep.f = addrr(addrr(fprep1.f,fprep2.f),mulir(powii(gen_2,negi(fprep.p)),mulrr(fprep1.f,fprep2.f)));
        gerepileall(ltop,4,&fprep.f,&fprep.b,&fprep.d,&fprep.k);
    }
    else
    {
        fprep.f = NULL;
        gerepileall(ltop,3,&fprep.b,&fprep.d,&fprep.k);
    }
    return fprep;
}