#include "real_quadratic_orders.h"

GEN
rqoinit(GEN d)
{
    if (typ(d) != t_INT) pari_err_TYPE("rqoinit",d);
    if (signe(d) < 0) pari_err_DOMAIN("rqoinit","sign(d)","<",gen_0,d);
    if (Z_issquare(d)) pari_err_DOMAIN("rqoinit","issquare(d)","=",gen_1,d);
    GEN c, res = cgetg(5,t_VEC);
    pari_sp ltop;
    long r;
    gel(res,1) = gcopy(d);
    c = core2(gel(res,1));
    gel(res,2) = c;
    r = mod4(gel(c,1)); 
    if (r == 1)
    {
        gel(res,3) = gen_2;
        gel(res,4) = gel(res,1); 
    }
    else
    {
        gel(res,3) = gen_1;
        ltop = avma; gel(res,4) = gerepileupto(ltop,mulii(gel(res,1),sqri(gen_2)));
    }
    return res;
}

GEN
rqiinit(GEN S, GEN Q, GEN P)
{
    if (typ(S) != t_INT) pari_err_TYPE("rqiinit0",S);
    if (typ(Q) != t_INT) pari_err_TYPE("rqiinit0",Q);
    if (typ(P) != t_INT) pari_err_TYPE("rqiinit0",P);
    GEN res, res2;
    res = cgetg(3,t_VEC);
    res2 = cgetg(3,t_VEC);
    gel(res2,1) = gcopy(Q);
    gel(res2,2) = gcopy(P);
    gel(res,1) = gcopy(S);
    gel(res,2) = res2;
    return res;
}

void
checkrqi(GEN O, GEN a, const char *f)
{
    pari_sp av = avma;
    GEN d;
    if (signe(gmael(a,2,1)) <= 0) pari_err_DOMAIN(f,itostr(gmael(a,2,1)),"<=",gen_0,a); //could be fixed instead of throwing error
    d = gcdii(gmael(a,2,1),gel(O,3));
    if (cmpii(d,gel(O,3))) pari_err_DOMAIN(f,itostr(gel(O,3)),"does not divide",gel(gel(a, 2),1),a);
    d = gcdii(mulii(gel(O,3),gmael(a,2,1)),subii(gel(O,1),sqri(gmael(a,2,2))));
    if (gcmp(mulii(gel(O,3),gmael(a,2,1)),d)) pari_err_DOMAIN(f,itostr(mulii(gel(O,3),gmael(a,2,1))),"does not divide",subii(gel(O,1),sqri(gmael(a,2,2))),a);
    set_avma(av); return;
}

GEN 
inucomp(GEN O, GEN a, GEN b, long flag) 
{
    //Which GENs can be turned into longs?
    //Control output of (A, B, C) via flag?
    if (flag)
    {
        checkrqi(O,a,"inucomp");
        checkrqi(O,b,"inucomp");
        if (cmpii(gen_1,gel(a,1))) pari_err_DOMAIN("inucomp","S","!=",gen_1,a);
        if (cmpii(gen_1,gel(b,1))) pari_err_DOMAIN("inucomp","S","!=",gen_1,b);
    }
    GEN swap, tmp1, tmp2, G, X, S, Y, Z, R_, R_1, R_0, C_0 = gen_0, C_1 = gen_m1, Q, P, B_0, B_1, q, M_1, M_2, Q_, k, P_, Q_old, Q_old_, A, B, C, res = cgetg(3,t_VEC);
    pari_sp ltop = avma, lbot, av, av2;
    long i = 0, j = 1;
    if (cmpii(gmael(a,2,1),gmael(b,2,1)) < 0)
    {
        swap = a;
        a = b;
        b = swap;
    }
    tmp1 = diviiexact(gmael(a,2,1),gel(O,3));
    G = bezout(tmp1,diviiexact(gmael(b,2,1),gel(O,3)),&tmp2,&X);
    if (cmpii(X,gen_0) < 0) X = addii(X,tmp1); 
    S = bezout(diviiexact(addii(gmael(a,2,2),gmael(b,2,2)),gel(O,3)),G,&Y,&Z);
    av = avma; R_ = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(gmael(b,2,2))),gmael(b,2,1)));
    av = avma; R_1 = gerepileupto(av,lift(gadd(gmul(X,gmul(Z,gsub(mkintmod(gmael(a,2,2),tmp1),gmael(b,2,2)))),gmul(Y,mkintmod(R_,tmp1)))));
    gerepileall(ltop,3,&R_1,&R_,&S);
    av2 = avma;  
    R_0 = diviiexact(gmael(a,2,1),S);
    av = avma; tmp1 = gerepileupto(av,floorr(gmul(gsqrt(gmul(gen_2,gel(O,3)),DEFAULTPREC),gsqrtn(gel(O,1),sqri(gen_2),NULL,DEFAULTPREC))));
    if (cmpii(R_0,tmp1) < 0)
    {
        j = 0;
        av = avma; Q = gerepileupto(av,diviiexact(mulii(gmael(a,2,1),gmael(b,2,1)),mulii(gel(O,3),sqri(S))));
        av = avma; P = gerepileupto(av,lift(gadd(mkintmod(gmael(b,2,2),Q),diviiexact(mulii(R_1,gmael(b,2,1)),mulii(gel(O,3),S)))));
        B_0 = gen_1;
        B_1 = gen_0;
    }
    else
    {
        while (cmpii(R_1,tmp1) > 0)
        {
            (i == 1) ? (i = 0) : (i = 1);
            q = gfloor(gdiv(R_0,R_1));
            swap = C_1;
            C_1 = subii(C_0,mulii(q,C_1));
            C_0 = swap;
            swap = R_1;
            R_1 = subii(R_0,mulii(q,R_1));
            R_0 = swap;
            gerepileall(av2,4,&R_0,&R_1,&C_0,&C_1);
        }
        av = avma; M_1 = gerepileupto(av,diviiexact(addii(mulii(divii(gmael(b,2,1),mulii(gel(O,3),S)),R_1),mulii(subii(gmael(a,2,2),gmael(b,2,2)),C_1)),divii(gmael(a,2,1),S)));
        av = avma; M_2 = gerepileupto(av,diviiexact(addii(mulii(addii(gmael(a,2,2),gmael(b,2,2)),R_1),mulii(mulii(mulii(gel(O,3),S),R_),C_1)),divii(gmael(a,2,1),S)));
        av = avma; Q = (i == 1) ? gerepileupto(av,subii(mulii(R_1,M_1),mulii(C_1,M_1))) : gerepileupto(av,negi(subii(mulii(R_1,M_1),mulii(C_1,M_2))));
        av = avma; P = gerepileupto(av,subii(diviiexact(addii(mulii(diviiexact(gmael(b,2,1),mulii(gel(O,3),S)),R_1),mulii(Q,C_0)),C_1),gmael(b,2,2)));
    }
    Q_ = absi(Q);
    av = avma; k = gerepileupto(av,gfloor(gdiv(subii(sqrti(gel(O,1)),P),Q_)));
    av = avma; P_ = gerepileupto(av,addii(mulii(k,Q_),P));
    if (j != 0)
    {
        B_0 = absi(C_0);
        av = avma; B_1 = gerepileupto(av, mulsi(signe(Q),absi(C_1)));
    }
    av = avma; tmp1 = gerepileupto(av,addii(P_,sqrti(gel(O,1))));
    if (cmpii(tmp1,Q_) < 0)
    {
        lbot = avma;
        av = avma; A = gerepileupto(av,mulii(S,addii(mulii(Q,B_0),mulii(P,B_1))));
        av = avma; B = gerepileupto(av,mulii(negi(S),B_1));
        C = gcopy(Q);
        return res = gerepile(ltop,lbot,mkvec2(rqiinit(gen_1,Q_,P_),mkvec3(A,B,C)));

        q = gfloor(gdiv(tmp1,Q_));
        Q_old = Q;
        av = avma; P = gerepileupto(av,subii(mulii(q, Q_),P));
        av = avma; Q = gerepileupto(av,diviiexact(gsub(gel(O,1),gsqr(P)),Q_));
        Q_ = absi(Q);
        av = avma; k = gerepileupto(av,gfloor(gdiv(subii(sqrti(gel(O,1)),P),Q_)));
        av = avma; P_ = gerepileupto(av,addii(mulii(k,Q_),P));
        swap = B_1;
        av = avma; B_1 = gerepileupto(av,addii(mulii(q,B_1),B_0));
        B_0 = swap;
        if (cmpii(tmp1,Q_) < 0)
        {
            Q_old_ = Q;
            av = avma; Q_ = gerepileupto(av,addii(subii(Q_old, Q_old_),mulii(gen_2,P)));
            Q = Q_;
            P = subii(Q_old_,P);
            P_ = subii(P,Q);
            swap = B_1;
            B_1 = addii(B_1, B_0);
            B_0 = swap;
        }
    }
    lbot = avma;
    av = avma; A = gerepileupto(av,mulii(S,addii(mulii(Q,B_0),mulii(P,B_1))));
    av = avma; B = gerepileupto(av,mulii(negi(S),B_1));
    C = gcopy(Q);
    return res = gerepile(ltop,lbot,mkvec2(rqiinit(gen_1,Q_,P_),mkvec3(A,B,C)));
}

GEN 
regulatorcf(GEN O, long prec, long flag)
{
    pari_sp ltop = avma, av, av2;
    GEN sqrtd_, sqrtd, a, Q_, Q, Q_0, P, q, B_0, B_1, G_0, G_1, swap, psi, theta;
    long n = 0;
    sqrtd_ = sqrti(gel(O,1));
    sqrtd = gsqrt(gel(O,1),prec);
    av = avma; a = gerepileupto(av,rqiinit(gen_1,gel(O,3),mulii(gmael(O,2,2),subii(gel(O,3),gen_1)))); //this is the identity in O
    av2 = avma;
    Q_0 = gmael(a,2,1);
    Q = Q_0;
    P = gmael(a,2,2);
    av = avma; q = gerepileupto(av,divii(addii(P,sqrtd_),Q));
    B_0 = gen_0;
    B_1 = gen_1;
    G_0 = Q;
    av = avma; G_1 = gerepileupto(av,subii(mulii(Q,q),P));
    if (flag == 2) theta = itor(gen_1,prec);
    do
    {
        n += 1;
        if (flag == 1) pari_printf("%ld: (%Ps, %Ps)\n",n,Q,P);
        else if (flag == 2)
        {
            av = avma; pari_printf("%ld: (%Ps, %Ps) %Ps\n",n,Q,P,gerepileupto(av,mplog(theta)));
        }
        Q_ = Q;
        av = avma; P = gerepileupto(av,subii(mulii(q,Q),P));
        av = avma; Q = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(P)),Q));
        av = avma; q = gerepileupto(av,divii(addii(P,sqrtd_),Q));
        swap = G_1;
        av = avma; G_1 = gerepileupto(av,addii(mulii(q,G_1),G_0));
        G_0 = swap;
        swap = B_1;
        av = avma; B_1 = gerepileupto(av,addii(mulii(q,B_1),B_0));
        B_0 = swap;
        if (flag == 2)
        {
            av = avma; psi = gerepileupto(av,divri(addir(P,sqrtd),Q_));
            av = avma; theta = gerepileupto(av,mulrr(theta,psi));
            gerepileall(av2,8,&q,&Q,&P,&G_0,&G_1,&B_0,&B_1,&theta);
        }
        else gerepileall(av2,7,&q,&Q,&P,&G_0,&G_1,&B_0,&B_1);
    } while (cmpii(Q,Q_0));
    n += 1; 
    if (flag == 1) pari_printf("%ld: (%Ps, %Ps)\n",n,Q,P);
    else if (flag == 2) 
    {
        av = avma; pari_printf("%ld: (%Ps, %Ps) %Ps\n",n,Q,P,gerepileupto(av,mplog(theta)));
    }
    return gerepileupto(ltop,mplog(divri(addir(G_0,mulri(sqrtd,B_0)),Q_0)));
}