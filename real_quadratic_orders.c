#include "real_quadratic_orders.h"
#include "utility.h"

GEN
rqoinit(GEN d)
{
    if (typ(d) != t_INT) pari_err_TYPE("rqoinit",d);
    if (signe(d) < 0) pari_err_DOMAIN("rqoinit","sign(d)","<",gen_0,d);
    if (Z_issquare(d)) pari_err_DOMAIN("rqoinit","issquare(d)","=",gen_1,d);
    GEN c, res = cgetg(5,t_VEC);
    pari_sp ltop;
    ulong r;
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
    if (typ(S) != t_INT) pari_err_TYPE("rqiinit",S);
    if (typ(Q) != t_INT) pari_err_TYPE("rqiinit",Q);
    if (typ(P) != t_INT) pari_err_TYPE("rqiinit",P);
    GEN res;
    res = cgetg(3,t_VEC);
    gel(res,1) =  (S == gen_1) ? gen_1 : gcopy(S); //common special case, no need to store S
    gel(res,2) = mkvec2copy(Q,P);
    return res;
}

void
checkrqi(GEN O, GEN a, const char *f)
{
    pari_sp ltop = avma;
    GEN d;
    if (signe(gmael(a,2,1)) <= 0) pari_err_DOMAIN(f,itostr(gmael(a,2,1)),"<=",gen_0,a); //could be fixed instead of throwing error
    d = gcdii(gmael(a,2,1),gel(O,3));
    if (cmpii(d,gel(O,3))) pari_err_DOMAIN(f,itostr(gel(O,3)),"does not divide",gel(gel(a, 2),1),a);
    d = gcdii(mulii(gel(O,3),gmael(a,2,1)),subii(gel(O,1),sqri(gmael(a,2,2))));
    if (cmpii(mulii(gel(O,3),gmael(a,2,1)),d)) pari_err_DOMAIN(f,itostr(mulii(gel(O,3),gmael(a,2,1))),"does not divide",subii(gel(O,1),sqri(gmael(a,2,2))),a);
    set_avma(ltop); return;
}

GEN
pci(GEN O)
{
    pari_sp ltop = avma;
    GEN sqrtd_;
    sqrtd_ = sqrti(gel(O,1));
    if (equali1(gel(O,3))) return gerepileupto(ltop,rqiinit(gen_1,gel(O,3),sqrtd_));
    else 
    {
        if (mod2(gmael(O,2,2)))
        {
            if (mod2(sqrtd_)) return gerepileupto(ltop,rqiinit(gen_1,gel(O,3),sqrtd_));
            else return gerepileupto(ltop,rqiinit(gen_1,gel(O,3),subii(sqrtd_,gen_1)));
        }
        else
        {
            if (mod2(sqrtd_)) return gerepileupto(ltop,rqiinit(gen_1,gel(O,3),subii(sqrtd_,gen_1)));
            else return gerepileupto(ltop,rqiinit(gen_1,gel(O,3),sqrtd_));
        }
    }
}

GEN 
imultiply(GEN O, GEN a, GEN b) // DEPRECATED
{
    checkrqi(O,a,"imultiply");
    checkrqi(O,b,"imultiply");
    if (cmpii(gen_1,gel(a,1))) pari_err_DOMAIN("imultiply","S","!=",gen_1,a);
    if (cmpii(gen_1,gel(b,1))) pari_err_DOMAIN("imultiply","S","!=",gen_1,b);
    pari_sp ltop = avma, av;
    GEN S, Q, P, U, V, W, X, Y, tmp;
    S = bezout(bezout(diviiexact(gmael(a,2,1),gel(O,3)),diviiexact(gmael(b,2,1),gel(O,3)),&V,&W),diviiexact(addii(gmael(a,2,2),gmael(b,2,2)),gel(O,3)),&X,&Y);
    V = mulii(X,V);
    W = mulii(X,W);
    av = avma; Q = gerepileupto(av,diviiexact(mulii(gmael(a,2,1),gmael(b,2,1)),mulii(sqri(S),gel(O,3))));
    tmp = diviiexact(gmael(a,2,1),S);
    av = avma; U = gerepileupto(av,lift(gadd(gmul(mkintmod(modii(W,tmp),tmp),subii(gmael(a,2,2),gmael(b,2,2))),gmul(mkintmod(modii(Y,tmp),tmp),diviiexact(subii(gel(O,1),sqri(gmael(b,2,2))),gmael(b,2,1))))));
    av = avma; P = gerepileupto(av,lift(gadd(gmael(b,2,2),gmul(mkintmod(U,Q),diviiexact(gmael(b,2,1),mulii(gel(O,3),S))))));
    return gerepileupto(ltop,rqiinit(S,Q,P));
}

GEN 
qiimultiply(GEN O, GEN qi, GEN i) // DEPRECATED
{
    checkrqi(O,i,"qiimultiply");
    if (cmpii(gen_1,gel(i,1))) pari_err_DOMAIN("qiimultiply","S","!=",gen_1,i);
    pari_sp ltop = avma, av;
    GEN a_1, a_2, b_1, b_2, p, q, a, b, c, tmp1, tmp2;
    av = avma; a_1 = gerepileupto(av,diviiexact(mulii(gel(qi,1),gmael(i,2,1)),mulii(gel(O,3),gel(qi,3))));
    av = avma; b_1 = gerepileupto(av,diviiexact(mulii(gel(qi,2),gmael(i,2,1)),mulii(gel(O,3),gel(qi,3))));
    av = avma; a_2 = gerepileupto(av,diviiexact(addii(mulii(gel(qi,1),gmael(i,2,2)),mulii(gel(qi,2),gel(O,1))),mulii(gel(O,3),gel(qi,3))));
    av = avma; b_2 = gerepileupto(av,diviiexact(addii(mulii(gel(qi,2),gmael(i,2,2)),gel(qi,1)),mulii(gel(O,3),gel(qi,3))));
    c = bezout(b_1,b_2,&p,&q);
    av = avma, a = gerepileupto(av,diviiexact(subii(mulii(a_1,b_2),mulii(b_1,a_2)),c));
    av = avma, b = gerepileupto(av,addii(mulii(p,a_1),mulii(q,a_2)));
    if(!equalii(c,gel(O,3))) pari_err_BUG("qiimultiply");
    tmp1 = diviiexact(b,gel(O,3));
    tmp2 = absi(diviiexact(a,gel(O,3)));
    while (cmpii(tmp1,gen_0) < 0)
    {
        av = avma; tmp1 = gerepileupto(av,addii(tmp1,tmp2));
    }
    return gerepileupto(ltop,rqiinit(gen_1,tmp2,diviiexact(b,gel(O,3))));
}

GEN 
inucomp(GEN O, GEN a, GEN b, long flag) 
{ //Which GENs can be turned into longs? Arithmetic on longs is about 5 times faster than on t_INTs
    if (flag)
    {
        checkrqi(O,a,"inucomp");
        checkrqi(O,b,"inucomp");
        if (cmpii(gen_1,gel(a,1))) pari_err_DOMAIN("inucomp","S","!=",gen_1,a);
        if (cmpii(gen_1,gel(b,1))) pari_err_DOMAIN("inucomp","S","!=",gen_1,b);
    }
    GEN swap, tmp1, tmp2, G, X, S, Y, Z, R_, R_1, R_0, C_0 = gen_0, C_1 = gen_m1, Q, P, B_0, B_1, q, M_1, M_2, Q_, k, P_, Q_old, Q_old_, res, resi;
    pari_sp ltop = avma, av, av2;
    ulong i = 0;
    if (cmpii(gmael(a,2,1),gmael(b,2,1)) < 0)
    {
        swap = a;
        a = b;
        b = swap;
    }
    G = bezout(diviiexact(gmael(a,2,1),gel(O,3)),diviiexact(gmael(b,2,1),gel(O,3)),&tmp2,&X);
    if (!cmpii(gmael(a,2,1),gmael(b,2,1))) X = gen_0;
    else if (cmpii(X,gen_0) < 0) X = addii(X,diviiexact(gmael(a,2,1),gel(O,3)));
    av2 = avma;
    S = bezout(diviiexact(addii(gmael(a,2,2),gmael(b,2,2)),gel(O,3)),G,&Y,&Z);
    av = avma; R_ = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(gmael(b,2,2))),gmael(b,2,1)));
    R_0 = diviiexact(gmael(a,2,1),S);
    av = avma; R_1 = gerepileupto(av,lift(gadd(gmul(X,gmul(Z,gsub(mkintmod(modii(gmael(a,2,2),R_0),R_0),gmael(b,2,2)))),gmul(Y,mkintmod(modii(R_,R_0),R_0)))));
    gerepileall(ltop,4,&R_0,&R_1,&R_,&S);
    av = avma; tmp1 = gerepileupto(av,floor_safe(gsqrtn(mulii(mulsi(4,sqri(gel(O,3))),gel(O,1)),stoi(4),NULL,maxss(DEFAULTPREC,(sigbits(gel(O,1))/(4*BITS_IN_LONG))+3)))); //flooring might cause precision loss
    if (cmpii(R_0,tmp1) < 0)
    {
        av = avma; Q = gerepileupto(av,diviiexact(mulii(gmael(a,2,1),gmael(b,2,1)),mulii(gel(O,3),sqri(S))));
        av = avma; P = gerepileupto(av,lift(gadd(mkintmod(modii(gmael(b,2,2),Q),Q),diviiexact(mulii(R_1,gmael(b,2,1)),mulii(gel(O,3),S)))));
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
            gerepileall(av2,6,&R_0,&R_1,&C_0,&C_1,&tmp1,&S);
        }
        av = avma; M_1 = gerepileupto(av,diviiexact(addii(mulii(divii(gmael(b,2,1),mulii(gel(O,3),S)),R_1),mulii(subii(gmael(a,2,2),gmael(b,2,2)),C_1)),divii(gmael(a,2,1),S)));
        av = avma; M_2 = gerepileupto(av,diviiexact(addii(mulii(addii(gmael(a,2,2),gmael(b,2,2)),R_1),mulii(mulii(mulii(gel(O,3),S),R_),C_1)),divii(gmael(a,2,1),S)));
        av = avma; Q = (i == 1) ? gerepileupto(av,subii(mulii(R_1,M_1),mulii(C_1,M_2))) : gerepileupto(av,negi(subii(mulii(R_1,M_1),mulii(C_1,M_2))));
        av = avma; P = gerepileupto(av,subii(diviiexact(addii(mulii(diviiexact(gmael(b,2,1),mulii(gel(O,3),S)),R_1),mulii(Q,C_0)),C_1),gmael(b,2,2)));
        av = avma; B_0 = gerepileupto(av,absi(C_0));
        av = avma; B_1 = gerepileupto(av,absi(C_1)); 
    }
    Q_ = absi(Q);
    av = avma; k = gerepileupto(av,gfloor(gdiv(subii(sqrti(gel(O,1)),P),Q_)));
    av = avma; P_ = gerepileupto(av,addii(mulii(k,Q_),P));
    if (cmpii(addii(P_,sqrti(gel(O,1))),Q_) < 0)
    {
        av = avma; q = gerepileupto(av,gfloor(gdiv(addii(P,sqrti(gel(O,1))),Q_)));
        Q_old = Q;
        av = avma; P = gerepileupto(av,subii(mulii(q, Q_),P));
        av = avma; Q = gerepileupto(av,diviiexact(gsub(gel(O,1),gsqr(P)),Q_));
        Q_ = absi(Q);
        av = avma; k = gerepileupto(av,gfloor(gdiv(subii(sqrti(gel(O,1)),P),Q_)));
        av = avma; P_ = gerepileupto(av,addii(mulii(k,Q_),P));
        swap = B_1;
        av = avma; B_1 = gerepileupto(av,addii(mulii(q,B_1),mulsi(signe(Q_old),B_0))); 
        B_0 = swap;
        if (cmpii(addii(P_,sqrti(gel(O,1))),Q_) < 0)
        {  
            Q_old_ = Q;
            av = avma; Q_ = gerepileupto(av,addii(subii(absi(Q_old),Q),mulii(gen_2,P)));
            Q = Q_;
            P = subii(Q_old_,P);
            P_ = subii(P,Q);
            swap = B_1;
            B_1 = addii(B_1,B_0); 
            B_0 = swap;
        }
    }
    res = cgetg(3,t_VEC);
    gel(res,1) = rqiinit(gen_1,Q_,P_);
    resi = cgetg(4,t_VEC);
    av = avma; gel(resi,1) = gerepileupto(av,mulii(S,addii(mulii(Q,B_0),mulii(P,B_1))));
    av = avma; gel(resi,2) = gerepileupto(av,mulii(negi(S),B_1));
    gel(resi,3) = gcopy(Q);
    gel(res,2) = resi;
    return gerepileupto(ltop,res);
}

GEN 
regulatorcf(GEN O, ulong prec, long flag)
{
    pari_sp ltop = avma, av, av2;
    GEN sqrtd_, sqrtd, a, Q_, Q, Q_0, P, q, B_0, B_1, G_0, G_1, swap, psi, theta, res, resi;
    ulong n = 0;
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
        if (flag == 1) pari_printf("%ld: [%Ps, %Ps]\n",n,Q,P);
        else if (flag == 2)
        {
            av = avma; pari_printf("%ld: [%Ps, %Ps] %Ps\n",n,Q,P,gerepileupto(av,mplog(theta)));
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
    if (flag == 1) pari_printf("%ld: [%Ps, %Ps]\n",n,Q,P);
    else if (flag == 2) 
    {
        av = avma; pari_printf("%ld: [%Ps, %Ps] %Ps\n",n,Q,P,gerepileupto(av,mplog(theta)));
    }
    res = cgetg(3,t_VEC);
    av = avma; gel(res,1) = gerepileupto(av,mplog(divri(addir(G_0,mulri(sqrtd,B_0)),Q_0)));
    resi = cgetg(4,t_VEC);
    av = avma; gel(resi,1) = gerepileupto(av,diviiexact(mulii(gen_2,G_0),Q_0));
    av = avma; gel(resi,2) = gerepileupto(av,diviiexact(mulii(gen_2,B_0),Q_0));
    gel(resi,3) = (n%2) ? gen_1 : gen_m1;
    gel(res,2) = resi;
    return gerepileupto(ltop,res);
}

GEN 
regulatorshanks(GEN O, ulong prec, long flag)
{ // turn L into hashtable
    pari_sp ltop = avma, av, av2;
    GEN sqrtd_, sqrtd, sqrt4d, a, b, L, Q_, Q_0, Q, P_, P, q, B_0, B_1, G_0, G_1, swap, psi, theta, lt, el, el2;
    ulong n = 0, i;
    sqrtd_ = sqrti(gel(O,1));
    sqrtd = gsqrt(gel(O,1),prec);
    sqrt4d = sqrtr(sqrtd);
    av = avma; a = gerepileupto(av,rqiinit(gen_1,gel(O,3),mulii(gmael(O,2,2),subii(gel(O,3),gen_1)))); //this is the identity in O
    el = cgetg(3,t_VEC);
    el2 = cgetg(3,t_VEC);
    gel(el,1) = el2;
    av2 = avma;
    L = cgetg(1,t_VEC);
    Q_0 = gmael(a,2,1);
    Q = Q_0;
    P = gmael(a,2,2);
    av = avma; q = gerepileupto(av,divii(addii(P,sqrtd_),Q));
    B_0 = gen_0;
    B_1 = gen_1;
    G_0 = Q;
    av = avma; G_1 = gerepileupto(av,subii(mulii(Q,q),P));
    theta = itor(gen_1,prec);
    lt = mplog(theta);
    do // baby steps
    {
        gel(el2,1) = Q;
        gel(el2,2) = P;
        gel(el,2) = lt;
        av = avma; L = gerepileupto(av,vec_append(L,el));
        n += 1;
        Q_ = Q; P_ = P;
        av = avma; P = gerepileupto(av,subii(mulii(q,Q),P));
        av = avma; Q = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(P)),Q));
        av = avma; q = gerepileupto(av,divii(addii(P,sqrtd_),Q));
        swap = G_1;
        av = avma; G_1 = gerepileupto(av,addii(mulii(q,G_1),G_0));
        G_0 = swap;
        swap = B_1;
        av = avma; B_1 = gerepileupto(av,addii(mulii(q,B_1),B_0));
        B_0 = swap;
        av = avma; psi = gerepileupto(av,divri(addir(P,sqrtd),Q_));
        av = avma; theta = gerepileupto(av,mulrr(theta,psi));
        lt = mplog(theta);
        if (equalii(Q,Q_0)) return gerepilecopy(ltop,lt);
        if (equalii(P_,P)) return gerepileupto(ltop,addrr(mulir(gen_2,subrr(lt,mplog(psi))),mplog(rdivii(Q_0,Q_,prec))));
        if (equalii(Q_,Q)) return gerepileupto(ltop,addrr(mulir(gen_2,subrr(lt,mplog(psi))),mplog(mulir(Q_0,divri(psi,Q_)))));
        gerepileall(av2,10,&q,&Q,&P,&G_0,&G_1,&B_0,&B_1,&theta,&lt,&L);
    } while (cmprr(lt,sqrt4d) <= 0);
    b = mkvec3(rqiinit(gen_1,Q,P),lt,theta);
    for (i = 0; i < 3; i++) // 3 more baby steps
    { // code duplication...
        gel(el2,1) = Q;
        gel(el2,2) = P;
        gel(el,2) = lt;
        av = avma; L = gerepileupto(av,vec_append(L,el));
        n += 1;
        Q_ = Q; P_ = P;
        av = avma; P = gerepileupto(av,subii(mulii(q,Q),P));
        av = avma; Q = gerepileupto(av,diviiexact(subii(gel(O,1),sqri(P)),Q));
        av = avma; q = gerepileupto(av,divii(addii(P,sqrtd_),Q));
        swap = G_1;
        av = avma; G_1 = gerepileupto(av,addii(mulii(q,G_1),G_0));
        G_0 = swap;
        swap = B_1;
        av = avma; B_1 = gerepileupto(av,addii(mulii(q,B_1),B_0));
        B_0 = swap;
        av = avma; psi = gerepileupto(av,divri(addir(P,sqrtd),Q_));
        av = avma; theta = gerepileupto(av,mulrr(theta,psi));
        lt = mplog(theta);
        if (equalii(Q,Q_0)) return gerepilecopy(ltop, lt);
        if (equalii(P_,P)) return gerepileupto(ltop,addrr(mulir(gen_2,subrr(lt,mplog(psi))),mplog(rdivii(Q_0,Q_,prec))));
        if (equalii(Q_,Q)) return gerepileupto(ltop,addrr(mulir(gen_2,subrr(lt,mplog(psi))),mplog(mulir(Q_0,divri(psi,Q_)))));
        gerepileall(av2,11,&q,&Q,&P,&G_0,&G_1,&B_0,&B_1,&theta,&lt,&L,&b);
    }
    gerepileall(ltop,2,&L,&b);
    sqrtd = gsqrt(gel(O,1),prec);
    lt = gel(b,2);
    av2 = avma;
    if (flag) pari_printf("Baby steps: %Ps\nGiant steps:\n%Ps %Ps\n",L,gmael(b,1,2),gel(b,2));
    a = b;
    do
    {
        a = inucomp(O,gel(a,1),gel(b,1),0);
        av = avma; lt = gerepileupto(av,addrr(addrr(lt,gel(b,2)),mplog(absr(divir(gmael(a,2,3),addir(gmael(a,2,1),mulir(gmael(a,2,2),sqrtd)))))));
        if (flag) pari_printf("%Ps %Ps\n",gmael(a,1,2),lt);
        for (i = 1; i <= n; i++)
        {
            if ((!cmpii(gmael3(a,1,2,1),gmael3(L,i,1,1))) && (!cmpii(gmael3(a,1,2,2),gmael3(L,i,1,2)))) return gerepileupto(ltop,subrr(lt,gmael(L,i,2)));
        }
        gerepileall(av2,2,&lt,&a);
    } while (1);
}