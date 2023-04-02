#include "fp_representations.h"

struct fp_rep 
fprepinit (GEN f, GEN p, GEN b, GEN d, GEN k)
{
    struct fp_rep a;
    pari_sp av;
    if (f != NULL && typ(f) != t_REAL && typ(f) != t_FRAC && typ(f) != t_INT && typ(f) != t_QUAD) pari_err_TYPE("fprepinit",f);
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