#include "utility.h"

long
sigbits (GEN x)
{
    if (typ(x) != t_INT) pari_err_TYPE("sigbits",x);
    unsigned long i = 0, x_;
    x_ = *int_MSW(x);
    while (x_ >>= 1) i++;
    return (lgefint(x)-3)*BITS_IN_LONG+i+1;
}