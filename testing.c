#include <stdio.h>
#include <pari/pari.h>
#include "real_quadratic_orders.h"
#include "testing.h"

int main() {
   pari_init(1000000,0);
   long j = 10, n = 100;
   printf("regulatorshanks: %ld\n",timedtest(regulatorshanks,j,n));
   printf("regulatorcf: %ld\n",timedtest(regulatorcf,j,n));
   pari_close();
   return 0;
}

long timedtest(GEN (*f)(GEN O, long prec, long flag), int i, int n)
{
   pari_timer timer;
   pari_sp ltop = avma, av;
   GEN x, tmp;
   int j; 
   long t;
   tmp = powis(strtoi("10"),i);
   timer_start(&timer);
   av = avma;
   for (j = 0; j < n; j++)
   {
      x = gerepileupto(av,addii(randomi(tmp),tmp));
      if (Z_issquare(x)) continue;
      f(rqoinit(x),DEFAULTPREC,0);
      set_avma(av);
   }
   t = timer_delay(&timer);
   set_avma(ltop);
   return t;
}
