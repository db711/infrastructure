#include <stdio.h>
#include "testing.h"
#include"utility.h"

int main() {
   pari_init(1000000000,0);
   return 0;
}

long 
timedtest(GEN (*f)(GEN O, long prec, long flag), int i, int n)
{
   pari_timer timer;
   pari_sp ltop = avma, av;
   GEN x, tmp; 
   long j, t;
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

long 
testcr(int i, int n)
{
   pari_sp ltop = avma, av, av2;
   GEN O, tmp, x, A, B, C, C_, y;
   long j, count = 0, size, size2;
   tmp = powis(strtoi("10"),i);
   av2 = avma;
   for (j = 0; j < n; j++)
   {
      av = avma; x = gerepileupto(av,addii(randomi(tmp),tmp));
      if (Z_issquare(x)) continue;
      pari_printf("d = %Ps\n",x);
      O = rqoinit(x);
      size = avma;
      A = regulatorcf(O,DEFAULTPREC,0);
      size -= avma;
      if (cmpri(gel(A,1),strtoi("10")) < 0) continue; // in this case cr might fail
      av = avma; y = gerepileupto(av,roundr(gdiv(gel(A,1),mplog2(DEFAULTPREC))));
      size2 = avma;
      B = cr(O,pci(O),y,ghalf);
      size2 -= avma;
      if (size2 < size) count++;
      av = avma; 
      C = expandcr(O,gel(B,2));
      if(cmpii(gel(O,3),gen_2))
      {
         C_ = cgetg(3,t_VEC);
         gel(C_,1) = mulii(gen_2,gel(C,1));
         gel(C_,2) = mulii(gen_2,gel(C,2));
         C = C_; C_ = NULL;
      }
      if (abscmpii(gel(C,1),gmael(A,2,1)) || abscmpii(gel(C,2),gmael(A,2,2))) pari_printf("Error found for d = %Ps\n%Ps\nvs\n%Ps\n\n",x,C,gel(A,2));
      else
      set_avma(av2);
   }
   set_avma(ltop);
   return count;
}