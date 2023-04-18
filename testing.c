#include <stdio.h>
#include "testing.h"
#include"utility.h"

int main() {
   pari_init(1000000000,0);
   pari_close();
   return 0;
}

ulong 
timedtest(GEN (*f)(GEN O, ulong prec, long flag), ulong i, ulong n)
{
   pari_timer timer;
   pari_sp ltop = avma, av;
   GEN x, tmp; 
   ulong j;
   tmp = powis(stoi(10),i);
   timer_start(&timer);
   av = avma;
   for (j = 0; j < n; j++)
   {
      x = gerepileupto(av,addii(randomi(tmp),tmp));
      if (Z_issquare(x)) continue;
      f(rqoinit(x),DEFAULTPREC,0);
      set_avma(av);
   }
   return gc_long(ltop,timer_delay(&timer));
}

ulong 
testcr(ulong i, ulong n, GEN m)
{
   pari_sp ltop = avma, av, av2;
   GEN O, tmp, x, A, B, C, C_, y;
   ulong j, count = 0, size, size2;
   tmp = powis(stoi(10),i);
   av2 = avma;
   for (j = 0; j < n; j++)
   {
      av = avma; x = gerepileupto(av,addii(randomi(tmp),tmp));
      if (Z_issquare(x)) continue;
      O = rqoinit(x);
      size = avma;
      A = regulatorcf(O,DEFAULTPREC,0);
      size -= avma;
      if (cmpri(gel(A,1),stoi(10)) < 0) continue; // in this case cr might fail
      av = avma; y = gerepileupto(av,roundr(gdiv(gel(A,1),mplog2(DEFAULTPREC))));
      size2 = avma;
      B = cr(O,pci(O),y,ghalf,m);
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
      set_avma(av2);
   }
   set_avma(ltop);
   return count;
}

ulong 
testcrsmoothpart(ulong i, ulong n, ulong B, ulong c, ulong sw)
{
   pari_timer timer;
   pari_sp ltop = avma, av, av2;
   GEN O, R, x, tmp, y, prod = gen_1, b;
   ulong j, p = 2;
   forprime_t S;
   u_forprime_init(&S,2,B);
   if (sw == 1)
   {
      while ((p = u_forprime_next(&S))) prod = gerepileupto(ltop,mulis(prod,p));
   }
   tmp = powis(stoi(10),i);
   timer_start(&timer);
   av2 = avma;
   for (j = 0; j < n; j++)
   {
      av = avma; x = gerepileupto(av,addii(randomi(tmp),tmp));
      if (Z_issquare(x)) continue;
      O = rqoinit(x);
      b = pci(O);
      if (Mod4(x) != 1) x = mulsi(4,x);
      if (i >= 10) R = gel(quadclassunit0(x,0,NULL,DEFAULTPREC),4);
      else R = quadregulator(x,DEFAULTPREC);
      if (cmpri(R,stoi(10)) < 0) continue; // in this case cr might fail
      av = avma; y = gerepileupto(av,roundr(gdiv(R,mplog2(DEFAULTPREC))));
      if (sw == 1) crsmoothpart(O,gel(cr(O,b,y,ghalf,prod),2),B,c);
      else if (sw == 2) crsmoothpart2(O,b,y,ghalf,B,c);
      else crsmoothpart_alt(O,gel(cr(O,b,y,ghalf,gen_0),2),B,c);
      set_avma(av2);
   }
   return gc_long(ltop,timer_delay(&timer));
}