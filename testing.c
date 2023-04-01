#include <stdio.h>
#include <pari/pari.h>
#include "real_quadratic_orders.h"

int main() {
   pari_init(1000000,0);
   GEN O, d;
   printf("D = "); d = gp_read_stream(stdin);
   O = rqoinit(d);
   pari_printf("%Ps\n", O);
   GEN Q, P, Q_, P_, I, J, K;
   printf("Q = "); Q = gp_read_stream(stdin);
   printf("P = "); P = gp_read_stream(stdin);
   I = rqiinit(gen_1, Q, P);
   printf("Q_ = "); Q_ = gp_read_stream(stdin);
   printf("P_ = "); P_ = gp_read_stream(stdin);
   J = rqiinit(gen_1, Q_, P_);
   pari_printf("%Ps\n", I);
   pari_printf("%Ps\n", J);
   K = inucomp(O,I,J,0);
   pari_printf("%Ps\n", K);
   pari_close();
   return 0;
}
