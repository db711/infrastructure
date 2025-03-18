#include "twin_smooths.h"
#include "stormer.h"
#include <pthread.h> 

#define NUM_THREADS 6 // number of threads to run (at once)
#define SECURITY_PARAM 128 // half the bit size we are looking for
#define LOWER_BOUND 250 // lower bound (in bits) for twin smooths
#define UPPER_BOUND 260 // upper bound (in bits) for twin smooths
#define ITER 8 // how many numbers in the Pell sequence are checked at most
#define SMOOTHNESS_BOUND 65536// 2^16
#define STARTING_D "2"
#define NUM_DISC 8192 // number of discriminants per thread
#define OUTPUT_FILE "output.txt"
#define STATUS_FILE "status.txt"

static inline GEN
bvtodisc(GEN bv, GEN start)
{
  pari_sp ltop = avma;
  ulong i;
  for (i = 1; i < lg(bv); i++) if (gel(bv,i)) start = gerepileupto(ltop,mulis(start,uprime(lg(bv)-i)));
  return gerepileupto(ltop,start);
}

static inline GEN
disctobv(GEN d, GEN start, long length)
{
   GEN bv;
   long i;
   pari_sp av, ltop = avma;
   bv = gtovecsmall0(gen_0,length);
   d = diviiexact(d,start);
   for (i = 1; i <= length; i++)
   {
      av = avma;
      if (!cmpii(gen_0,modii(d,prime(i))))
      {
         gel(bv,length-i+1) = 1;
         d = gerepileupto(av,diviiexact(d,prime(i)));
         if (!cmpii(gen_1,d)) break;
      }
      else set_avma(av);
   }
   return gerepileupto(ltop,bv);
}

void *
twin_smooth_range_d_small_bulk(void *arg)
{
  GEN ds, tmp, ret;
  long i;
  ds = pari_thread_start((struct pari_thread*) arg);
  ret = cgetg(1, t_VEC);
  for (i = 1; i < lg(ds); i++)
  {
    tmp = twin_smooth_range_d_small(SMOOTHNESS_BOUND, gel(ds, i), LOWER_BOUND, UPPER_BOUND, ITER);
    if (lg(tmp) > 1) ret = vec_append(ret, mkvec2(gel(ds, i), tmp));
  }
  pari_thread_close();
  return (void*)ret;
}

int
main(void)
{
  long i, j, np;
  pari_init(1048576000,SMOOTHNESS_BOUND+1);
  pthread_t th[NUM_THREADS];
  struct pari_thread pth[NUM_THREADS];

  GEN d_start, ub, bv, stormer, in, out;
  pari_sp ltop, av;
  pari_timer timer;
  FILE* output = NULL;
  FILE* status = NULL;
  d_start = strtoi(STARTING_D);
  av = avma; ub = gerepileupto(av, powis(stoi(2), 2*UPPER_BOUND));
  av = avma; np = itos(primepi(stoi(SMOOTHNESS_BOUND))); set_avma(av);
  av = avma; bv = gerepileupto(av,disctobv(strtoi("254"),d_start,np)); // 254 = 127*2
  stormer = stormer_gen(np, d_start, ub, bv);
  set_avma(avma - SECURITY_PARAM); // kinda hacky

  if (NULL == (output = fopen(OUTPUT_FILE, "w")))
  {
    pari_printf("Error opening output file %s. Aborting.\n", OUTPUT_FILE);
    pari_close();
    return 1;
  }
  if (NULL == (status = fopen(STATUS_FILE, "w")))
  {
    pari_printf("Error opening status file %s. Aborting.\n", STATUS_FILE);
    pari_close();
    return 1;
  }

  timer_start(&timer);
  ltop = avma;
  while (NULL != stormer)
  {
    in = cgetg(NUM_THREADS+1, t_VEC);
    for (i = 1; i < lg(in); i++)
    {
      gel(in, i) = cgetg(NUM_DISC+1, t_VEC);
      for (j = 1; j < lg(gel(in, i)); j++)
      {
        if (NULL == stormer) break;
        gmael2(in, i, j) = gcopy(gel(stormer, 2));
        stormer = stormer_next(stormer, np, ub);
      }
    }
    pari_fprintf(status, "%d: in created\n", timer_get(&timer));
    for (i = 0; i < NUM_THREADS; i++) pari_thread_alloc(&pth[i], 1048576000, gel(in, i+1));
    for (i = 0; i < NUM_THREADS; i++) pthread_create(&th[i], NULL, &twin_smooth_range_d_small_bulk, (void*)&pth[i]);
    out = cgetg(NUM_THREADS+1, t_VEC);
    for (i = 0; i < NUM_THREADS; i++) pthread_join(th[i],(void*)&gel(out, i+1));
    for (i = 1; i < lg(out); i++) for (j = 1; j < lg(gel(out, i)); j++) pari_fprintf(output, "%Ps: %Ps\n", gmael3(out, i, j, 1), gmael3(out, i, j, 2));
    for (i = 1; i < NUM_THREADS; i++) pari_thread_free(&pth[i]);
    pari_fprintf(status, "%d: out written, latest disc %Ps\n", timer_get(&timer), gel(gel(in, lg(in)-1), lg(gel(in, lg(in)-1))-1));
    set_avma(ltop);
  }
  fclose(output); fclose(status);
  pari_close();
  return 0;
}