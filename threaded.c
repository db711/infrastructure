#include "twin_smooths.h"
#include "stormer.h"
#include <pthread.h> 

/*#define NUM_THREADS 30 // number of threads to run (at once)
#define SMOOTHNESS_BOUND 65536// 2^16
#define LOWER_BOUND "1427247692705959881058285969449495136382746624" // 2^150
#define UPPER_BOUND "1606938044258990275541962092341162602522202993782792835301376" // 2^200
#define MIN_FACTOR 15 // minimum number of factors in Stormer discriminant
#define MAX_FACTOR 20 // maximum number of factors in Stormer discriminant
#define STARTING_D "36893488147419103232" // 2^63
#define CONDUCTOR "2147483648" // 2^31 = sqrt(STARTING_D/2)*/

#define NUM_THREADS 10 // number of threads to run (at once)
#define SMOOTHNESS_BOUND 65536// 2^16
#define LOWER_BOUND "4294967296" // 2^32
#define UPPER_BOUND "340282366920938463463374607431768211456" // 2^128
#define MIN_FACTOR 4 // minimum number of factors in Stormer discriminant
#define MAX_FACTOR 16 // maximum number of factors in Stormer discriminant
#define STARTING_D "131072" // 2^15
#define CONDUCTOR "128" // 2^7 = sqrt(STARTING_D/2)

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
regulator_cryptographic_(void *arg)
{
  GEN stormer, ret, in; // [k, [np, l, m, lb, ub, d, f]]
  long i, k, np, l, m, h;
  pari_sp av;
  char output_s[100], status_s[100];
  FILE* output = NULL;
  FILE* status = NULL;
  pari_timer timer;
  in = pari_thread_start((struct pari_thread*) arg);
  k = itos(gel(in,1));
  np = itos(gmael(in,2,1));
  l = itos(gmael2(in,2,2));
  m = itos(gmael2(in,2,3));
  h = 0;
  sprintf(output_s, "output_%ld", k);
  sprintf(status_s, "status_%ld", k);
  if (NULL == (output = fopen(output_s, "w")))
  {
    pari_printf("Error opening output file %s. Aborting.\n", output_s);
    pari_close();
    return NULL;
  }
  if (NULL == (status = fopen(status_s, "w")))
  {
    pari_printf("Error opening status file %s. Aborting.\n", status_s);
    pari_close();
    return NULL;
  }
  stormer = stormer_gen(np, gmael2(in,2,6), gmael2(in,2,4), gmael2(in,2,5), NULL, &h, l);
  for (i = 0; i < k; i++) stormer = stormer_next(stormer, np, gmael2(in,2,5), &h, l, m);
  timer_start(&timer);
  while (NULL != stormer)
  {
    av = avma;
    ret = regulator_cryptographic(gel(stormer,2), gmael2(in,2,7));
    if (lg(ret) > 1) pari_fprintf(output, "%Ps: %Ps\n", gel(stormer,2), ret);
    if (timer_get(&timer) > 3600000)
    {
      pari_fprintf(status, "%Ps\n", gel(stormer,2));
      timer_delay(&timer);
    }
    set_avma(av);
    for (i = 0; i < NUM_THREADS; i++) stormer = stormer_next(stormer, np, gmael2(in,2,5), &h, l, m);
  }
  pari_thread_close();
  return NULL;
}

void *
pell_and_boost_(void *arg)
{
  GEN stormer, ret, in; // [k, [np, l, m, lb, ub, d, f]]
  long i, k, np, l, m, h;
  pari_sp av;
  char output_s[100], status_s[100];
  FILE* output = NULL;
  FILE* status = NULL;
  pari_timer timer;
  in = pari_thread_start((struct pari_thread*) arg);
  k = itos(gel(in,1));
  np = itos(gmael(in,2,1));
  l = itos(gmael2(in,2,2));
  m = itos(gmael2(in,2,3));
  h = 0;
  sprintf(output_s, "output_%ld", k);
  sprintf(status_s, "status_%ld", k);
  if (NULL == (output = fopen(output_s, "w")))
  {
    pari_printf("Error opening output file %s. Aborting.\n", output_s);
    pari_close();
    return NULL;
  }
  if (NULL == (status = fopen(status_s, "w")))
  {
    pari_printf("Error opening status file %s. Aborting.\n", status_s);
    pari_close();
    return NULL;
  }
  stormer = stormer_gen(np, gmael2(in,2,6), gmael2(in,2,4), gmael2(in,2,5), NULL, &h, l);
  for (i = 0; i < k; i++) stormer = stormer_next(stormer, np, gmael2(in,2,5), &h, l, m);
  timer_start(&timer);
  while (NULL != stormer)
  {
    av = avma;
    ret = pell_and_boost(gel(stormer,2), gmael2(in,2,7));
    if (lg(ret) > 1) pari_fprintf(output, "%Ps: %Ps\n", gel(stormer,2), ret);
    if (timer_get(&timer) > 3600000)
    {
      pari_fprintf(status, "%Ps\n", gel(stormer,2));
      timer_delay(&timer);
    }
    set_avma(av);
    for (i = 0; i < NUM_THREADS; i++) stormer = stormer_next(stormer, np, gmael2(in,2,5), &h, l, m);
  }
  pari_thread_close();
  return NULL;
}

/*
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
*/

int
main(void)
{
  long np, i;
  pari_init(1048576000,SMOOTHNESS_BOUND+1);
  pthread_t th[NUM_THREADS];
  struct pari_thread pth[NUM_THREADS];
  GEN d_start, f, lb, ub, in;
  pari_sp av;
  d_start = strtoi(STARTING_D);
  f = strtoi(CONDUCTOR);
  lb = strtoi(LOWER_BOUND);
  ub = strtoi(UPPER_BOUND);
  av = avma; np = itos(primepi(stoi(SMOOTHNESS_BOUND))); set_avma(av);

  in = cgetg(8, t_VEC);
  gel(in, 1) = stoi(np);
  gel(in, 2) = stoi(MIN_FACTOR);
  gel(in, 3) = stoi(MAX_FACTOR);
  gel(in, 4) = lb;
  gel(in, 5) = ub;
  gel(in, 6) = d_start;
  gel(in, 7) = f;
  for (i = 0; i < NUM_THREADS; i++) pari_thread_alloc(&pth[i], 1048576000, mkvec2(stoi(i),in));
  for (i = 0; i < NUM_THREADS; i++) pthread_create(&th[i], NULL, &regulator_cryptographic_, (void*)&pth[i]);
  for (i = 0; i < NUM_THREADS; i++) pthread_join(th[i],NULL);
  for (i = 0; i < NUM_THREADS; i++) pari_thread_free(&pth[i]);

  pari_close();
  return 0;
}