#include "twin_smooths.h"
#include <pthread.h> 

#define NUM_THREADS 4 
#define LOWER_BOUND 245
#define UPPER_BOUND 275
#define SMOOTHNESS_BOUND 100

void *
find_smooths(void *arg) // UNUSED
{
    GEN F, M;
    M = pari_thread_start((struct pari_thread*) arg); // = [S, P, d]
    F = regulator_range(rqoinit(gel(M,2)),LOWER_BOUND,UPPER_BOUND) == NULL ? NULL : twin_smooth_d(gel(M,1),gel(M,2),0,gel(M,3));
    pari_thread_close();
    return (void*)F;
}

GEN 
twin_smooth_range_d_small_wrapper(GEN i, ulong B, GEN S, ulong bot, ulong top, ulong m)
{
    GEN d, D;
    ulong j;
    pari_sp ltop = avma, av;
    d = gen_2;
    D = shallowextract(S,i); 
    av = avma;
    for (j = 1; j < lg(D); j++) d = gerepileupto(av,mulis(d,D[j]));
    return gerepileupto(ltop,twin_smooth_range_d_small(B,d,bot,top,m));
}

int
main(void)
{
    pari_init(1073741824,0); //1GB
    default0("nbthreads","NUM_THREADS");
    GEN end, e, S, fun;
    parfor_t T;
    ulong B = SMOOTHNESS_BOUND;
    pari_sp ltop = avma, av;
    S = primes_interval_zv(2,B); 
    //install((void *)&twin_smooth_range_d_small_wrapper,"TwinSmoothRangeDSmallWrapper","GUGUUU");
    //fun = strtoclosure("TwinSmoothRangeDSmallWrapper",5,stoi(B),S,stoi(245),stoi(275),stoi(4));
    fun = snm_closure(install((void *)&twin_smooth_range_d_small_wrapper,"TwinSmoothRangeDSmallWrapper","GUGUUU"),mkvec5(stoi(B),S,stoi(LOWER_BOUND),stoi(UPPER_BOUND),stoi(4)));
    end = subii(powis(gen_2,lg(S)-1),gen_1);
    parfor_init(&T,gen_2,end,fun);
    av = avma;
    while ((e = parfor_next(&T)))
    {
        if (lg(gel(e,2)) > 1) pari_printf("%Ps\n",gel(e,2));
        set_avma(av);
    }
    set_avma(ltop);
    pari_close();
    return 0;
    /*GEN S, s, d, P, M;
    ulong B = SMOOTHNESS_BOUND, i;
    forsubset_t T;
    pari_sp av;
    pthread_t ths[NUM_THREADS];
    struct pari_thread pths[NUM_THREADS];
    pari_init(1073741824,0); //1GB
    S = primes_upto_zv(B);
    P = gen_1;
    av = avma;
    for (i = 1; i < lg(S); i++) P = gerepileupto(av,mulis(P,S[i]));
    forallsubset_init(&T,lg(S)-1);
    av = avma;
    forsubset_next(&T); forsubset_next(&T);
    //if(regulator_range(rqoinit(gen_2),245,275)) pari_printf("%Ps\n",gen_2);
    for (i = 0; i < NUM_THREADS; i++)
    {
        pari_thread_alloc(&pths[i],1073741824,M); //1GB
    }*/
}