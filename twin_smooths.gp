{
    install(twin_smooth_d, "GGUG", TwinSmoothD, "./infrastructure.so");
    addhelp(TwinSmoothD, "TwinSmoothD(S, d, m, P):
        Compute the twin B-smooth integers corresponding to d (even, squarefree and B-smooth); S is the t_VECSMALL of primes smaller than B, m is a technical upperbound, not greater than (B - 1)/2, P is the product of all numbers in S");
    install(twin_smooth_range_d_small, "UGUUU", TwinSmoothRangeDSmall, "./infrastructure.so");
    addhelp(TwinSmoothRangeDSmall, "TwinSmoothRangeDSmall(B, d, bot, top, m):
        Compute twin B-smooth integers corresponding to d(even, squarefree and B-smooth); only check solutions that could lead to x having bits in [bot, top], m is an upper bound for solutions checked (choose m <= 4)");
}

twin_smooths(B) =
{
    my(res = []);
    S = Vecsmall(primes(primepi(B)));export(S);
    P = multiply(S);export(P);
    m = (S[length(S)]+1)/2;export(m);
    export(multiply);
    start = getwalltime;
    parfor(i=0,2^(length(S))-3,TwinSmoothD(S,2*multiply(vecextract(S,2^(length(S))-1-i)),m,P),r,if(r,write(output,i,": ",r)));
    write(output,2^(length(S))-3,": ",TwinSmoothD(S,2,m,P));
    write(output,"Computation took ",getwalltime-start,"ms");
}

twin_smooths_range_small(B) = 
{
    my(res = []);
    S = Vecsmall(primes(primepi(B)));export(S);
    m = 4;export(m);
    export(multiply);
    start = getwalltime;
    write(output,0,": ",TwinSmoothRangeDSmall(S[length(S)],2,245,275,m));
    parfor(i=2,2^(length(S))-1,iferr(TwinSmoothRangeDSmall(S[length(S)],2*multiply(vecextract(S,i)),245,275,m),E,print(i,": ",component(E,1))),r,if(r,write(output,i,": ",r)));
    write(output,"Computation took ",getwalltime-start,"ms");
}

multiply(S) =
{
    my(p = 1);
    for(i=1,length(S),p*=S[i]);
    p;
}

default(parisize,1073741824);
default(nbthreads,48);
