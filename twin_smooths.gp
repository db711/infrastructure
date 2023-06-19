{
    install(twin_smooth_d, "GGUG", TwinSmoothD, "./infrastructure.so");
    addhelp(TwinSmoothD, "TwinSmoothD(S, d, m, P):
        Compute the twin B-smooth integers corresponding to d (even, squarefree and B-smooth); S is the t_VECSMALL of primes smaller than B, m is a technical upperbound, not greater than (B - 1)/2, P is the product of all numbers in S");
}

twin_smooths(B) =
{
    my(res = []);
    S = Vecsmall(primes(primepi(B)));export(S);
    P = multiply(S);export(P);
    m = (S[length(S)]-1)/2;export(m);
    export(multiply);
    start = getwalltime;
    parfor(i=0,2^(length(S))-3,TwinSmoothD(S,2*multiply(vecextract(S,2^(length(S))-1-i)),m,P),r,if(r,write(output,i,": ",r)));
    write(output,2^(length(S))-3,": ",TwinSmoothD(S,2,m,P));
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