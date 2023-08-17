# Introduction
Those are algorithms related to computations in the infrastructure of a real quadratic number field, based on so-called (f, p)-representations of ideals ([[1]](#1), primarily chapters 11 and 12).

The goal of this implementation was to compute twin smooth integers of cryptographic size following [[2]](#2) and [[3]](#3); for more background information on this problem see [[4]](#4).

This is an improved version of work I did for my Master Thesis, which can be found [here](https://github.com/db711/twin_smooth).

# How to use
All of the code is written for use with [PARI/GP](https://pari.math.u-bordeaux.fr/). If you plan on using this, you can compile a shared library (for use in GP) with the included Makefile via `make shared`; to use specific functions, they have to be installed in GP first, see the PARI/GP documentation on how to do this or possibly also the example in `twin_smooths.gp`.
All the functions are well documented in the respective `*.h` files.

# Data
The `data` directory contains some files that were created with the algorithms, including a comparison of the algorithms to test compact representations for smoothness as well as some lists of twin smooths (that is a list of numbers $m$, such that $m(m+1)$ is $B$-smooth) for various values of $B$.

# References
<a id="1">[1]</a> \
M. J. Jacobson, H. C. Williams (2009) \
Solving the Pell Equation \
Springer \
https://link.springer.com/book/10.1007/978-0-387-84923-2

<a id="2">[2]</a> \
C. Størmer (1897) \
Quelques théorèmes sur l’équation de Pell $x^2 −Dy^2 = \pm 1$ et leurs applications \
Skrift. Vidensk. Christiania I. Math.-naturv. Klasse, no. 2, p. 48

<a id="3">[3]</a> \
D. H. Lehmer (1964) \
On a problem of Störmer \
Illinois Journal of Mathematics, vol. 8, no. 1, pp. 57–79 \
https://projecteuclid.org/journals/illinois-journal-of-mathematics/volume-8/issue-1/On-a-problem-of-St%C3%B8rmer/10.1215/ijm/1256067456.pdf

<a id="4">[4]</a> \
Giacomo Bruno, Maria Corte-Real Santos, Craig Costello, Jonathan Komada Eriksen, Michael Meyer, Michael Naehrig, Bruno Sterner \
Cryptographic Smooth Neighbors \
Preprint, Cryptology ePrint Archive: Report 2022/1439 \
https://eprint.iacr.org/2022/1439.pdf
