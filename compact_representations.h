#ifndef COMPACT_REPRESENTATIONS_H
#define COMPACT_REPRESENTATIONS_H
#include "fp_representations.h"

GEN expandcr(GEN O, GEN cr);
/* Expand compact representation.
 * Input:   Real quadratic order O (as output by rqoinit);
            list [[m_1, n_1], L_1], ..., [[m_l, n_l], L_l], [[m_{l+1}, n_{l+1}], L_{l+1}]] (as output for example by cr), such that
                theta = lambda_{l+1} * \prod_{i=0}^l (lambda_i/L_{i+1})^{2^(l-i)}, where lambda_i = (m_i + n_i*sqrtd(d))/s
            is a compact representation of some theta in O.
 * Output:  [a, b], such that theta = (a + b*sqrt(d))/s.
*/

#endif