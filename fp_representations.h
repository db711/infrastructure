#ifndef FP_REPRESENTATIONS_H
#define FP_REPRESENTATIONS_H
#include "real_quadratic_orders.h"

struct fp_rep 
{
    GEN f;
    GEN p;
    GEN b;
    GEN d;
    GEN k;
};

struct fp_rep fprepinit (GEN f, GEN p, GEN b, GEN d, GEN k);
/*
 * (f, p) representation initialization.
 * Input: t_REAL f >= 1 (or NULL for unknown), t_INT p, d and k, real quadratic ideal b (as output by rqoinit).
 *        Here p is a precision value, that is already on the stack and only stored as a pointer.
 * Output: struct fp_rep with the input data.
 * (Representation is checked for validity, ideal b is not.)
 */

struct fp_rep fpremove (struct fp_rep fprep, GEN T, GEN C, GEN s);
/*
 * Input: An (f, p) representation fprep (as output by fprepinit), t_INT T, C and s.
 * Output: An (f + 9/8, p) representation.
 */

struct fp_rep numult (GEN O, struct fp_rep fprep1, struct fp_rep fprep2); 
/*
 * Input: Reduced (f', p) (resp. (f'', p)) representations fprep1 of ideal a' and fprep2 of ideal a''.
 * Output: Reduced (f, p) representation of the product a = a'*a''.
 */

struct fp_rep wnear (struct fp_rep fprep, GEN w);
/*
 * Input: A reduced (f, p) representation fprep of ideal a, t_INT w.
 * Output: A w-near (f + 9/8, p) representation of ideal a.
 */

#endif