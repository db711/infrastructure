#ifndef FP_REPRESENTATIONS_H
#define FP_REPRESENTATIONS_H
#include "real_quadratic_orders.h"

GEN fprepinit (GEN f, GEN p, GEN b, GEN d, GEN k);
/*
 * (f, p) representation initialization.
 * Input: t_REAL f >= 1 (or NULL for unknown), t_INT p, d and k, real quadratic ideal b (as output by rqoinit).
 *        Here p is a precision value, that is already on the stack and only stored as a pointer.
 * Output: [[f, p], [S, [Q, P]], d, k].
 * (Does not perform any checks for the validity of this representation.)
*/

GEN remove (GEN fprep, GEN T, GEN C, GEN s);
/*
 * Input: (f, p) representation fprep (as output by fprepinit), integers C and s.
 * Output: An (f + 9/8, p) representation.
*/

GEN numult (GEN fprep1, GEN fprep2);
/*
 * Input: Reduced (f', p) (resp. (f'', p)) representation fprep1 of ideal a' (resp. fprep2 of ideal a'') (as output by fprepinit).
 * Output: Reduced (f, p) representation of the product a = a'*a''.
*/

#endif