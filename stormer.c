#include "stormer.h"

static inline GEN
createnode(GEN bv, GEN d, long size, GEN prev)
/* Create node.
 * Input:   bit vector bv \in {0, 1}^*;
            t_INT d;
            size (in longwords) of the largest d to be stored in the tree;
            other node prev or NULL.
 * Output:  [bv, d, prev].
            (bv and d are copied, prev is stored as a pointer only)
*/
{
    if (typ(bv) != t_VECSMALL) pari_err_TYPE("createnode",bv);
    if (typ(d) != t_INT) pari_err_TYPE("createnode",d);
    if (prev != NULL && typ(prev) != t_VEC) pari_err_TYPE("createnode",prev);
    GEN res = cgetg(4,t_VEC);
    gel(res,1) = gcopy(bv);
    gel(res,2) = cgeti(size); gaffect(d,gel(res,2));
    gel(res,3) = prev;
    return res;
}

static inline int
isleaf (GEN node, long length) { return (lg(gel(node,1)) > length); }

static inline int
isleftchild (GEN node) { return (lg(gel(node,1)) > 1 && (long)gel(gel(node,1),lg(gel(node,1))-1) == 0); }

static inline GEN
leftchild(GEN node)
/* Left child.
 * Input:   node [bv, d, prev] (as returned by createnode or *child);
 * Output:  Node [vec_append(bv,0), d, node].
*/
{
    GEN res = cgetg(4,t_VEC);
    gel(res,1) = vecsmall_append(gel(node,1),0);
    gel(res,2) = cgeti(gsizeword(gel(node,2))); gaffect(gel(node,2),gel(res,2));
    gel(res,3) = node;
    return res;
}

static inline GEN
rightchild(GEN node, long length)
/* Right child.
 * Input:   node [bv, d, prev] (as returned by createnode or *child);
            length (of the considered list of primes).
 * Output:  Node [vec_append(bv,1), d*prime(lg(bv)), node];
*/
{
    GEN res = cgetg(4,t_VEC);
    pari_sp ltop;
    gel(res,1) = vecsmall_append(gel(node,1),1);
    gel(res,2) = cgeti(gsizeword(gel(node,2)));
    ltop = avma; gaffect(mulii(gel(node,2),prime(length-lg(gel(node,1))+1)),gel(res,2)); set_avma(ltop);
    gel(res,3) = node;
    return res;
}

GEN
stormer_gen(long length, GEN d, GEN ub, GEN bv)
{ 
    if (typ(d) != t_INT) pari_err_TYPE("stormer_gen",d);
    if (bv != NULL && typ(bv) != t_VECSMALL) pari_err_TYPE("stormer_gen",bv);
    GEN node;
    pari_sp av = avma;
    long size, i;
    if (cmpii(d,ub) > 0) return NULL;
    size = gsizeword(ub);
    node = gerepileupto(av,createnode(cgetg(1,t_VECSMALL),d,size,NULL));
    if (bv == NULL) while (!isleaf(node,length)) node = leftchild(node);
    else
    {
        for (i = 1; i < lg(bv); i++)
        {
            if ((long)gel(bv,i) == 1) node = rightchild(node,length);
            else node = leftchild(node);
        }
    }
    return node;
}

GEN 
stormer_next(GEN node, long length, GEN ub, long* h, long l, long m)
{ 
    GEN prev, rc;
    pari_sp lbot = avma;
    do
    {
        prev = node;
        node = gel(node,3);
        if (node == NULL) return NULL;
        if (isleftchild(prev))
        {
            set_avma((pari_sp)gel(node,2));
            rc = rightchild(node,length); 
            (*h)++;
            if ( *h > m || cmpii(ub,gel(rc,2)) < 0) 
            {
                (*h)--;
                set_avma(lbot);
            }
            else
            {
                node = rc;
                while(!isleaf(node,length) && length-lg(gel(node,1)) >= l-*h) node = leftchild(node);
                while(!isleaf(node,length) && *h < l) 
                {
                    node = rightchild(node,length);
                    (*h)++;
                }
                while(!isleaf(node,length)) node = leftchild(node);
                set_avma(lbot);
                return node;
            }
        }
        else (*h)--;
    } while (1);
}