#include "stormer.h"

static inline GEN
createnode(GEN bv, GEN d, long size, GEN p, GEN prev)
/* Create node.
 * Input:   bit vector bv \in {0, 1}^*;
            t_INT d;
            size (in longwords) of the largest d to be stored in the tree;
            p (potential);
            other node prev or NULL.
 * Output:  [bv, d, prev].
            (bv and d are copied, prev is stored as a pointer only)
*/
{
    if (typ(bv) != t_VECSMALL) pari_err_TYPE("createnode",bv);
    if (typ(d) != t_INT) pari_err_TYPE("createnode",d);
    if (prev != NULL && typ(prev) != t_VEC) pari_err_TYPE("createnode",prev);
    GEN res = cgetg(5,t_VEC);
    gel(res,1) = gcopy(bv);
    gel(res,2) = cgeti(size); gaffect(d,gel(res,2));
    gel(res,3) = gcopy(p);
    gel(res,4) = prev;
    return res;
}

static inline int
isleaf (GEN node, long length) { return (lg(gel(node,1)) > length); }

static inline int
isleftchild (GEN node) { return (lg(gel(node,1)) > 1 && (long)gel(gel(node,1),lg(gel(node,1))-1) == 0); }

static inline GEN
leftchild(GEN node, long length)
/* Left child.
 * Input:   node [bv, d, prev] (as returned by createnode or *child);
 * Output:  Node [vec_append(bv,0), d, gel(node,3)-log(prime(lg(bv))), node].
*/
{
    GEN res = cgetg(5,t_VEC);
    pari_sp av;
    gel(res,1) = vecsmall_append(gel(node,1),0);
    gel(res,2) = cgeti(gsizeword(gel(node,2))); gaffect(gel(node,2),gel(res,2)); //remove gaffect?
    av = avma; gel(res,3) = gerepileupto(av,subrr(gel(node,3),glog(prime(length-lg(gel(node,1))+1),DEFAULTPREC)));
    gel(res,4) = node;
    return res;
}

static inline GEN
rightchild(GEN node, long length)
/* Right child.
 * Input:   node [bv, d, prev] (as returned by createnode or *child);
            length (of the considered list of primes).
 * Output:  Node [vec_append(bv,1), d*prime(lg(bv)), gel(node,3), node];
*/
{
    GEN res = cgetg(5,t_VEC);
    pari_sp av;
    gel(res,1) = vecsmall_append(gel(node,1),1);
    gel(res,2) = cgeti(gsizeword(gel(node,2)));
    av = avma; gaffect(mulii(gel(node,2),prime(length-lg(gel(node,1))+1)),gel(res,2)); set_avma(av);
    gel(res,3) = gcopy(gel(node,3)); //remove gcopy?
    gel(res,4) = node;
    return res;
}

GEN
stormer_gen(long length, GEN d, GEN lb, GEN ub, GEN bv, long* h, long l)
{ 
    if (typ(d) != t_INT) pari_err_TYPE("stormer_gen",d);
    if (bv != NULL && typ(bv) != t_VECSMALL) pari_err_TYPE("stormer_gen",bv);
    GEN node, p, lc;
    pari_sp av = avma;
    long size, i;
    if (cmpii(d,ub) > 0 || length < l-*h) return NULL;
    size = gsizeword(ub);
    p = glog(d,DEFAULTPREC);
    for (i=1; i <= length; i++) p = addrr(p,glog(prime(i),DEFAULTPREC));
    p = subrr(p,glog(lb,DEFAULTPREC));
    if (cmpri(p,gen_0) < 0) return NULL;
    node = gerepileupto(av,createnode(cgetg(1,t_VECSMALL),d,size,p,NULL));
    if (bv == NULL)
    {
        while (!isleaf(node,length)) 
        {
            if (length-lg(gel(node,1)) < l-*h)
            {
                node = rightchild(node,length);
                (*h)++;
            }
            else
            {
                av = avma;
                lc = leftchild(node,length);
                if (cmpri(gel(lc,3),gen_0) < 0)
                {
                    set_avma(av);
                    node = rightchild(node,length);
                    (*h)++;
                }
                else node = lc;
            }
        }
    }
    else
    {
        for (i = 1; i < lg(bv); i++)
        {
            if ((long)gel(bv,i) == 1) node = rightchild(node,length);
            else node = leftchild(node,length);
        }
    }
    return node;
}

GEN 
stormer_next(GEN node, long length, GEN ub, long* h, long l, long m)
{ 
    GEN prev, rc, lc;
    pari_sp lbot = avma;
    long cont;
    do
    {
        prev = node;
        node = gel(node,4);
        if (isleftchild(prev))
        {
            if (*h >= m) continue;
            cont = 0;
            set_avma((pari_sp)gel(node,3)); //we start appending new nodes here
            rc = rightchild(node,length);
            (*h)++;

            if (cmpii(ub,gel(rc,2)) < 0)
            {
                (*h)--;
                set_avma(lbot);
                continue;
            }
            else
            {
                node = rc;
                while (!isleaf(node,length))
                {
                    lbot = avma;
                    lc = leftchild(node,length);
                    if ((length-lg(gel(node,1)) < l-*h) || (cmpri(gel(lc,3),gen_0) < 0))
                    {
                        rc = rightchild(node,length);
                        (*h)++;
                        if (cmpii(ub,gel(rc,2)) < 0) 
                        {
                            (*h)--;
                            set_avma(lbot);
                            cont = 1;
                        }
                        else
                        {
                            node = rc;
                            continue;
                        }
                    }
                    if (cont) continue;
                    else node = lc;
                }
            }
            return node;
        }
        else (*h)--;
    } while (node != NULL);
    return NULL;
}