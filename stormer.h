#ifndef STORMER_H
#define STORMER_H
#include <stdlib.h>
#include <pari/pari.h>

typedef enum
{
   zero,
   one,
   undefined
} tl; //three-valued logic

typedef struct 
{
   tl* data;
   const ulong* size;
} tlv; //three-valued logic vector

typedef struct 
{
    tlv* vector;
    GEN sol; //sum of logs
} node;

GEN logsofprimes(ulong B);
/* Logarithms of primes.
 * Input:   Smoothness bound B.
 * Output:  A GEN holding a t_VEC containing a sorted list of the logarithms of the primes <= B (which are themselves t_REAL of precision DEFAULTPREC).
*/

char* tlv2str(const tlv* vector, char* str);
/* three-valued logic vector to string.
 * Input:   tlv vector; char* str.
            User must ensure that str can store *(tlv.size)+1 many characters.
 * Ouput:   tlv.data converted to a word over the alphabet {0, 1, x}^* and stored in str.
*/

tlv* str2tlv(const char* str, const ulong* size);
/* string to three-valued logic vector.
 * Input:   char* str, holding a word over the alphabet {0, 1, x}^*;
            ulong* size, pointing to a global length parameter for tlvs.
            User must ensure that the string is (at least) (*size) characters long.
 * Output:  A pointer to a tlv that is allocated in memory.
*/

node* rightchild(const node* parent, node* child, GEN lop);
/* right child.
 * Input:   node* parent, node* child;
            GEN lop (logarithms of primes) as returned by logsofprimes (defined globally).
 * Output:  The right child of node parent is stored in child and returned.
*/

node* leftchild(const node* parent, node* child, GEN lop);
/* left child.
 * Input:   node* parent, node* child;
            GEN lop (logarithms of primes) as returned by logsofprimes (defined globally).
 * Output:  The left child of node parent is stored in child and returned.
*/

#endif