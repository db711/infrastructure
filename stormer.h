#ifndef STORMER_H
#define STORMER_H
#include <stdlib.h>
#include <pari/pari.h>

//#define LOG_2_B 64 //an upper bound on the number of bits in log_2(B) for the smoothness bound B

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

char* tlv2str(const tlv* vector, char* str);
/* three-valued logic vector to string
 * Input:   tlv vector; char* str.
            User must ensure that str can store *(tlv.size)+1 many characters.
 * Ouput:   tlv.data converted to a word over the alphabet {0, 1, x}^*.
*/

tlv* str2tlv(const char* str, const ulong* size);
/* string to three-valued logic vector
 * Input:   char* str, holding a word over the alphabet {0, 1, x}^*;
            ulong* size, pointing to a global length parameter for tlvs.
            User must ensure that the string is (at least) (*size) characters long.
 * Output:  A pointer to a tlv that is allocated in memory.
*/

node* rightchild(const node* parent, node* child);
node* leftchild(const node* parent, node* child);

#endif