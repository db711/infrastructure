#include "stormer.h"

char*
tlv2str(const tlv* vector, char* str)
{
    ulong i;
    i = *(vector->size);
    while (i > 0)
    {
        i--;
        if (vector->data[i] == zero) str[i] = '0';
        else if (vector->data[i] == one) str[i] = '1';
        else str[i] = 'x'; //fallthrough?
    }
    str[*(vector->size)] =  '\0';
    return str;
}

tlv*
str2tlv(const char* str, const ulong* size)
{
    ulong i;
    tl* data = malloc((*size)*sizeof(tl));
    tlv* ret = malloc(sizeof(tlv));
    i = (*size);
    while (i > 0)
    {
        i--;
        if (str[i] == '0') data[i] = zero;
        else if (str[i] == '1')  data[i] = one;
        else data[i] = undefined; //fallthrough?
    }
    ret->data = data;
    ret->size = size;
    return ret;
}

node* 
rightchild(const node* parent, node* child)
{
    if (((parent->vector)->data)[0] != undefined) return NULL; //parent is a leaf
    ulong i;
    i = 1;
    while (i < *(parent->vector->size) && (parent->vector->data)[i] == undefined)
    {
        (child->vector->data)[i-1] = undefined;
        i++;
    }
    (child->vector->data)[i-1] = one;
    for (;i < *(parent->vector->size); i++) (child->vector->data)[i] = (parent->vector->data)[i];
    return child;
}

node* 
leftchild(const node* parent, node* child)
{
    if (((parent->vector)->data)[0] != undefined) return NULL; //parent is a leaf
    ulong i;
    i = 1;
    while (i < *(parent->vector->size) && (parent->vector->data)[i] == undefined)
    {
        (child->vector->data)[i-1] = undefined;
        i++;
    }
    (child->vector->data)[i-1] = zero;
    for (;i < *(parent->vector->size); i++) (child->vector->data)[i] = (parent->vector->data)[i];
    return child;
}