#ifndef _list_h
#define _list_h
/**
 * @file list.h
 * @author Mark Hoemmen
 * @since 11 July 2006
 * @date Time-stamp: <2008-07-16 10:18:01 mhoemmen>
 *
 * Definition of a general singly-linked list datatype that looks like
 * the lists found in Lisp.
 *
 * Copyright (c) 2008, Regents of the University of California 
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright 
 *   notice, this list of conditions and the following disclaimer in 
 *   the documentation and/or other materials provided with the 
 *   distribution.
 *
 * * Neither the name of the University of California, Berkeley, nor
 *   the names of its contributors may be used to endorse or promote
 *   products derived from this software without specific prior
 *   written permission.  
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <stdio.h>

/**
 * Don't use this one.  We need it so we can define the recursive structure.
 */
struct __list_node_t
{
  void* item;
  struct __list_node_t* next;
};

/**
 * Use this one.
 */
typedef struct __list_node_t list_node_t;

/**
 * Returns nonzero if there is circular structure in the list whose
 * head is HEAD; else returns zero.
 */
int
circular_structure_p (list_node_t* head);

/**
 * A list datatype that includes the tail (last element) of the list.  
 * Handy for appending elements in O(1) time, but takes more space and
 * costs a little runtime.  As a result, list_t member functions do 
 * more safety checking than the alist_t datatype below (which see).
 */
typedef struct {
  list_node_t* head;
  list_node_t* tail;
} list_t;

/**
 * Returns an empty list.
 */
list_t
list_create ();

/**
 * Returns a new list node.
 */
list_node_t*
list_node_create (void* item, list_node_t* next);


/**
 * Pops off the head of the list L and stores it in *pResult.
 * Be sure to assign the return value back to L, otherwise
 * you'll lose the rest of the list!
 */
list_t
list_pop (list_node_t** pResult, list_t L);

/**
 * Appends the given node to the end of the given list
 * and returns the resulting new list.
 */
list_t
list_append_node (list_t L, list_node_t* node);


/**
 * Creates a new node with "item" as its datum and returns
 * the list which is the original list appended by the new
 * node.
 */
list_t
list_append_item (list_t L, void* item);

/**
 * Frees the given list and its data and returns the empty list.
 */
list_t
list_destroy (list_t L);

/**
 * Like list_destroy(), except calls the given destructor on each node.
 */
list_t
list_destroy_custom (list_t L, void (*destructor) (void*));

/**
 * Returns zero if the given list is empty and nonzero otherwise.
 */
int
list_empty_p (list_t L);

/**
 * Returns (a pointer to) the n-th element in the given list, 
 * in which n is zero-based.
 */
void*
list_nth (list_t L, const int n);

/**
 * Returns the CDR of the list (the list consisting of every node in L
 * but the first node).  Useful for writing recursive functions, e.g. 
 * 
 * void foo (L) { 
 *   if (! list_empty_p (L)) { 
 *     dosomething (list_car (L)); 
 *     foo (list_cdr (L)) 
 *   } 
 * }
 */
list_t
list_cdr (list_t L);

/**
 * Returns a pointer to the value stored at the head of the given list.
 * Useful for writing recursive functions; see list_cdr.
 */
void*
list_car (list_t L);

void*
list_cadr (list_t L);

void*
list_caddr (list_t L);

void*
list_cadddr (list_t L);

/**
 * (setf (car L) item)
 */
void
list_set_car (list_t L, void* item);

/**
 * Returns the length (number of elements) of the given list.
 */
int
list_length (list_t L);

/**
 * Searches for the given ITEM in the given list L, using the predicate
 * EQUAL to test for equality.  If no items are in the list L that are 
 * EQUAL to ITEM, returns NULL, else returns a pointer to the first node
 * containing ITEM.
 */
list_node_t* 
list_find (list_t L, void* item, int (*equal) (void*, void*));

/**
 * Sorts the given list in place, using the given predicate PRED
 * (which works like strcmp) to compare values.  Be sure to assign the
 * return value back to L to recover the correct head and tail.
 *
 * @warn The sort used may not be stable.
 */
list_t 
list_sort (list_t L, int (*pred) (void*, void*));

/**
 * Assuming that the given node contains char* data, prints the node
 * to the given output stream.
 */
void
string_list_node_print (FILE* out, list_node_t* node);

/**
 * Assuming that the given list contains char* data, prints the list
 * to the given output stream.
 */
void 
string_list_print (FILE* out, list_t L);


/* ========================================================================= */
/* alist_t                                                                   */
/* ========================================================================= */

typedef struct {
  list_node_t* head;
} alist_t;

alist_t 
alist_create ();

alist_t
alist_push_node (alist_t L, list_node_t* node);

alist_t
alist_push_item (alist_t L, void* item);

alist_t
alist_pop (list_node_t** pResult, alist_t L);

alist_t
alist_append_node (alist_t L, list_node_t* node);

alist_t
alist_append_item (alist_t L, void* item);

alist_t
alist_destroy (alist_t L);

alist_t
alist_destroy_custom (alist_t L, void (*destructor) (void*));

int
alist_empty_p (alist_t L);

void*
alist_nth (alist_t L, const int n);

alist_t
alist_cdr (alist_t L);

void*
alist_car (alist_t L);

int
alist_length (alist_t L);

list_node_t*
alist_find (alist_t L, void* item, int (*equal) (void*, void*));

alist_t
alist_sort (alist_t L, int (*pred) (void*, void*));



#endif /* _list_h */
