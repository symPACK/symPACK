/**
 * @author Mark Hoemmen
 * @file list.c
 * @date Time-stamp: <2008-07-16 10:10:12 mhoemmen>
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
#include <bebop/util/list.h>
#include <bebop/util/log.h>
#include <bebop/util/util.h>

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>


list_t
list_create ()
{
  list_t s;
  s.head = NULL;
  s.tail = NULL;
  return s;
}

alist_t
alist_create ()
{
  alist_t L;
  L.head = NULL;
  return L;
}

list_node_t*
list_node_create (void* item, list_node_t* next)
{
  list_node_t* node = malloc (sizeof (list_node_t));
  assert (node != NULL);
  node->item = item;
  node->next = next;
  return node;
}


int
circular_structure_p (list_node_t* head)
{
  /*
   * We detect circular structure in the list of nodes dangling off of
   * HEAD by stepping simultaneously by ones and twos through the list.
   * If these ever equal one another, we've found circular structure.
   */
  list_node_t* p1 = head;
  list_node_t* p2 = head;

  if (head == NULL)
    return 0;

  while (1)
    {
      p1 = p1->next;
      if (p1 == NULL)
	return 0;

      p2 = p2->next;
      if (p2 == NULL)
	return 0;
      p2 = p2->next;
      if (p2 == NULL)
	return 0;

      if (p1 == p2)
	return 1;
    }
}

/**
 * Returns the tail (last non-NULL node in the list).
 *
 * @warn Does not check for circular structure.
 */
static list_node_t*
find_tail_helper (list_node_t* node, list_node_t* prev)
{
  if (node == NULL)
    return prev;
  else
    return find_tail_helper (node->next, node);
}

/**
 * Returns the tail (last non-NULL node in the list).
 * 
 * @warn Does not check for circular structure.
 */
static list_node_t*
find_tail (list_node_t* node)
{
  return find_tail_helper (node, node);
}


list_t
list_append_node (list_t L, list_node_t* node)
{
  if (node != NULL)
    {
      list_node_t* tail = NULL;

      if (L.head == NULL)
	L.head = node;
      else
	L.tail->next = node;

      /* 
       * Be sure to append whatever structure is dangling off
       * NODE.  In order to get L.tail right, we have to find the
       * end of the string of nodes dangling off NODE.  If there
       * is circular structure, we'll never be able to get the
       * tail right, so we need to check for this.
       */
      assert (! circular_structure_p (node));
      tail = find_tail (node);
      assert (tail != NULL);
      L.tail = tail;
    }

  return L;
}

alist_t
alist_append_node (alist_t L, list_node_t* node)
{
  list_node_t* tail = find_tail (L.head);
  if (tail == NULL) 
    L.head = node; /* the list was empty */
  else
    tail->next = node;

  return L;
}

alist_t
alist_push_node (alist_t L, list_node_t* node)
{
  if (L.head == NULL)
    L.head = node;
  else
    {
      /* WARNING: any structure dangling from NODE will be lost */
      node->next = L.head;
      L.head = node;
    }
  return L;
}

alist_t
alist_push_item (alist_t L, void* item)
{
  return alist_push_node (L, list_node_create (item, NULL));
}


static list_node_t* 
__pop (list_node_t** pResult, list_node_t* L)
{
  if (L == NULL)
    *pResult = NULL;
  else
    {
      list_node_t* head = L;
      L = head->next;
      head->next = NULL; /* ground the free node */
      *pResult = head;
    }
  return L;
}

list_t
list_pop (list_node_t** pResult, list_t L)
{
  L.head = __pop (pResult, L.head);
  if (L.head == NULL)
    L.tail = NULL;

  return L;
}

alist_t
alist_pop (list_node_t** pResult, alist_t L)
{
  L.head = __pop (pResult, L.head);
  return L;
}

list_t
list_append_item (list_t L, void* item)
{
  return list_append_node (L, list_node_create (item, NULL));
}

alist_t
alist_append_item (alist_t L, void* item)
{ 
  return alist_append_node (L, list_node_create (item, NULL));
}

static void
__destroy (list_node_t* L, void (*destructor) (void*))
{
  if (L == NULL)
    return;
  else
    {
      if (L->item != NULL) 
	destructor (L->item);

      L->item = NULL;
      __destroy (L->next, destructor);
      free (L);
    }
}

list_t
list_destroy_custom (list_t L, void (*destructor) (void*))
{
  /* We don't want any infinite loops. */
  assert (! circular_structure_p (L.head));
  __destroy (L.head, destructor);
  return list_create ();
}

alist_t
alist_destroy_custom (alist_t L, void (*destructor) (void*))
{
  __destroy (L.head, destructor);
  return alist_create ();
}

list_t
list_destroy (list_t L)
{
  return list_destroy_custom (L, &free);
}

alist_t
alist_destroy (alist_t L)
{ 
  return alist_destroy_custom (L, &free);
}

void
string_list_node_print (FILE* out, list_node_t* node)
{
  if (node != NULL)
    fprintf (out, "Item: %s\n", (char*) (node->item));
}

void 
string_list_print (FILE* out, list_t L)
{
  list_node_t* node = L.head;

  assert (! circular_structure_p (L.head));
  fprintf (out, "String list:\n");
  while (node != NULL)
    {
      string_list_node_print (out, node);
      node = node->next;
    }
  fprintf (out, "\n");
}


int
list_empty_p (list_t L)
{
  return (L.head == NULL);
}

int
alist_empty_p (alist_t L)
{
  return (L.head == NULL);
}

/**
 * Returns the n-th item (zero-based) in the given list whose head is L.
 */
static void*
__nth (list_node_t* L, const int n)
{
  list_node_t* node = L;
  int i;

  for (i = 0; i < n; i++)
    node = node->next;

  return node->item;
}

void*
list_nth (list_t L, const int n)
{
  assert (n >= 0);
  assert (! list_empty_p (L));

  return __nth (L.head, n);
}

void*
alist_nth (alist_t L, const int n)
{
  assert (n >= 0);
  return __nth (L.head, n);
}

list_t
list_cdr (list_t L)
{
  list_t cdr;

  assert (L.head != NULL);
  if (L.head == L.tail)
    {
      /* The input list contains only one element; 
       * return the empty list. */
      cdr.head = NULL;
      cdr.tail = NULL;
      return cdr;
    }

  cdr.head = (L.head)->next;
  cdr.tail = L.tail;

  return cdr;
}

alist_t
alist_cdr (alist_t L)
{
  alist_t cdr;

  assert (L.head != NULL);
  cdr.head = (L.head)->next;
  return cdr;
}

static inline void*
__car (list_node_t* head)
{
  assert (head != NULL);
  return head->item;
}

void*
list_car (list_t L)
{
  return __car (L.head);
}

void*
alist_car (alist_t L)
{
  return __car (L.head);
}

/**
 * Returns the length of the list L (which is given by its head).
 *
 * @warn This function doesn't check for circular structure.
 */
static int
__length (list_node_t* L)
{
  int len = 0;
  while (L != NULL)
    {
      L = L->next;
      len++;
    }
  return len;
}


int
list_length (list_t L)
{
  assert (! circular_structure_p (L.head));
  return __length (L.head);
}

int
alist_length (alist_t L)
{
  return __length (L.head);
}


/**
 * Merges two lists, each given by their head node.  The predicate
 * function PRED works like strcmp (returns -1 if first < second, 0 if
 * equal, +1 if first > second).
 * 
 * @note We use list_node_t instead of list_t here to avoid the overhead
 * of searching for tails (which would be necessary to maintain the list_t
 * invariants).
 * 
 * @warn This function doesn't check for circular structure.
 */
static list_node_t* 
__merge (list_node_t* L1, list_node_t* L2,
       int (*pred) (void*, void*))
{
  if (L1 == NULL)
    return L2;
  else if (L2 == NULL)
    return L1;
  else
    {
      list_node_t* head;

      if (pred (L1->item, L2->item) < 0)
	{
	  head = L1;
	  head->next = __merge (L1->next, L2, pred);
	}
      else 
	{
	  head = L2;
	  head->next = __merge (L1, L2->next, pred);
	}
      return head;
    }
}


/**
 * Sorts the list L (given by its head) using the predicate PRED
 * (which works like strcmp) and returns the head of the resulting
 * list.
 *
 * @note We use list_node_t instead of list_t here to avoid the overhead
 * of searching for tails (which would be necessary to maintain the list_t
 * invariants).
 *
 * @warn This function doesn't check for circular structure.
 */
static list_node_t*
__merge_sort (list_node_t* L, int (*pred) (void*, void*))
{
  const int len = __length (L);
  if (len < 2)
    return L; /* A list of zero or one element is already sorted */
  else
    {
      list_node_t* L1 = L;
      list_node_t* L2 = L;
      list_node_t* prev = NULL;
      int i;

      /* Set L2 to the start of the midpoint of the list */
      for (i = 0; i < len/2; i++)
	{
	  prev = L2;
	  L2 = L2->next;
	}

      /* Cut the list L in two parts, L1 and L2.  
       * prev is the tail of L1. */
      prev->next = NULL;

      /* Sort the sublists L1 and L2 */
      __merge_sort (L1, pred);
      __merge_sort (L2, pred);

      /* Merge the sublists and return the result */
      return __merge (L1, L2, pred);
    }
}



list_t
list_sort (list_t L, int (*pred) (void*, void*))
{
  assert (! circular_structure_p (L.head));
  L.head = __merge_sort (L.head, pred);
  L.tail = find_tail (L.head);
  return L;
}

alist_t
alist_sort (alist_t L, int (*pred) (void*, void*))
{
  L.head = __merge_sort (L.head, pred);
  return L;
}


/**
 * Searches for the given ITEM in the given list (which is represented
 * by its head L), using the predicate EQUAL to test for equality.  If
 * no items are in the list L that are EQUAL to ITEM, returns NULL,
 * else returns a pointer to the first node containing ITEM.
 */
static list_node_t*
__find (list_node_t* L, void* item, int (*equal) (void*, void*))
{
  if (L == NULL)
    return NULL;
  else
    {
      if (equal (L->item, item))
	return L;
      else
	return __find (L->next, item, equal);
    }
}

list_node_t* 
list_find (list_t L, void* item, int (*equal) (void*, void*))
{
  assert (! circular_structure_p (L.head));
  return __find (L.head, item, equal);
}

list_node_t*
alist_find (alist_t L, void* item, int (*equal) (void*, void*))
{
  return __find (L.head, item, equal);
}


void
list_set_car (list_t L, void* item)
{
  if (list_empty_p (L))
    {
      bebop_error ("list:empty", "Attempt to set CAR of an empty list");
      bebop_exit (EXIT_FAILURE);
    }
  else
    {
      list_node_t* node = L.head;
      node->item = item;
    }
}

void*
list_cadr (list_t L)
{
  return (L.head)->next->item;
}

void*
list_caddr (list_t L)
{
  return (L.head)->next->next->item;
}

void*
list_cadddr (list_t L)
{
  return (L.head)->next->next->next->item;
}
