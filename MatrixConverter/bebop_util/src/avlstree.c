/**
 * @file avlstree.c
 * @author Mark Hoemmen
 * @since 25 Aug 2006 
 * @date Time-stamp: <2008-07-16 10:08:56 mhoemmen>
 *
 * Implementation of an AVL tree with strings as keys and arbitrary
 * values.
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
#include <bebop/util/config.h>
#include <bebop/util/avlstree.h>
#include <bebop/util/malloc.h>
#include <bebop/util/string.h>
#include <bebop/util/util.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>


static inline int
max (const int a, const int b)
{
  return (a > b) ? a : b;
}

static inline int
min (const int a, const int b)
{
  return (a < b) ? a : b;
}

static inline const char* 
KEY (avlspos_t p)
{ 
  return p->key;
}

static inline void*
VALUE (avlspos_t p)
{
  return p->value;
}


char*
avlspos_key (avlspos_t p)
{
  return p->value;
}

void*
avlspos_value (avlspos_t p)
{
  return p->value;
}



static inline void
avlsnode_delete (avlsnode_t* node)
{
  bebop_free (node->key);
  bebop_free (node->value);
  bebop_free (node);
}

void
avlstree_delete (avlstree_t T)
{
  if (T != NULL)
    {
      avlstree_delete (T->left);
      avlstree_delete (T->right);
      avlsnode_delete ((avlsnode_t*) T);
    }
}

avlstree_t
make_avlstree ()
{
  return (avlstree_t) NULL;
}


avlspos_t
avlstree_find (const char* key, avlstree_t T)
{
  int c;
 
  /* C compilers aren't always so clever at unrolling recursions
     automatically, so we unrolled it ourselves. */
  while (T != NULL)
    {
      c = strcmp (key, KEY (T));
      if (c < 0)
	T = T->left;
      else if (c > 0)
	T = T->right;
      else
	return (avlspos_t) T;
    }
  return (avlspos_t) NULL;
}

avlspos_t
avlstree_findmin (avlstree_t T)
{
  if (T != NULL)
    while (T->left != NULL)
      T = T->left;

  return (avlspos_t) T;
}

avlspos_t
avlstree_findmax (avlstree_t T)
{
  if (T != NULL)
    while (T->right != NULL)
      T = T->right;
  
  return (avlspos_t) T;
}

static inline int
height (avlspos_t P)
{
  if (P == NULL)
    return -1;
  else
    return P->height;
}


/**
 * Should only be called if K2 has a left child.  Do a rotation
 * between K2 and its left child.  Update heights, and return the new
 * root.
 */
static avlspos_t
single_rotate_with_left (avlspos_t K2)
{
  avlspos_t K1;
  
  K1 = K2->left;
  K2->left = K1->right;
  K1->right = K2;
  
  K2->height = max (height (K2->left), height (K2->right)) + 1;
  K1->height = max (height (K1->left), height (K2)) + 1;
  
  return K1;  /* New root */
}


/**
 * Should be called only if K1 has a right child.  Do a rotate between
 * K1 and its right child.  Update heights and return new root.
 */
static avlspos_t
single_rotate_with_right (avlspos_t K1)
{
  avlspos_t K2;
  
  K2 = K1->right;
  K1->right = K2->left;
  K2->left = K1;
  
  K1->height = max (height (K1->left), height (K1->right)) + 1;
  K2->height = max (height (K2->right), height (K1)) + 1;
  
  return K2;  /* New root */
}



/**
 * Should be called only if K3 has a left child and K3's left child
 * has a right child.  Do the left-right double rotation.  Update 
 * heights, and return new root.
 */
static avlspos_t 
double_rotate_with_left (avlspos_t K3)
{
  /* Rotate between K1 and K2 */
  K3->left = single_rotate_with_right (K3->left);

  /* Rotate between K3 and K2 */
  return single_rotate_with_left (K3);
}


/**
 * Should be called only if K1 has a right child, and K1's right child
 * has a left child.  Do the right-left double rotation.  Update
 * heights, and return new root.
 */
static avlspos_t
double_rotate_with_right (avlspos_t K1)
{
  /* Rotate between K3 and K2 */
  K1->right = single_rotate_with_left (K1->right);

  /* Rotate between K1 and K2 */
  return single_rotate_with_right (K1);
}


static inline avlsnode_t* 
make_avlsnode (const char* key, void* value, int height, 
	       avlsnode_t* left, avlsnode_t* right)
{
  avlsnode_t* T = bebop_malloc (sizeof (avlsnode_t));
 
  T->key = bebop_strdup (key);
  T->value = value;
  T->height = height;
  T->left = left;
  T->right = right;

  return T;
}


avlstree_t
avlstree_insert (const char* key, void* value, avlstree_t T)
{
  if (T == NULL)
    {
      /* Create and return a one-node tree */
      T = make_avlsnode (key, value, 0, NULL, NULL);
    }
  else
    {
      const int c = strcmp (key, KEY (T));
   
      /* if (key < T->key) */
      if (c < 0)
	{
	  T->left = avlstree_insert (key, value, T->left);
	  if (height (T->left) - height (T->right) == 2)
	    {
	      /* if (key < (T->left)->key) */
	      if (strcmp (key, KEY (T->left)) < 0)
		T = single_rotate_with_left (T);
	      else
		T = double_rotate_with_left (T);
	    }
	  }
      else if (c > 0) /* if (key > T->key) */
	{
	  T->right = avlstree_insert (key, value, T->right);
	  if (height (T->right) - height (T->left) == 2)
	    {
	      /* if (key > KEY (T->right)) */
	      if (strcmp (key, KEY (T->right)) > 0)
		T = single_rotate_with_right (T);
	      else
		T = double_rotate_with_right (T);
	    }
	}
    }
  /* Else X is in the tree already; we'll do nothing */

  T->height = max (height (T->left), height (T->right)) + 1;
  return T;
}


void 
avlstree_testvisit (avlstree_t T, void* data)
{
  assert (T != NULL);
  printf ("%s\n", KEY (T));
}



void
avlstree_traverse_in_order (avlstree_t T, 
			    avlstree_visitfn_t visit, 
			    void* data)
{
  if (T == NULL)
    return;  /* base case */

  avlstree_traverse_in_order (T->left, visit, data);
  visit (T, data);
  avlstree_traverse_in_order (T->right, visit, data);
}
