/**
 * @file avltree_intpair.c
 * @author Mark Hoemmen
 * @since 2004
 * @date Time-stamp: <2008-07-16 10:09:10 mhoemmen>
 *
 * Implementation of an AVL tree holding pairs of integers.
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
#include <bebop/util/avltree_intpair.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>

#include <stdlib.h>


/**
 * Lexicographical comparison function for a pair of ints.  
 * Follows the name return convention as (e.g.) strcmp.
 *
 * @return  -1 if p1 < p2, 0 if p1 == p2, and +1 if p1 > p2.
 */
static int
compare (struct int_pair p1, struct int_pair p2)
{
  if (p1.first < p2.first)
    return -1;
  else if (p1.first == p2.first)
    {
      if (p1.second < p2.second)
	return -1;
      else if (p1.second == p2.second)
	return 0;
      else  /* p1.second > p2.second */
	return 1;
    }
  /* else */  /* Compilers might whine about "no return statement" */
    return 1;
}


/**
 * Uses the above lexicographical comparison to find the max of two int pairs.
 */
/*
static struct int_pair
max (struct int_pair lhs, struct int_pair rhs)
{
  if (compare (lhs,rhs) < 0)
    return rhs;
  else
    return lhs;
}
*/
/* The standard int max. */
static int
max (int lhs, int rhs)
{
  return lhs < rhs ? rhs : lhs;
}

/**
 * Convenient shorthand. 
 */
static element_type
element (avl_position_intpair p)
{
  return p->element;
}




/*******************************************************
 * Now follows the actual implementation of the tree.  *
 *******************************************************/

avl_tree_intpair
make_empty_avl_tree_intpair (avl_tree_intpair T)
{
  if (T != NULL)
    {
      make_empty_avl_tree_intpair (T->left);
      make_empty_avl_tree_intpair (T->right);
      bebop_free (T);
    }
  return NULL;
}

avl_position_intpair
find_intpair (element_type X, avl_tree_intpair T)
{
  if (T == NULL)
	return NULL;

  /* if (X < element (T)) */
  if (compare (X, element (T)) < 0)
	return find_intpair (X, T->left);
  else
	{
	  /* if (X > element (T)) */
	  if (compare (X, element (T)) > 0)
		return find_intpair (X, T->right);
	  else
		return T;
	}
}

avl_position_intpair
find_min_intpair (avl_tree_intpair T)
{
  if (T == NULL)
	return NULL;
  else
	{
	  if (T->left == NULL)
		return T;
	  else
		return find_min_intpair (T->left);
	}
}

avl_position_intpair
find_max_intpair (avl_tree_intpair T)
{
  if (T != NULL)
    while (T->right != NULL)
      T = T->right;
  
  return T;
}

static int
height (avl_position_intpair P)
{
  if (P == NULL)
    return -1;
  else
    return P->height;
}


/**
 * Should only be called if K2 has a left child.  Do a rotation
 * between K2 and its left child.  Update heights, and return 
 * the new root.
 */
static avl_position_intpair
single_rotate_with_left (avl_position_intpair K2)
{
  avl_position_intpair K1;
  
  K1 = K2->left;
  K2->left = K1->right;
  K1->right = K2;
  
  K2->height = max (height (K2->left), height (K2->right)) + 1;
  K1->height = max (height (K1->left), height (K2)) + 1;
  
  return K1;  /* New root */
}


/**
 * Should be called only if K1 has a right child.  Do a rotate
 * between K1 and its right child.  Update heights and return 
 * new root.
 */
static avl_position_intpair
single_rotate_with_right (avl_position_intpair K1)
{
  avl_position_intpair K2;
  
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
static avl_position_intpair
double_rotate_with_left (avl_position_intpair K3)
{
  /* Rotate between K1 and K2 */
  K3->left = single_rotate_with_right (K3->left);

  /* Rotate between K3 and K2 */
  return single_rotate_with_left (K3);
}


/**
 * Should be called only if K1 has a right child, and K1's right child
 * has a left child.  Do the right-left double rotation.  Update heights,
 * and return new root.
 */
static avl_position_intpair
double_rotate_with_right (avl_position_intpair K1)
{
  /* Rotate between K3 and K2 */
  K1->right = single_rotate_with_left (K1->right);

  /* Rotate between K1 and K2 */
  return single_rotate_with_right (K1);
}


avl_tree_intpair
insert_intpair (element_type X, avl_tree_intpair T)
{
  if (T == NULL)
	{
	  /* Create and return a one-node tree */
	  T = bebop_malloc (sizeof (struct avl_node_intpair));
	  T->element = X; 
	  T->height  = 0;
	  T->left = T->right = NULL;
	}
  else
	/* if (X < T->element) */
	if (compare (X, element (T)) < 0)
	  {
		T->left = insert_intpair (X, T->left);
		if (height (T->left) - height (T->right) == 2)
		  {
			/* if (X < element (T->left)) */
			if (compare (X, element (T->left)) < 0)
			  T = single_rotate_with_left (T);
			else
			  T = double_rotate_with_left (T);
		  }
	  }
	else
	  /* if (X > element (T)) */
	  if (compare (X, element (T)) > 0)
		{
		  T->right = insert_intpair (X, T->right);
		  if (height (T->right) - height (T->left) == 2)
			{
			  /* if (X > element (T->right)) */
			  if (compare (X, element (T->right)) > 0)
				T = single_rotate_with_right (T);
			  else
				T = double_rotate_with_right (T);
			}
		}
  /* Else X is in the tree already; we'll do nothing */

  T->height = max (height (T->left), height (T->right)) + 1;
  return T;
}



element_type
retrieve_intpair (avl_position_intpair P)
{
  return P->element;
}


void 
test_visit_avl_tree_intpair (avl_tree_intpair T, void* data)
{
  if (T == NULL)
    return;  /* Sanity check (traversal fn should not visit NULL nodes) */

  printf ("(%d,%d)\n", (T->element).first, (T->element).second);
}



void
traverse_in_order_avl_tree_intpair (avl_tree_intpair T, 
				    visit_function visit, void* data)
{
  if (T == NULL)
    return;  /* base case */

  traverse_in_order_avl_tree_intpair (T->left, visit, data);
  visit (T, data);
  traverse_in_order_avl_tree_intpair (T->right, visit, data);
}
