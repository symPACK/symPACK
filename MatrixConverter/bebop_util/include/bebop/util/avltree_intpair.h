#ifndef _avl_tree_intpair_h
#define _avl_tree_intpair_h
/**
 * @file avltree_intpair.h
 * @author Mark Hoemmen
 * @since 23 Mar 2004
 * @date Time-stamp: <2008-07-16 10:14:14 mhoemmen>
 *
 * @brief An AVL tree holding pairs of integers.
 *
 * An AVL tree holding pairs of integers.  Supplied with a dictionary
 * interface.  Borrowed from the Weiss data structures book web site (which
 * see) and modified to store pairs of integers with a lexicographic
 * comparison.
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
 ****************************************************************************/

/** A pair of integers. */
struct int_pair
{
  /** First element in the pair. */
  int first;  
  /** second element in the pair. */
  int second; 
};

/** 
 * The element in each avl_node_intpair is an integer pair.  
 *
 * @see avl_node_intpair 
 */
typedef struct int_pair element_type;

/**
 * @struct avl_node_intpair
 *
 * Node of an AVL tree containing pairs of integers.
 */
struct avl_node_intpair
{
  /** The left subtree. */
  struct avl_node_intpair* left;  

  /** The right subtree. */
  struct avl_node_intpair* right; 

  /** The datum stored in a node. */
  element_type element;   

  /** Tree height of this node.   */
  int height;             
};


typedef struct avl_node_intpair *avl_position_intpair;
typedef struct avl_node_intpair *avl_tree_intpair;

/**
 * Clears the contents of T (if there are any), and returns NULL, 
 * the empty tree.
 */
avl_tree_intpair
make_empty_avl_tree_intpair (avl_tree_intpair T);

/**
 * Returns NULL if X is not in the tree T.  Otherwise, returns
 * the position of X in T.
 */
avl_position_intpair
find_intpair (element_type X, avl_tree_intpair T);

/**
 * Returns NULL if T is empty, otherwise returns the position of the
 * minimum element in the tree.
 */
avl_position_intpair
find_min_intpair (avl_tree_intpair T);

/**
 * Returns NULL if T is empty, otherwise returns the position of the
 * maximum element in the tree.
 */
avl_position_intpair
find_max_intpair (avl_tree_intpair T);

/**
 * `T = Insert (X,T)' is the correct usage.  Inserts X into T.
 */
avl_tree_intpair
insert_intpair (element_type X, avl_tree_intpair T);

/**
 * Returns the value stored at position p in an avl_tree_intpair.
 */
element_type 
retrieve_intpair (avl_position_intpair p);

/**
 * Type of a function used to implement visiting each node, and doing
 * something at that node.  The function takes two arguments: first, the
 * current subtree, and second, a void pointer to a working set of data.  Your
 * function is responsible for interpreting the data correctly (this is an
 * unfortunate consequence of using a language like C, which is not
 * type-safe).  The function returns nothing.
 * 
 * @note For programmers not used to function pointers: A typedef for a
 * function pointer looks different than a usual typedef: the type's name is
 * `visit_function'.
 * 
 * @warn The structure of the tree must not be modified during traversal, or 
 * results are undefined.
 */
typedef  void(*visit_function)  (avl_tree_intpair, void*);


/**
 * A simple visiting function that just prints out the value at the node.
 * Useful for testing.  Pass in NULL (or anything) for data (no reads or
 * writes to data).
 */
void 
test_visit_avl_tree_intpair (avl_tree_intpair T, void* data);


/**
 * Traverses the tree in order, visiting each node using the given visit
 * function and using the given data during the visit.  Your visit function is
 * responsible for interpreting the data correctly.
 * 
 * @warn The structure of the tree must not be modified during traversal, or 
 * results are undefined.
 *
 * @param T       Current subtree.
 * @param visit   Function that `visits' (does operations on) the node.
 * @param data    Data to pass into the visit function.
 * 
 * @see visit_function
 */
void
traverse_in_order_avl_tree_intpair (avl_tree_intpair T, visit_function visit, void* data);


#endif  /* _avl_tree_intpair_h */

