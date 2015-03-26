#ifndef _avlstree_h
#define _avlstree_h
/**
 * @file avlstree.h
 * @author Mark Hoemmen
 * @since 25 Aug 2006 
 * @date Time-stamp: <2008-07-16 10:14:03 mhoemmen>
 *
 * @brief An AVL tree with strings as keys and arbitrary values.
 *
 * An AVL tree (sorted balanced binary tree) holding strings (stored
 * as char*) as keys and arbitrary values associated with the keys.
 * Supplied with a dictionary interface.  Borrowed from the Weiss data
 * structures book web site (which see) and modified to store string
 * keys instead of integers.
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


struct __avlsnode 
{
  struct __avlsnode* left;
  struct __avlsnode* right;
  char* key;
  void* value;
  int height;
};
typedef struct __avlsnode avlsnode_t;

typedef avlsnode_t* avlspos_t;
typedef avlsnode_t* avlstree_t;

void
avlstree_delete (avlstree_t T);

/**
 * Clears the contents of T (if there are any), and returns NULL, 
 * the empty tree.
 */
avlstree_t
make_avlstree ();

/**
 * Returns NULL if the given key is not in the tree T.  Otherwise,
 * returns the position of the key in T.
 */
avlspos_t
avlstree_find (const char* key, avlstree_t T);

/**
 * Returns NULL if T is empty, otherwise returns the position of the
 * minimum element in the tree.
 */
avlspos_t
avlstree_findmin (avlstree_t T);

/**
 * Returns NULL if T is empty, otherwise returns the position of the
 * maximum element in the tree.
 */
avlspos_t
avlstree_findmax (avlstree_t T);

/**
 * `T = Insert (X,T)' is the correct usage.  Inserts X into T.
 */
avlstree_t
avlstree_insert (const char* key, void* value, avlstree_t T);

/**
 * Returns the value stored at position p in an avl_tree_string.
 */
void*
avlspos_value (avlspos_t p);

char*
avlspos_key (avlspos_t p);

/**
 * Type of a function used to implement visiting each node, and doing
 * something at that node.  The function takes two arguments: first, the
 * current subtree, and second, a void pointer to a working set of data.  Your
 * function is responsible for interpreting the data correctly (this is an
 * unfortunate consequence of using a language like C, which is not
 * type-safe).  The function returns nothing.
 * 
 * @warn The structure of the tree must not be modified during traversal, or 
 * results are undefined.
 */
typedef void (*avlstree_visitfn_t) (avlstree_t, void*);


/**
 * A simple visiting function that just prints out the value at the node.
 * Useful for testing.  Pass in NULL (or anything) for data (no reads or
 * writes to data).
 */
void 
avlstree_testvisit (avlstree_t T, void* data);

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
avlstree_traverse_in_order (avlstree_t T, avlstree_visitfn_t visit, void* data);


#endif  /* _avl_tree_string_h */

