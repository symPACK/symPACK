/*
	 Copyright (c) 2016 The Regents of the University of California,
	 through Lawrence Berkeley National Laboratory.  

   Author: Mathias Jacquelin
	 
   This file is part of symPACK. All rights reserved.

	 Redistribution and use in source and binary forms, with or without
	 modification, are permitted provided that the following conditions are met:

	 (1) Redistributions of source code must retain the above copyright notice, this
	 list of conditions and the following disclaimer.
	 (2) Redistributions in binary form must reproduce the above copyright notice,
	 this list of conditions and the following disclaimer in the documentation
	 and/or other materials provided with the distribution.
	 (3) Neither the name of the University of California, Lawrence Berkeley
	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
	 be used to endorse or promote products derived from this software without
	 specific prior written permission.

	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	 You are under no obligation whatsoever to provide any bug fixes, patches, or
	 upgrades to the features, functionality or performance of the source code
	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
	 available either publicly, or directly to Lawrence Berkeley National
	 Laboratory, without imposing a separate written license agreement for such
	 Enhancements, then you hereby grant the following license: a non-exclusive,
	 royalty-free perpetual license to install, use, modify, prepare derivative
	 works, incorporate into other computer software, distribute, and sublicense
	 such enhancements or derivative works thereof, in binary and source code form.
*/


#ifndef _INTTREE_IMPL_HPP_
#define _INTTREE_IMPL_HPP_




// A utility function to create a new Interval Search Tree Node
inline ITree::ITNode * ITree::newNode_(ITree::Interval & i)
{
  ITree::ITNode *temp = new ITree::ITNode;
  temp->i = new ITree::Interval(i);
  temp->max = i.high;
  temp->left = temp->right = NULL;
  temp->height = 1;

  temp->min = i.low;

  return temp;
};


inline   Int ITree::height_(ITree::ITNode *N)
{
  if (N == NULL)
    return 0;
  return N->height;
}

inline   Int ITree::max_(ITree::ITNode *N)
{
  if (N == NULL)
    return 0;
  return N->max;
}

inline   Int ITree::min_(ITree::ITNode *N)
{
  if (N == NULL)
    return 0;
  return N->min;
}





// A utility function to insert a new ITree::Interval Search Tree Node
// This is similar to BST Insert.  Here the low value of interval
// is used tomaintain BST property
inline   ITree::ITNode * ITree::insert_(ITree::ITNode *root, ITree::Interval & i)
{
  assert(i.low<=i.high);

  // Base case: Tree is empty, new node becomes root
  if (root == NULL)
    return newNode_(i);

  // Get low value of interval at root
  Int l = root->i->low;

  // If root's low value is smaller, then new interval goes to
  // left subtree
  if (i.low < l)
    root->left = insert_(root->left, i);
  // Else, new node goes to right subtree.
  else
    root->right = insert_(root->right, i);


  // Update the max value of this ancestor if needed
  if (root->max < i.high)
    root->max = i.high;

  // Update the min value of this ancestor if needed
  if (root->min > i.low)
    root->min = i.low;

  /* 2. Update height of this ancestor node */
  root->height = std::max(height_(root->left), height_(root->right)) + 1;

  return root;
}

inline   Int ITree::recomputeMax_(ITree::ITNode * root)
{
  root->max = root->i->high;

  if(root->left != NULL){
    recomputeMax_(root->left);
    if(root->max < max_(root->left)){
      root->max = max_(root->left);
    } 
  }

  if(root->right != NULL){
    recomputeMax_(root->right);
    if(root->max < max_(root->right)){
      root->max = max_(root->right);
    } 
  }
  return root->max;
}


inline   Int ITree::recomputeMinMax_(ITree::ITNode * root)
{
  root->max = root->i->high;
  root->min = root->i->low;

  if(root->left != NULL){
    recomputeMinMax_(root->left);
    Int maxL = max_(root->left);
    if(root->max < maxL){
      root->max = maxL;
    } 

    Int minL = min_(root->left);
    if(root->min > minL){
      root->min = minL;
    } 
  }

  if(root->right != NULL){
    recomputeMinMax_(root->right);
    Int maxL = max_(root->right);
    if(root->max < maxL){
      root->max = maxL;
    } 

    Int minL = min_(root->right);
    if(root->min > minL){
      root->min = minL;
    } 

  }
  return root->max;
}





// A utility function to check if given two intervals overlap
inline   bool ITree::doOVerlap_(const ITree::Interval &i1, const ITree::Interval &i2)
{
  if (i1.low <= i2.high && i2.low <= i1.high)
    return true;
  return false;
}

inline   bool ITree::doOVerlap_(const ITree::Interval &i1, const Int & low, const Int & high)
{
  if (i1.low <= high && low <= i1.high)
    return true;
  return false;
}



// The main function that searches a given interval i in a given
// Interval Tree.
inline   ITree::Interval * ITree::intervalSearch_(ITree::ITNode *root,const Int & begin,const Int &end)
//  ITree::Interval * ITree::intervalSearch_(ITree::ITNode *root, ITree::Interval &i)
{
  // Base Case, tree is empty
  if (root == NULL){
    return NULL;
  }

  // If given interval overlaps with root
  //if (doOVerlap_(*(root->i), begin,end)){
  //  return root->i;
  //}
  //without the function call
  if (root->i->low <= end && begin <= root->i->high){
    return root->i;
  }

  // If left child of root is present and max of left child is
  // greater than or equal to given interval, then i may
  // overlap with an interval is left subtree
  //if (root->left != NULL && root->left->max >= i.low)
  if (root->left != NULL && root->left->max >= begin){
    return intervalSearch_(root->left, begin,end);
  }
  //return intervalSearch_(root->left, i);

  // Else interval can only overlap with right subtree
  return intervalSearch_(root->right, begin,end);
  //return intervalSearch_(root->right, i);
}


inline   ITree::Interval * ITree::intervalSearch_(ITree::ITNode *root,const Int & begin,const Int &end, Interval * & closestR, Interval * & closestL)
//  ITree::Interval * ITree::intervalSearch_(ITree::ITNode *root, ITree::Interval &i)
{
  // Base Case, tree is empty
  if (root == NULL){
    return NULL;
  }

  // If given interval overlaps with root
  //if (doOVerlap_(*(root->i), begin,end)){
  //  return root->i;
  //}
  //without the function call
  if (root->i->low <= end && begin <= root->i->high){
    closestL = root->i;
    closestR = root->i;
    return root->i;
  }

  if(root->left!=NULL){
    if(closestL!=NULL){
      if( std::min(0,end - root->left->i->high) < std::min(0,end - closestL->high)){
        closestL = root->left->i;
      }
    }
    else{
      closestL = root->left->i;
    }

    if(closestR!=NULL){
      if( std::min(0,root->left->i->low - begin) < std::min(0,closestL->low - begin)){
        closestR = root->left->i;
      }
    }
    else{
      closestR = root->left->i;
    }
  }

  if(root->right!=NULL){
    if(closestL!=NULL){
      if( std::min(0,end - root->right->i->high) < std::min(0,end - closestL->high)){
        closestL = root->right->i;
      }
    }
    else{
      closestL = root->right->i;
    }

    if(closestR!=NULL){
      if( std::min(0,root->right->i->low - begin) < std::min(0,closestL->low - begin)){
        closestR = root->right->i;
      }
    }
    else{
      closestR = root->right->i;
    }
  }

  // If left child of root is present and max of left child is
  // greater than or equal to given interval, then i may
  // overlap with an interval in left subtree
  //if (root->left != NULL && root->left->max >= i.low)
  if (root->left != NULL && root->left->max >= begin && root->left->min <= end){
    return intervalSearch_(root->left, begin,end,closestR,closestL);
  }
  //return intervalSearch_(root->left, i);

  // Else interval can only overlap with right subtree
  if(root->right != NULL && root->right->max >= begin && root->right->min <= end){
    return intervalSearch_(root->right, begin,end,closestR,closestL);
  }

  //interval is outside the range
  return NULL;
  //return intervalSearch_(root->right, i);
}





inline   void ITree::inorder_(ITree::ITNode *root)
{
  if (root == NULL) return;

  logfileptr->OFS()<< " LEFT of "<< "[" << root->i->low << ", " << root->i->high << "]"<<": "<<std::endl;
  inorder_(root->left);

  logfileptr->OFS()<< " MIDDLE: "<<std::endl;
  logfileptr->OFS()<< "[" << root->i->low << ", " << root->i->high << "] on "<<root->i->block_idx
    << " max = " << root->max << std::endl;

  logfileptr->OFS()<< " RIGHT of "<< "[" << root->i->low << ", " << root->i->high << "]"<<": "<<std::endl;
  inorder_(root->right);
}


inline   Int ITree::getSize_(ITree::ITNode *root)
{
  if (root == NULL) return 0;

  Int size = getSize_(root->left);
  size += sizeof(*root);
  size += getSize_(root->right);
  return size;
}


//AVLITree



// A utility function to insert a new ITree::Interval Search Tree Node
// This is similar to BST Insert.  Here the low value of interval
// is used tomaintain BST property
inline   ITree::ITNode * AVLITree::insert_(ITree::ITNode *root, ITree::Interval & i)
{
  root = ITree::insert_(root,i);

  SYMPACK_TIMER_START(ITREE_BALANCE_AVL);
  /* 3. Get the balance factor of this ancestor node to check whether
     this node became unbalanced */
  Int balance = getBalance_(root);

  // If this node becomes unbalanced, then there are 4 cases

  // Left Left Case
  if (balance > 1 && i.low < root->left->i->low)
  {
    //logfileptr->OFS()<<" ROOT "; root->Dump(); logfileptr->OFS()<<std::endl;
    //logfileptr->OFS()<<" LEFT LEFT ROTATION "<<std::endl;
    root = rightRotate_(root);
  }
  // Right Right Case
  else if (balance < -1 && i.low > root->right->i->low)
  {
    //logfileptr->OFS()<<" ROOT "; root->Dump(); logfileptr->OFS()<<std::endl;
    //logfileptr->OFS()<<" RIGHT RIGHT ROTATION "<<std::endl;
    root = leftRotate_(root);
  }
  // Left Right Case
  else if (balance > 1 && i.low > root->left->i->low)
  {
    //logfileptr->OFS()<<" ROOT "; root->Dump(); logfileptr->OFS()<<std::endl;
    //logfileptr->OFS()<<" LEFT RIGHT ROTATION "<<std::endl;
    root->left =  leftRotate_(root->left);
    root = rightRotate_(root);
  }
  // Right Left Case
  else if (balance < -1 && i.low < root->right->i->low)
  {
    //logfileptr->OFS()<<" ROOT "; root->Dump(); logfileptr->OFS()<<std::endl;
    //logfileptr->OFS()<<" RIGHT LEFT ROTATION "<<std::endl;
    root->right = rightRotate_(root->right);
    root = leftRotate_(root);
  }
  //recompute max
  //recomputeMax_(root);
  recomputeMinMax_(root);
  SYMPACK_TIMER_STOP(ITREE_BALANCE_AVL);

  return root;
}


// A utility function to right rotate subtree rooted with y
// See the diagram given above.
inline   ITree::ITNode * AVLITree::rightRotate_(ITree::ITNode *y)
{
  ITree::ITNode *x = y->left;
  ITree::ITNode *T2 = x->right;

  // Perform rotation
  x->right = y;
  y->left = T2;

  // Update heights
  y->height = std::max(height_(y->left), height_(y->right))+1;
  x->height = std::max(height_(x->left), height_(x->right))+1;

  // Return new root
  return x;
}
inline   ITree::ITNode * AVLITree::leftRotate_(ITree::ITNode *x)
{
  ITree::ITNode *y = x->right;
  ITree::ITNode *T2 = y->left;

  // Perform rotation
  y->left = x;
  x->right = T2;

  //  Update heights
  x->height = std::max(height_(x->left), height_(x->right))+1;
  y->height = std::max(height_(y->left), height_(y->right))+1;

  // Return new root
  return y;
}

// Get Balance factor of node N
inline   Int AVLITree::getBalance_(ITree::ITNode *N)
{
  if (N == NULL)
    return 0;
  return height_(N->left) - height_(N->right);
}


//DSWITree

// Tree to Vine algorithm:  a "pseudo-root" is passed ---
// comparable with a dummy header for a linked list.
inline void DSWITree::tree_to_vine_( ITree::ITNode* root, Int &size )
{  
  ITree::ITNode* vineTail, *remainder, *tempPtr;

  vineTail = root;
  remainder = vineTail->right;
  size = 0;

  while ( remainder != NULL )
  {//If no leftward subtree, move rightward
    if ( remainder->left == NULL )
    {  vineTail = remainder;
      remainder = remainder->right;
      size++;
    }
    //    else eliminate the leftward subtree by rotations
    else  /* Rightward rotation */
    {  tempPtr = remainder->left;
      remainder->left = tempPtr->right;
      tempPtr->right = remainder;
      remainder = tempPtr;
      vineTail->right = tempPtr;
    }
  }
}

inline void DSWITree::compression_( ITree::ITNode* root, Int count )
{  
  ITree::ITNode * scanner, *child;
  Int     j;

  scanner = root;
  for ( j = 0; j < count; j++ )
  {//Leftward rotation
    child = scanner->right;
    scanner->right = child->right;
    scanner = scanner->right;
    child->right = scanner->left;
    scanner->left = child;
  }  // end for
}  // end compression

// Code added by Tim Rolfe:  Expands on Warren & Stout's
// notation involving powers, floors, and base-2 logs
inline Int DSWITree::fullSize_( Int size )    // Full portion complete tree
{  Int Rtn = 1;
  while ( Rtn <= size )      // Drive one step PAST FULL
    Rtn = Rtn + Rtn + 1;   // next pow(2,k)-1
  return Rtn/2;
}

inline void DSWITree::vine_to_tree_ ( ITree::ITNode * root, Int size )
{
  Int full_count = fullSize_(size);
  compression_(root, size - full_count);
  for ( size = full_count ; size > 1 ; size /= 2 )
    compression_( root, size / 2 );
}

// Traverse entire tree, correcting heights and parents
  inline void DSWITree::correctTree_( ITree::ITNode* node )
{  if ( node != NULL )
  {  Int LtHt, RtHt;

    correctTree_ (node->left);
    correctTree_ (node->right);
    LtHt = node->left  ? node->left->height  : 0;
    RtHt = node->right ? node->right->height : 0;
    node->height = 1 + std::max( LtHt, RtHt );
  }
}

inline void DSWITree::Rebalance()
  // Public member function:  Do the DSW algorithm to balance the tree
{//Declare as automatic variable; remember to pass as pointer
  SYMPACK_TIMER_START(ITREE_BALANCE_DSW);
  //   BaseCell pseudo_root( -1, NULL, NULL, Root );

  Int size;

  // Stout/Warren transformation of tree to vine
  tree_to_vine_ (root_, size);

  vine_to_tree_ (root_, size);

  correctTree_ (root_->right);
  //recomputeMax_(root_);
  recomputeMinMax_(root_);
  SYMPACK_TIMER_STOP(ITREE_BALANCE_DSW);
}


#endif // _INTTREE_IMPL_HPP_
