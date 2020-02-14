#ifndef _INTTREE_IMPL_HPP_
#define _INTTREE_IMPL_HPP_




// A utility function to create a new Interval Search Tree Node
  template<typename F>
inline typename ITree<F>::template ITNode<F> * ITree<F>::newNode_(typename ITree<F>::template Interval<F> & i)
{
  ITree<F>::ITNode<F> *temp = new typename ITree<F>::template ITNode<F>;
  temp->i = new typename ITree<F>::template Interval<F>(i);
  temp->max = i.high;
  temp->left = temp->right = nullptr;
  temp->height = 1;

  temp->min = i.low;

  return temp;
}


  template<typename F>
inline   Int ITree<F>::height_(ITree<F>::ITNode<F> *N)
{
  if (N == nullptr)
    return 0;
  return N->height;
}

  template<typename F>
inline   Int ITree<F>::max_(ITree<F>::ITNode<F> *N)
{
  if (N == nullptr)
    return 0;
  return N->max;
}

  template<typename F>
inline   Int ITree<F>::min_(ITree<F>::ITNode<F> *N)
{
  if (N == nullptr)
    return 0;
  return N->min;
}





// A utility function to insert a new ITree::Interval Search Tree Node
// This is similar to BST Insert.  Here the low value of interval
// is used tomaintain BST property
  template<typename F>
inline typename  ITree<F>::template ITNode<F> * ITree<F>::insert_(typename ITree<F>::template ITNode<F> *root, typename ITree<F>::template Interval<F> & i)
{
  assert(i.low<=i.high);

  // Base case: Tree is empty, new node becomes root
  if (root == nullptr)
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
  root->max = std::max( root->max, i.high);

  // Update the min value of this ancestor if needed
  root->min = std::min( root->min, i.low);

  /* 2. Update height of this ancestor node */
  root->height = std::max(height_(root->left), height_(root->right)) + 1;

  return root;
}

  template<typename F>
inline   Int ITree<F>::recomputeMax_(ITree<F>::ITNode<F> * root)
{
  root->max = root->i->high;

  if(root->left != nullptr){
    recomputeMax_(root->left);
    root->max = std::max( root->max, max_(root->left));
  }

  if(root->right != nullptr){
    recomputeMax_(root->right);
    root->max = std::max( root->max, max_(root->right));
  }
  return root->max;
}


  template<typename F>
inline   Int ITree<F>::recomputeMinMax_(ITree<F>::ITNode<F> * root)
{
  root->max = root->i->high;
  root->min = root->i->low;

  if(root->left != nullptr){
    recomputeMinMax_(root->left);
    Int maxL = max_(root->left);
    root->max = std::max( root->max, maxL);

    Int minL = min_(root->left);
    root->min = std::min( root->min, minL);
  }

  if(root->right != nullptr){
    recomputeMinMax_(root->right);
    Int maxL = max_(root->right);
    root->max = std::max( root->max, maxL);

    Int minL = min_(root->right);
    root->min = std::min( root->min, minL);
  }
  return root->max;
}





// A utility function to check if given two intervals overlap
  template<typename F>
inline   bool ITree<F>::doOVerlap_(const ITree<F>::Interval<F> &i1, const ITree<F>::Interval<F> &i2)
{
  if (i1.low <= i2.high && i2.low <= i1.high)
    return true;
  return false;
}

  template<typename F>
inline   bool ITree<F>::doOVerlap_(const ITree<F>::Interval<F> &i1, const Int & low, const Int & high)
{
  if (i1.low <= high && low <= i1.high)
    return true;
  return false;
}



// The main function that searches a given interval i in a given
// Interval Tree.
  template<typename F>
inline typename ITree<F>::template Interval<F> * ITree<F>::intervalSearch_(typename ITree<F>::template ITNode<F> *root,const Int & begin,const Int &end)
{
  // Base Case, tree is empty
  if (root == nullptr){
    return nullptr;
  }

  //without the function call
  if (root->i->low <= end && begin <= root->i->high){
    return root->i;
  }

  // If left child of root is present and max of left child is
  // greater than or equal to given interval, then i may
  // overlap with an interval is left subtree
  if (root->left != nullptr && root->left->max >= begin){
    return intervalSearch_(root->left, begin,end);
  }

  // Else interval can only overlap with right subtree
  return intervalSearch_(root->right, begin,end);
}


  template<typename F>
inline typename ITree<F>::template Interval<F> * ITree<F>::intervalSearch_(typename ITree<F>::template ITNode<F> *root,const Int & begin,const Int &end, Interval<F> * & closestR, Interval<F> * & closestL)
{
  // Base Case, tree is empty
  if (root == nullptr){
    return nullptr;
  }

  //without the function call
  if (root->i->low <= end && begin <= root->i->high){
    closestL = root->i;
    closestR = root->i;
    return root->i;
  }

  if(root->left!=nullptr){
    if(closestL!=nullptr){
      if( std::min((Int)0,end - root->left->i->high) < std::min((Int)0,end - closestL->high)){
        closestL = root->left->i;
      }
    }
    else{
      closestL = root->left->i;
    }

    if(closestR!=nullptr){
      if( std::min((Int)0,root->left->i->low - begin) < std::min((Int)0,closestL->low - begin)){
        closestR = root->left->i;
      }
    }
    else{
      closestR = root->left->i;
    }
  }

  if(root->right!=nullptr){
    if(closestL!=nullptr){
      if( std::min((Int)0,end - root->right->i->high) < std::min((Int)0,end - closestL->high)){
        closestL = root->right->i;
      }
    }
    else{
      closestL = root->right->i;
    }

    if(closestR!=nullptr){
      if( std::min((Int)0,root->right->i->low - begin) < std::min((Int)0,closestL->low - begin)){
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
  if (root->left != nullptr && root->left->max >= begin && root->left->min <= end){
    return intervalSearch_(root->left, begin,end,closestR,closestL);
  }

  // Else interval can only overlap with right subtree
  if(root->right != nullptr && root->right->max >= begin && root->right->min <= end){
    return intervalSearch_(root->right, begin,end,closestR,closestL);
  }

  //interval is outside the range
  return nullptr;
}





  template<typename F>
inline   void ITree<F>::inorder_(ITree<F>::ITNode<F> *root)
{
  if (root == nullptr) return;

  logfileptr->OFS()<< " LEFT of "<< "[" << root->i->low << ", " << root->i->high << "]"<<": "<<std::endl;
  inorder_(root->left);

  logfileptr->OFS()<< " MIDDLE: "<<std::endl;
  logfileptr->OFS()<< "[" << root->i->low << ", " << root->i->high << "] on "<<root->i->data
    << " max = " << root->max << std::endl;

  logfileptr->OFS()<< " RIGHT of "<< "[" << root->i->low << ", " << root->i->high << "]"<<": "<<std::endl;
  inorder_(root->right);
}


  template<typename F>
inline   Int ITree<F>::getSize_(ITree<F>::ITNode<F> *root)
{
  if (root == nullptr) return 0;

  Int size = getSize_(root->left);
  size += sizeof(*root);
  size += getSize_(root->right);
  return size;
}


//AVLITree



// A utility function to insert a new ITree::Interval Search Tree Node
// This is similar to BST Insert.  Here the low value of interval
// is used tomaintain BST property
  template<typename F>
inline   typename ITree<F>::template ITNode<F> * AVLITree<F>::insert_(typename ITree<F>::template ITNode<F> *root, typename ITree<F>::template Interval<F> & i)
{
  root = ITree<F>::insert_(root,i);

  SYMPACK_TIMER_START(ITREE_BALANCE_AVL);
  /* 3. Get the balance factor of this ancestor node to check whether
     this node became unbalanced */
  Int balance = getBalance_(root);

  // If this node becomes unbalanced, then there are 4 cases

  // Left Left Case
  if (balance > 1 && i.low < root->left->i->low)
  {
    root = rightRotate_(root);
  }
  // Right Right Case
  else if (balance < -1 && i.low > root->right->i->low)
  {
    root = leftRotate_(root);
  }
  // Left Right Case
  else if (balance > 1 && i.low > root->left->i->low)
  {
    root->left =  leftRotate_(root->left);
    root = rightRotate_(root);
  }
  // Right Left Case
  else if (balance < -1 && i.low < root->right->i->low)
  {
    root->right = rightRotate_(root->right);
    root = leftRotate_(root);
  }
  recomputeMinMax_(root);
  SYMPACK_TIMER_STOP(ITREE_BALANCE_AVL);

  return root;
}


// A utility function to right rotate subtree rooted with y
// See the diagram given above.
  template<typename F>
inline   typename ITree<F>::template ITNode<F> * AVLITree<F>::rightRotate_(typename ITree<F>::template ITNode<F> *y)
{
  typename ITree<F>::template ITNode<F> *x = y->left;
  typename ITree<F>::template ITNode<F> *T2 = x->right;

  // Perform rotation
  x->right = y;
  y->left = T2;

  // Update heights
  y->height = std::max(height_(y->left), height_(y->right))+1;
  x->height = std::max(height_(x->left), height_(x->right))+1;

  // Return new root
  return x;
}

  template<typename F>
inline   typename ITree<F>::template ITNode<F> * AVLITree<F>::leftRotate_(typename ITree<F>::template ITNode<F> *x)
{
  typename ITree<F>::template ITNode<F> *y = x->right;
  typename ITree<F>::template ITNode<F> *T2 = y->left;

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
  template<typename F>
inline   Int AVLITree<F>::getBalance_(typename ITree<F>::template ITNode<F> *N)
{
  if (N == nullptr)
    return 0;
  return height_(N->left) - height_(N->right);
}


//DSWITree

// Tree to Vine algorithm:  a "pseudo-root" is passed ---
// comparable with a dummy header for a linked list.
  template<typename F>
inline void DSWITree<F>::tree_to_vine_( typename ITree<F>::template ITNode<F>* root, Int &size )
{  
  typename ITree<F>::template ITNode<F>* vineTail, *remainder, *tempPtr;

  vineTail = root;
  remainder = vineTail->right;
  size = 0;

  while ( remainder != nullptr )
  {//If no leftward subtree, move rightward
    if ( remainder->left == nullptr )
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

  template<typename F>
inline void DSWITree<F>::compression_( typename ITree<F>::template ITNode<F>* root, Int count )
{  
  typename ITree<F>::template ITNode<F> * scanner, *child;
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
  template<typename F>
inline Int DSWITree<F>::fullSize_( Int size )    // Full portion complete tree
{  Int Rtn = 1;
  while ( Rtn <= size )      // Drive one step PAST FULL
    Rtn = Rtn + Rtn + 1;   // next pow(2,k)-1
  return Rtn/2;
}

  template<typename F>
inline void DSWITree<F>::vine_to_tree_ ( typename ITree<F>::template ITNode<F> * root, Int size )
{
  Int full_count = fullSize_(size);
  compression_(root, size - full_count);
  for ( size = full_count ; size > 1 ; size /= 2 )
    compression_( root, size / 2 );
}

// Traverse entire tree, correcting heights and parents
  template<typename F>
  inline void DSWITree<F>::correctTree_( typename ITree<F>::template ITNode<F>* node )
{  if ( node != nullptr )
  {  Int LtHt, RtHt;

    correctTree_ (node->left);
    correctTree_ (node->right);
    LtHt = node->left  ? node->left->height  : 0;
    RtHt = node->right ? node->right->height : 0;
    node->height = 1 + std::max( LtHt, RtHt );
  }
}

  template<typename F>
inline void DSWITree<F>::Rebalance()
  // Public member function:  Do the DSW algorithm to balance the tree
{//Declare as automatic variable; remember to pass as pointer
  SYMPACK_TIMER_START(ITREE_BALANCE_DSW);

  Int size;

  // Stout/Warren transformation of tree to vine
  this->tree_to_vine_ (this->root_, size);

  this->vine_to_tree_ (this->root_, size);

  this->correctTree_ (this->root_->right);
  this->recomputeMinMax_(this->root_);
  SYMPACK_TIMER_STOP(ITREE_BALANCE_DSW);
}


#endif // _INTTREE_IMPL_HPP_
