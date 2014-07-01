#include "IntervalTree.hpp"



















namespace LIBCHOLESKY{


  // A utility function to create a new Interval Search Tree Node
  ITree::ITNode * ITree::newNode_(ITree::Interval & i)
  {
    ITree::ITNode *temp = new ITree::ITNode;
    temp->i = new ITree::Interval(i);
    temp->max = i.high;
    temp->left = temp->right = NULL;
    temp->height = 1;
    return temp;
  };


  Int ITree::height_(ITree::ITNode *N)
  {
    if (N == NULL)
      return 0;
    return N->height;
  }

  Int ITree::max_(ITree::ITNode *N)
  {
    if (N == NULL)
      return 0;
    return N->max;
  }


  // A utility function to right rotate subtree rooted with y
  // See the diagram given above.
  ITree::ITNode * ITree::rightRotate_(ITree::ITNode *y)
  {
    ITree::ITNode *x = y->left;
    ITree::ITNode *T2 = x->right;

    // Perform rotation
    x->right = y;
    y->left = T2;

    // Update heights
    y->height = max(height_(y->left), height_(y->right))+1;
    x->height = max(height_(x->left), height_(x->right))+1;

    // Return new root
    return x;
  }
  ITree::ITNode * ITree::leftRotate_(ITree::ITNode *x)
  {
    ITree::ITNode *y = x->right;
    ITree::ITNode *T2 = y->left;

    // Perform rotation
    y->left = x;
    x->right = T2;

    //  Update heights
    x->height = max(height_(x->left), height_(x->right))+1;
    y->height = max(height_(y->left), height_(y->right))+1;

    // Return new root
    return y;
  }

  // Get Balance factor of node N
  Int ITree::getBalance_(ITree::ITNode *N)
  {
    if (N == NULL)
      return 0;
    return height_(N->left) - height_(N->right);
  }


  // A utility function to insert a new ITree::Interval Search Tree Node
  // This is similar to BST Insert.  Here the low value of interval
  // is used tomaintain BST property
  ITree::ITNode * ITree::insert_(ITree::ITNode *root, ITree::Interval & i)
  {
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




    /* 2. Update height of this ancestor node */
    root->height = max(height_(root->left), height_(root->right)) + 1;

    /* 3. Get the balance factor of this ancestor node to check whether
       this node became unbalanced */
    Int balance = getBalance_(root);

    // If this node becomes unbalanced, then there are 4 cases

    // Left Left Case
    if (balance > 1 && i.low < root->left->i->low)
      root = rightRotate_(root);

    // Right Right Case
    if (balance < -1 && i.low > root->right->i->low)
      root = leftRotate_(root);

    // Left Right Case
    if (balance > 1 && i.low > root->left->i->low)
    {
      root->left =  leftRotate_(root->left);
      root = rightRotate_(root);
    }

    // Right Left Case
    if (balance < -1 && i.low < root->right->i->low)
    {
      root->right = rightRotate_(root->right);
      root = leftRotate_(root);
    }
    //recompute max
    recomputeMax_(root);

    return root;
  }

  Int ITree::recomputeMax_(ITree::ITNode * root)
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


  // A utility function to check if given two intervals overlap
  bool ITree::doOVerlap_(ITree::Interval &i1, ITree::Interval &i2)
  {
    if (i1.low <= i2.high && i2.low <= i1.high)
      return true;
    return false;
  }

  bool ITree::doOVerlap_(ITree::Interval &i1, Int & low, Int & high)
  {
    if (i1.low <= high && low <= i1.high)
      return true;
    return false;
  }



  // The main function that searches a given interval i in a given
  // Interval Tree.
  ITree::Interval * ITree::intervalSearch_(ITree::ITNode *root, Int & begin, Int &end)
//  ITree::Interval * ITree::intervalSearch_(ITree::ITNode *root, ITree::Interval &i)
  {
    // Base Case, tree is empty
    if (root == NULL) return NULL;

    // If given interval overlaps with root
    //if (doOVerlap_(*(root->i), i))
    if (doOVerlap_(*(root->i), begin,end))
      return root->i;

    // If left child of root is present and max of left child is
    // greater than or equal to given interval, then i may
    // overlap with an interval is left subtree
    //if (root->left != NULL && root->left->max >= i.low)
    if (root->left != NULL && root->left->max >= begin)
      return intervalSearch_(root->left, begin,end);
      //return intervalSearch_(root->left, i);

    // Else interval can only overlap with right subtree
    return intervalSearch_(root->right, begin,end);
    //return intervalSearch_(root->right, i);
  }


  void ITree::inorder_(ITree::ITNode *root)
  {
    if (root == NULL) return;

    inorder_(root->left);

    logfileptr->OFS()<< "[" << root->i->low << ", " << root->i->high << "] on "<<root->i->block_idx
      << " max = " << root->max << endl;

    inorder_(root->right);
  }


  Int ITree::getSize_(ITree::ITNode *root)
  {
    if (root == NULL) return 0;

    Int size = getSize_(root->left);
    size += sizeof(*root);
    size += getSize_(root->right);
    return size;
  }

}
