/// @file IntervalTree.hpp
/// @brief Interval Tree.
/// @author Mathias Jacquelin
/// @date 2010-09-27
#ifndef _INTTREE_DECL_HPP_
#define _INTTREE_DECL_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/LogFile.hpp"

using namespace std;

namespace LIBCHOLESKY{

  class ITree{
    public:
      // Structure to represent an interval
      struct Interval
      {
        Int low, high;
        Int block_idx;
      };


protected:
      // Structure to represent a node in Interval Search Tree
      class ITNode
      {
        public:
        Interval *i;  // 'i' could also be a normal variable
        Int max;
        Int height;
        ITNode *left, *right;
        ~ITNode(){
          delete i;

          if(left!=NULL){
            delete left;
          }

          if(right!=NULL){
            delete right;
          }
        }
        void Dump(){
          logfileptr->OFS()<< "[" << i->low << ", " << i->high << "] on "<<i->block_idx;
        }
      };


      ITNode *root_;



      // A utility function to create a new Interval Search Tree Node
      ITNode * newNode_(Interval & i);

      Int height_(ITNode *N);
      Int max_(ITNode *N);

      // A utility function to right rotate subtree rooted with y
      // See the diagram given above.
      ITNode * rightRotate_(ITNode *y);

      ITNode * leftRotate_(ITNode *x);
      // Get Balance factor of node N
      Int getBalance_(ITNode *N);
      // A utility function to insert a new Interval Search Tree Node
      // This is similar to BST Insert.  Here the low value of interval
      // is used tomaintain BST property
      ITNode *insert_(ITNode *root, Interval & i);
      Int recomputeMax_(ITNode * root);

      // A utility function to check if given two intervals overlap
      bool doOVerlap_(Interval &i1, Interval &i2);
      bool doOVerlap_(ITree::Interval &i1, Int & low, Int & high);
      //Interval *intervalSearch_(ITNode *root, Interval &i);
      Interval *intervalSearch_(ITNode *root, Int & begin, Int & end);

      void inorder_(ITNode *root);

      Int getSize_(ITNode *root);




///
///
///      // A utility function to create a new Interval Search Tree Node
///      ITNode * newNode_(Interval & i)
///      {
///        ITNode *temp = new ITNode;
///        temp->i = new Interval(i);
///        temp->max = i.high;
///        temp->left = temp->right = NULL;
///        temp->height = 1;
///        return temp;
///      };
///
///
///Int height_(ITNode *N)
///{
///    if (N == NULL)
///        return 0;
///    return N->height;
///}
///
///Int max_(ITNode *N)
///{
///    if (N == NULL)
///        return 0;
///    return N->max;
///}
///
///
///// A utility function to right rotate subtree rooted with y
///// See the diagram given above.
///ITNode * rightRotate_(ITNode *y)
///{
///    ITNode *x = y->left;
///    ITNode *T2 = x->right;
/// 
///    // Perform rotation
///    x->right = y;
///    y->left = T2;
/// 
///    // Update heights
///    y->height = max(height_(y->left), height_(y->right))+1;
///    x->height = max(height_(x->left), height_(x->right))+1;
/// 
///    // Return new root
///    return x;
///}
///ITNode * leftRotate_(ITNode *x)
///{
///    ITNode *y = x->right;
///    ITNode *T2 = y->left;
/// 
///    // Perform rotation
///    y->left = x;
///    x->right = T2;
/// 
///    //  Update heights
///    x->height = max(height_(x->left), height_(x->right))+1;
///    y->height = max(height_(y->left), height_(y->right))+1;
/// 
///    // Return new root
///    return y;
///}
///
///// Get Balance factor of node N
///Int getBalance_(ITNode *N)
///{
///    if (N == NULL)
///        return 0;
///    return height_(N->left) - height_(N->right);
///}
///
///
///      // A utility function to insert a new Interval Search Tree Node
///      // This is similar to BST Insert.  Here the low value of interval
///      // is used tomaintain BST property
///      ITNode *insert_(ITNode *root, Interval & i)
///      {
///        // Base case: Tree is empty, new node becomes root
///        if (root == NULL)
///          return newNode_(i);
///
///        // Get low value of interval at root
///        Int l = root->i->low;
///
///        // If root's low value is smaller, then new interval goes to
///        // left subtree
///        if (i.low < l)
///          root->left = insert_(root->left, i);
///
///        // Else, new node goes to right subtree.
///        else
///          root->right = insert_(root->right, i);
///
///        // Update the max value of this ancestor if needed
///        if (root->max < i.high)
///          root->max = i.high;
///
///
///
///
///    /* 2. Update height of this ancestor node */
///    root->height = max(height_(root->left), height_(root->right)) + 1;
/// 
///    /* 3. Get the balance factor of this ancestor node to check whether
///       this node became unbalanced */
///    Int balance = getBalance_(root);
/// 
///    // If this node becomes unbalanced, then there are 4 cases
/// 
///    // Left Left Case
///    if (balance > 1 && i.low < root->left->i->low)
///        root = rightRotate_(root);
/// 
///    // Right Right Case
///    if (balance < -1 && i.low > root->right->i->low)
///        root = leftRotate_(root);
/// 
///    // Left Right Case
///    if (balance > 1 && i.low > root->left->i->low)
///    {
///       root->left =  leftRotate_(root->left);
///       root = rightRotate_(root);
///    }
/// 
///    // Right Left Case
///    if (balance < -1 && i.low < root->right->i->low)
///    {
///        root->right = rightRotate_(root->right);
///        root = leftRotate_(root);
///    }
///    //recompute max
///    recomputeMax_(root);
///
///        return root;
///      }
///
///      Int recomputeMax_(ITNode * root){
///        root->max = root->i->high;
///
///        if(root->left != NULL){
///          recomputeMax_(root->left);
///          if(root->max < max_(root->left)){
///            root->max = max_(root->left);
///          } 
///        }
///
///        if(root->right != NULL){
///          recomputeMax_(root->right);
///          if(root->max < max_(root->right)){
///            root->max = max_(root->right);
///          } 
///        }
///        return root->max;
///      }
///
///
///      // A utility function to check if given two intervals overlap
///      bool doOVerlap_(Interval &i1, Interval &i2)
///      {
///        if (i1.low <= i2.high && i2.low <= i1.high)
///          return true;
///        return false;
///      }
///
///
///      // The main function that searches a given interval i in a given
///      // Interval Tree.
///      Interval *intervalSearch_(ITNode *root, Interval &i)
///      {
///        // Base Case, tree is empty
///        if (root == NULL) return NULL;
///
///        // If given interval overlaps with root
///        if (doOVerlap_(*(root->i), i))
///          return root->i;
///
///        // If left child of root is present and max of left child is
///        // greater than or equal to given interval, then i may
///        // overlap with an interval is left subtree
///        if (root->left != NULL && root->left->max >= i.low)
///          return intervalSearch_(root->left, i);
///
///        // Else interval can only overlap with right subtree
///        return intervalSearch_(root->right, i);
///      }
///
///
///      void inorder_(ITNode *root)
///      {
///        if (root == NULL) return;
///
///        inorder_(root->left);
///
///        logfileptr->OFS()<< "[" << root->i->low << ", " << root->i->high << "] on "<<root->i->block_idx
///          << " max = " << root->max << endl;
///
///        inorder_(root->right);
///      }
///
///
///      Int getSize_(ITNode *root)
///      {
///        if (root == NULL) return 0;
///
///        Int size = getSize_(root->left);
///        size += sizeof(*root);
///        size += getSize_(root->right);
///        return size;
///      }
///

public:


      ITree():root_(NULL)
      {
      }


      ~ITree()
      {
        if(root_!=NULL){
          delete root_;
        }
      }


      void Dump()
      {
        logfileptr->OFS()<<"Height of tree is "<<height_(root_)<<endl;
        inorder_(root_);
      }


      void Insert(Interval & i)
      {
        //logfileptr->OFS()<<"Interval ["<<i.low<<" -- "<<i.high<<"] is inserted"<<endl;
        root_ = insert_(root_,i);
        
      }

      void Insert(Int i)
      {
        Interval it;  
        it.low = i;
        it.high = i;
        Insert(it);
      }


      Interval * IntervalSearch(Interval & i){
        return intervalSearch_(root_,i.low,i.high);
      }

      Interval * IntervalSearch(Int & low, Int & high){
        return intervalSearch_(root_,low,high);
      }

      Int StorageSize(){
        Int size = 0;
        size = getSize_(root_);
        return size;
      }


  };

  namespace UnitTest{
    bool ITree_Test();
  }
}

#endif
