/// @file IntervalTree.hpp
/// @brief Interval Tree.
/// @author Mathias Jacquelin
/// @date 2010-09-27
#ifndef _INTTREE_DECL_HPP_
#define _INTTREE_DECL_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/LogFile.hpp"

#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) 
#define TIMER_STOP(a) 
#endif
#endif


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


      // A utility function to insert a new Interval Search Tree Node
      // This is similar to BST Insert.  Here the low value of interval
      // is used tomaintain BST property
      virtual ITNode *insert_(ITNode *root, Interval & i);
      Int recomputeMax_(ITNode * root);

      // A utility function to check if given two intervals overlap
      bool doOVerlap_(Interval &i1, Interval &i2);
      bool doOVerlap_(ITree::Interval &i1, Int & low, Int & high);
      //Interval *intervalSearch_(ITNode *root, Interval &i);
      Interval *intervalSearch_(ITNode *root, Int & begin, Int & end);

      void inorder_(ITNode *root);

      Int getSize_(ITNode *root);


public:


      ITree():root_(NULL)
      {
      }


      virtual ~ITree()
      {
        if(root_!=NULL){
          delete root_;
        }
      }


      virtual void Rebalance(){};

      void Dump()
      {
        logfileptr->OFS()<<"Height of tree is "<<height_(root_)<<endl;
        inorder_(root_);
      }


      void Insert(Interval & i)
      {
        TIMER_START(ITREE_INSERT);
        //logfileptr->OFS()<<"Interval ["<<i.low<<" -- "<<i.high<<"] is inserted"<<endl;
        root_ = insert_(root_,i);
        
        TIMER_STOP(ITREE_INSERT);
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


  class AVLITree: public ITree{
    public:

protected:
      // A utility function to right rotate subtree rooted with y
      // See the diagram given above.
      ITNode * rightRotate_(ITNode *y);

      ITNode * leftRotate_(ITNode *x);
      // Get Balance factor of node N
      Int getBalance_(ITNode *N);

      // A utility function to insert a new Interval Search Tree Node
      // This is similar to BST Insert.  Here the low value of interval
      // is used tomaintain BST property
      virtual ITNode *insert_(ITNode *root, Interval & i);


public:


      AVLITree():ITree()
      {
      }


      virtual ~AVLITree()
      {
      }
  };

  class DSWITree: public ITree{
    public:

protected:
      inline Int fullSize_ ( Int size );
      void tree_to_vine_ (ITNode * root, Int &size);
      void compression_  (ITNode * root, Int count);
      void vine_to_tree_ (ITNode * root, Int n);
      void correctTree_  (ITNode * root);


public:


      DSWITree():ITree()
      {
      }


      virtual ~DSWITree()
      {
      }
      
      virtual void Rebalance();
  };







  namespace UnitTest{
    bool ITree_Test();
  }














////
////struct TreeNode
////{
////    enum Kind {RED, BLUE};
////
////    TreeNode(Kind kind_, TreeNode* left_ = NULL, TreeNode* right_ = NULL)
////        : kind(kind_), left(left_), right(right_)
////    {}
////
////    Kind kind;
////    TreeNode *left, *right;
////};
////
////template <typename Derived>
////class GenericVisitor
////{
////public:
////    void visit_preorder(TreeNode* node)
////    {
////        if (node) {
////            dispatch_node(node);
////            visit_preorder(node->left);
////            visit_preorder(node->right);
////        }
////    }
////
////    void visit_inorder(TreeNode* node)
////    {
////        if (node) {
////            visit_inorder(node->left);
////            dispatch_node(node);
////            visit_inorder(node->right);
////        }
////    }
////
////    void visit_postorder(TreeNode* node)
////    {
////        if (node) {
////            visit_postorder(node->left);
////            visit_postorder(node->right);
////            dispatch_node(node);
////        }
////    }
////
////    void handle_RED(TreeNode* node)
////    {
////        cerr << "Generic handle RED\n";
////    }
////
////    void handle_BLUE(TreeNode* node)
////    {
////        cerr << "Generic handle BLUE\n";
////    }
////
////private:
////    // Convenience method for CRTP
////    //
////    Derived& derived()
////    {
////        return *static_cast<Derived*>(this);
////    }
////
////    void dispatch_node(TreeNode* node)
////    {
////        switch (node->kind) {
////            case TreeNode::RED:
////                derived().handle_RED(node);
////                break;
////            case TreeNode::BLUE:
////                derived().handle_BLUE(node);
////                break;
////            default:
////                assert(0);
////        }
////    }
////};
////
////class SpecialVisitor : public GenericVisitor<SpecialVisitor>
////{
////public:
////    void handle_RED(TreeNode* node)
////    {
////        cerr << "RED is special\n";
////    }
////};
////








}
#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#endif
#endif


#endif
