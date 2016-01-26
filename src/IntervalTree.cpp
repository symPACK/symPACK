#include "sympack/IntervalTree.hpp"



















namespace SYMPACK{






  namespace UnitTest{

    bool ITree_Test(){


      ITree * tree;
      bool success = true;



      ITree::Interval * resIt;
      ITree::Interval it;
      it.low = 0;
      it.high = 1;

      //Empty tree
      // Interval search should return NULL
////      tree = new ITree();
////      resIt = tree->IntervalSearch(it);
////      success &= (resIt==NULL);
////      delete tree;
////
////      //Insertion
////      //Insert intervals in increasing order
////      tree = new ITree();
////      it.low = 2;
////      it.high = 2;
////      tree->Insert(it);
////
////      it.low = 3;
////      it.high = 3;
////      tree->Insert(it);
////
////      it.low = 4;
////      it.high = 4;
////      tree->Insert(it);
////
////
//////      tree->Dump();
////
////      delete tree;
////      //Insert two intervals with the second going to the left
////      tree = new ITree();
////      it.low = 4;
////      it.high = 4;
////      tree->Insert(it);
////
////      it.low = 3;
////      it.high = 3;
////      tree->Insert(it);
////
////      it.low = 2;
////      it.high = 2;
////      tree->Insert(it);
////
//////      tree->Dump();
////
////      delete tree;


if(iam==0){


      tree = new ITree();
      tree->Insert(14);
//      tree->Dump();
      tree->Insert(17);
//     tree->Dump();
      tree->Insert(46);
//      tree->Dump();
      tree->Insert(15);
//      tree->Dump();
      tree->Insert(16);
//      tree->Dump();
      tree->Insert(35);
//      tree->Dump();
      tree->Insert(36);
//      tree->Dump();
      tree->Insert(40);
//      tree->Dump();
      tree->Insert(42);
//      tree->Dump();
      tree->Insert(43);
//      tree->Dump();
      tree->Insert(45);
//      tree->Dump();
      tree->Insert(47);
      tree->Dump();

logfileptr->OFS()<<"Rebalancing the tree"<<endl;
      tree->Rebalance();
      tree->Dump();

      delete tree;

}


      return success;

    }
  }

}
