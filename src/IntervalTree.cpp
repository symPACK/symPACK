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
#include "sympack/IntervalTree.hpp"

namespace symPACK{

//
//
//
//
//
//  namespace UnitTest{
//
//    bool ITree_Test(){
//int iam = upcxx::myrank();
////    int iam =0;
////    int np =1;
////    MPI_Comm_rank(comm,&iam);
////    MPI_Comm_size(comm,&np);
//
//
//
//      ITree * tree;
//      bool success = true;
//
//
//
//      ITree::Interval * resIt;
//      ITree::Interval it;
//      it.low = 0;
//      it.high = 1;
//
//      //Empty tree
//      // Interval search should return NULL
//////      tree = new ITree();
//////      resIt = tree->IntervalSearch(it);
//////      success &= (resIt==NULL);
//////      delete tree;
//////
//////      //Insertion
//////      //Insert intervals in increasing order
//////      tree = new ITree();
//////      it.low = 2;
//////      it.high = 2;
//////      tree->Insert(it);
//////
//////      it.low = 3;
//////      it.high = 3;
//////      tree->Insert(it);
//////
//////      it.low = 4;
//////      it.high = 4;
//////      tree->Insert(it);
//////
//////
////////      tree->Dump();
//////
//////      delete tree;
//////      //Insert two intervals with the second going to the left
//////      tree = new ITree();
//////      it.low = 4;
//////      it.high = 4;
//////      tree->Insert(it);
//////
//////      it.low = 3;
//////      it.high = 3;
//////      tree->Insert(it);
//////
//////      it.low = 2;
//////      it.high = 2;
//////      tree->Insert(it);
//////
////////      tree->Dump();
//////
//////      delete tree;
//
//
//if(iam==0){
//
//
//      tree = new ITree();
//      tree->Insert(14);
////      tree->Dump();
//      tree->Insert(17);
////     tree->Dump();
//      tree->Insert(46);
////      tree->Dump();
//      tree->Insert(15);
////      tree->Dump();
//      tree->Insert(16);
////      tree->Dump();
//      tree->Insert(35);
////      tree->Dump();
//      tree->Insert(36);
////      tree->Dump();
//      tree->Insert(40);
////      tree->Dump();
//      tree->Insert(42);
////      tree->Dump();
//      tree->Insert(43);
////      tree->Dump();
//      tree->Insert(45);
////      tree->Dump();
//      tree->Insert(47);
//      tree->Dump();
//
//logfileptr->OFS()<<"Rebalancing the tree"<<std::endl;
//      tree->Rebalance();
//      tree->Dump();
//
//      delete tree;
//
//}
//
//
//      return success;
//
//    }
//  }
//
}
