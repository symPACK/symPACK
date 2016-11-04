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
/// @file ETree.hpp
/// @brief Various elimination tree related subroutines.
/// @author Mathias Jacquelin
/// @date 2013-08-31
#ifndef _ETREE_HPP_ 
#define _ETREE_HPP_

#include  <stdlib.h>

#include "sympack/Environment.hpp"
#include "sympack/DistSparseMatrixGraph.hpp"
#include "sympack/Ordering.hpp"


namespace symPACK{
  class DisjointSet{
    protected:
      std::vector<Int> pp_;
      std::vector<Int> root_;
    public:
      void Initialize(Int n);
      void Finalize();
      inline Int makeSet(Int i);
      inline Int link(Int s, Int t);
      inline void Union(Int s, Int t, Int root = -1);
      inline Int find(Int i);
      inline Int & Root(Int i){return root_[i];};
  };

class ETree{

  friend std::ostream& operator<<( std::ostream& os, const ETree& tree);
  template <class F> friend class DistSparseMatrix;

protected:

  void BTreeToPO(std::vector<Int> & fson, std::vector<Int> & brother, std::vector<Int> & perm);

public:
  ETree();

  //not working
  void ConstructETree(DistSparseMatrixGraph & aDistExp, Ordering & aOrder);
  
  void ConstructETree(SparseMatrixGraph & sgraph, Ordering & aOrder, MPI_Comm & aComm);
  void PostOrderTree(Ordering & aOrder);
  void DeepestFirst(Ordering & aOrder);

  ETree ToSupernodalETree(std::vector<Int> & aXsuper,std::vector<Int> & aSupMembership,Ordering & aOrder) const;

  inline bool IsPostOrdered() const { return bIsPostOrdered_;};
  inline Int n() const { return n_; };
  inline Int PostParent(Int i) const { return poparent_[i]; }
  inline Int Parent(Int i) const { return parent_[i]; };
  inline Int Size() const { return parent_.size(); };
  void SortChildren(std::vector<Int> & cc, Ordering & aOrder);

protected:
  Int n_;
  bool bIsPostOrdered_;

  std::vector<Int> parent_;
  std::vector<Int> poparent_;
};






  inline Int DisjointSet::makeSet(Int i){
    pp_[i-1]=i;
    return i;
  }

  inline Int DisjointSet::link(Int s, Int t){
    pp_[s-1]=t;

    return t;
  }

  inline Int DisjointSet::find(Int i){
    Int p, gp;

    p = pp_[i-1];
    gp=pp_[p-1];

    while(gp!=p){
      i = makeSet(gp);
      p = pp_[i-1];
      gp = pp_[p-1];
    }

    return p;
  }


  inline void DisjointSet::Union(Int s, Int t, Int root){
    Int tSet= find(t);
    Int sSet= find(s);
    sSet = link(sSet, tSet );
    root_[sSet-1] = root==-1?t:root;
  }














}
#endif //_ETREE_HPP_
