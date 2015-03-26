;;;; File: conversion.lisp
;;;; Author: mfh
;;;; Since: 09 May 2007
;;;; Time-stamp: <2008-07-16 11:01:27 mhoemmen>
;;;;
;;;; Copyright (c) 2008, Regents of the University of California 
;;;; All rights reserved.
;;;; Redistribution and use in source and binary forms, with or
;;;; without modification, are permitted provided that the
;;;; following conditions are met:
;;;; 
;;;; * Redistributions of source code must retain the above copyright
;;;;   notice, this list of conditions and the following disclaimer.
;;;;
;;;; * Redistributions in binary form must reproduce the above copyright 
;;;;   notice, this list of conditions and the following disclaimer in 
;;;;   the documentation and/or other materials provided with the 
;;;;   distribution.
;;;;
;;;; * Neither the name of the University of California, Berkeley, nor
;;;;   the names of its contributors may be used to endorse or promote
;;;;   products derived from this software without specific prior
;;;;   written permission.  
;;;; 
;;;; THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
;;;; "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
;;;; LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
;;;; FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
;;;; COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
;;;; INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
;;;; (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
;;;; SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
;;;; HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
;;;; STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
;;;; ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
;;;; OF THE POSSIBILITY OF SUCH DAMAGE.


;; A type conversion graph is a hash table from sparse matrix type, to
;; a list of types to which the first type can be converted (the
;; opposite conversion direction may not hold; this is a directed
;; graph).
(defvar *typegraph* (make-hash-table :test #'eq))
;; key is ordered pair (parent child)
(defvar *costs* (make-hash-table :test #'equal))

;; Builds the typegraph from a neighbor list
(defun typegraph (csr)
  (loop for neighbor-list in csr do
	(setf (gethash *typegraph* (car neighbor-list))
	      (cdr neighbor-list))))

;; Adds a node and one of its children to the typegraph.
(defun register-conversion! (parent child &optional (cost 1 cost-provided-p))
  (progn
    (setf (gethash *typegraph* parent)
	  (cons child (gethash *typegraph* parent)))
    (if cost-provided-p
	(setf (gethash *costs* (list parent child)) cost))))
	      
;; Returns the children of NODE
(defun children (node)
  (gethash *typegraph* node))

(defun cost (parent child)
  (gethash *costs* (list parent child)))

;; What you really want here is BFS, ideally with weighted edges for
;; cost estimates.  When we do BFS, as soon as we reach the target
;; node, we follow the backwards links (which are unique -- they
;; encode the distance information) to get to the start node.  The
;; path from start to target encodes the conversions needed.
;;
;; We can make BFS! a "greedy best-first search" by making the queue a
;; priority queue, prioritized on edge cost.

(defvar *parent* (make-hash-table :test #'eq))

(defun visitedp (node)
  (gethash *parent* node))
(defun visit! (node parent)
  (setf (gethash *parent* node) parent))
(defun parent (node)
  (gethash *parent* node))

(defun bfs! (start target)
  (let ((Q (make-queue)))
    (progn
      (setf (parent start) nil)
      (enqueue! Q start)
      (loop while (not (queue-emptyp Q))
	   (let ((node (dequeue! Q)))
	     (if (eq target node)
		 (return-from bfs!)
		 (mapc #'(lambda (child)
			   (if (not (visitedp child))
			       (visit! child node)))
		       (children node))))))))

(defun path (start target)
  "Returns the path from start to target, represented as a list
   of ordered pairs which identify the individual conversion steps."
  (labels ((path-helper (start target)
	     (if (eq start target) 
		 nil
		 (let ((target-parent (parent target)))
		   (if (not target-parent)
		       (error "No path from ~A to ~A" start target)
		       (cons (list target-parent target)
			     (path start target-parent)))))))
    (reverse (path-helper start target))))

