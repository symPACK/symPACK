;;;; File: interface.lisp
;;;; Author: mfh
;;;; Since: 09 May 2007
;;;; Time-stamp: <2008-07-16 11:01:49 mhoemmen>
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


(require 'asdf)
(pushnew "/Users/mhoemmen/pkg/cffi/" asdf:*central-registry* :test 'equal)
(asdf:oos 'asdf:load-op :cffi)

(defpackage #:sp
  (:use #:common-lisp #:cffi #:cffi-utils)
  (:export load! save! format? convert! mult triple-product destroy! validp print!))

(in-package #:sp)

;;; The Sparse Matrix Converter library
(define-foreign-library lib-sparse-matrix-converter
  (t (:default "libsparse_matrix_converter")))

;;; The BeBOP Utility Library
(define-foreign-library lib-bebop-util
  (t (:default "libbebop_util")))

(use-foreign-library lib-bebop-util)
(use-foreign-library lib-sparse-matrix-converter)


(defcfun ("sp_load" %sp-load) :pointer
  (path :string)
  (fmt :string))
(defun load! (path &optional fmt)
  (let ((the-fmt (cond ((null fmt) "MATRIX_MARKET")
		       (t fmt))))
    (%sp-load path the-fmt)))

(defcfun ("sp_save" %sp-save) :int
  (A :pointer)
  (path :string)
  (fmt :string))
(defun save! (A path &optional fmt)
  (let ((the-fmt (cond ((null fmt) "MATRIX_MARKET")
		       (t fmt))))
    (%sp-save A path the-fmt)))

(defcfun ("sp_format" format?) :void 
  (A :pointer))

(defcfun ("sp_convert" convert!) :int
  (A :pointer)
  (type :string))

(defcfun ("sp_mult" mult) :pointer
  (B :pointer)
  (A :pointer))

;;; Returns R*A*P, where R is given as R^T.
(defcfun ("sp_triprod" triple-product) :pointer
  (RT :pointer)
  (A :pointer)
  (P :pointer))

(defcfun ("sp_destroy" %destroy!) :void
  (A :pointer))

(defun destroy! (A)
  (progn
    (assert (and (not (null A)) (not (cffi:null-pointer-p A))))
    (%destroy! A)))

(defcfun ("valid_sparse_matrix" %valid-sparse-matrix) :int
  (A :pointer))
(defun validp (A)
  (/= 0 (%valid-sparse-matrix A)))

(defcfun ("sp_print" %sp-print) :void
  (A :pointer))
(defun print! (A)
  (%sp-print A))
