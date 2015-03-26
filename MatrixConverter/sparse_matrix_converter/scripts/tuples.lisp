;;; -*- Mode: lisp -*-
;;; template-parser (BeBOP C "templates" parser and code generator)
;;; Copyright (C) 2008 Mark Hoemmen <mhoemmen@cs.berkeley.edu>
;;;
;;; This file is part of the template-parser library, which parses C
;;; "templates" and generates code for different data types.
;;; template-parser is free software: you can redistribute it and/or
;;; modify it under the terms of the GNU General Public License as
;;; published by the Free Software Foundation, either version 3 of the
;;; License, or (at your option) any later version.
;;;
;;; template-parser is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;; General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with the source code of this library.  If not, see
;;; <http://www.gnu.org/licenses/>.
;;;

(in-package #:template-parser)

(defun singletonp (x)
  (and (listp x)
       (not (endp x))
       (not (atom x))
       (not (cdr x))))

(defun map-tuples-helper (ls f curtupl)
  (if (endp ls)
      nil
      (if (singletonp ls)
	  (dolist (x (car ls))
	    (let ((tuple (append curtupl (list x))))
	      (funcall f tuple)))
	  (dolist (x (car ls))
	    (map-tuples-helper (cdr ls) f (append curtupl (list x)))))))

(defun map-tuples (ls f)
  (map-tuples-helper ls f nil))

;; Conses a lot, and does a lot of slow appending.
(defun map-alist-tuples-helper (alists f curtupl)
  (if (endp alists)
      nil
      (destructuring-bind (key &rest values) (car alists)
	(if (singletonp alists)
	    (dolist (val values)
	      (funcall f (append curtupl (list (cons key val)))))
	    (dolist (val values)
	      (map-alist-tuples-helper (cdr alists) f (append curtupl (list (cons key val)))))))))


;; Conses less than the above version, and doesn't call APPEND.
(defun map-alist-tuples-helper2 (alists f curtupl)
  (if (endp alists)
      nil
      (destructuring-bind (key &rest values) (car alists)
	;; We avoid calling O(n) APPEND at each iteration of DOLIST by
	;; copying curtupl, keeping a ref to the end of the copy, and
	;; destructively appending the last tuple cell at each
	;; iteration.  So there's only one O(n) copy, instead of n of
	;; the O(n) APPENDs (each of which would CONS O(n) times
	;; anyway).
	(let* ((the-curtupl (copy-list curtupl))
	       (end-the-curtupl (last the-curtupl)))
	  (if (endp the-curtupl)
	      (if (singletonp alists)
		  (dolist (val values)
		    (setf the-curtupl (list (cons key val)))
		    (funcall f the-curtupl))
		  (dolist (val values)
		    (setf the-curtupl (list (cons key val)))
		    (map-alist-tuples-helper2 (cdr alists) f the-curtupl)))
	      (if (singletonp alists)
		  (dolist (val values)
		    (setf (cdr end-the-curtupl) (list (cons key val)))
		    (funcall f the-curtupl))
		  (dolist (val values)
		    (setf (cdr end-the-curtupl) (list (cons key val)))
		    (map-alist-tuples-helper2 (cdr alists) f the-curtupl))))))))

;; Recycles the "current pair," but calls APPEND.  Conses less than
;; #1, but more than #2.
(defun map-alist-tuples-helper3 (alists f curtupl)
  (if (endp alists)
      nil
      (destructuring-bind (key &rest values) (car alists)
	(let ((pair (cons nil nil))) ; to be recycled to avoid consing
	  (if (singletonp alists)
	      (dolist (val values)
		(setf (car pair) key)
		(setf (cdr pair) val)
		(funcall f (append curtupl (list pair))))
	      (dolist (val values)
		(setf (car pair) key)
		(setf (cdr pair) val)
		(map-alist-tuples-helper3 (cdr alists) f (append curtupl (list pair)))))))))

;; Conses less than #1, #2, or #3 above.
(defun map-alist-tuples-helper4 (alists f curtupl)
  (if (endp alists)
      nil
      (destructuring-bind (key &rest values) (car alists)
	(let ((end-curtupl (last curtupl)))
	  (if (endp curtupl)
	      (if (singletonp alists)
		  (dolist (val values)
		    (setf curtupl (list (cons key val)))
		    (funcall f curtupl))
		  (dolist (val values)
		    (setf curtupl (list (cons key val)))
		    (map-alist-tuples-helper4 (cdr alists) f curtupl)))
	      (if (singletonp alists)
		  (dolist (val values)
		    (setf (cdr end-curtupl) (list (cons key val)))
		    (funcall f curtupl))
		  (dolist (val values)
		    (setf (cdr end-curtupl) (list (cons key val)))
		    (map-alist-tuples-helper4 (cdr alists) f curtupl))))))))

;; The way to improve on the consing is to recycle space for the
;; current tuple.  Keep a junk list of pairs around and append
;; destructively as necessary, rather than consing new (key . val)
;; pairs.
;;
;; This one seems to do the best, in terms of speed and consing.
;; Constant space per key, and constant time per (key . val) pair.
;; You have to supply scratch space: curtupl is an alist with the same
;; number of entries as the number of keys in ALISTS, and curspot
;; should be EQ to curtupl on first call.
(defun mapath5 (alists f curtupl curspot)
  "Helper function for MAP-ALIST-TUPLES"
  (if (endp alists)
      nil
      (destructuring-bind (key &rest values) (car alists)
	(if (endp values)
	    ;; Key but no values; skip it because it doesn't assume
	    ;; any values at all.
	    (mapath5 (cdr alists) f curtupl curspot)
	    (if (singletonp alists)
		(dolist (val values)
		  (setf (caar curspot) key)
		  (setf (cdar curspot) val)
		  (funcall f curtupl))
		(dolist (val values)
		  (setf (caar curspot) key)
		  (setf (cdar curspot) val)
		  (mapath5 (cdr alists) f curtupl (cdr curspot))))))))

;; Make scratch space for MAPATH5.
(defun make-empty-alist-n (n)
  (if (<= n 0)
      nil
      (cons (cons nil nil) (make-empty-alist-n (1- n)))))

(defun map-alist-tuples (alists f)
  "ALISTS is a list of (key . vals) pairs, in which VALS consists of
  one or more values.  KEY is the name of a variable, and VALS the
  possible values which that variable can take on.  The function F
  operates on a single KEY, VALUE pair.  MAP-ALIST-TUPLES maps F over
  all possible alists which can be generated from the given
  combinations of key, value pairs.

  Here is an example run:

  (defvar x (list (cons \"foo\" (list 1 2 3)) 
                  (cons \"bar\" (list 11 12)) 
                  (cons \"baz\" (list 21 22 23 24))))      
  X
  * (map-alist-tuples x #'pprint)

  ((\"foo\" . 1) (\"bar\" . 11) (\"baz\" . 21))
  ((\"foo\" . 1) (\"bar\" . 11) (\"baz\" . 22))
  ((\"foo\" . 1) (\"bar\" . 11) (\"baz\" . 23))
  ((\"foo\" . 1) (\"bar\" . 11) (\"baz\" . 24))
  ((\"foo\" . 1) (\"bar\" . 12) (\"baz\" . 21))
  ((\"foo\" . 1) (\"bar\" . 12) (\"baz\" . 22))
  ((\"foo\" . 1) (\"bar\" . 12) (\"baz\" . 23))
  ((\"foo\" . 1) (\"bar\" . 12) (\"baz\" . 24))
  ((\"foo\" . 2) (\"bar\" . 11) (\"baz\" . 21))
  ((\"foo\" . 2) (\"bar\" . 11) (\"baz\" . 22))
  ((\"foo\" . 2) (\"bar\" . 11) (\"baz\" . 23))
  ((\"foo\" . 2) (\"bar\" . 11) (\"baz\" . 24))
  ((\"foo\" . 2) (\"bar\" . 12) (\"baz\" . 21))
  ((\"foo\" . 2) (\"bar\" . 12) (\"baz\" . 22))
  ((\"foo\" . 2) (\"bar\" . 12) (\"baz\" . 23))
  ((\"foo\" . 2) (\"bar\" . 12) (\"baz\" . 24))
  ((\"foo\" . 3) (\"bar\" . 11) (\"baz\" . 21))
  ((\"foo\" . 3) (\"bar\" . 11) (\"baz\" . 22))
  ((\"foo\" . 3) (\"bar\" . 11) (\"baz\" . 23))
  ((\"foo\" . 3) (\"bar\" . 11) (\"baz\" . 24))
  ((\"foo\" . 3) (\"bar\" . 12) (\"baz\" . 21))
  ((\"foo\" . 3) (\"bar\" . 12) (\"baz\" . 22))
  ((\"foo\" . 3) (\"bar\" . 12) (\"baz\" . 23))
  ((\"foo\" . 3) (\"bar\" . 12) (\"baz\" . 24))
  "
  (let ((current-tuple (make-empty-alist-n (length alists))))
    (mapath5 alists f current-tuple current-tuple)))

;;
;; SERIES:SCAN-LIST-OF-LISTS or SERIES:SCAN-LIST-OF-LISTS-FRINGE could
;; be helpful, but you'll still have to strip off the keys, and you'll
;; also have to know the correspondence between key and value list.  
;;
;; A better approach might be to use SERIES:SCAN-LIST to get (key
;; . value-list) pairs, then call SCAN-LIST iteratively (once for each
;; key).  You'll still want scratch space in order to avoid consing,
;; so I'm not sure how much it would help.
;;
		 
