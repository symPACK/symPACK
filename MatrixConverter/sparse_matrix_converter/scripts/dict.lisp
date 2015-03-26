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

(in-package :template-parser)

(defun listify (thing)
  (cond ((consp thing)
	 thing)
	((atom thing)
	 (list thing))
	(t 
	 (error "Thing must be an atom or a list"))))

(let ((*associated-variables* (make-hash-table :test #'equal))
      (*var-vals*             (make-hash-table :test #'equal)))
  (defun associate-variables (sets)
    "SETS is a list of lists.  Each member list contains all the
variables that are to be associated together into a group.  Within
each group, the variables take on parameter values in lockstep, rather
than independently."
    (labels ((add-set (cur-set)
	       ;; CUR-SET contains at least one variable name.  Union
	       ;; up the sets associated with those variable(s) in the
	       ;; set that may already have entries in the hash table.
	       (let ((the-union (keys-union *associated-variables* cur-set nil)))
		 (mapc (lambda (name) (setf (gethash name *associated-variables*) the-union)) cur-set)))
	     (current-set (sets) (listify (car sets))))
      (if (endp sets)
	  nil
	  (if (not (consp sets))
	      (error "SETS must be a list")
	      (progn
		(add-set (current-set sets))
		(associate-variables (cdr sets)))))))

  (defun associatedp (name1 name2)
    "Return true iff variable NAME1 is associated with variable NAME2.
ASSOCIATE-VARIABLES ensures that this relation is symmetric."
    (multiple-value-bind (lst win) (gethash name1 *hash*)
      (if (not win)
	  (error "Variable ~A is not in the table" name1)
	  (if (atom lst) 
	      (equal name2 lst)
	      (member name2 lst :test #'equal)))))

  (defun read-dict (dict)
    (labels ((first-pass (dict)
	       (cond ((endp dict)
		      nil)
		     ((not (consp dict))
		      (error "DICT must be a list of lists"))
		     (t
		      (let ((cur-entry (car dict)))
			(cond ((not (consp dict))
			       (error "Current entry ~A must be a cons" cur-entry))
			      ((endp (cdr cur-entry))
			       (error "Current entry ~A has no values" cur-entry))
			      (t
			       (let ((name (car cur-entry))
				     (vals (listify (cdr cur-entry))))
				 (setf (gethash name *var-vals*) (union vals (gethash name *var-vals*)))
				 (first-pass (cdr dict)))))))))
	     (second-pass (dict) ;; should go backwards
	       (cond ((endp dict)
		      nil)
		     (t
		      (let ((name (car cur-entry))
			    (vals (listify (cdr cur-entry))))
			(setf (gethash name *var-vals*) (union vals (gethash name *var-vals*)))
			(second-pass (cdr dict))))))
	     (third-pass (dict) ;; should go forwards; check that
				;; bound variables have the same # of
				;; values.
	       (if (endp dict)
		   t
		   (let ((name (car cur-entry))
			 (vals (listify (cdr cur-entry)))
			 (bound-names (gethash name *associated-variables*)))
		     (and (reduce (lambda (x cur-name) (and x (equal vals (gethash cur-name *var-vals*))))
				  bound-names)
			  (third-pass (cdr dict)))))))
      (first-pass dict)
      (second-pass (reverse dict))
      (third-pass dict))))
			     

				 
		      

		      
		    


	



;; A dictionary is stored in a file as an alist.
;; Keys should all be strings.
(defun valid-dict-p (dict)
  ;; Note that (length (cons 1 2)) doesn't work.  Checking (cdr x) for
  ;; non-NIL excludes the possibility of an empty set of
  ;; substitutions.  The substitution set may be either an atom (and
  ;; therefore a single value), or a list of values.
  (flet ((dict-entry-p (x)
	   (and (consp x) 
		(stringp (car x))
		(cdr x))))
    (or (endp dict)
	(and (dict-entry-p (cdar dict))
	     (valid-dict-p (cdr dict))))))

(defun read-dict (dictfile)
  (let ((dict (read dictfile)))
    (if (valid-dict-p dict)
	dict
	(error "Parse error in dictionary file"))))


