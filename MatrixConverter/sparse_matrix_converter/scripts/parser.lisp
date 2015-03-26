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


;;;;
;;;; Parser and code generator for doing "templates" without the pain
;;;; of C++ name mangling (which is platform-dependent, so it's hard
;;;; to generate bindings).
;;;;
(in-package :template-parser)

(defconstant +pragma-template-start+ 
  (cl-ppcre:create-scanner "^\s*#pragma \s*bebop \s*template \s*start\s*$"))
(defconstant +pragma-template-end+ 
  (cl-ppcre:create-scanner "^\s*#pragma \s*bebop \s*template \s*end\s*$"))
(defconstant +variable+ 
  (cl-ppcre:create-scanner "\$\{(.+)\}"))

(defun matchesp (string scanner)
  (not (null (cl-ppcre:scan scanner string))))

(defun matches-variable-p (var-name)
  (matchesp var-name +variable+))

(defun var-subst (line tuple)
  "Tuple is a key,value pair; key is a variable name and value is the
  value to substitute."
  (destructuring-bind (var-name . value) tuple
    ;; Have to add extra backlash for FORMAT's sake
    (let ((var-pattern (format nil "\\$\\{~A\\}" var-name)))
      ;(pprint var-pattern)
      (cl-ppcre:regex-replace-all var-pattern line value :preserve-case t))))

(defun parse-file (infile outfile dict)
  (let ((line (read-line infile nil nil t)))
    (when line
      (cond ((matchesp line +pragma-macro-start+)
	     (parse-macro infile outfile)
	     (parse-file infile outfile dict nil))
	    ((matchesp line +pragma-template-start+)
	     (parse-template infile outfile dict)
	     (parse-file infile outfile dict nil))
	    ((matchesp line +pragma-template-end+)
	     (error "\"#pragma bebop template end\" without \"start\""))
	    (t 
	     (format outfile line)))))) ; emit line verbatim

(defun parse-template (infile outfile dict)
  (write-variants outfile
		  (generate-variants (collect-template-lines infile) 
				     dict)))

;; OK
(defun collect-template-lines (infile)
  (labels ((collect-helper (lines lines-end)
	     (let ((line (read-line infile nil nil t)))
	       (when line
		 (if (matchesp line +pragma-template-end+)
		     lines
		     (if (endp lines)
			 (let ((lines (list line)))
			   (collect-helper lines lines))
			 (let ((last-elt (list line)))
			   (rplacd lines-end last-elt)
			   (collect-helper lines last-elt))))))))
    (collect-helper nil nil)))

(defun gen-line (cur-alist line)
  (if (endp cur-alist)
      line
      (gen-line (cdr cur-alist)
		(var-subst line (car cur-alist)))))
;; FIXME
(defun generate-variants (code-lines dict)
  (map-alist-tuples dict 
		    #'(lambda (cur-alist)
			(mapcar #'(lambda (line) (gen-line cur-alist line)) 
				code-lines))))

(defun write-variants (outfile variants)
  (mapc #'(lambda (variant) (mapc #'(lambda (line) (format outfile "~A" line)) variant))))


(defun generate-code (infilename outfilename)
  (let ((dict (with-open-file (dictfile dictfilename) (read-dict dictfile))))
    (with-open-file (infile infilename)
      (with-open-file (outfile outfilename)
	(parse-file infile outfile dict))))


(defun keys-union (hash keys initial-set)
  "Return the union of lists which are the values of the given set of
KEYS in the given HASH table, unioned with INITIAL-SET."
  (if (endp keys)
      initial-set
      (let ((cur-key (car keys)))
	(multiple-value-bind (cur-set win) (gethash cur-key hash)
	  (if win
	      (keys-union hash (cdr keys) (union initial-set cur-set))
	      (keys-union hash (cdr keys) initial-set))))))

(defmacro if-let (name test &body body)
  (let ((test-val (gensym)))
    `(let ((,test-val ,test))
       (if ,test-val
	   (let ((,name ,test-val))
	     ,@body)
	   (progn
	     ,@body)))))

    
  
