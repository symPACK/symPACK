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

;; List of macros
(defvar *macros* nil)

(defun parse-macro-line (line)
  (destructuring-bind (name value)
      (cl-ppcre:split "\\s+" line)
    (push (cons name value) *macros*)))

(defun parse-macro (infile)
  (let ((line (read-line infile nil nil t)))
    (when line
      (parse-macro-line line))))

(defconstant +pragma-macro-start+ 
  (cl-ppcre:create-scanner "^#\s*pragma\s+bebop\s+macro\s+$"))
