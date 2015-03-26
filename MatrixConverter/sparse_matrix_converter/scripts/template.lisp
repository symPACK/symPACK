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

(defmacro over-key-tuples (alists &body body)
  (if (endp alists)
      ,@body
      `(loop :for (k . v) :in ,(car alists) :do 
	  ,(over-key-tuples (cdr alists) body))))

(defmacro over-tuples ((&rest lists) &body body)
  (progn
    (print lists)
    (print body)
    'nil))


  (if (endp lists)
      'nil
      (if (endp (car lists))
	  'nil
	  (let ((loopvar (gensym)))
	    `(dolist (,loopvar ,(car lists))
	       ,@body)))))


;      `(loop :for (k v) :in ,(car alists) :do ,(over-key-tuples (cdr alists) body))))
