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

(defpackage #:edu.berkeley.cs.bebop.template-parser-system
  (:use #:common-lisp #:asdf)
  (:nicknames #:template-parser-system))

(in-package #:edu.berkeley.cs.bebop.template-parser-system)

(asdf:defsystem :template-parser
    :description "BeBOP template parser and code generator"
    :author "Mark Hoemmen <mhoemmen@cs.berkeley.edu>"
    :maintainer "Mark Hoemmen <mhoemmen@cs.berkeley.edu>"
    :licence "GPL v. 3"
    :version "1.0.0"
    :depends-on () ; (:split-sequence :dir)
    :serial t
    :components ((:file "package") 
		 (:file "util") 
		 (:file "dict")
		 (:file "tuples")
		 (:file "parser")
		 ))

		 
