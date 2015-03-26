(defparameter +data+ '(("Var"
			("Name"
			 "r32"
			 "r64"
			 "c64"
			 "c128")
			("Type"
			 "float"
			 "double"
			 "float _complex"
			 "double _complex"))
		       ("Ind"
			("Name"
			 "i32"
			 "sizet")
			("Type"
			 "int32_t"
			 "size_t"))
		       ("Size"
			("Name"
			 "i32"
			 "sizet")
			("Type"
			 "int32_t"
			 "size_t"))))

(defun gen-pair (parent child val-spot)
  (if (endp val-spot)
      nil
    (cons (concatenate 'string parent ":" child)
	  (car val-spot))))

;;; DONE above here

(defun gen-group (grouphead)
  (let ((parentname (car grouphead))
	(childnames (heads (cdr grouphead)))
	(val-spots  (mapcar #'cdr (cdr grouphead))))
    (make-alist :var-names
		(mapcar #'(lambda (s) (concatenate 'string parentname ":" s))
			childnames)
		:cur-val-spots 
		val-spots
		:orig-val-spots
		val-spots)))

(defun gen-groups (lst)
  (mapcar #'gen-group lst))

(defun cur-pairs (groups)
  (labels ((pairs-list (grp)
	     (with-alist (var-names cur-val-spots) grp
	       (mapcar #'(lambda (v c) (cons v c)) var-names cur-val-spots))))
    (if (endp groups)
	nil
      (append (pairs-list (car groups)) (cur-pairs (cdr groups))))))

(defun advance-group (group &key (cyclicp t))
  (with-alist (var-names cur-val-spots orig-val-spots) group
    (if (endp var-names) 
	nil
      (if (endp (car cur-val-spots)) ; ran out of val-spots; return to start
	  (if cyclicp
	      (make-alist :var-names var-names
			  :cur-val-spots orig-val-spots
			  :orig-val-spots orig-val-spots)
	    (make-alist :var-names var-names
			:cur-val-spots 'end
			:orig-val-spots orig-val-spots))
	(make-alist :var-names var-names
		    :cur-val-spots (mapcar #'cdr cur-val-spots)
		    :orig-val-spots orig-val-spots)))))

(defun at-end-p (group)
  (with-alist (cur-val-spots) group
    (eq 'end cur-val-spots)))

(defun advance-groups (groups)
  (cond ((endp groups)
	 nil)
	((endp (cdr groups))
	 (advance-group (first groups)))
	((at-end-p (second groups))
	 (list (advance-group (first groups))
	       (advance-group (second groups))))
	(t
	 (cons (first groups) (advance-groups (cdr groups))))))

(defun map-combs (fn! data)
  (loop for groups = (gen-groups data)
	then (advance-groups groups)
	until (endp groups)
	do (fn! (cur-pairs groups))))
	   
(defun make-sed-subst (var-name var-value)
  (concatenate 'string "-e 's/${" var-name "}/" var-value "/'"))

(defun join (sep &rest args)
  (cond ((endp args)
	 "")
	((endp (cdr args))
	 (concatenate 'string (first args) sep (second args)))
	(t
	 (concatenate 'string (first args) sep (join sep (cdr args))))))

(defun invoke-sed! (inpath outpath &rest substitutions)
  (let* ((sedargs (join " " substitutions))
	 (command (format t "cat ~A | sed ~A > ~A" inpath sedargs outpath)))
    (if (zerop (ext:system command))
	t
	nil)))

(defun strcat (&rest args)
  (cond ((endp args)
	 "")
	(t
	 (concatenate 'string (car args) (strcat (cdr args))))))

(defun string-subst-at (in-string find-string repl-string start find-string-len)
  (let ((spos (search find-string in-string :test #'string= :start2 start)))
    (if (endp pos)
	(subseq in-string spos)
	(let ((epos (+ spos find-string-len)))
	  (strcat (subseq in-string 0 spos)
		  repl-string
		  (string-subst-at (subseq in-string epos)
				   find-string
				   epos
				   find-string-len))))))

(defun string-subst (in-string find-string repl-string)
  (string-subst-at in-string find-string repl-string 0 (length find-string)))

(defun make-outpath (path-tmpl pairs)
  (if (endp pairs)
      template
    (let ((var (caar pairs))
	  (val (cdar pairs)))
      ;; FIXME: if PATH-TMPL is a path rather than a string,
      ;; convert it to a string first.
      (make-outpath (string-subst path-tmpl (strcat "${" var "}") val)
		    (cdr pairs)))))
  

(defun generate-code! (inpath outpath-template data)
  (map-combs data
	     (lambda (pairs)
	       (invoke-sed! inpath 
			    (make-outpath outpath-template pairs)
			    (mapcar (lambda (pair) 
				      (make-sed-subst (car pair) 
						      (cdr pair))) 
				    pairs)))))
