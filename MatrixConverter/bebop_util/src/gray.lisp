(defclass c-input-stram (gray:fundamental-character-input-stream)
  ((file-ptr :initarg :file-ptr :reader :file-ptr)))




/**
 * Dynamically create a Gray stream class that wraps a C output stream
 * (FILE*).  The idea and framework for the code come from Dustin
 * Long's "EclGui", a project that uses ECL as the back end of a
 * Windows Lisp REPL widget:
 * 
 * http://www.progmatism.com/hacks/eclgui/index.php
 */
void
ecl_create_gray_stream_class ()
{
#define SYMBOL_READ "%TERM-IO-READ-CHAR%"
#define SYMBOL_FLUSH "%TERM-IO-FLUSH%"

  /* 
   * Define ECL callbacks that implement Gray stream functionality 
   */
  cl_def_c_function (cl_intern(1, make_simple_base_string(SYMBOL_WRITE)),
		     (void*) &fgetc, 1)
  cl_def_c_function (cl_intern(1, make_simple_base_string(SYMBOL_FLUSH)),
		     (void*) &fflush, 0);
  ecl_run ("(defclass c-input-stream "
	   "  (gray:fundamental-character-input-stream) "
	   "    (last-char))", 1);
  ecl_run ("(defmethod gray:stream-write-char "
           "  ((stream %term-io%) c) "
           "    (" SYMBOL_WRITE " c))", 1);
  ecl_run ("(defmethod gray:stream-force-output "
	    "  ((stream %term-io%)) "
	    "    (" SYMBOL_FLUSH ") t)", 1);
  ecl_run ("(setq *terminal-io* "
           "  (make-two-way-stream "
           "    (make-string-input-stream \"\") "
           "    (make-instance '%term-io%)))", 1);
  ecl_run ("(setq *error-output* "
           "  (two-way-stream-output-stream *terminal-io*))", 1);
#undef SYMBOL_FLUSH
#undef SYMBOL_WRITE
}



(defmethod gray:stream-unread-char ((stream c-input-stream) character)
  (progn 
    ;; We ignore the return value (ungetc() may fail, in which case it
    ;; returns EOF).
    (%ungetc% character (file-ptr stream))
    nil))

(defmethod gray:stream-peek-char ((stream c-input-stream))
  (with-slots ((ptr file-ptr))
	      (let ((c (%fgetc% ptr)))
		(if (c-eof-p c)
		    :eof
		  (if (c-eof-p (%ungetc% c ptr))
		      :eof
		    c)))))

(defmethod gray:stream-read-char ((stream c-input-stream))
  (let ((c (%fgetc% (file-ptr stream))))
    (if (c-eof-p c)
	:eof
      c)))

;; We assume for now that the stream is noninteractive; see the CL
;; Hyperspec entry at
;;
;; http://www.lisp.org/HyperSpec/Body/fun_listen.html
;;
(defmethod gray:stream-listen ((stream c-input-stream))
  (if (%feof (file-ptr stream))
      nil
    t))

