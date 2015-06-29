;;; apt-mode.el -- SPT code editing
;;; 
;;; based on awk-mode.el
;;;
;;; drop the code into a place where emacs can load the lisp
;;;  do "M-x byte-compile"
;;;
;;; add the following two line to the .emacs file
;;; (autoload 'spt-mode "~/.my_emacs_addons/spt-mode.el" "Enter Spt mode." t)
;;; (setq auto-mode-alist (cons '("\\.spt\\'" . spt-mode) auto-mode-alist))
;;; 
;;; make sure emacs knows where to load customerized .el files by adding
;;; the following line in the $HOME/.emacs file
;;; (setq load-path
;;;      (append (list nil 
;;;                    "~/.my_emacs_addons")
;;;              load-path))
;;;
;;; Limitations:
;;;   1) only accepts "#" for comments, don't know how to add "%" and "*"
;;;   2) the keywords are case sensitive, so only all lower case and 
;;;      upper-case for the first letter is implemented. don't know how to 
;;;      code keywords to case insensitive
;;;   3) indention doesn't work all the time

;;; Code:

(defvar spt-mode-syntax-table nil
  "Syntax table in use in Spt-mode buffers.")

(if spt-mode-syntax-table
    ()
  (setq spt-mode-syntax-table (make-syntax-table))
  (modify-syntax-entry ?\\ "\\" spt-mode-syntax-table)
  (modify-syntax-entry ?\n ">   " spt-mode-syntax-table)
  (modify-syntax-entry ?\f ">   " spt-mode-syntax-table)
  (modify-syntax-entry ?\# "<   " spt-mode-syntax-table)
  (modify-syntax-entry ?/ "." spt-mode-syntax-table)
  (modify-syntax-entry ?* "." spt-mode-syntax-table)
  (modify-syntax-entry ?+ "." spt-mode-syntax-table)
  (modify-syntax-entry ?- "." spt-mode-syntax-table)
  (modify-syntax-entry ?= "." spt-mode-syntax-table)
  (modify-syntax-entry ?% "." spt-mode-syntax-table)
  (modify-syntax-entry ?< "." spt-mode-syntax-table)
  (modify-syntax-entry ?> "." spt-mode-syntax-table)
  (modify-syntax-entry ?& "." spt-mode-syntax-table)
  (modify-syntax-entry ?| "." spt-mode-syntax-table)
  (modify-syntax-entry ?_ "_" spt-mode-syntax-table)
  (modify-syntax-entry ?\' "\"" spt-mode-syntax-table))


;; Regexps written with help from Peter Galbraith <galbraith@mixing.qc.dfo.ca>.
(defconst spt-font-lock-keywords
  (eval-when-compile
    (list
     ;;
     ;; Function names.
     '("^[ \t]*\\(function\\)\\>[ \t]*\\(\\sw+\\)?"
       (1 font-lock-keyword-face) (2 font-lock-function-name-face nil t))
     ;;
     ;; Variable names.
     (cons (regexp-opt
	    '("Second" "second" "Minute" "minute" "Hour" "hour"
	      "area" "Area" "depth" "Depth"
	      ) 'words)
	   'font-lock-builtin-face)
     ;;
     ;; Keywords.
     (regexp-opt
      '("DEF" "END" "def" "end" "Def" "End" "options" "node" "segment" "junction"
	"qsource" "lateralsource" "boundarycondition" "Options" "Node" "Segment" "Junction" "QSource"
	"LateralSource" "BoundaryCondition" "TimeSeries" "timeseries"
	"Trapezoidal" "trapezoidal" "Rectangular" "rectangular" "XY" "xy" "INTRINSIC" "intrinsic"
	) 'words)
     ;;
     ;; Builtins.
     (list (regexp-opt
	    '("metric" "Metric" "id" "Id" "BottomWidth" "bottomwidth"
	      "verbose" "Verbose" "epoch" "Epoch"
	      "Slope" "slope" "StopTime"
	      "stoptime" "StopTimeUnit" "stoptimeunit" "TimeStep" "timestep" "TimeStepUnit"
	      "timestepunit" "ssfile" "TimeUnit" "timeunit" 
	      "PrtInterval" "prtinterval" "PrtIntervalUnit" "prtintervalunit"
	      "PrtStart" "prtstart" "PrtStartUnit" "prtstartunit"
	      "sr" "SR" "n" "N" "zr" "ZR" "hr" "HR" "T" "t" "V" "v"
	      "Location" "location" "up" "Up" "Down" "down" "Length" "length"
	      "Type" "type" "x" "X" "y" "Y" "A" "a" "P" "p" "W" "w"
	      "up1" "Up1" "up2" "Up2" "coeff1" "Coeff1" "coeff2" "Coeff2"
	      "xcoord" "XCoord" "ycoord" "YCoord"
	      "lmax" "LMax" "lmin" "LMin"
	      "prtdepth" "PrtDepth" "prtsurfelev" "PrtSurfElev" "prtq" "PrtQ" "prta" "PrtA" "prtcoord" "PrtCoord"
	      "prtxy" "Prtxy" "PrtXY"
	      "checkonly" "CheckOnly"
	      ) 'words)
	   1 'font-lock-variable-name-face)
     ;;
     ;; Operators.  Is this too much?
     (cons (regexp-opt '("&&" "||" "<=" "<" ">=" ">" "==" "!=" "!~" "~"))
	   'font-lock-constant-face)
     ))
 "Default expressions to highlight in SPT mode.")

;;;###autoload
(define-derived-mode spt-mode sh-mode "SPT"
  "Major mode for editing SPT code.
This is much like C mode except for the syntax of comments.  Its keymap
inherits from C mode's and it has the same variables for customizing
indentation.  It has its own abbrev table and its own syntax table.

Turning on SPT mode runs `spt-mode-hook'."
  (set (make-local-variable 'paragraph-start) (concat "$\\|" page-delimiter))
  (set (make-local-variable 'paragraph-separate) paragraph-start)
  (set (make-local-variable 'comment-start) "# ")
  (set (make-local-variable 'comment-end) "")
  (set (make-local-variable 'comment-start-skip) "#+ *")
  (setq font-lock-defaults '(spt-font-lock-keywords nil nil ((?_ . "w"))))
  (setq font-lock-keywords-case-fold-search t)
)
(provide 'spt-mode)

;;; spt-mode.el ends here
