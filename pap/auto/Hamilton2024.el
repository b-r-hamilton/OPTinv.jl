(TeX-add-style-hook
 "Hamilton2024"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("agujournal2019" "draft")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("trackchanges" "inline")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "agujournal2019"
    "agujournal201910"
    "url"
    "lineno"
    "trackchanges"
    "wasysym"
    "soul"
    "siunitx"
    "bm"
    "xr")
   (LaTeX-add-labels
    "sec:spatmode"
    "fig:svdmodes"
    "sec:prop"
    "eq:diffeq"
    "sec:agemodel"
    "eq:cibs"
    "fig:offset"
    "sec:E"
    "eq:impresp"
    "eq:inverse"
    "sec:solution"
    "eq:costfunc"
    "eq:twls"
    "eq:twlsCuu"
    "sec:solchar"
    "fig:utilde"
    "tab:cost"
    "sec:reconcore"
    "eq:solcov"
    "fig:recon"
    "sec:surfrecon"
    "eq:surfrecon"
    "eq:surfreconunc"
    "fig:meants"
    "fig:sst"
    "fig:plankticstack"
    "fig:ocean2kcomp"
    "sec:singlemode"
    "fig:reconalt"
    "fig:meantsgm2"
    "fig:meantsgh19")
   (LaTeX-add-bibliographies))
 :latex)

