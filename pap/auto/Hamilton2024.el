(TeX-add-style-hook
 "Hamilton2024"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("agujournal2019" "draft")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("trackchanges" "inline")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (TeX-run-style-hooks
    "latex2e"
    "agujournal2019"
    "agujournal201910"
    "url"
    "lineno"
    "trackchanges"
    "wasysym"
    "soul")
   (LaTeX-add-labels
    "sec:ludata"
    "sec:overview"
    "eq:inverse"
    "eq:costfunc"
    "eq:twls"
    "eq:twlsCuu"
    "sec:agemodel"
    "eq:covariance"
    "sec:tanddelta"
    "eq:cibs"
    "fig:offset"
    "sec:spatmode"
    "eq:swinv"
    "fig:svdmodes"
    "sec:Cuu"
    "eq:Cuustats"
    "eq:offdiags"
    "sec:E"
    "eq:diffeq"
    "eq:impresp"
    "eq:solcov"
    "eq:surfrecon"
    "fig:mags"
    "fig:recon"
    "fig:labcomp"
    "fig:sst"
    "fig:meants"
    "fig:EGWGSI"
    "sec:oc2k")
   (LaTeX-add-bibliographies))
 :latex)

