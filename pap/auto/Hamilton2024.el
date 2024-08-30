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
    "soul"
    "siunitx"
    "bm"
    "xr")
   (LaTeX-add-labels
    "sec:ludata"
    "sec:overview"
    "eq:inverse"
    "sec:agemodel"
    "eq:covariance"
    "eq:cibs"
    "sec:tanddelta"
    "fig:offset"
    "sec:spatmode"
    "eq:swinv"
    "eq:svd"
    "eq:pervar"
    "fig:svdmodes"
    "sec:Cuu"
    "eq:Cuustats"
    "eq:corr"
    "sec:E"
    "eq:diffeq"
    "eq:2"
    "eq:impresp"
    "sec:solution"
    "eq:costfunc"
    "eq:twls"
    "eq:twlsCuu"
    "eq:solcov"
    "eq:surfrecon"
    "eq:surfreconunc"
    "tab:cost"
    "fig:utilde"
    "fig:recon"
    "sec:surfrecon"
    "fig:sst"
    "fig:plankticstack"
    "fig:meants"
    "fig:ocean2kcomp"
    "fig:EGWGSI")
   (LaTeX-add-bibliographies))
 :latex)

