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
    "amsmath"
    "url"
    "lineno"
    "trackchanges"
    "wasysym"
    "soul"
    "siunitx"
    "bm"
    "xr")
   (LaTeX-add-labels
    "sec:y"
    "eq:cibs"
    "fig:offset"
    "sec:prop"
    "eq:diffeq"
    "eq:diffeqspec"
    "eq:yeff"
    "eq:expand"
    "eq:modeleq"
    "sec:spatmode"
    "eq:surfbc"
    "eq:modeleqfin"
    "fig:svdmodes"
    "sec:solution"
    "eq:costfunc"
    "eq:twls"
    "eq:twlsCuu"
    "sec:som"
    "sec:solchar"
    "tab:cost"
    "fig:utilde"
    "sec:reconcore"
    "eq:solcov"
    "fig:recon"
    "sec:surfrecon"
    "eq:surfrecon"
    "eq:surfreconunc"
    "fig:meants"
    "fig:sst"
    "fig:ocean2kcomp"
    "fig:plankticstack"
    "tab:costalt"
    "fig:reconalt"
    "fig:meantsgm2"
    "sec:gammasens"
    "eq:2"
    "sec:agemodel"
    "eq:covariance"
    "eq:sumC"
    "sec:corr"
    "fig:td18O_cesm"
    "fig:td18O")
   (LaTeX-add-environments
    '("aligned" LaTeX-env-args ["argument"] 0)
    '("MT_gathered_env" LaTeX-env-args ["argument"] 0)
    '("#1*" LaTeX-env-args ["argument"] 0)
    '("smallmatrix*" LaTeX-env-args ["argument"] 0)
    '("Vmatrix*" LaTeX-env-args ["argument"] 0)
    '("vmatrix*" LaTeX-env-args ["argument"] 0)
    '("Bmatrix*" LaTeX-env-args ["argument"] 0)
    '("bmatrix*" LaTeX-env-args ["argument"] 0)
    '("pmatrix*" LaTeX-env-args ["argument"] 0)
    '("matrix*" LaTeX-env-args ["argument"] 0)
    '("multlined" LaTeX-env-args ["argument"] 0))
   (LaTeX-add-bibliographies))
 :latex)

