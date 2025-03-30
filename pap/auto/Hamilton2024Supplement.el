(TeX-add-style-hook
 "Hamilton2024Supplement"
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
    "xr"
    "mathtools"
    "array"
    "longtable")
   (LaTeX-add-labels
    "sec:offset"
    "sec:ttm"
    "eq:mtau"
    "eq:tau"
    "eq:2"
    "eq:1"
    "eq:3"
    "sec:spatmode"
    "eq:ss"
    "eq:svd"
    "eq:pervar"
    "fig:U"
    "sec:irf"
    "fig:modemags"
    "sec:Cuu"
    "sec:tempvar"
    "eq:Cuustats"
    "fig:umodel"
    "sec:oxygenvar"
    "fig:cesmts"
    "fig:td18O"
    "fig:u0utilde"
    "sec:instrumental"
    "fig:labcomp"
    "sec:smest"
    "fig:impulse"
    "sec:gammasens"
    "eq:4"
    "fig:gamma"
    "var")
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
   (LaTeX-add-bibliographies
    "Hamilton2024"))
 :latex)

