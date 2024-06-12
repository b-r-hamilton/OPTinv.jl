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
    "eq:covariance"
    "eq:cibs"
    "fig:offset"
    "eq:inverse"
    "sec:spatmode"
    "fig:svdmodes"
    "fig:recon"
    "fig:mags"
    "fig:labcomp"
    "fig:sst"
    "fig:meants"
    "fig:EGWGSI")
   (LaTeX-add-bibliographies))
 :latex)

