(TeX-add-style-hook
 "OPTinv"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("geometry" "portrait" "margin=0.5in")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "inputenc"
    "indentfirst"
    "geometry"
    "array"
    "graphicx"
    "amsmath"
    "caption"
    "subcaption"
    "natbib"
    "textcomp"
    "siunitx"
    "hyperref"
    "wasysym"
    "setspace")
   (TeX-add-symbols
    "numberthis")
   (LaTeX-add-labels
    "fig:svd"
    "fig:nmf"
    "fig:recon"
    "fig:mags"
    "fig:sst"
    "fig:meants"
    "fig:labcomp"
    "fig:EGWGSI"
    "fig:post1900")
   (LaTeX-add-bibliographies))
 :latex)

