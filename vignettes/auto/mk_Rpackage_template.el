(TeX-add-style-hook
 "mk_Rpackage_template"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "$if(fontsize)$$fontsize$" "$endif$$if(lang)$$babel-lang$" "$endif$$if(papersize)$$papersize$" "$endif$$for(classoption)$$classoption$$sep$" "$endfor$")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "$geometry$") ("$fontfamily$" "$fontfamilyoptions$") ("newtxmath" "cmintegrals" "cmbraces" "libertine") ("fontenc" "$if(fontenc)$$fontenc$$else$T1$endif$") ("inputenc" "utf8") ("underscore" "strings") ("hyperref" "setpagesize=false" "unicode=false" "xetex" "unicode=true")))
   (add-to-list 'LaTeX-verbatim-environments-local "code")
   (add-to-list 'LaTeX-verbatim-environments-local "VerbatimOut")
   (add-to-list 'LaTeX-verbatim-environments-local "SaveVerbatim")
   (add-to-list 'LaTeX-verbatim-environments-local "LVerbatim")
   (add-to-list 'LaTeX-verbatim-environments-local "BVerbatim")
   (add-to-list 'LaTeX-verbatim-environments-local "Verbatim")
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "Verb")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "Verb")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "geometry"
    "$fontfamily$"
    "libertine"
    "newtxmath"
    "fontenc"
    "inputenc"
    "listings"
    "fancyvrb"
    "longtable"
    "booktabs"
    "graphicx"
    "grffile"
    "titlesec"
    "natbib"
    "underscore"
    "setspace"
    "hyperref"
    "xcolor"
    "endnotes")
   (TeX-add-symbols
    "tightlist"
    "maxwidth"
    "maxheight"
    "footnote")
   (LaTeX-add-environments
    '("code"))
   (LaTeX-add-bibliographies
    "$for(bibliography)$$bibliography$$sep$"
    "$endfor$")
   (LaTeX-add-xcolor-definecolors
    "webbrown"))
 :latex)

