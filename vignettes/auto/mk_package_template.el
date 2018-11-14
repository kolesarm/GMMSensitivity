(TeX-add-style-hook
 "mk_package_template"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "$if(fontsize)$$fontsize$" "$endif$$if(lang)$$babel-lang$" "$endif$$if(papersize)$$papersize$" "$endif$$for(classoption)$$classoption$$sep$" "$endfor$")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "left=1in" "top=1in" "right=1in" "bottom=1in") ("$fontfamily$" "$fontfamilyoptions$") ("fontenc" "$if(fontenc)$$fontenc$$else$T1$endif$") ("inputenc" "utf8") ("underscore" "strings") ("hyperref" "setpagesize=false" "unicode=false" "xetex" "unicode=true") ("color" "usenames" "dvipsnames")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "geometry"
    "$fontfamily$"
    "lmodern"
    "fontenc"
    "inputenc"
    "abstract"
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
    "color"
    "endnotes")
   (TeX-add-symbols
    "authorfont"
    "tightlist"
    "maxwidth"
    "maxheight"
    "footnote")
   (LaTeX-add-environments
    "hypothesis")
   (LaTeX-add-bibliographies
    "$for(bibliography)$$bibliography$$sep$"
    "$endfor$"))
 :latex)

