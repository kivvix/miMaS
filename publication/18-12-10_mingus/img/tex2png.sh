#! /bin/bash

pdflatex -shell-escape <<EOF
\documentclass[preview]{standalone}
\usepackage{graphicx}
\begin{document}
\input{ordre_weno.tex}
\end{document}
EOF

#pdflatex -shell-escape "$input"

