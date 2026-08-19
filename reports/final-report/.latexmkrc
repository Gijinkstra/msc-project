$pdf_mode = 1;          # use pdflatex (use 4 for lualatex, 5 for xelatex)
$out_dir  = 'build';
$aux_dir = 'build/auxiliary';
$synctex = 1;
# Should be either absent (latexmk auto-detects) or:
@default_files = ('main.tex');