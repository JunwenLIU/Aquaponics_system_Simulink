DATE:=$(shell date '+%B %d, %Y')
TITLE:=Sample Title
AUTHOR:=Sean Lane
default: *.md
	pandoc -o output.pdf *.md -V geometry:margin=1in --highlight-style=tango  -V author='$(AUTHOR)' -V date='$(DATE)' -V title='$(TITLE)' # -V classoption=twocolumn

linux: *.md
	pandoc -o output.pdf *.md -V geometry:margin=1in --highlight-style=tango -V author='$(AUTHOR)' -V date='$(DATE)' -V title='$(DATE)' # --variable classoption=twocolumn
	evince output.pdf

latex: *.md
	pandoc -o output.tex *.md -V geometry:margin=1in  --highlight-style=tango --template default.latex -V author='$(AUTHOR)' -V date='$(DATE)' -V title='$(DATE)'
