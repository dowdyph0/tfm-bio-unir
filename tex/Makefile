all:
	xelatex --shell-escape work.tex
	biber work
	xelatex --shell-escape work.tex
	xelatex --shell-escape work.tex

clean:
	rm -f *.aux *.bbl *.bcf *.blg *.log *.run.xml *.fdb_latexmk *.fls
