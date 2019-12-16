all: README.md

biblio.png: biblio.bib
	./authorsgraph.py $^ | neato -T png -o $@

README.md: biblio.bib pdf biblio.png
	./bibdisplay.py > $@

clean:
	rm README.md biblio.png
