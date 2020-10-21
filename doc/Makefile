# Makefile to work with JupyterBook

SRC=$(wildcard *.ipynb) $(wildcard *.yml) references.bib
TARGET=_build/html/index.html

all: $(TARGET)

$(TARGET): $(SRC)
	$(RM) -R _build
	jupyter-book build .

publish: $(TARGET)
	ghp-import --no-jekyll --push --force 
