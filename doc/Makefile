# Makefile to work with JupyterBook

SRC=$(wildcard *.ipynb) $(wildcard *.yml)

all: $(SRC) | jupyterbook
	( \
	  source jupyterbook/bin/activate; \
	  jupyter-book build .; \
	)
	#xdg-open _build/html/index.html

jupyterbook:
	python3 -m venv jupyterbook
	( \
	  source jupyterbook/bin/activate; \
	  pip install --upgrade pip jupyter-book ghp-import; \
	)
