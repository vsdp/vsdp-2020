# Conversion from Jupyter Notebook format to Markdown (for the GitHub website)
# and to LaTeX for the published manual.

SRC_FILES       := $(wildcard *.ipynb)
LATEX_FILES     := $(patsubst %.ipynb,latex/%.tex,$(SRC_FILES))
WEBSITE_FILES   := $(patsubst %.ipynb,html/%.html,$(SRC_FILES)) html/README.md

all: $(LATEX_FILES) $(WEBSITE_FILES)

latex/%.tex: %.ipynb
	jupyter nbconvert --output-dir=latex --to latex    $^

html/%.html: %.ipynb
	jupyter nbconvert --output-dir=html  --to html     $^
	# Fix URLs for displayed in browser.
	sed -i -e 's/\.ipynb/.html/g' $@

html/%.md: %.ipynb
	jupyter nbconvert --output-dir=html  --to markdown $^
	# Fix URLs for displayed in browser.
	sed -i -e 's/\.ipynb/.html/g' $@
