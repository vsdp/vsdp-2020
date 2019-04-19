# Conversion from Jupyter Notebook format to Markdown (for the GitHub website)
# and to LaTeX for the published manual.

SRC_FILES       := $(wildcard *.ipynb)
LATEX_FILES     := $(patsubst %.ipynb,latex/%.tex,$(SRC_FILES))
MARKDOWN_FILES  := $(patsubst %.ipynb,markdown/%.md,$(SRC_FILES))

all: $(LATEX_FILES) $(MARKDOWN_FILES)

latex/%.tex: %.ipynb
	jupyter nbconvert --output-dir=latex    --to latex    $^

markdown/%.md: %.ipynb
	jupyter nbconvert --output-dir=markdown --to markdown $^
	# GitHub pages does not support Octave syntax highlighting.
	sed -i -e 's/^```octave/```matlab/g' $@
