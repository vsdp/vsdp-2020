SRC_FILES       := $(wildcard *.ipynb)
LATEX_FILES     := $(patsubst %.ipynb,latex/%.tex,$(SRC_FILES))
MARKDOWN_FILES  := $(patsubst %.ipynb,markdown/%.md,$(SRC_FILES))

all: $(LATEX_FILES) $(MARKDOWN_FILES)

latex/%.tex: %.ipynb
	jupyter nbconvert --output-dir=latex    --to latex    $^

markdown/%.md: %.ipynb
	jupyter nbconvert --output-dir=markdown --to markdown $^
