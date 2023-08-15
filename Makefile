SLIDE_QMD_FILES := $(wildcard static/slides/*.qmd)
SLIDE_HTML_FILES  := $(subst qmd,html,$(SLIDE_QMD_FILES))
SLIDE_PDF_FILES  := $(subst qmd,pdf,$(SLIDE_QMD_FILES))

.PHONY: clean push build all pdf clean_docs clean_html clean_pdf

build: $(SLIDE_HTML_FILES) $(SLIDE_PDF_FILES)
	hugo
	rm -rf docs/slides/prev
	rm -rf docs/slides/*.rds
	rm -rf docs/slides/*.Rdata
	rm -rf docs/*_cache
	rm -rf docs/wip
	
all: pdf build

html: $(SLIDE_HTML_FILES)
	echo $(SLIDE_HTML_FILES)

pdf: $(SLIDE_PDF_FILES)
	echo $(SLIDE_PDF_FILES)

open: build
	open docs/index.html

clean_docs:
	rm -rf docs/

clean_pdf:
	rm -f static/slides/*.html
	
clean_html:
	rm -f static/slides/*.pdf

clean: clean_docs clean_pdf clean_html
	
static/slides/%.html: static/slides/%.qmd
	quarto render $<
	
static/slides/%.pdf: static/slides/%.html
	Rscript -e "renderthis::to_pdf('$<')"

push: build
	git pull
	git add .
	git commit -m "Make update"
	git push