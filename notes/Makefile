SRC_MD = $(wildcard sources/*.md)
BUILD_dir = journal
JOURNAL_file = index.html
PWD_file ?= $(shell pwd $(BUILD_dir)/$(JOURNAL_file))

SRC_CSS = $(wildcard template/style/*.css)
DES_CSS = $(SRC_CSS:template/style/%.css=$(BUILD_dir)/style/%.css)
SRC_JS = $(wildcard template/js/*.js)
DES_JS = $(SRC_JS:template/js/%.js=$(BUILD_dir)/js/%.js)

SRC_IMG = $(wildcard sources/img/*.png)
DES_IMG = $(SRC_IMG:sources/img/%.png=$(BUILD_dir)/img/%.png)

PANDOC_format = markdown yaml_metadata_block implicit_figures definition_lists example_lists pipe_tables table_captions smart inline_notes task_lists

all:
	@echo "make something"

# make directories ==========
$(BUILD_dir):
	@echo "\033[95m+ $@/\033[0m"
	mkdir -p $@

$(BUILD_dir)/style:
	@echo "\033[95m+ $@/\033[0m"
	mkdir -p $@

$(BUILD_dir)/js:
	@echo "\033[95m+ $@/\033[0m"
	mkdir -p $@

$(BUILD_dir)/img:
	@echo "\033[95m+ $@/\033[0m"
	mkdir -p $@

# build rules ===============
html: $(BUILD_dir) $(DES_CSS) $(DES_JS) $(SRC_MD) $(DES_IMG) $(BUILD_dir)/$(JOURNAL_file)
	@echo "\033[36m> Convert mardown files into a HTML file\033[0m"

$(BUILD_dir)/$(JOURNAL_file): $(SRC_MD) template/after.html
	@echo "\033[35m+ $@\033[0m"
	pandoc --wrap=preserve \
		--metadata-file=metadata.yaml \
		--css=style/boucdoc.css --css=style/pandoc.css \
		-f $(subst $(eval) ,+,$(PANDOC_format)) \
		--standalone -t html5 \
		--include-after-body=template/after.html \
		--katex=https://cdn.jsdelivr.net/npm/katex@0.10.0-beta/dist/ \
		$(filter %.md,$^) \
		-o $@

# copy css files from template directory
$(BUILD_dir)/style/%.css: template/style/%.css $(BUILD_dir)/style
	@echo "\033[35m+ $@\033[0m"
	cp $< $@

# copy js files from template directory
$(BUILD_dir)/js/%.js: template/js/%.js $(BUILD_dir)/js
	@echo "\033[35m+ $@\033[0m"
	cp $< $@

# copy png images from source/img directory
# TODO: make a rule to run all *.plt files before (gnuplot)
$(BUILD_dir)/img/%.png: sources/img/%.png $(BUILD_dir)/img
	@echo "\033[35m+ $@\033[0m"
	cp $< $@

# make file which contains inclusion of all js files
template/after.html: $(DES_JS)
	@echo "\033[35m+ $@\033[0m"
	touch $@
	echo "$(foreach js,$(DES_JS),<script type='text/javascript' src='js/$(notdir js/$(js))' ></script>)" >> $@

# misc rules ================
# make a new day if it doesn't exist
new: sources/$(shell date "+%Y-%m-%d").md
	@echo "\033[36m> New file generated\033[0m"

# make a new file with date as a title
sources/%.md:
	@echo "\033[35m+ $@\033[0m"
	echo "# $(basename $(notdir $@))\n\n" >> $@

# edit the newest file (but doesn't make a new file if it's a new day)
edit:
	vim $(shell ls -t sources/*.md | head -1)

# remove all $(BUILD_dir)
clear:
	@echo "\033[31m- $@\033[0m"
	rm -Rf $(BUILD_dir) template/after.html

# a test to send automatically to ftp server
USER = josmassot@perso.univ-rennes1.fr/josselin.massot
FILE = $(BUILD_dir)/$(JOURNAL_file)

