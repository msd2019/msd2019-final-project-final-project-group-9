DATA_DIRS = $(shell find data/ -type d)
DATA_FILES = $(shell find data/ -type f -name '*')

all: 05_final_report.pdf

download_data: 01_get_original_data.sh
	bash 01_get_original_data.sh

05_final_report.pdf: 02_final_report.Rmd 01_get_original_data.sh $(DATA_DIRS) $(DATA_FILES)
	Rscript -e "rmarkdown::render('02_final_report.Rmd')"

clean:
	rm -f data/original/*.txt
	rm -f data/original/*.zip
	rm -rf data/original/__MACOSX

