all: 05_final_report.pdf

original_data:
	bash 01_get_original_data.sh

05_final_report.pdf:
	Rscript -e "rmarkdown::{05_final_report.Rmd}"

clean:
	rm -f data/original/*.txt
	rm -f data/original/*.zip
	rm -rf data/original/__MACOSX

