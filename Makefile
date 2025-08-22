.PHONY: pull push

MSG ?= Auto-commit from makefile

pull:
	git pull
	Rscript -e "renv::restore()"
	R CMD INSTALL .

push:
	Rscript -e "renv::snapshot()"
	Rscript -e "roxygen2::roxygenize()"
	git add .
	git commit -m "$(MSG)"
	git push
