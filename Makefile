.PHONY: pull push

pull:
	git pull
	Rscript -e "renv::restore()"

push:
	Rscript -e "renv::snapshot()"
	git add .
	git commit -m "Update project dependencies"
	git push