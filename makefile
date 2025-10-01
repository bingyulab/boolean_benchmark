install:
	pip install --upgrade pip
# 	conda install -c potassco clingo
	pip install --force-reinstall git+https://github.com/hklarner/pyboolnet
	pip install --no-cache-dir numpy==1.25.2 pandas==1.4.0
	pip install -r requirements.txt
	Rscript Install/install.R

test:
	@if [ -f network_analysis.log ]; then rm network_analysis.log; fi
	python Step_03_Performance.py 

toy: 
	@if [ -f network_analysis.log ]; then rm network_analysis.log; fi
	python Step_03_Performance.py -p -d toy
	python Step_04_Plot.py -d toy

dream:
	@if [ -f network_analysis.log ]; then rm network_analysis.log; fi
	python Step_03_Performance.py -d dream -p
	python Step_04_Plot.py -d dream

TCell:
	@if [ -f network_analysis.log ]; then rm network_analysis.log; fi
	python Step_03_Performance.py -d TCell -p
	python Step_04_Plot.py -d TCell
.PHONY: install