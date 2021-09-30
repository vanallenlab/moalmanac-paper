# Molecular Oncology Almanac Paper
Analyses for the Molecular Oncology Almanac publication. Some common code used by multiple notebooks, mostly style preferences for figures, can be found in `common/`. All figures produced in this repository can be found in `figures/` and tables under `tables/`. 

## Installation
Code in this repository uses Python 3.7 and there is one R script in `analyses/knowledge-bases/`. The R script contains the package installation. For Python,  we recommend using a [virtual environment](https://docs.python.org/3/tutorial/venv.html) and running Python with either [Anaconda](https://www.anaconda.com/download/) or  [Miniconda](https://conda.io/miniconda.html). After installing Anaconda or Miniconda, you can set up by running

```bash
conda create -y -n moalmanac-paper python=3.7
conda activate moalmanac-paper
pip install -r requirements.txt
ipython kernel install --user --name=moalmanac-paper
```
## Setting font to arial
Several notebooks in this repository change the default font to Arial. The following command is used to install Arial as a font option for matplotlib,
```
conda install -n moalmanac-paper -c conda-forge mscorefonts
```
Afterwards, you will have to edit your `matplotlibrc` file for your jupyter notebook to uncomment line 207 and change Arial to the first item. For me on a macbook pro, this file was located here: `/Users/brendan/opt/miniconda3/envs/moalmanac-paper/lib/python3.7/site-packages/matplotlib/mpl-data/matplotlibrc`. [This guide from the Fowler lab](http://fowlerlab.org/2019/01/03/changing-the-sans-serif-font-to-helvetica/) was used to change font preferences with matplotlib. 


## Citation
Please cite our paper if using any information or code from this repository  
> [Reardon, B., Moore, N.D., Moore, N.S., *et al*. Integrating molecular profiles into clinical frameworks through the Molecular Oncology Almanac to prospectively guide precision oncology. *Nat Cancer* (2021). https://doi.org/10.1038/s43018-021-00243-3](https://www.nature.com/articles/s43018-021-00243-3)

You can also see prior iterations of this work from AACR abstracts over the years, 
- AACR 2019, [Abstract 2470: A molecular oncology almanac for integrative clinical interpretation of molecular profiles to guide precision cancer medicine](https://cancerres.aacrjournals.org/content/79/13_Supplement/2470.short)
- AACR 2018, [Abstract 2286: Feature-based clinical interpretation of whole exome and transcriptome data for precision cancer medicine](https://cancerres.aacrjournals.org/content/78/13_Supplement/2286.short)
- AACR 2017, [Abstract 558: Computational analysis of clinically actionable genomic features: precision heuristics for interpreting the alteration landscape (PHIAL)](https://cancerres.aacrjournals.org/content/77/13_Supplement/558.short) with the associated [github repo](https://github.com/vanallenlab/2017-aacr_phial2). 

