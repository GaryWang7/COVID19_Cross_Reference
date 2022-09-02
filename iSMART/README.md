# Readme
This folder contains adapted versions of iSMART and GIANA (GIANA is a newer algorithm from the same group that runs faster with comparable clustering performance). 

Both algorithms has been modified so it can be run in a Windows 10 conda environment with Python 3. We thank the authors for open access to their algorithms. Please refer to their publications for more details.

iSMART: https://aacrjournals.org/clincancerres/article/26/6/1359/268472/Investigation-of-Antigen-Specific-T-Cell-Receptor

GIANA: https://www.nature.com/articles/s41467-021-25006-7

# Links to original packages
To run iSMART, follow instructions from the authors (https://github.com/s175573/iSMART) to format your TCR file. 

The original package of GIANA can be found at https://github.com/s175573/GIANA.

# Running the codes:
## iSMART
Download the Imgt_Human_TRBV.fasta and iSMARTv3_Rev2.py files into the same folder, and run iSMART according to the authors' instructions.

One cmd example: 
```
python iSMARTv3_Rev2.py -f tcr_info_iSMART.tsv -t 3 -o./iSMART/
```
## GIANA
Download all files in /GIANA folder into the same folder, and run the code according to the authors' instructions.
