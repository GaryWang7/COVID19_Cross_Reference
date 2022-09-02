This folder contains adapted versions of iSMART and GIANA (a newer algorithm from the same group). Both algorithms has been modified so it can be run in a Windows 10 conda environment with Python 3.

To run iSMART, follow instructions from the authors (https://github.com/s175573/iSMART) to format your TCR file. The GIANA original package can be found at https://github.com/s175573/GIANA.

One cmd example: 
python iSMARTv3_Rev2.py -f tcr_info_iSMART.tsv -t 3 -o./iSMART/
