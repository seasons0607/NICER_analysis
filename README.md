# NICER automated extraction tool
## Overview
- This tool is a simple script to extract light curve and Obs-ID averaged and GTI-divided spectra of the NICER data.
- When you have use this script in your paper, please cite Inoue et al., 2024, MNRAS, ??, ??.
- If you have any bugs, please let me know!
--  📧 inoue *at* cr.scphys.kyoto-u.ac.jp

## How to use
1. Download the NICER data from the HEASARC archive
2. Write the Obs-ID and the path to the data in line 662 and 663 of <code>nicer.py</code>, respectively.
3. Run nicer.py as <code>python nicer.py</code> 
4. Input the some parameters (Source Type, Source Name, RA, DEC) and push Enter. You can input any string for Source Type and Source Name.

## Extracted products
- All products are outputed in the <code>/ObsID/analysis/</code>.
- <code>lc_ObsID.pdf</code> shows the band-sliced light curve and each GTI time interval.
- All light curve fits file are in <code>/ObsID/analysis/lc/</code>. The files with <code>_bary</code> mean that they are barycentric-corrected.
- The ObsID-averaged spectrum and response files are in <code>/ObsID/analysis/spec/block_all</code>.
- The GTI-divided spectra and response files are in <code>/ObsID/analysis/spec/block_Number</code>. If there is only one GTI in the ObsID, the GTI-divided spectra are not extracted.

## Extra Option
- If you change the criteria of the screening of <code>nicerl2</code>, you can add them in line 61. Please also see the help page of nicerl2.

## Enviroment
This code can be runned with Python 3 and Heasoft newer than 6.32.

## History
2024-10-04 new version 0.0 is created.
