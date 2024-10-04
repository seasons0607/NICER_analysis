# -*- coding: utf-8 -*-
# Last modified 2024-07-03
from astropy.io import fits
from matplotlib import gridspec
from astropy.stats import bayesian_blocks
import os
import sys
import numpy as np
import pandas as pd 
import random
import subprocess
from astropy.time import Time
from datetime import datetime, timedelta, timezone
import matplotlib.pyplot as plt 
import matplotlib.dates as dates

def day_to_sec(d):
    return d * 24.0 * 60.0 * 60.0
def sec_to_min(s):
    return s / 60.0
def sec_to_day(s):
    return s / (60.0 * 60.0 * 24.0)
def print_red(text):
    print("\033[91m" + text + "\033[0m")


class RawData():
    """
    Step 0. make analysis directory (make_dir)
    Step 1. screening data (nicerl2)
    Step 2. barycentric correction (barycorr)
    Step 3. making band-sliced light curves with xselect (make_band_sliced_lc, make_band_sliced_lc_bary) 
    Step 4. plotting light curves and making the GTI file (plot_qlcurves)
    Step 5. making a Obs-ID-averaged spectrum (source/background), arf, and rmf (make_all_GTI_spectrum)
    Step 6. making GTI-divided GTI files, event files, spectra (source/background), arf, and rmf (make_all_GTI_spectrum)
    """
    def __init__(self, stellar_type, stellar_name, ra, dec, obsid, obsid_dir):	
        self.stellar_type = stellar_type		
        self.stellar_name = stellar_name
        self.ra = ra
        self.dec = dec
        self.obsid = obsid
        self.obsid_dir = obsid_dir

    def make_dir(self):
        print("-----[START:"+str(sys._getframe().f_code.co_name)+"]----")

        if os.path.isfile(self.obsid_dir+"/analysis") == False:
            cmd = "mkdir "+str(self.obsid_dir)+"/analysis\n"
            cmd += "mkdir "+str(self.obsid_dir)+"/analysis/lc\n"
            cmd += "mkdir "+str(self.obsid_dir)+"/analysis/spec\n"
            print(cmd)
            os.system(cmd)


        print("-----[END:"+str(sys._getframe().f_code.co_name)+"]----")

    def nicerl2(self):
        print("-----[START:"+str(sys._getframe().f_code.co_name)+"]----")

        cmd = "nicerl2 indir="+str(self.obsid_dir)+" overonly_range=0-5 cor_range='1.5-*' clobber=YES niprefilter2_coltypes=base,3c50 detlist=launch,-14,-34\n"
        # cmd = "nicerl2 indir="+str(self.obsid_dir)+" clobber=YES niprefilter2_coltypes=base,3c50 detlist=launch,-14,-34\n"
        print(cmd)
        os.system(cmd)

        print("-----[END:"+str(sys._getframe().f_code.co_name)+"]----")

    def barycorr(self):
        print("-----[START:"+str(sys._getframe().f_code.co_name)+"]----")

        cmd = "barycorr infile="+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"_0mpu7_cl.evt outfile="+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"_0mpu7_cl_bary.evt orbitfiles="+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".orb.gz ra="+str(self.ra)+" dec="+str(self.dec)+" refframe=ICRS ephem=JPLEPH.430 barytime=no clobber=yes\n"
        cmd += "barycorr infile="+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"_0mpu7_ufa.evt outfile="+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"_0mpu7_ufa_bary.evt orbitfiles="+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".orb.gz ra="+str(self.ra)+" dec="+str(self.dec)+" refframe=ICRS ephem=JPLEPH.430 barytime=no clobber=yes\n"
        print(cmd)
        os.system(cmd)

        print("-----[END:"+str(sys._getframe().f_code.co_name)+"]----")

    def make_band_sliced_lc(self):
        print("-----[START:"+str(sys._getframe().f_code.co_name)+"]----")

        cmd = "cd "+str(self.obsid_dir)+"/xti/event_cl/\n"
        if os.path.isfile(str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"_0mpu7_cl.evt.gz") == True and os.path.isfile(str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"_0mpu7_cl.evt") == False:
                cmd += "gunzip -k "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"_0mpu7_cl.evt.gz\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "30\n"
        cmd += "400\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_0340.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_ufa.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "30\n"
        cmd += "400\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_0340_ufa.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "30\n"
        cmd += "60\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_0306.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "60\n"
        cmd += "120\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_0612.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "120\n"
        cmd += "200\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_1220.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "200\n"
        cmd += "400\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_2040.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "30\n"
        cmd += "120\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_0312.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "120\n"
        cmd += "400\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_1240.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "mv lc_*.fits "+str(self.obsid_dir)+"/analysis/lc"

        print(cmd)
        os.system(cmd)

        print("-----[END:"+str(sys._getframe().f_code.co_name)+"]----")
    
    def make_band_sliced_lc_bary(self):
        print("-----[START:"+str(sys._getframe().f_code.co_name)+"]----")

        cmd = "cd "+str(self.obsid_dir)+"/xti/event_cl/\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl_bary.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "30\n"
        cmd += "400\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_0340_bary.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl_bary.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "30\n"
        cmd += "60\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_0306_bary.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl_bary.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "60\n"
        cmd += "120\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_0612_bary.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl_bary.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "120\n"
        cmd += "200\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_1220_bary.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl_bary.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "200\n"
        cmd += "400\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_2040_bary.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl_bary.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "30\n"
        cmd += "120\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_0312_bary.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "xselect <<EOF\n"
        cmd += str(random.randint(1, 1000000))+"select\n"
        cmd += "read event\n"
        cmd += "./\n"
        cmd += "ni"+str(self.obsid)+"_0mpu7_cl_bary.evt\n"
        cmd += "yes\n"
        cmd += "filter pha_cutoff\n"
        cmd += "120\n"
        cmd += "400\n"
        cmd += "extract curve binsize=64.0\n"
        cmd += "save curve\n"
        cmd += "lc_"+str(self.obsid)+"_1240_bary.fits\n"
        cmd += "exit\n"
        cmd += "no\n"
        cmd += "EOF\n"

        cmd += "mv lc_*_bary.fits "+str(self.obsid_dir)+"/analysis/lc"

        print(cmd)
        os.system(cmd)

        print("-----[END:"+str(sys._getframe().f_code.co_name)+"]----")

    def plot_qlcurves(self):
        print("-----[START:"+str(sys._getframe().f_code.co_name)+"]----")

        #Reading files
        filename_0340 = str(self.obsid_dir)+"/analysis/lc/lc_"+str(self.obsid)+"_0340.fits"
        fitsFile_0340 = fits.open(filename_0340)
        header_0340 = fitsFile_0340[0].header
        data_0340 = fitsFile_0340[1].data
        mjd_start_0340 = header_0340["MJD-OBS"]
        t_start_0340 = header_0340["TSTART"]
        sec_0340 = data_0340.field("time")
        rate_0340 = data_0340.field("RATE")
        rate_err_0340 = data_0340.field("ERROR")

        filename_0306 = str(self.obsid_dir)+"/analysis/lc/lc_"+str(self.obsid)+"_0306.fits"
        fitsFile_0306 = fits.open(filename_0306)
        header_0306 = fitsFile_0306[0].header
        data_0306 = fitsFile_0306[1].data
        mjd_start_0306 = header_0306["MJD-OBS"]
        sec_0306 = data_0306.field("time")
        rate_0306 = data_0306.field("RATE")
        rate_err_0306 = data_0306.field("ERROR")

        filename_0612 = str(self.obsid_dir)+"/analysis/lc/lc_"+str(self.obsid)+"_0612.fits"
        fitsFile_0612 = fits.open(filename_0612)
        header_0612 = fitsFile_0612[0].header
        data_0612 = fitsFile_0612[1].data
        mjd_start_0612 = header_0612["MJD-OBS"]
        sec_0612 = data_0612.field("time")
        rate_0612 = data_0612.field("RATE")
        rate_err_0612 = data_0612.field("ERROR")

        filename_1220 = str(self.obsid_dir)+"/analysis/lc/lc_"+str(self.obsid)+"_1220.fits"
        fitsFile_1220 = fits.open(filename_1220)
        header_1220 = fitsFile_1220[0].header
        data_1220 = fitsFile_1220[1].data
        mjd_start_1220 = header_1220["MJD-OBS"]
        sec_1220 = data_1220.field("time")
        rate_1220 = data_1220.field("RATE")
        rate_err_1220 = data_1220.field("ERROR")

        filename_2040 = str(self.obsid_dir)+"/analysis/lc/lc_"+str(self.obsid)+"_2040.fits"
        fitsFile_2040 = fits.open(filename_2040)
        header_2040 = fitsFile_2040[0].header
        data_2040 = fitsFile_2040[1].data
        mjd_start_2040 = header_2040["MJD-OBS"]
        sec_2040 = data_2040.field("time")
        rate_2040 = data_2040.field("RATE")
        rate_err_2040 = data_2040.field("ERROR")

        filename_0312 = str(self.obsid_dir)+"/analysis/lc/lc_"+str(self.obsid)+"_0312.fits"
        fitsFile_0312 = fits.open(filename_0312)
        header_0312 = fitsFile_0312[0].header
        data_0312 = fitsFile_0312[1].data
        mjd_start_0312 = header_0312["MJD-OBS"]
        sec_0312 = data_0312.field("time")
        rate_0312 = data_0312.field("RATE")
        rate_err_0312 = data_0312.field("ERROR")

        filename_1240 = str(self.obsid_dir)+"/analysis/lc/lc_"+str(self.obsid)+"_1240.fits"
        fitsFile_1240 = fits.open(filename_1240)
        header_1240 = fitsFile_1240[0].header
        data_1240 = fitsFile_1240[1].data
        mjd_start_1240 = header_1240["MJD-OBS"]
        sec_1240 = data_1240.field("time")
        rate_1240 = data_1240.field("RATE")
        rate_err_1240 = data_1240.field("ERROR")

        ##reading bary file for getting BJD
        filename_0340_bary = str(self.obsid_dir)+"/analysis/lc/lc_"+str(self.obsid)+"_0340_bary.fits"
        fitsFile_0340_bary = fits.open(filename_0340_bary)
        header_0340_bary = fitsFile_0340_bary[0].header
        mjd_start_0340_bary = header_0340_bary["MJD-OBS"]

        #Calculating hardness-ratio
        hr = [0 for i in range(len(rate_0312))]
        hr_err = [0 for i in range(len(rate_0312))]

        for i in range(len(hr)):
            numerator = rate_1240[i] - rate_0312[i]
            denominator = rate_1240[i] + rate_0312[i]
            result = numerator / denominator

            partial_derivative_numerator_1240 = 1
            partial_derivative_numerator_0312 = -1
            error_numerator = (partial_derivative_numerator_1240 * rate_err_1240[i]) + (partial_derivative_numerator_0312 * rate_err_0312[i])
            partial_derivative_denominator_1240 = 1
            partial_derivative_denominator_0312 = 1
            error_denominator = (partial_derivative_denominator_1240 * rate_err_1240[i]) + (partial_derivative_denominator_0312 * rate_err_0312[i])
            error_result =  result * ((error_numerator / numerator) - (error_denominator / denominator))

            hr[i] = result
            hr_err[i] = error_result
        
        #dividing by GTI (Note that nicerl2 may leave the incomplete mkf file)
        if os.path.exists(str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".att.gz") == True and os.path.exists(str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".att") == False:
                os.system("gunzip -k "+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".att.gz")
        if os.path.exists(str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".mkf.gz") == True and os.path.exists(str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".mkf") == False:
                os.system("niprefilter2 indir="+str(self.obsid_dir)+" infile="+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".mkf.gz outfile="+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".mkf clobber=YES")
        if os.path.exists(str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".att.gz") == True and os.path.exists(str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".att") == True:
                os.system("rm "+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".att")
                os.system("gunzip -k "+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".att.gz")
        if os.path.exists(str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".mkf.gz") == True and os.path.exists(str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".mkf") == True:
                os.system("rm "+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".mkf")
                os.system("niprefilter2 indir="+str(self.obsid_dir)+" infile="+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".mkf.gz outfile="+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".mkf clobber=YES")
        
        os.system("nimaketime "+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".mkf "+str(self.obsid_dir)+"/auxil/"+str(self.obsid)+"_all.gti mingti=100.0")
        if os.path.isfile(str(self.obsid_dir)+"/auxil/"+str(self.obsid)+"_all.gti") == True:
            filename_all_gti = str(self.obsid_dir)+"/auxil/"+str(self.obsid)+"_all.gti"
            fitsFile_all_gti = fits.open(filename_all_gti)
            header_all_gti = fitsFile_all_gti[0].header
            data_all_gti = fitsFile_all_gti[1].data
            mjd_start_all_gti = header_all_gti["MJDREFI"] + header_all_gti["MJDREFF"]
            t_start_all_gti = header_all_gti["TSTART"]
            gti_start = data_all_gti.field("START")
            gti_end = data_all_gti.field("STOP")
            gti_interval = [0 for i in range(len(gti_start))]
            for i in range(len(gti_start)):
                gti_interval[i] = gti_end[i] - gti_start[i]
                gti_start[i] = sec_to_min(gti_start[i] - t_start_0340)
                gti_end[i] = sec_to_min(gti_end[i] - t_start_0340)
        


        #Setting time (min)
        time_binsize_min = 64.0 / 60.0
        bin_number = len(sec_0340)

        min_0340 = [0 for i in range(bin_number)]
        min_0306 = [0 for i in range(bin_number)]
        min_0612 = [0 for i in range(bin_number)]
        min_1220 = [0 for i in range(bin_number)]
        min_2040 = [0 for i in range(bin_number)]
        min_0312 = [0 for i in range(bin_number)]

        for i in range(bin_number):
            min_0340[i] = sec_to_min(sec_0340[i])
            min_0306[i] = sec_to_min(sec_0306[i])
            min_0612[i] = sec_to_min(sec_0612[i])
            min_1220[i] = sec_to_min(sec_1220[i])
            min_2040[i] = sec_to_min(sec_2040[i])
            min_0312[i] = sec_to_min(sec_0312[i])

        #Plotting
        title = "Stellar Type: "+str(self.stellar_type)+" / Stellar Name: "+str(self.stellar_name)+"\n"
        title += "(ra, dec) = ("+str(self.ra)+", "+str(self.dec)+")\n"
        title += "Obs-ID: "+str(self.obsid)+" / Exposure: "+str(sum(gti_interval))+" (sec)\n"

        fig=plt.figure(figsize=(8.8, 10.0))
        gs = gridspec.GridSpec(6, 1, height_ratios=(1, 1, 1, 1, 1, 1))
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        axs = []
        for i in range(6):
            axs.append(fig.add_subplot(gs[i, 0])) 

        time_min = min([min(min_0340), min(gti_start)])
        time_max = max([max(min_0340), max(gti_end)])

        xstart = time_min - ((time_max-time_min) * 0.1)
        xend = time_max + ((time_max-time_min) * 0.1)
        label_fontsize = 14.5
        markersize = 6.1
        markeredgewidth = 1.7
        xy_fontsize = 13

        for ax in axs:
            ax.set_xlim(xstart, xend)
            ax.set_ylabel(r"cps", fontsize = label_fontsize)
            ax.tick_params(direction = "in", right = True, labelright = False, top = True, labeltop = False, labelleft=True, labelbottom = False, labelsize=xy_fontsize, which="both")
            ax.yaxis.set_label_coords(-0.1, 0.5)
            ax.minorticks_on()
            ax.grid(False, which="minor")
            ax.grid(True, which="major", color="lightgray", linewidth=1.0, linestyle="dashed", axis='x')
            if os.path.isfile(str(self.obsid_dir)+"/auxil/"+str(self.obsid)+"_all.gti") == True:
                for i in range(len(gti_start)):
                    ax.axvspan(gti_start[i]-(32.0/60.0), gti_end[i]+(32.0/60.0), color="cornflowerblue", alpha=0.1)
                    axs[0].text(((gti_start[i]+gti_end[i])/2.0 - xstart) / (xend - xstart), 1.1, str(i), fontsize=10, bbox=dict(boxstyle="round",facecolor="white",edgecolor="black",linewidth=1.8,linestyle="-"), ha="center", transform=axs[0].transAxes)
        
        axs[0].set_title(title, fontsize = 17, y=1.05)
        axs[0].errorbar(min_0340, rate_0340, xerr = time_binsize_min / 2.0, yerr = [abs(i) for i in rate_err_0340], capsize=0, fmt='s', markersize=markersize, ecolor='black', markeredgecolor = "black", markerfacecolor="palegreen", linewidth=markeredgewidth, markeredgewidth=markeredgewidth)
        axs[0].text(0.99, 0.89, r"0.3-4 keV", fontsize=label_fontsize, bbox=dict(boxstyle="round",facecolor="None",edgecolor="None",linewidth=1.3,linestyle="-"), transform=axs[0].transAxes, fontweight = "bold", ha="right", va="top")
        axs[1].errorbar(min_0306, rate_0306, xerr = time_binsize_min / 2.0, yerr = [abs(i) for i in rate_err_0306], capsize=0, fmt='s', markersize=markersize, ecolor='black', markeredgecolor = "black", markerfacecolor="aqua", linewidth=markeredgewidth, markeredgewidth=markeredgewidth)
        axs[1].text(0.99, 0.89, r"0.3-0.6 keV", fontsize=label_fontsize, bbox=dict(boxstyle="round",facecolor="None",edgecolor="None",linewidth=1.3,linestyle="-"), transform=axs[1].transAxes, fontweight = "bold", ha="right", va="top")
        axs[2].errorbar(min_0612, rate_0612, xerr = time_binsize_min / 2.0, yerr = [abs(i) for i in rate_err_0612], capsize=0, fmt='s', markersize=markersize, ecolor='black', markeredgecolor = "black", markerfacecolor="lightskyblue", linewidth=markeredgewidth, markeredgewidth=markeredgewidth)
        axs[2].text(0.99, 0.89, r"0.6-1.2 keV", fontsize=label_fontsize, bbox=dict(boxstyle="round",facecolor="None",edgecolor="None",linewidth=1.3,linestyle="-"), transform=axs[2].transAxes, fontweight = "bold", ha="right", va="top")
        axs[3].errorbar(min_1220, rate_1220, xerr = time_binsize_min / 2.0, yerr = [abs(i) for i in rate_err_1220], capsize=0, fmt='s', markersize=markersize, ecolor='black', markeredgecolor = "black", markerfacecolor="lightseagreen", linewidth=markeredgewidth, markeredgewidth=markeredgewidth)
        axs[3].text(0.99, 0.89, r"1.2-2.0 keV", fontsize=label_fontsize, bbox=dict(boxstyle="round",facecolor="None",edgecolor="None",linewidth=1.3,linestyle="-"), transform=axs[3].transAxes, fontweight = "bold", ha="right", va="top")
        axs[4].errorbar(min_2040, rate_2040, xerr = time_binsize_min / 2.0, yerr = [abs(i) for i in rate_err_2040], capsize=0, fmt='s', markersize=markersize, ecolor='black', markeredgecolor = "black", markerfacecolor="royalblue", linewidth=markeredgewidth, markeredgewidth=markeredgewidth)
        axs[4].text(0.99, 0.89, r"2.0-4.0 keV", fontsize=label_fontsize, bbox=dict(boxstyle="round",facecolor="None",edgecolor="None",linewidth=1.3,linestyle="-"), transform=axs[4].transAxes, fontweight = "bold", ha="right", va="top")
        axs[5].errorbar(min_0312, hr, xerr = time_binsize_min / 2.0, yerr = [abs(i) for i in hr_err], capsize=0, fmt='s', markersize=markersize, ecolor='black', markeredgecolor = "black", markerfacecolor="orangered", linewidth=markeredgewidth, markeredgewidth=markeredgewidth)
        axs[5].text(0.99, 0.89, r"$\mathbf{(H-S)/(H+S)}$", fontsize=label_fontsize, bbox=dict(boxstyle="round",facecolor="None",edgecolor="None",linewidth=1.3,linestyle="-"), transform=axs[5].transAxes, fontweight = "bold", ha="right", va="top")
        axs[5].tick_params(direction = "in", right = True, labelright = False, top = True, labeltop = False, labelleft=True, labelbottom = True, labelsize=xy_fontsize, which="both")
        axs[5].set_xlabel(r"Time from BJD = "+str(mjd_start_0340_bary+2400000.5)+" (min)", fontsize = label_fontsize)
        axs[5].set_ylabel(r"HR", fontsize = label_fontsize)

        plt.subplots_adjust(hspace=0, left=0.11, right=1-0.11, bottom=0.06, top=0.88)
        plt.savefig(str(self.obsid_dir)+"/analysis/lc_"+str(self.obsid)+".pdf", format="pdf", dpi=1000)


        print("-----[END:"+str(sys._getframe().f_code.co_name)+"]----")


    def make_all_GTI_spectrum(self):

        print("-----[START:"+str(sys._getframe().f_code.co_name)+"]----")

        cmd = "\n"


        cmd += "mkdir "+str(self.obsid_dir)+"/analysis/spec/block_all\n"
        print(cmd)
        os.system(cmd)

        try:
            subprocess.run("niextract-events "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"_0mpu7_cl.evt "+str(self.obsid_dir)+"/analysis/spec/block_all/"+str(self.obsid)+"_block_all_cl.evt timefile="+str(self.obsid_dir)+"/auxil/"+str(self.obsid)+"_all.gti", shell=True, check=True, timeout=60)
        except subprocess.TimeoutExpired:
            print_red("Timeout!!: niextract-events for cl file")
            os.system("rm -r "+str(self.obsid_dir)+"/analysis/spec/block_all")
        
        if os.path.isfile(str(self.obsid_dir)+"/analysis/spec/block_all/"+str(self.obsid)+"_block_all_cl.evt") == True:
            try:
                subprocess.run("niextract-events "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"_0mpu7_ufa.evt "+str(self.obsid_dir)+"/analysis/spec/block_all/"+str(self.obsid)+"_block_all_ufa.evt timefile="+str(self.obsid_dir)+"/auxil/"+str(self.obsid)+"_all.gti", shell=True, check=True, timeout=60)
            except subprocess.TimeoutExpired:
                print_red("Timeout!!: niextract-events for ufa file")
                os.system("rm -r "+str(self.obsid_dir)+"/analysis/spec/block_all")

        cmd = "\n"
        cmd += "cd "+str(self.obsid_dir)+"/analysis/spec/block_all\n"
        cmd += "cp "+str(self.obsid_dir)+"/auxil/"+str(self.obsid)+"_all.gti "+str(self.obsid_dir)+"/analysis/spec/block_all/\n"
        cmd += "cd "+str(self.obsid_dir)+"\n"
        cmd += "cd ..\n"
        cmd += "nicerl3-spect "+str(self.obsid)+" bkgmodeltype=3c50 clobber=YES grouptype=NONE\n"
        cmd += "mv "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"mpu7_bg.pha "+str(self.obsid_dir)+"/analysis/spec/block_all/\n"
        cmd += "mv "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"mpu7_bg.pha_bkg_day.pi "+str(self.obsid_dir)+"/analysis/spec/block_all/\n"
        cmd += "mv "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"mpu7_bg.pha_bkg_ngt.pi "+str(self.obsid_dir)+"/analysis/spec/block_all/\n"
        cmd += "mv "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"mpu7.arf "+str(self.obsid_dir)+"/analysis/spec/block_all/\n"
        cmd += "mv "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"mpu7_sk.arf "+str(self.obsid_dir)+"/analysis/spec/block_all/\n"
        cmd += "mv "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"mpu7_load.xcm "+str(self.obsid_dir)+"/analysis/spec/block_all/\n"
        cmd += "mv "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"mpu7_sr.pha "+str(self.obsid_dir)+"/analysis/spec/block_all/\n"
        cmd += "mv "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"mpu7.rmf "+str(self.obsid_dir)+"/analysis/spec/block_all/\n"
        cmd += "cd "+str(self.obsid_dir)+"/analysis/spec/block_all/\n"
        cmd += "mv ni"+str(self.obsid)+"mpu7_bg.pha "+str(self.obsid)+"_block_all_bkg.pha\n"
        cmd += "mv ni"+str(self.obsid)+"mpu7_bg.pha_bkg_day.pi "+str(self.obsid)+"_block_all_bkg.pha_bkg_day.pi\n"
        cmd += "mv ni"+str(self.obsid)+"mpu7_bg.pha_bkg_ngt.pi "+str(self.obsid)+"_block_all_bkg.pha_bkg_ngt.pi\n"
        cmd += "mv ni"+str(self.obsid)+"mpu7.arf "+str(self.obsid)+"_block_all.arf\n"
        cmd += "mv ni"+str(self.obsid)+"mpu7_sk.arf "+str(self.obsid)+"_block_all_ak.arf\n"
        cmd += "mv ni"+str(self.obsid)+"mpu7_load.xcm "+str(self.obsid)+"_block_all_load.xcm\n"
        cmd += "mv ni"+str(self.obsid)+"mpu7_sr.pha "+str(self.obsid)+"_block_all_tot.pha\n"
        cmd += "mv ni"+str(self.obsid)+"mpu7.rmf "+str(self.obsid)+"_block_all.rmf\n"

        if os.path.isfile(str(self.obsid_dir)+"/analysis/spec/block_all/"+str(self.obsid)+"_block_all_cl.evt") == True and os.path.isfile(str(self.obsid_dir)+"/analysis/spec/block_all/"+str(self.obsid)+"_block_all_ufa.evt") == True:
            print(cmd)
            os.system(cmd)

        
        print("-----[END:"+str(sys._getframe().f_code.co_name)+"]----")

    def make_GTI_divided_spectrum(self):

        print("-----[START:"+str(sys._getframe().f_code.co_name)+"]----")

        filename_all_gti = str(self.obsid_dir)+"/auxil/"+str(self.obsid)+"_all.gti"
        fitsFile_all_gti = fits.open(filename_all_gti)
        header_all_gti = fitsFile_all_gti[0].header
        data_all_gti = fitsFile_all_gti[1].data
        gti_start = data_all_gti.field("START")
        gti_end = data_all_gti.field("STOP")
        gti_interval = [0 for i in range(len(gti_start))]

        if len(gti_start) == 1:
            print("The number of GTI is one!!")
        
        if len(gti_start) > 1:
            for i in range(len(gti_start)):
                gti_interval[i] = gti_end[i] - gti_start[i]
                os.system("nimaketime "+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".mkf "+str(self.obsid_dir)+"/auxil/"+str(self.obsid)+"_block"+str(i)+".gti expr='((TIME >="+str(gti_start[i])+").and.(TIME <="+str(gti_end[i])+"))' history = YES")
            
            for i in range(len(gti_start)):
                cmd = "\n"
                cmd += "mkdir "+str(self.obsid_dir)+"/analysis/spec/block"+str(i)+"\n"

                print(cmd)
                os.system(cmd)

                try:
                    subprocess.run("niextract-events "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"_0mpu7_cl.evt "+str(self.obsid_dir)+"/analysis/spec/block"+str(i)+"/"+str(self.obsid)+"_block"+str(i)+"_cl.evt timefile="+str(self.obsid_dir)+"/auxil/"+str(self.obsid)+"_block"+str(i)+".gti", shell=True, check=True, timeout=60)
                except subprocess.TimeoutExpired:
                    print_red("Timeout!!: niextract-events for cl file")
                    os.system("rm -r "+str(self.obsid_dir)+"/analysis/spec/block"+str(i))
                
                if os.path.isfile(str(self.obsid_dir)+"/analysis/spec/block"+str(i)+"/"+str(self.obsid)+"_block"+str(i)+"_cl.evt") == True:
                    try:
                        subprocess.run("niextract-events "+str(self.obsid_dir)+"/xti/event_cl/ni"+str(self.obsid)+"_0mpu7_ufa.evt "+str(self.obsid_dir)+"/analysis/spec/block"+str(i)+"/"+str(self.obsid)+"_block"+str(i)+"_ufa.evt timefile="+str(self.obsid_dir)+"/auxil/"+str(self.obsid)+"_block"+str(i)+".gti", shell=True, check=True, timeout=60)
                    except subprocess.TimeoutExpired:
                        print_red("Timeout!!: niextract-events for cl file")
                        os.system("rm -r "+str(self.obsid_dir)+"/analysis/spec/block"+str(i))

                cmd = "\n"				
                cmd += "cd "+str(self.obsid_dir)+"/analysis/spec/block"+str(i)+"\n"
                cmd += "mv "+str(self.obsid_dir)+"/auxil/"+str(self.obsid)+"_block"+str(i)+".gti "+str(self.obsid_dir)+"/analysis/spec/block"+str(i)+"/\n"
                cmd += "nicerl3-spect indir='.' clfile="+str(self.obsid)+"_block"+str(i)+"_cl.evt ufafile="+str(self.obsid)+"_block"+str(i)+"_ufa.evt mkfile="+str(self.obsid_dir)+"/auxil/ni"+str(self.obsid)+".mkf phafile="+str(self.obsid)+"_block"+str(i)+"_tot.pi bkgfile="+str(self.obsid)+"_block"+str(i)+"_bkg.pi loadfile="+str(self.obsid)+"_block"+str(i)+"_load.xcm arffile="+str(self.obsid)+"_block"+str(i)+".arf skyarffile="+str(self.obsid)+"_block"+str(i)+"_sk.arf rmffile="+str(self.obsid)+"_block"+str(i)+".rmf bkgrmffile="+str(self.obsid)+"_block"+str(i)+"_bkg.rmf bkgmodeltype=3c50 clobber=yes grouptype=NONE detlist=launch,-14,-34\n"
                
                if os.path.isfile(str(self.obsid_dir)+"/analysis/spec/block"+str(i)+"/"+str(self.obsid)+"_block"+str(i)+"_cl.evt") == True and os.path.isfile(str(self.obsid_dir)+"/analysis/spec/block"+str(i)+"/"+str(self.obsid)+"_block"+str(i)+"_ufa.evt") == True:
                    print(cmd)
                    os.system(cmd)

        
        print("-----[END:"+str(sys._getframe().f_code.co_name)+"]----")


def main():
    source_type = input("Source Type: ")
    source_name = input("Source Name: ")
    ra = input("RA: ")
    dec = input("DEC: ")
    
    Obs_id = "7555020203" #input Obs-ID
    star_path = "/Users/Desktop/CR_Dra/7555020203" #input the path to the directory of the Obs-ID
    Star = RawData(str(source_type), str(source_name), ra, dec, Obs_id, star_path)
    Star.make_dir()
    Star.nicerl2()
    Star.barycorr()
    Star.make_band_sliced_lc()
    Star.make_band_sliced_lc_bary()
    Star.plot_qlcurves()
    Star.make_all_GTI_spectrum()
    Star.make_GTI_divided_spectrum()  



if __name__ == "__main__":
    main()
