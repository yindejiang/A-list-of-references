import numpy as np
import astropy.io.fits as fits
import os
import datetime
import sys
from array import array
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import *

##############################################################
# 20161008 adapted from cut_FASTpsrfits_freq_time.py
#          output 2 pol and pol averaged data
# 20161009 dimension of DAT_OFFS changed from chnum*2 to chnum
#          format of DAT_OFFS changed from dataformat3 to dataformat2
#          size(float_data)/nline/nchan/npol=nsblk
##############################################################
print('Last updated 20161219')

#mpl.rcParams['image.interpolation']='none'
if len(sys.argv) < 5:
    print('too few inputs!')
    print('example:')
    print('python cut_FASTpsrfits.py startchan endchan startn endn FAST.fits')
    sys.exit()

starttime = datetime.datetime.now()

startfreq = int(sys.argv[1])
endfreq = int(sys.argv[2])
startn = int(sys.argv[3])
endn = int(sys.argv[4])
filename = sys.argv[5]

fileroot = filename[:-5]
print(fileroot)

hdulist = fits.open(filename)

hdu0 = hdulist[0]
data0 = hdu0.data
header0 = hdu0.header
print(data0)

hdu1 = hdulist[1]
data1 = hdu1.data
header1 = hdu1.header

nchan = header1['NCHAN']
nsblk = header1['NSBLK']
npol = header1['NPOL']
tbin = header1['TBIN']
chan_bw = header1['CHAN_BW']
print('check1:', nchan, nsblk, npol)

print(header0['OBSBW'])
print(header0['OBSNCHAN'])

chnum = endfreq - startfreq + 1
linenum = endn - startn + 1
freq = hdu0.header['OBSFREQ']
print('hehe', hdu0.header['OBSFREQ'])
hdu0.header['OBSFREQ'] = ((startfreq + endfreq) * 1.0 / 2 + 1.0) / ((nchan - 1.0) * 1.0 / 2 + 1.0) * freq
print('hehe', hdu0.header['OBSFREQ'])
hdu0.header['OBSBW'] = chnum * 1.0
hdu0.header['OBSNCHAN'] = chnum

print(hdu0.header['OBSBW'])
print(hdu0.header['OBSNCHAN'])

float_tsubint = np.array(data1['TSUBINT'])[startn:endn+1]
float_offs_sub = np.array(data1['OFFS_SUB'])[startn:endn+1]
float_lst_sub = np.array(data1['LST_SUB'])[startn:endn+1]
float_ra_sub = np.array(data1['RA_SUB'])[startn:endn+1]
float_dec_sub = np.array(data1['DEC_SUB'])[startn:endn+1]
float_glon_sub = np.array(data1['GLON_SUB'])[startn:endn+1]
float_glat_sub = np.array(data1['GLAT_SUB'])[startn:endn+1]
float_fd_ang = np.array(data1['FD_ANG'])[startn:endn+1]
float_pos_ang = np.array(data1['POS_ANG'])[startn:endn+1]
float_par_ang = np.array(data1['PAR_ANG'])[startn:endn+1]
float_tel_az = np.array(data1['TEL_AZ'])[startn:endn+1]
float_tel_zen = np.array(data1['TEL_ZEN'])[startn:endn+1]

float_data = np.array(data1['DATA'])
temp_float_dat_scl = np.array(data1['DAT_SCL'])
print(np.size(float_data))
print(np.size(temp_float_dat_scl) / npol / nchan)
nline = header1['NAXIS2']

print(np.size(float_data) / nline)
print('nsblk ', nsblk, np.size(float_data) / nline / npol / nchan, np.size(float_data) / nline / nchan / npol)

float_dat_freq = np.array(data1['DAT_FREQ'])[startn:endn+1, startfreq:endfreq+1]
float_dat_wts = np.array(data1['DAT_WTS'])[startn:endn+1, startfreq:endfreq+1]
float_dat_offs = np.array(data1['DAT_OFFS'])[startn:endn+1, startfreq:endfreq+1]
float_dat_scl = np.array(data1['DAT_SCL'])[startn:endn+1, startfreq:endfreq+1]

print(np.size(float_dat_freq), np.size(np.array(data1['DAT_FREQ'])))

float_data2 = np.zeros([linenum, nsblk * chnum])
float_data3 = np.zeros([linenum, nsblk * chnum])
float_data_tot = np.zeros([linenum, nsblk * chnum])

dataformat = str(nsblk * chnum) + 'B'
print(dataformat, np.size(float_data2), linenum, nline)

for i in range(linenum):
    temp_data = float_data[i + startn, :].reshape([np.size(float_data[i + startn, :]) // nchan // npol, npol * nchan])
    print(nchan, np.size(temp_data), np.size(float_data[i + startn, :]) // nchan // npol)
    temp_data2 = temp_data[:, startfreq:endfreq+1].reshape(np.size(float_data[i + startn, :]) // nchan // npol * chnum)
    temp_data3 = temp_data[:, nchan + startfreq:nchan + endfreq+1].reshape(np.size(float_data[i + startn, :]) // nchan // npol * chnum)
    temp_data_tot = (temp_data2 + temp_data3) / 2
    float_data2[i, :] = temp_data2
    float_data3[i, :] = temp_data3
    float_data_tot[i, :] = temp_data_tot

dataformat2 = str(chnum) + 'E'
print(dataformat, dataformat2)

column2_data = fits.Column(name='TSUBINT', format='1D', array=float_tsubint, unit='s')
column3_data = fits.Column(name='OFFS_SUB', format='1D', array=float_offs_sub, unit='s')
column4_data = fits.Column(name='LST_SUB', format='1D', array=float_lst_sub, unit='s')
column5_data = fits.Column(name='RA_SUB', format='1D', array=float_ra_sub, unit='deg')
column6_data = fits.Column(name='DEC_SUB', format='1D', array=float_dec_sub, unit='deg')
column7_data = fits.Column(name='GLON_SUB', format='1D', array=float_glon_sub, unit='deg')
column8_data = fits.Column(name='GLAT_SUB', format='1D', array=float_glat_sub, unit='deg')
column9_data = fits.Column(name='FD_ANG', format='1E', array=float_fd_ang, unit='deg')
column10_data = fits.Column(name='POS_ANG', format='1E', array=float_pos_ang, unit='deg')
column11_data = fits.Column(name='PAR_ANG', format='1E', array=float_par_ang, unit='deg')
column12_data = fits.Column(name='TEL_AZ', format='1E', array=float_tel_az, unit='deg')
column13_data = fits.Column(name='TEL_ZEN', format='1E', array=float_tel_zen, unit='deg')
column16_data = fits.Column(name='DAT_FREQ', format=dataformat2, array=float_dat_freq, unit='MHz')
column17_data = fits.Column(name='DAT_WTS', format=dataformat2, array=float_dat_wts, unit='MHz')
column18_data = fits.Column(name='DAT_OFFS', format=dataformat2, array=float_dat_offs, unit='MHz')
column19_data = fits.Column(name='DAT_SCL', format=dataformat2, array=float_dat_scl, unit='MHz')
column20_data = fits.Column(name='DATA', format=dataformat, array=float_data2, unit='Jy')
column20_data_2 = fits.Column(name='DATA', format=dataformat, array=float_data3, unit='Jy')
column20_data_tot = fits.Column(name='DATA', format=dataformat, array=float_data_tot, unit='Jy')

table_hdu = fits.BinTableHDU.from_columns([column2_data, column3_data, column4_data, column5_data, column6_data, column7_data, column8_data, column9_data, column10_data, column11_data, column12_data, column13_data, column16_data, column17_data, column18_data, column19_data, column20_data])

table_hdu.header.append(('INT_TYPE', 'TIME', 'Time axis (TIME, BINPHSPERI, BINLNGASC, etc)'))
table_hdu.header.append(('INT_UNIT', 'SEC', 'Unit of time axis (SEC, PHS (0-1),DEG)'))
table_hdu.header.append(('SCALE', 'FluxDec', 'Intensiy units (FluxDec/RefFlux/Jansky)'))
table_hdu.header.append(('NPOL', 1, 'Nr of polarisations'))
table_hdu.header.append(('POL_TYPE', 'AABB', 'Polarisation identifier (e.g., AABBCRCI, AA+BB)'))
table_hdu.header.append(('TBIN', tbin, '[s] Time per bin or sample'))
table_hdu.header.append(('NBIN', 1, 'Nr of bins (PSR/CAL mode; else 1)'))
table_hdu.header.append(('NBIN_PRD', 0, 'Nr of bins/pulse period (for gated data)'))
table_hdu.header.append(('PHS_OFFS', 0.0, 'Phase offset of bin 0 for gated data'))
table_hdu.header.append(('NBITS', 8, 'Nr of bits/datum (SEARCH mode "X" data, else 1)'))
table_hdu.header.append(('NSUBOFFS', 0, 'Subint offset (Contiguous SEARCH-mode files)'))
table_hdu.header.append(('NCHAN', chnum, 'Number of channels/sub-bands in this file'))
table_hdu.header.append(('CHAN_BW', chan_bw, '[MHz] Channel/sub-band width'))
table_hdu.header.append(('NCHNOFFS', 0, 'Channel/sub-band offset for split files'))
table_hdu.header.append(('NSBLK', nsblk, 'Samples/row (SEARCH mode, else 1)'))
table_hdu.header.append(('EXTNAME', 'SUBINT  ', 'name of this binary table extension'))

table_hdu2 = fits.BinTableHDU.from_columns([column2_data, column3_data, column4_data, column5_data, column6_data, column7_data, column8_data, column9_data, column10_data, column11_data, column12_data, column13_data, column16_data, column17_data, column18_data, column19_data, column20_data_2])

table_hdu2.header.append(('INT_TYPE', 'TIME', 'Time axis (TIME, BINPHSPERI, BINLNGASC, etc)'))
table_hdu2.header.append(('INT_UNIT', 'SEC', 'Unit of time axis (SEC, PHS (0-1),DEG)'))
table_hdu2.header.append(('SCALE', 'FluxDec', 'Intensiy units (FluxDec/RefFlux/Jansky)'))
table_hdu2.header.append(('NPOL', 1, 'Nr of polarisations'))
table_hdu2.header.append(('POL_TYPE', 'AABB', 'Polarisation identifier (e.g., AABBCRCI, AA+BB)'))
table_hdu2.header.append(('TBIN', tbin, '[s] Time per bin or sample'))
table_hdu2.header.append(('NBIN', 1, 'Nr of bins (PSR/CAL mode; else 1)'))
table_hdu2.header.append(('NBIN_PRD', 0, 'Nr of bins/pulse period (for gated data)'))
table_hdu2.header.append(('PHS_OFFS', 0.0, 'Phase offset of bin 0 for gated data'))
table_hdu2.header.append(('NBITS', 8, 'Nr of bits/datum (SEARCH mode "X" data, else 1)'))
table_hdu2.header.append(('NSUBOFFS', 0, 'Subint offset (Contiguous SEARCH-mode files)'))
table_hdu2.header.append(('NCHAN', chnum, 'Number of channels/sub-bands in this file'))
table_hdu2.header.append(('CHAN_BW', chan_bw, '[MHz] Channel/sub-band width'))
table_hdu2.header.append(('NCHNOFFS', 0, 'Channel/sub-band offset for split files'))
table_hdu2.header.append(('NSBLK', nsblk, 'Samples/row (SEARCH mode, else 1)'))
table_hdu2.header.append(('EXTNAME', 'SUBINT  ', 'name of this binary table extension'))

table_hdu3 = fits.BinTableHDU.from_columns([column2_data, column3_data, column4_data, column5_data, column6_data, column7_data, column8_data, column9_data, column10_data, column11_data, column12_data, column13_data, column16_data, column17_data, column18_data, column19_data, column20_data_tot])

table_hdu3.header.append(('INT_TYPE', 'TIME', 'Time axis (TIME, BINPHSPERI, BINLNGASC, etc)'))
table_hdu3.header.append(('INT_UNIT', 'SEC', 'Unit of time axis (SEC, PHS (0-1),DEG)'))
table_hdu3.header.append(('SCALE', 'FluxDec', 'Intensiy units (FluxDec/RefFlux/Jansky)'))
table_hdu3.header.append(('NPOL', 1, 'Nr of polarisations'))
table_hdu3.header.append(('POL_TYPE', 'AABB', 'Polarisation identifier (e.g., AABBCRCI, AA+BB)'))
table_hdu3.header.append(('TBIN', tbin, '[s] Time per bin or sample'))
table_hdu3.header.append(('NBIN', 1, 'Nr of bins (PSR/CAL mode; else 1)'))
table_hdu3.header.append(('NBIN_PRD', 0, 'Nr of bins/pulse period (for gated data)'))
table_hdu3.header.append(('PHS_OFFS', 0.0, 'Phase offset of bin 0 for gated data'))
table_hdu3.header.append(('NBITS', 8, 'Nr of bits/datum (SEARCH mode "X" data, else 1)'))
table_hdu3.header.append(('NSUBOFFS', 0, 'Subint offset (Contiguous SEARCH-mode files)'))
table_hdu3.header.append(('NCHAN', chnum, 'Number of channels/sub-bands in this file'))
table_hdu3.header.append(('CHAN_BW', chan_bw, '[MHz] Channel/sub-band width'))
table_hdu3.header.append(('NCHNOFFS', 0, 'Channel/sub-band offset for split files'))
table_hdu3.header.append(('NSBLK', nsblk, 'Samples/row (SEARCH mode, else 1)'))
table_hdu3.header.append(('EXTNAME', 'SUBINT  ', 'name of this binary table extension'))

hdulist4 = fits.HDUList([hdu0, table_hdu3])
outname3 = fileroot + '_tot_' + sys.argv[1] + '_' + sys.argv[2] + '_' + sys.argv[3] + '_' + sys.argv[4] + '.fits'
rmcomm3 = 'rm -f ' + outname3
os.system(rmcomm3)
hdulist4.writeto(outname3)

print('--------------------------------------------')
print('             Finished!                      ')

endtime = datetime.datetime.now()
print('START:', starttime)
print('END:', endtime)
duration = endtime - starttime
print('DURATION:', duration.seconds, ' sec')
