#!/usr/bin/env python
# coding: utf-8

# ### Listen to the merger converted to audio ###

from gwdatafind import (
    find_urls,
    Session,
)

# search data with gwdatafind
with Session() as sess:
    # connection info
    urls = find_urls(
        "K",                        # required (observatory)
        "K1_R",                     # required (frametype)
        1402299104,                 # required (start)
        1402299136,                 # required (end)
        host="http://127.0.0.1:80", # required (server-url)
        session=sess,               # required (session)
    )
    target = urls[0].replace("file://localhost", "")

import os
#Set up the libraries we need and download the data
os.system(f'cp {target} .') 

from pycbc.frame import read_frame
from pycbc.filter import highpass_fir, lowpass_fir
from pycbc.psd import welch, interpolate
from pycbc.types import TimeSeries

# Read data and remove low frequency content
fname = os.path.basename(target)
h1 = highpass_fir(read_frame(fname, 'K1:CAL-CS_PROC_DARM_STRAIN_DBL_DQ'), 15.0, 8)

# estimate the noise spectrum and whiten
psd = interpolate(welch(h1), 1.0 / 32)
white_strain = (h1.to_frequencyseries() / psd ** 0.5 * psd.delta_f).to_timeseries()

# remove some of the high and low frequencies
smooth = highpass_fir(white_strain, 25, 8)
smooth = lowpass_fir(white_strain, 250, 8)

# slow the data down by a factor of 4
smooth = TimeSeries(smooth, delta_t=smooth.delta_t*4)

#strech out and shift the frequency upwards by 300 Hz to aid human hearing
frequency_shift = 300.0
fdata = smooth.to_frequencyseries()
fdata.roll(int(frequency_shift / fdata.delta_f))
smooth = fdata.to_timeseries()

#Take slice around signal
smooth = smooth[len(smooth)//2 - 1500:len(smooth)//2 + 3000]
smooth.save_to_wav('gw140229_k1_chirp.wav')

# ### Downlod your audio file ###
#  1. If using JupyterLab, your file should show up in the sidebar
#  2. Select and download gw150914_h1_chirp.wav
