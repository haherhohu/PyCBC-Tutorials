#!/usr/bin/env python
# coding: utf-8

# ## See the signal in the data ##
# 
# Because GW150914 was relatively loud, it is readily visible in the data. 

# ### See GW150914 in the time series ###

from pycbc.filter import highpass_fir, lowpass_fir
from pycbc.psd import welch, interpolate
from matplotlib import pyplot

from pycbc.catalog import Merger

for ifo in ['H1', 'L1']:
    # Read data: The data is available to the public through losc.ligo.org!
    ts = Merger("GW150914").strain(ifo)
    
    # Estimate the noise spectrum and use it to whiten the data
    psd = interpolate(welch(ts), 1.0 / ts.duration)
    white_strain = (ts.to_frequencyseries() / psd ** 0.5).to_timeseries()

    # remove frequencies below and above where the main portion of the signal lies
    smooth = highpass_fir(white_strain, 35, 8)
    smooth = lowpass_fir(smooth, 300, 8)

    # time shift and flip L1 to match the time of arrival and phase of Hanford
    if ifo == 'L1':
        smooth *= -1
        smooth.roll(int(.007 / smooth.delta_t))

    pyplot.plot(smooth.sample_times.numpy(), smooth.numpy(), label=ifo)
    pyplot.savefig(f"gw150914_stain_{ifo}.png")
    pyplot.clf()

# Plot the region around the signal (time is in gps seconds)
pyplot.legend()
pyplot.xlim(1126259462.26, 1126259462.48)
pyplot.ylabel('Smoothed-Whitened Strain')
pyplot.grid()
pyplot.ylim(-100, 100)
pyplot.xlabel('GPS Time (s)')
pyplot.savefig("gw150914_smoothed_whitened_stain.png")
pyplot.clf()


# ### See the time frequency evolution of GW150914 in Hanford ###
ts = Merger("GW150914").strain("H1")

ts = ts.whiten(4, 4)
zoom = ts.time_slice(1126259462.4 - 1, 1126259462.4 + 1)
times, freqs, power = zoom.qtransform(.001, 1, frange=(20, 512), qrange=(4, 64))

pyplot.figure(figsize=(18, 3))
pyplot.pcolormesh(times, freqs, power)
pyplot.ylim(20, 512)
pyplot.xlabel('Time (s)')
pyplot.ylabel('Frequency (Hz)')
pyplot.xlim(times.min(), times.max())
pyplot.yscale('log')
pyplot.savefig("gw150914_time_frequency_evolution_h1.png")
pyplot.clf()


# ### See the time frequency evolution of GW150914 in Livingston ###
from pycbc.frame import read_frame
ts = Merger("GW150914").strain("L1")

ts = ts.whiten(4, 4)
zoom = ts.time_slice(1126259462.4 - 1, 1126259462.4 + 1)
times, freqs, power = zoom.qtransform(.001, 1, frange=(20, 512), qrange=(4, 64))

pyplot.figure(figsize=(18, 3))
pyplot.pcolormesh(times, freqs, power)
pyplot.ylim(20, 512)
pyplot.xlabel('Time (s)')
pyplot.ylabel('Frequency (Hz)')
pyplot.xlim(times.min(), times.max())
pyplot.yscale('log')

pyplot.savefig("gw150914_time_frequency_evolution_l1.png")
pyplot.clf()


