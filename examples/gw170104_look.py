#!/usr/bin/env python
# coding: utf-8

# Get the data for H1 from the LOSC site
import os
os.system('wget https://losc.ligo.org/s/events/GW170104/H-H1_LOSC_4_V1-1167559920-32.gwf')
os.system('wget https://losc.ligo.org/s/events/GW170104/L-L1_LOSC_4_V1-1167559920-32.gwf')


# ### See the track of GW170104 in both the Hanford and Livingston data ###

from pycbc.catalog import Merger
from matplotlib import pyplot

for ifo in ['H1', 'L1']:
    pyplot.figure()
    ts = Merger("GW170104").strain(ifo)

    ts = ts.whiten(4, 4)
    zoom = ts.time_slice(1167559936.6 - .75, 1167559936.6 + .75)
    times, freqs, power = zoom.qtransform(.001, 1, frange=(20, 512), qrange=(4, 64))

    pyplot.figure(figsize=(18, 3))
    pyplot.pcolormesh(times, freqs, power)
    
    pyplot.ylim(20, 512)
    pyplot.xlabel('Time (s)')
    pyplot.ylabel('Frequency (Hz)')
    pyplot.xlim(times.min(), times.max())
    pyplot.yscale('log')

    pyplot.savefig(f"gw170104_{ifo}_look.png")
    pyplot.clf()
