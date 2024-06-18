#!/usr/bin/env python
# coding: utf-8

# ### Comparing Gravitational waveforms to each other###

#get_ipython().run_line_magic('matplotlib', 'inline')
# We learn about the potential parameters of a source by comparing it to many different waveforms
# each of which represents a possible source with different properties. 
from matplotlib import pyplot
from pycbc.waveform import get_td_waveform

# We can directly compare how similar waveforms are to each other using an inner product between then called 
# a 'match'. This maximizes over the possible time of arrival and phase. We'll generate a reference waveform
# which we'll compare to.
m1 = m2 = 20
f_lower = 20
approximant = "SEOBNRv4"
delta_t = 1.0 / 2048
hp, _ = get_td_waveform(approximant=approximant,
                         mass1=m1, mass2=m2,
                         delta_t=delta_t, f_lower=f_lower)
pyplot.plot(hp.sample_times, hp)
pyplot.xlabel('Time (s)')
pyplot.ylabel('Strain')
pyplot.savefig('waveform_sampled.png')
pyplot.clf()

# How similar waveforms are to each other depends on how important we consider different frequencies, we 
# can account for this by weighting with an estimated power spectral density. We'll use here 
# the predicted final Advanced LIGO final design sensitivity
from pycbc.psd import aLIGOZeroDetHighPower
psd = aLIGOZeroDetHighPower(len(hp) // 2 + 1, 1.0 / hp.duration, f_lower)

pyplot.loglog(psd.sample_frequencies, psd)
pyplot.xlabel('Frequency (Hz)')
pyplot.ylabel('Strain**2 / Hz')
pyplot.xlim(20, 1000)

pyplot.savefig('waveform_sample_frequency.png')
pyplot.clf()


# We can now compare how similar our waveform is to others with different masses
from pycbc.filter import match
import numpy

masses = numpy.arange(19, 21, .2)
matches = []
for m2 in masses:
    hp2, _ = get_td_waveform(approximant=approximant,
                         mass1=m1, mass2=m2,
                         delta_t=delta_t, f_lower=f_lower)
    hp2 = hp2[:len(hp)] if len(hp) < len(hp2) else hp2
    hp2.resize(len(hp))
    
    m, idx = match(hp, hp2, psd=psd, low_frequency_cutoff=f_lower)
    matches.append(m)
    pyplot.plot(hp2.sample_times, hp2)
pyplot.xlim(-.05, .02)
pyplot.savefig('waveform_compare_other_masses.png')
pyplot.clf()


pyplot.plot(masses, matches)
pyplot.ylabel('Match')
pyplot.xlabel('Mass of second object (Solar Masses)')
pyplot.savefig('waveform_match_masses.png')
pyplot.clf()

# You can think of the match also as the fraction of signal-to-noise that you could recover with a template that 
# doesn't *exactly* look like your source


