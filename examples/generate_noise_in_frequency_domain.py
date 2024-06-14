#!/usr/bin/env python
# coding: utf-8

# ### Generate Noise Directly in the Frequency Domain ###
# 
# It is sometimes useful to generate noise directly in the frequency domain. PyCBC provides functionality to generate noise that has the same spectrum as a PSD you provide. Note that due to how this is made, the equivelant time series will have circular boundary conditions.

import pycbc.noise
import pycbc.psd
from matplotlib import pyplot

# Generate a PSD using an analytic expression for 
# the full design Advanced LIGO noise curve
f_lower = 10
duration = 128
sample_rate = 4096
tsamples = sample_rate * duration
fsamples = tsamples // 2 + 1
df = 1.0 / duration

psd = pycbc.psd.from_string('aLIGOZeroDetHighPower', fsamples, df, f_lower)

# Let's take a look at the spectrum
pyplot.loglog(psd.sample_frequencies, psd)
pyplot.xlim(20, 1024)
pyplot.ylim(1e-48, 1e-45)
pyplot.xlabel('Frequency (Hz)')
pyplot.ylabel('Strain^2 / Hz')
pyplot.grid()
pyplot.savefig('generate_sample_frequencies.png')
pyplot.clf()

# Now, let's generate noise that has the same spectrum
htilde = pycbc.noise.frequency_noise_from_psd(psd, seed=857)

# Equivelantly in the time domain
hoft = htilde.to_timeseries()

# Well zoom in around a short time
hoft_zoom = hoft.time_slice(43.5, 44)
pyplot.plot(hoft_zoom.sample_times, hoft_zoom)
pyplot.savefig('generate_noise.png')
pyplot.clf()


# Notice the extreme difference in amplitudes of noise at different frequencies.
# You can see clearly the loud low frequency noise.

# Now we can verify that the PSD of this noise is as expected. 

# 여기에서 이거 출력되는데 뭔지는 모르겠음.
# lal.MSUN_SI != Msun
# __name__ = gwsurrogate.new.spline_evaluation
# __package__= gwsurrogate.new
# Use Welch's method with 4s segments
psd_estimated = hoft.psd(4)

pyplot.loglog(psd_estimated.sample_frequencies, psd_estimated, label='Estimated')
pyplot.loglog(psd.sample_frequencies, psd, label='Original', linewidth=3)
pyplot.legend()
pyplot.xlim(20, 1024)
pyplot.ylim(1e-48, 1e-45)
pyplot.savefig('generate_sample_frequencies_psd_estimated.png')
