# pycbc-tutorial for kagra

## New Features
 - Test code verified in the KAGRA-igwn environment
 - Added *.py files that can be executed directly in the KAGRA terminal
 - Provided examples of using KAGRA data
 - Provided examples of integrating with gwdatafind

## Major Changes
 - pylab -> pyplot: pylab is no longer used in matplotlib (deprecated). Therefore, it is strongly recommended not to use it, so everything has been changed to pyplot
 - Directory structure changes due to added examples and updated examples Differentiating examples, tutorial, inference, and separating by extensions such as ipynb, py
 - Organized naming rules, and accordingly, cleaned up names, code naming/output file names, etc. (planned) <- gwosc probably won't like this much