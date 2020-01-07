0.1.3 (2020-01-07)
-------------------
- Fixed issue 10 whereby `FitOverlaps` would return illegitimate overlap fits
that violated one-to-one correspondance between a red side emission feature and
its 'duplicate' on the blue side.
- The BackgroundSubtractSpectrum nor the BackgroundSubtract stages no longer effect the statistical uncertainty
of the image.
- RuntimeContext raises a custom attribute error if an attribute is missing, this
error suggests to the user that the attribute name is likely missing from the configuration
file.

0.1.2 (2019-12-06)
-------------------
- Fixed xwavecal.fibers.IdentifyFibers to work for instruments with any number
of wavelength calibration fibers or lit fibers. This should resolve issues #6 and #7.
- Not giving a --frame-type to the xwavecal_reduce command will default to 'any',
which will reduce any frame type, even frame types with custom names that are defined
in the config.ini file and have associated custom sets of reduction stages.

0.1.1 (2019-10-23)
-------------------
- Minor refactor. Removed IRD config file from the examples.
- Every logger message now displays the stage name (if applicable)
- The wavelength solution aborts if there is no main spectrum on the image.
- FileNotFoundError if an incorrect config file path is provided. 

0.1.0 (2019-10-20)
-------------------
- Initial release.