0.1.11 (2020-05-21)
-------------------
- Limited the number of lines to try to pair during overlap fitting (to 20 from blue and 15 from the red side). This addresses a bug whereby overlap fitting would stall indefinitely if there were too many pairs of lines to try.

0.1.10 (2020-04-23)
-------------------
- Removed peakutils requirement. Updated astropy and scipy dependencies.

0.1.9 (2020-04-23)
-------------------
- Removed requirements.txt. Placed requirements directly into install_requires in setup.py

0.1.8 (2020-04-21)
-------------------
- Hotfix for a bug in ``find_feature_wavelengths`` where the function would work only 
if the features were marked as belonging to fiber 1. Changed ``find_feature_wavelengths``
so that it assumes all the features come from a single fiber.

0.1.7 (2020-04-20)
-------------------
- Added ``min_num_matches``, ``global_scale_spacing``, ``initial_mad_clip``, ``final_mad_clip``
as possible options for the config files.

0.1.6 (2020-04-20)
-------------------
- Added a convienience function which wavelength calibrates from a list
of spectral feature pixel and order positions, and a provided reference line list.

0.1.5 (2020-03-22)
-------------------
- Feature identification stage now fetches the pixel error in the centroid, and the flux error
of the feature.
- lines can be weighted. This is not easily implemented yet, but implementing this in the future is easy.
We now need to only set the 'weight' key of measured_lines in any stage prior to refine_wcs().

0.1.4 (2020-03-22)
-------------------
- Fixed bug where two identical diffraction orders had different reference ID assigned by fibers.IdentifyFibers. This 
bug caused the e2e tests to yield a poor wavelength solution on one fiber of the test data, because
the principle order number was wrong for that fiber (because the reference id for each order was wrong). 
- changed required numpy version to >=1.16 .

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

