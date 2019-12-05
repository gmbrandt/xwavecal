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