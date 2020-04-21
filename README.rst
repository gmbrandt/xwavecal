xwavecal
========
.. image:: https://coveralls.io/repos/github/gmbrandt/xwavecal/badge.svg?branch=master
    :target: https://coveralls.io/github/gmbrandt/xwavecal?branch=master

.. image:: https://travis-ci.org/gmbrandt/xwavecal.svg?branch=master
    :target: https://travis-ci.org/gmbrandt/xwavecal

.. image:: https://zenodo.org/badge/207922849.svg
   :target: https://zenodo.org/badge/latestdoi/207922849

A library of routines for automatically wavelength calibrating echelle
spectrographs for high precision radial velocity work. ``xwavecal`` is designed to operate on data with
extracted 1D spectra.

:Authors: Mirek Brandt, Curtis McCully, and Timothy D. Brandt.

If you use this code, please cite **Brandt, G.M. et al. (2019)** which can be found
on ArXiv here: https://arxiv.org/abs/1910.08079

As well, please cite the Zenodo DOI above.

At best, using ``xwavecal`` only requires editing a config.ini file for your data.
I cover how to do that in this readme.

Installation
============
``xwavecal`` is installed via pip by running

.. code-block:: bash

    pip install .

While in the root directory of this repository. It can also be installed by running

.. code-block:: bash

    pip install git+https://github.com/gmbrandt/xwavecal

Wavelength Calibrating Spectrum
===============================
This section covers how to wavelength calibrate data which already have a spectrum and a blaze
corrected version of the same spectrum. Using ``xwavecal`` with spectra is preferred.

If you have raw data only and extracting a spectrum is difficult, you may try the experimental data
reduction pipeline included with ``xwavecal``, see section "Configuring for full data reduction".
However, I highly recommend extracting the spectrum first and running ``xwavecal`` in the preferred way.

Information about the input data products are to
be set by users via a :code:`config.ini` file. See the file ``xwavecal/example_config/nres_config.ini``
for the configuration to reduce NRES data (starting from extracted spectrum).
I now cover how to wavelength calibrate data, using the Network of Robotic Echelle Spectrographs (NRES) from Las Cumbres Observatory
as an example.

To reduce the directory of NRES test data included
in this repo, you would run from the command line (after modifying the two paths ``line_list_path`` 
and ``database_path`` in the nres_config.ini file):

.. code-block:: bash

    xwavecal_reduce_dir --input-dir xwavecal/tests/data/nres_test_data/
     --output-dir ~/Downloads --config-file xwavecal/example_config/nres_config.ini

To do the same reduction by specifying paths, you would run:

.. code-block:: bash

    xwavecal_reduce --data-paths
     xwavecal/tests/data/nres_test_data/cptnrs03-fa13-20190405-0004-w00.fits.fz
      xwavecal/tests/data/nres_test_data/cptnrs03-fa13-20190405-0014-a00.fits.fz
       --output-dir ~/Downloads --config-file xwavecal/example_config/nres_config.ini

For NRES, ``a00.fits.fz`` files are ThAr wavelength calibrations, ``w00.fits.fz`` are lampflats.
Configuring this wavelength solution to work for your instrument should only involve
making a new config.ini file. The rest of this readme is devoted to setting the config
file for a new instrument where the input data are extracted 1D spectra. I use
NRES as an example.


Wavelength Calibrating from a list of spectral features
=======================================================
If you do not want to use xwavecal's line identification and spectral reduction utilities (which I will dive into shortly),
there is a convienience function which will return a list of wavelengths from just a list of
spectral feature coordinates (pixel and order) and a reference line list. The returned wavelengths are the
 wavelengths of the measured spectral features under the best fit wavelength model. This function is
``xwavecal.wavelength.wavelength_calibrate()``. It is a wrapper for all the xwavecal stages that
occur after spectral features have been centroided. This would be useful if your calibration lamp is e.g.
not a lamp, but an absorption cell, or if you only want to use xwavecal
as a fallback calibration in the existing pipeline that you use. See the docstring for
``xwavecal.wavelength.wavelength_calibrate()`` for more details.


Configuring for wavelength calibration
======================================
``xwavecal`` is designed in a modular fashion. Each step of the wavelength
calibration is a stage which can be disabled by removing the associated line
in the config.ini file. Wavelength calibrating data which already have spectra
means only using the wavelength calibration stages. Using the full experimental pipeline
means enabling the other data reduction stages (e.g. overscan subtraction etc.).

The completed config.ini file is "nres_config_wcs_only.ini", this contains
all the options and settings to reduce NRES data which already has a 1D spectrum
and a 1D blaze corrected spectrum. This repo includes raw NRES data, which has to
be reduced with nres_config.ini (which includes all the overscan subtraction, spectral extraction etc. stages).


Pointing to the database and line list
--------------------------------------

We start by telling the config.ini where the database for the reduced data should live.

Before reducing, copy the nres_config_wcs_only.ini file to a new location, rename it for your instrument, and
change :code:`database_path` under the [reduction] section to the path where you
want to the database to exist. The parent folder for the database must already exist. E.g. for myself,
this is :code:`"/home/gmbrandt/Downloads/pipeline.db"` . The surrounding :code:`" "` quotes must be there for
the config file to process properly.

The database will keep track of all your processed files. All processed calibration files are saved under the
table :code:`caldata` in the .db file specified.

In [reduction] change the ``line_list_path`` as well. We include the original ThAr (Thorium-Argon) atlas
from the European Southern Observatory (ESO). This was retrieved
from http://www.eso.org/sci/facilities/paranal/instruments/uves/tools/tharatlas.html in late
2019. This line list was designed for spectrographs with a resolving power (R) of 100,000, and thus
it may not be suited for your instrument if it has a lower or larger R. Moreover, the wavelengths are air wavelengths.
It is up to you to download a line list suitable for your instrument (if the ThAr atlas is not suitable)
and correct the line list for the index of refraction of air if necessary.


Data settings
-------------
Here we tell ``xwavecal`` via the config file where various information lies in the header of
your input data.

In section [data] we will need to edit:

- ``primary_data_extension``
- ``files_contain``
- ``header_keys``
- ``type_keys``

data_class is also editable, but most likely will not need to be changed. data_class is the
Python object used to load in your data. The default ``xwavecal.images.Image`` should be fine for your data.

I describe the four items above with examples of setting them. See the full config file
``xwavecal/example_config/nres_config_wcs_only.ini`` for an example of setting all the above.

- ``primary_data_extension`` is the number label of the
  fits extension (e.g. ``0``)  where all the relevant header data is stored such as
  the observation date, instrument name etc. These are used for writing out the file with an informative name.
- ``files_contain`` is a list of strings, where each string must be present in the input file types. The default is
  ['.fits'] in which case only files with '.fits' in the name are reduced. For example:

 * If I had two files in the directory I was about to reduce: 'IRDA003.fits' and 'IRDB002.fits',
   and I wanted to only process 'IRDA003.fits',
   I would set ``files_contain = ['.fits', 'IRDA']``


header_keys
~~~~~~~~~~~

``header_keys`` is a python dictionary. The *values* of the dictionary are the header keywords
in your raw data that give things like the read noise, the observation date, etc. The *keys*
are the standard keys understood by ``xwavecal``. Some of these keys are:

- ``type`` : the frame type e.g. lampflat
- ``gain`` : the gain in e-/ADU
- ``read_noise`` : the read noise in e-
- ``fiber_state`` : the string which gives which fibers are lit and with what. See fiber_state in its subsection.
- ``observation_date`` : observation date, see time_format in its subsection.
- ``instrument``: see below
- ``instrument2``
- ``site_name``
- ``unique_id``: Some running identifier for the input frames. If none, choose a stagnant one -- just
  beware of accidental overwrites if you do not choose a unique identifer for your data.

``instrument``, ``instrument2``, ``site_name`` are used to index the processed data in the
sqlite database. E.g. for NRES, I set:

.. code-block:: python

               ...
               'instrument': 'TELESCOP',
               'instrument2': 'INSTRUME',
               'site_name': 'SITEID',
               ...

This means that processed data will be stored in the database with telescope name, instrument name, and the
ID of our site. These data are stored in NRES frames under the header keys 'TELESCOP', 'INSTRUME', 'SITEID'.

``observation_date`` is the .fits header key which gives the observation date of the frame.
One must set time_format (see further down in this section) to agree with the format of the .fits value given
by the ``observation_date`` key.

For ``fiber_state``, the NRES and HARPS store this in a single string in 'OBJECTS' and 'ESO DPR TYPE', respectively.
For NRES the value of the header looks like ``thar&thar&none`` for a frame with Thorium-Argon (ThAr) lit on fibers 0,1 and
fiber 2 unlit. For HARPS, the same configuration (but no third fiber since it does not exist) would be
``WAVE,WAVE,THAR2``. We will convert ``WAVE,WAVE,THAR2`` to ``thar&thar&none`` with the type_keys next.

type_keys
~~~~~~~~~

``type_keys`` is by far the most confusing part of configuring an instrument. This may get easier in a future release.
``type_keys`` is a dictionary which takes the value of any .fits header value and converts it in place. Consider if the
``fiber_state`` key in my .fits header was ``ESO DPR TYPE``, and that portion of the header looked like:
``{'ESO DPR TYPE': 'WAVE,WAVE,THAR2'}``. I could set
``type_keys = {..., 'WAVE,WAVE,THAR2': 'thar&thar&none'}``, then any time ``xwavecal`` reads the ``fiber_state`` item
it will read 'thar&thar&none'.

fiber_state
~~~~~~~~~~~
A note on ``fiber_state``: One must convert whatever ``fiber_state`` value in your .fits file to be
of the string format interpretable by ``wavecal``. This format is always ``fiber0lamp&fiber1lamp&fiber2lamp``.
Where ``fiberxlamp`` is the type of light coming through that lamp. If your instrument
only has two fibers, leave the last entry as 'none'.

If I had a fictional instrument with two
lamps, quartz and thorium argon and only two fibers, then in type_keys I would have to add all expected permutations thereof:

.. code-block:: python

    type_keys = {...,
                'quartzANDquartz': 'other&other&none',
                'tharANDthar': 'thar&thar&none',
                'unlitANDthar': 'none&thar&none',
                 ...}

and so forth. It does not matter what you call lampflat or other lamps that are not calibration lamps. All
wavelength calibration lamp states must be called ``thar`` (regardless of whether the lamp is ThAr, or NeAr, or some other
gaseous mixture, although be sure to point ``xwavecal`` to an appropriate line list).

Important note
~~~~~~~~~~~~~~

Setting header_keys and type_keys
builds a translator which understands how to interpret your fits header, ``xwavecal`` does not modify existing header keys.
None of these translations will ever be saved onto the fits header of your output data product. The fits
header of your data will *not* have ``read_noise`` etc appended as extra headers.

time_format
~~~~~~~~~~~

In [reduction], ``time_format`` is the time format of the ``observation_date`` output from
the fits header. This must be a string contained in double quotes ``" "`` and understood by
``datetime.datetime.strftime``. Then replace single ``%`` with ``%%`` (to fix a quirk of using a config file).

Other parameters
~~~~~~~~~~~~~~~~
There are other type_keys and header_keys that need to be set only if you run the full data reduction pipeline. Because
I prefer one to run ``xwavecal`` with extracted spectra, I will cover and document these at a later date.

Wavelength calibration settings
-------------------------------
To wavelength calibrate your data, the following settings in config.ini may need to be changed:

- ``main_spectrum_name``
- ``blaze_corrected_spectrum_name``
- ``ref_id``
- ``template_trace_id``
- ``overlap_min_peak_snr``
- ``max_red_overlap``
- ``max_blue_overlap``
- ``global_scale_range``
- ``min_peak_snr``
- ``approx_detector_range_angstroms``
- ``approx_num_orders``
- ``principle_order_number``
- ``m0_range``
- ``min_num_overlaps``
- ``min_num_matches``
- ``flux_tol``
- ``global_scale_spacing``
- ``initial_mad_clip``
- ``final_mad_clip``


There are several other parameters you will most likely not need to change.
Let us go through the pertinent ones in the list above one-by-one:

- ``main_spectrum_name`` : this is the name of the .fits extension that contains
  the BinTableHDU of the spectrum that ``xwavecal`` will calibrate.
- ``blaze_corrected_spectrum_name`` : this is the name of the .fits extension that contains
  the BinTableHDU of the blaze corrected spectrum that ``xwavecal`` will use to find the overlaps.
  If you do not have a blaze corrected spectrum, set
  this to some string (that is not in the raw data) such as ``'None'``.
- ``template_trace_id`` : this is the trace id (id column in the input spectrum) for the
  diffraction order that you want to save as a template. This template will be used to identify this same
  diffraction order in all subsequent spectra you reduce. It will have a ref_id associated with it
  such that the diffraction order number understood by ``xwavecal`` is ``ref_id + m0`` where
  ``m0`` is the principle order number. I recommend setting the ``template_trace_id`` to some middle order on the detector.
- ``ref_id`` : this is the reference id you wish to assign the template spectrum (the order which has the ``id`` of
  ``template_trace_id``) such that the
  diffraction order number understood by ``xwavecal`` for the template spectrum is ``ref_id + m0`` where
  ``m0`` is the principle order number.
- ``overlap_min_peak_snr`` : the minimum signal to noise for an emission peak to be considered in the overlap algorithm.
  see Brandt et al. 2019 for a discussion of the overlap algorithm. I recommend this be set to something low like 5. In
  general, overlap fitting works better if more peaks are detected. For NRES we use 5 and detect ~4000 peaks.
- ``flux_tol`` : If two emission peaks from neighboring orders have flux f1 and f2, ``flux_tol`` is
  the maximum allowed value of abs(f1 - f2)/(mean(f1, f2)) for two peaks to be considered
  a matched pair in the overlap algorithm. For decent blaze correction, use 0.2.
  For bad, or an absence of, correction, use 0.5.
- ``min_num_overlaps`` : The minimum number of properly fit overlaps required for the wavelength solution to proceed.
  Default is 5, and should not need to be changed.
- ``min_num_matches`` : The minimum number of features that are marked as matched for an overlap to be counted as
  properly fit. Default is 6, and should not need to be changed.
- ``min_peak_snr`` : the minimum signal to noise for an emission peak to be used to constrain the wavelength
  solution after overlap detection. This should be something reasonable like 10 or 20 so
  as to detect between 1000 and 2000 emission lines. Weak lines are often contamination from trace elements
  (which are not in reference line lists and so would throw off our algorithm).
- ``max_red_overlap`` : The maximum allowed pixel coordinate for a red-side peak to be considered for our overlap algorithm.
- ``max_blue_overlap`` : The minimum allowed pixel coordinate for a blue-side peak to be considered for our overlap algorithm.

  * The overlap algorithm will try to match peaks from
    (0, ``max_red_overlap``) to (max_pixel, max_pixel - ``max_blue_overlap``). Where max_pixel is the width of
    your detector in x (i.e. the number of columns; e.g. 4096).

- ``approx_detector_range_angstroms`` : If the spectrograph covers the spectral range 3000A to 9000A, then set
  ``approx_detector_range_angstroms = 5000``. Note this value does not need to be precise.
- ``approx_num_orders`` : approximate number of distinct diffraction orders in the spectrum. E.g. 67 for NRES.
  Note this is not the number of traces (visible light streaks on the echelle detector) but the number of diffraction orders.
  I.e. num_of_traces/num_of_lit_fibers. This does not need to be precise.
- ``global_scale_range`` : See Brandt et al. 2019 for a discussion of the global scale.
  This is the range about the initial guess where ``xwavecal`` will search for the global scale. We
  recommend ``global_scale_range = (0.5, 1.5)``.

  * For example: if the guess generated by ``xwavecal`` is ``K`` and if ``global_scale_range = (0.8, 1.2)``
    then ``xwavecal`` will search for the global scale between ``0.8K`` and ``1.2K``.
- ``global_scale_spacing`` : The spacing in Angstroms used the global scale search. Default is 10 for an R 50,000 spectrograph.
  One should change it to 1 for an R 500,000 spectrograph.
- ``initial_mad_clip`` : Lines will be marked as outliers if their residuals exceed initial_mad_clip. Outliers are
  recomputed and ignored during each solving iteration. ``initial_mad_clip`` sets outlier rejection during initial
  refinement (SolutionRefineInitial)
- ``final_mad_clip`` : sets outlier rejection during final refinement (SolutionRefineFinal)

- ``principle_order_number``: This is an integer and needs to exactly correct. This is the true diffraction order
  number of the diffraction order with ref_id = 0. If you do not know this, insert the m0 identification stage
  (I will cover how to do this in a following section), and set ``m0_range`` to a reasonable range of values.
- ``m0_range`` : the range of possible ``m0`` (principle order number) values. This is only used if you
  are searching for ``m0`` (i.e. if you have included 'xwavecal.wavelength.IdentifyPrincipleOrderNumber' in
  the set of stages for wavecal frames).


The Wavelength Models
---------------------
The wavelength models at each of the three stages are set by the parameters ``initial_wavelength_model``,
``intermediate_wavelength_model``, and ``final_wavelength_model`` in the config.ini file. They are dictionaries,
the format is: {xpower: [ipower, ipower,...],..}. The default basis functions are legendre polynomials: ``P^m(x)`` and ``P^m(i)``, where subscript m denotes the
mth order basis function, so that ``P^m(x)`` is akin to ``x^m``.
The default ``initial_wavelength_model`` wavelength model is:

.. code-block:: python

    initial_wavelength_model = {1: [0, 1, 2],
                                2: [0, 1, 2]}


This wavelength solution model above then is:
:math:`\lambda(x, i) = \frac{1}{m0 + i} \left( P^1(x) * [P^0(i) + P^1(i) + P^2(i)] + P^2(x) *[P^0(i) + P^1(i) + P^2(i)] \right)`


The ``(1/(m_0 + i))`` prefactor is always included. If one instead made the model:

.. code-block:: python

    initial_wavelength_model = {1: [0, 1, 2],
                                2: [0, 1]}

That would set the initial wavelength solution model to


:math:`\lambda(x, i) =  \frac{1}{m0 + i} \left( P^1(x) * [P^0(i) + P^1(i) + P^2(i)] + P^2(x) *[P^0(i) + P^1(i)] \right)`

Formatting the input data
-------------------------
The input data should be a .fits file with three data extensions:

- A primary data extension (e.g. one that contains the raw 2d frame). Its header must contain all the necessary
  information like ``fiber_state`` etc. If this data is in extension 0, then set ``primary_data_extension=0``
- An extracted spectrum (e.g. box or optimally extracted) as a ``astropy.fits.BinTableHDU``.
  Set ``main_spectrum_name`` in the config.ini to the extension name of this spectrum.
- A blaze corrected version of the same above spectrum as a ``astropy.fits.BinTableHDU``.
  Set ``blaze_corrected_spectrum_name`` in the config.ini to the name of this spectrum.

For example, below is an exploration of an NRES frame with the spectra attached.

.. code-block:: python

    from astropy.io import fits
    from astropy.table import table

    im = fits.open('/some/example/image.fits.fz')
    im.info()
    >>> No.    Name      Ver    Type      Cards   Dimensions   Format
    >>> 0  SPECTRUM      1 PrimaryHDU     186   (4096, 4096)   float64
    >>> 1  SPECBOX       1 BinTableHDU     24   135R x 7C   [K, 4096D, 4096D, 4096K, K, K, 4096D]
    >>> 2  BLZCORR       1 BinTableHDU     24   135R x 7C   [K, 4096D, 4096D, 4096K, K, K, 4096D]

I have three extensions here. ``im[0].data`` would gives the 2d frame of raw data. ``im[0].header['OBSTYPE']`` would
give the observation type (remember your data does not have to have the key 'OBSTYPE', you set those in config.ini).
Ignore the fact that the raw 2d data is named ``SPECTRUM`` yet the 1D spectra have names ``SPECBOX`` and ``BLZCORR``.
In ``xwavecal/example_config/nres_config.ini`` or ``xwavecal/example_config/nres_config_wcs_only.ini``,
``blaze_corrected_spectrum_name`` and ``main_spectrum_name`` are set to ``BLZCORR`` and ``SPECBOX``, respectively.


Now let us look at the 1D spectra extension closely (the blaze corrected 1D spectrum im['BLZCORR'] has the same format).

.. code-block:: python

    type(im['SPECBOX'])
    >>> astropy.io.fits.hdu.table.BinTableHDU
    # The type must be a table, so that we can do the following.
    spec = Table(im['SPECBOX'].data)
    spec.info()
    >>> <Table length=135>
    >>>    name     dtype   shape
    >>> ---------- ------- -------
    >>>         id   int64
    >>>     ref_id   int64
    >>>       flux float64 (4096,)
    >>>     stderr float64 (4096,)
    >>>      pixel   int64 (4096,)
    >>>      fiber   int64
    >>> wavelength float64 (4096,)

Every spectrum attached to the image must have this format with these columns. Let N be the number of traces.
For NRES, N~135 for 2 lit fibers (so ~67 orders per fiber). ``id, ref_id`` and  ``fiber`` are
1d columns of length N.
``id`` is an arbitrary identification number for each trace. ``ref_id`` is the absolute identifcation number for that
trace. The ``id`` of a diffraction order may change, however the ``ref_id`` will not because that is found by cross
correlating the spectrum with a template (which ``xwavecal`` will create automatically). ``fiber`` is the fiber id
for each row of the spectrum. If you only have one fiber lit, this column can be all 0's or 1's as long as it is consistent
with your .fits header ``fiber_state``.

If you do not want to use ``xwavecal``'s order identification routine: comment out the ``fibers.IdentifyFibers``
stage in the configuration file. In this case, your input spectrum must have the reference_id (ref_id) column correctly
filled out with the reference id ``i``: each consecutive diffraction order must have a reference_id of 1 higher
than the previous. This is important because the grating equation prefactor in the wavelength solution is ``1/(m0 + i)``

Let the detector be X pixels wide, where the echelle grating has dispersed each order across the width. For NRES, X=4096,
where pixel 0 is bluer than pixel 1. ``flux`` are the counts as a function of ``pixel`` (Both shape (N, X) (rows, columns).
``stderr`` is the 1-sigma error for each point in ``flux``. ``wavelength`` is the wavelength of each pixel in ``pixel``.
Of course, ``wavelength`` can be set to 0's or ``np.nan`` or whatever you like -- ``xwavecal`` will populate ``wavelength``
for you.

The spectrum **have to be ordered** such that ``spec[0]`` is redder than ``spec[1]`` (on average) and such that
``spec[0]['flux'][0]`` is bluer than ``spec[0]['flux'][1]``. In other words, the spectrum get bluer on average as one
proceeds down the table, and within an order: pixels on the left are bluer than pixels on the right. If you have no
idea which way is which, make the four possible trial spectra which are flipped relative to each other and run ``xwavecal``
on all of them. The one where ``xwavecal`` succeeds has the correct orientation.

For perspective, here is a print of an NRES spectrum. It is wavelength calibrated so the ``wavelength`` column has meaningful
values here (in Angstroms).

.. code-block:: python

    spec = Table(im['SPECBOX'].data)
    print(spec)

    >>>  id               flux [4096]                            stderr [4096]              pixel [4096] fiber ref_id            wavelength [4096]
    >>> --- --------------------------------------- --------------------------------------- ------------ ----- ------ ----------------------------------------
    >>>   0                     1236.144 .. 567.132  46.16381699989722 .. 33.45280257317763    0 .. 4095     2      0   8875.365322050326 .. 9052.794682947573
    >>>   1            906.7319999999999 .. 455.064  46.49367698945739 .. 33.45280257317763    0 .. 4095     1      1    8707.754989719553 .. 8881.80763072762
    >>>   2                      1120.68 .. 652.032  48.00306240230929 .. 34.35430104077217    0 .. 4095     2      1    8707.822142311728 .. 8881.94973673945
    >>>   3            967.8600000000004 .. 736.932  45.83158299688109 .. 40.22812449021207    0 .. 4095     1      2     8546.46058531058 .. 8717.32420220928
    >>>   4          1161.4319999999998 .. 1124.076 48.285215128442786 .. 45.19736717995861    0 .. 4095     2      2    8546.478280151588 .. 8717.42523057298
    >>>   5                    1008.612 .. 1134.264 48.565728657150814 .. 50.31725350215371    0 .. 4095     1      3   8391.017900052297 .. 8558.812280103835
    >>>   6          1208.976 .. 1630.0800000000004  50.24971641711026 .. 54.74557516366048    0 .. 4095     2      3   8390.995629540508 .. 8558.876525008069
    >>> ...                                     ...                                     ...          ...   ...    ...                                      ...
    >>> 128          1008.612 .. 125.65200000000002  38.41445040606464 .. 33.45280257317763    0 .. 4095     2     64   3963.128098400572 .. 4046.554824188698
    >>> 129 910.1279999999998 .. 146.02800000000005  34.30483930876225 .. 33.45280257317763    0 .. 4095     1     65 3928.6597621432047 .. 4011.7277354354555
    >>> 130            937.2959999999999 .. 139.236  35.13622062772261 .. 33.45280257317763    0 .. 4095     2     65  3928.593357878421 .. 4011.4417999949746
    >>> 131                       47.544 .. 149.424  33.45280257317763 .. 33.45280257317763    0 .. 4095     1     66  3894.679458299859 .. 3977.1857184717724
    >>> 132               0.0 .. 203.75999999999993  33.45280257317763 .. 33.45280257317763    0 .. 4095     2     66  3894.623034269695 .. 3976.9033509112373
    >>> 133               0.0 .. 247.90799999999996  33.45280257317763 .. 33.45280257317763    0 .. 4095     1     67  3861.250017262523 .. 3943.2015758208286
    >>> 134                           0.0 .. 220.74  33.45280257317763 .. 33.45280257317763    0 .. 4095     2     67 3861.2025523440852 .. 3942.9243187156476


Note that if you do not have a blaze corrected spectrum (so your input data only has 2 extensions),
go into the config.ini file and set:
``flux_tol = 0.5`` (to account for bad blaze correction); and ``blaze_corrected_spectrum_name`` to 'None'
or 'empty', or some extension which does not exist.

If you want to look at the processed NRES file I used to make the above example, then process the NRES data contained
in ``xwavecal/tests/data`` with the config file ``xwavecal/data/nres_config.ini``. Note that this will run the full
data reduction pipeline.

Now that the input data is a .fits file with the appropriate data extensions, we proceed.

Setting the reduction stages
----------------------------
In nres_config_wcs_only.ini you will see the section [stages]. This section contains the ordered list of operations
to be done to each input image. You should only need to toggle on or off a few optional stages. The list
looks something like:

.. code-block:: python

    [stages]
    # Reduction stages for a wavelength calibration frame, in order.
    wavecal = [
              #'xwavecal.fibers.MakeFiberTemplate',
              'xwavecal.fibers.IdentifyFibers',
              ...
              'xwavecal.wavelength.IdentifyArcEmissionLines',
              #'xwavecal.wavelength.IdentifyPrincipleOrderNumber',
              ...
              'xwavecal.wavelength.IdentifyArcEmissionLinesLowSN',
              'xwavecal.wavelength.ApplyToSpectrum',
              'xwavecal.wavelength.TabulateArcEmissionLines']

I have shortened the list in places with ... to be brief. This is a list of xwavecal.stages.Stage objects from
``xwavecal``. In principle, they can come from any package you want that conforms to the xwavecal.stages.Stage template.

On your first reduction, you will want to uncomment ``'xwavecal.fibers.MakeFiberTemplate'``. This will make
and write out a few orders of your input spectra as templates. These templates are cross correlated with
later spectra so that the same diffraction order always has the same ``ref_id``. See Section Wavelength calibration settings
for how to change the settings in the config.ini file to select which diffraction orders are saved.

If you do not know the principle order number m0, then uncomment  ``'xwavecal.wavelength.IdentifyPrincipleOrderNumber'``.
This will iterate the entire ``xwavecal`` procedure over the range of trial m0 specified in the config.ini file.

If you do not want the low signal to noise lines saved with your spectrum, comment or delete the last
``'xwavecal.wavelength.IdentifyArcEmissionLinesLowSN'`` stage. Doing so will then save only the lines with a S/N higher
than ``min_peak_snr`` (instead of all those with S/N higher than ``min_overlap_peak_snr``).
See the discussion on the 'LINES' extension in Section 'Output files: Spectra' for more.

Now we can reduce data.

Reducing data
-------------
There are two ways to reduce data: reducing a directory or reducing select files. Both were covered
at the top of this readme for the case of the full reduction pipeline on the included test NRES data. The commands
are identical, except for reducing a directory we specify ``--frame-type wavecal`` so that we do not attempt to
process lampflat files (which is relevant only for the full pipeline).

To reduce a batch of example wavelength calibrations (wavecal types),
we would run:

.. code-block:: bash

    xwavecal_reduce_dir --input-dir data/path/
     --output-dir ~/Downloads --config-file path/to/config.ini --frame-type wavecal

.. code-block:: bash

    xwavecal_reduce --data-paths data/path/1.fits data/path/2.fits
       --output-dir ~/Downloads --config-file path/to/config.ini

where data/path/1.fits data/path/2.fits are wavecal frames.

A .db file will be created at the path specified in ``path/to/config.ini``. If you
re-reduce the same data, the entries in the .db will be updated appropriately. A fiber_template file
will be written out for each wavecal file (and it's path saved in the .db) if you have that stage enabled.

When reducing wavecals, ``xwavecal`` will automatically select the fiber_template
files created which have the nearest observation date.

If you want to fpack (.fz) the output files. You must first install ``libcfitsio``.
E.g. via :code:`sudo apt install libcfitsio-bin` on linux.
Then run the xwavecal reduction command with the added flag: ``--fpack``. The files
are fpacked with a quantization of 10^6 by default. This gives an average error of roughly 10^(-7) on a frame
consisting of gaussian noise only.


Output files
------------

If you are using ``xwavecal`` with 1D extracted spectra, the only output files will be
the wavelength calibrated spectrum and fiber template(s).

Spectra
~~~~~~~

the wavelength calibrated files will be written to the output directory specified in the command
line call. The output file will be almost exactly like that shown in Section 'Formatting the input data',
in that the wavelength column of the 'main' spectrum is now populated.
The blaze corrected spectrum will not have the wavelength column filled in.

the wavelength calibrated files will look like the following.

.. code-block:: python

    from astropy.io import fits
    from astropy.table import table

    im = fits.open('/some/example/image.fits.fz')
    im.info()
    >>> No.    Name      Ver    Type      Cards   Dimensions   Format
    >>> 0  SPECTRUM      1 PrimaryHDU     186   (4096, 4096)   float64
    >>> 1  SPECBOX       1 BinTableHDU     24   135R x 7C   [K, 4096D, 4096D, 4096K, K, K, 4096D]
    >>> 2  BLZCORR       1 BinTableHDU     24   135R x 7C   [K, 4096D, 4096D, 4096K, K, K, 4096D]
    >>> 3  OVERLAP       1 BinTableHDU     23   115R x 7C   [K, K, K, 1000D, 1000D, D, L]
    >>> 4  LINES         1 BinTableHDU     27   4875R x 8C   [K, E, E, D, E, K, D, D]

Notice the two new extensions 'OVERLAP' and 'LINES'. 'OVERLAP' gives the pixel positions of each peak from the red
side of an overlap, and the pixel positions of the matched peaks on the blue side. For example:

.. code-block:: python

    overlaps = Table(im['overlap'].data)
    overlaps.info()
    >>> <Table length=115>
    >>>      name       dtype   shape  n_bad
    >>> -------------- ------- ------- ------
    >>>         ref_id   int64              0
    >>>          fiber   int64              0
    >>> matched_ref_id   int64              0
    >>>          pixel float64 (1000,) 114624
    >>>  matched_pixel float64 (1000,) 114624
    >>>          peaks float64              0
    >>>           good    bool              0

'peaks' gives the number of matched peaks in the overlap between the orders 'ref_id' and 'matched_ref_id'. 'good' is
whether ``xwavecal`` used that overlap to constrain the wavelength solution. `pixel` and `matched_pixel` are best shown
by example:

.. code-block:: python

    print(overlaps[20:25])
    >>> ref_id fiber matched_ref_id        pixel [1000]         matched_pixel [1000]   peaks  good
    >>> ------ ----- -------------- ------------------------- ------------------------ ----- -----
    >>>     20     2             21 137.82643127441406 .. nan  2726.89306640625 .. nan   5.0 False
    >>>     21     2             22 156.71871948242188 .. nan 2711.098388671875 .. nan  13.0  True
    >>>     22     2             23 163.01547241210938 .. nan  2675.88037109375 .. nan   7.0  True
    >>>     23     2             24 25.796588897705078 .. nan  2431.62548828125 .. nan  14.0  True
    >>>     24     2             25 182.21432495117188 .. nan  2622.63330078125 .. nan  14.0  True

    print(overlaps[21]['pixel'][:5])
    print(overlaps[21]['matched_pixel'][:5])
    >>> [156.71871948 178.88464355 307.34054565 323.81674194 436.28128052]
    >>> [2711.09838867 2744.41796875 2939.70263672 2965.02099609 3139.48120117]

In this example, pixel 156.71871948 from the order labelled by ref_id=21 matches pixel 2711.09838867
from the order labelled by matched_ref_id=22. Same with 178.88464355 and 2744.41796875, and so forth. In that overlap
13 such peaks were matched and so ``overlaps[21]['pixel']`` will have 13 non ``np.nan`` elements. The rest will be
``np.nan``.

Now for the 'LINES' extension. This gives the table of pixel and order (ref_id) positions of emission lines, the errors on the line centroid in pixels (``pixel_err``), their wavelengths
under the final model fit by ``xwavecal`` (which you can change in config.ini), and the closest reference wavelength
in the reference line list. The ``pixel_err`` for any line should be close to the width of the line divided by the square root of the counts in the line.

.. code-block:: python

    lines = Table(im['lines'].data)
    lines.info()
    >>> <Table length=4875>
    >>>         name          dtype
    >>> -------------------- -------
    >>>                order   int64
    >>>                pixel float32
    >>>            pixel_err float32
    >>>                 flux float32
    >>>         normed_order float64
    >>>         normed_pixel float32
    >>>                fiber   int64
    >>>           wavelength float64
    >>> reference_wavelength float64

There are 4875 emission lines across both fibers, so roughly 2300 found in either. Note the number found depends directly
on what you set for the emission line signal to noise in config.ini. 'normed_order' and 'normed_pixel' are for calculation
purposes only. 'wavelength' is the wavelength of the line as calculated from the model, and the reference_wavelength is
the reference wavelength. Printing this table gives:

.. code-block:: python

    lines = Table(im['lines'].data)
    print(lines)

    >>> order   pixel      flux   normed_order normed_pixel fiber     wavelength     reference_wavelength
    >>> ----- ---------- -------- ------------ ------------ ----- ------------------ --------------------
    >>>     1  154.85875 3542.028         -1.0   -0.9243669     1  8716.591549446843             8713.654
    >>>     1  220.09575  736.932         -1.0   -0.8925051     1  8720.252997377796             8719.629
    >>>     1  317.38748  669.012         -1.0   -0.8449878     1  8725.647254580183             8724.376
    >>>     1   636.6035   832.02         -1.0   -0.6890825     1   8742.80120092068             8739.781
    >>>     1  736.34924   730.14         -1.0   -0.6403667     1  8747.994093409337             8748.031
    >>>     1 1006.75824 1283.688         -1.0  -0.50829875     1  8761.683388879328              8761.72
    >>>     1     2085.0 1253.124         -1.0  0.018315077     1   8810.92878975186             8810.254
    >>>   ...        ...      ...          ...          ...   ...                ...                  ...
    >>>    67     2591.0  1731.96          1.0    0.2654457     2  3919.178244503971             3919.023
    >>>    67   2927.275  988.236          1.0    0.4296825     2 3925.1112628448705             3925.093
    >>>    67  2963.9255 2822.076          1.0   0.44758272     2  3925.736424580697              3925.72
    >>>    67  3034.8652  7454.22          1.0    0.4822297     2 3926.9344446157725             3927.421
    >>>    67   3137.707 2142.876          1.0    0.5324576     2  3928.643034664589              3928.62
    >>>    67  3201.8215  685.992          1.0   0.56377125     2  3929.691309680273             3929.669
    >>>    67  3381.7493  692.784          1.0   0.65164804     2  3932.563496920231              3932.55
    >>> Length = 4875 rows


We imagine that one can use this table to initialize any other pipeline's wavelength solution.

Fiber templates
~~~~~~~~~~~~~~~

These output files will be a .fits file with one extension. This extension will contain 3 rows (three orders)
of the spectrum processed while ``'xwavecal.fibers.MakeFiberTemplate'`` was included in the ordered stages.
consequently, the fiber template data will be in the exact same format as the 'main' spectrum extension of the input data.

Notes on reduction
------------------

The ``xwavecal`` database handles instruments independently. You can safely reduce data from
separate instruments simultaneously, provided the .fits keywords in :code:`config.ini` are enough
to specify each input .fits file to a unique instrument. By default, ``xwavecal`` uses the instrument
name (nres03 for instance) and the site name (cpt for instance) and a third designator ``instrument2``. All three
identifiers are pulled from the header of the primary .fits extension of the raw data.

One sets in the :code:`config.ini` where to find these specifiers in a .fits header and under what keywords. See
Section 'Data settings'.


Configuring for full data reduction (experimental)
==================================================

One can use ``xwavecal`` to fully reduce their data by adding stages to the [stages] section, and
by adding options to the [reduction] section of the config.ini file. The pipeline is
automatic, however you have to change roughly twice the number of options in the config.ini file and so
errors are more likely to occur. Example configuration files for HARPS and NRES spectrographs
are in the ``xwavecal/example_config/``. The HARPS configuration files are meant to be examples only:
they were made on a limited set of HARPS data. The value of each configuration parameter in
those example files will change often as I tweak the files.

I may document the full data reduction pipeline a later release (perhaps much later). Or, I may move that functionality
to a new git repository.

End note
========
Please contact me if you have issues or find the documentation confusing.

Contributions
=============
We encourage and welcome contributions to ``xwavecal``. The master branch is protected
so the workflow for contributing is first to open a branch and then make a pull request.
One approving review from an administrator is required before the branch can be merged.

License
=======
MIT license, see LICENSE for more details.

