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

:Author: Mirek Brandt

If you use this code, please cite **Brandt, G.M. et al. (2019)** which can be found on ArXiv here: XXX
and cite please cite the Zendo DOI: https://zenodo.org/record/3494618

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
This section covers how to wavelength calibrate data which already have a spectrum, and a blaze
corrected spectrum. Using ``xwavecal`` with spectra is preferred.

If you have raw data only and extracting a spectrum is difficult, you may try the experimental data
reduction pipeline included with ``xwavecal``, see section "Configuring for full data reduction".
However, I highly recommend extracting the spectrum first and running ``xwavecal`` in the preferred way.

Information about the input data products are to
be set by users via a :code:`config.ini` file. See the files
:code:`xwavecal/tests/data/test_config.ini` or :code:`xwavecal/example_config/nres_config.ini`
for the configuration to reduce NRES data. I now cover how to wavelength calibrate data, using the
Network of Robotic Echelle Spectrographs (NRES) from Las Cumbres Observatory
as an example.

To reduce the directory of NRES test data included
in this repo, you would run from the command line (after modifying a couple paths
in the config.ini file):

.. code-block:: bash

    xwavecal_reduce_dir --input-dir xwavecal/tests/data/
     --output-dir ~/Downloads --config-file xwavecal/data/nres_config.ini

To do the same reduction by specifying paths, you would run:

.. code-block:: bash

    xwavecal_reduce --data-paths
     xwavecal/tests/data/nres_test_data/cptnrs03-fa13-20190405-0004-w00.fits.fz
      xwavecal/tests/data/nres_test_data/cptnrs03-fa13-20190405-0014-a00.fits.fz
       --output-dir ~/Downloads --config-file xwavecal/data/nres_config.ini

Configuring this wavelength solution to work for your instrument should only involve
making a new config.ini file. The rest of this readme is devoted to setting the config
file for a new instrument where the input data are extracted 1D spectra. I use
NRES as an example.

Configuring for wavelength calibration
======================================
``xwavecal`` is designed in a modular fashion. Each step of the wavelength
calibration is a stage which can be disabled by removing the associated line
in the config.ini file. Wavelength calibrating data which already has spectra simply
means only using the wavelength calibration stages. Using the full experimental pipeline
means enabling the other data reduction stages (e.g. overscan subtraction etc.).

The completed config.ini file is "nres_config_wcs_only.ini", this contains
all the options and settings to reduce NRES data which already has a 1D spectrum
and a 1D blaze corrected spectrum. This repo include raw NRES data, which has to
be reduced with nres_config.ini (which includes all the overscan subtraction, spectral extraction etc. stages).

We start by telling the config.ini where the database for the reduced data should live.

Pointing to the database and line list
--------------------------------------
Before reducing, copy the nres_config_wcs_only.ini file to a new location, rename it for your instrument, and
change :code:`database_path` under the [reduction] section to the path where you
want to the database to exist. The parent folder for the database must already exist. E.g. for myself,
this is :code:`"/home/gmbrandt/Downloads/pipeline.db"` . The surrounding :code:`" "` quotes must be there for
the config file to process properly.

The database will keep track of all your processed files. All processed calibration files are saved under the
table :code:`caldata` in the .db file specified.

In [reduction] change the ``line_list_path`` as well. If your instrument is optical,
the ESO ThAr atlas *may* work fine (included with this repo). However you must still
provide a valid path in the config.ini file.


Data settings
---------------
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

- ``primary_data_extension`` is the fits extension where all the relevant header data is stored such as
the observation date, instrument name etc. These are used for writing out the file with an informative name.
- ``files_contain`` is a list of strings, where each string must be present in the input file types. The default
is ['.fits'] in which case only files with '.fits' in the name are reduced. For example:
  * If I had two files: 'IRDA003.fits' and 'IRDB002.fits', and I wanted to only process IRDA and .fits files,
    I would set ``files_contain = ['.fits', 'IRDA']``

header_keys
~~~~~~~~~~~

``header_keys`` is a python dictionary. The values of the dictionary are the header keywords
in your raw data that give things like the read noise, the observation date, etc. The keys
are the standard keys understood by ``xwavecal``. Some of these keys are
  * 'type' (the frame type e.g. lampflat)
  * 'gain' (the gain in e-/ADU)
  * 'read_noise' (the read noise in e-)
  * 'fiber_state' (the string which gives which fibers are lit and with what. See more later)
  * 'observation_date' (observation date, see time_format later.)
  * 'instrument'
  * 'instrument2'
  * 'site_name'
  * 'unique_id'

``instrument``, ``instrument2``, ``site_name`` are designators which are how the data base would look up
processed data. E.g. for NRES, I set

.. code-block:: python

               ...
               'instrument': 'TELESCOP',
               'instrument2': 'INSTRUME',
               'site_name': 'SITEID',
               ...

This means that processed data will be stored in the database with telescope name, instrument name, and the
ID of our site. These data are stored in NRES frames under the header keys 'TELESCOP', 'INSTRUME' 'SITEID'.

``observation_date`` is the .fits header key which gives the observation date of the frame.
One must set time_format (see further down in this section) to agree with the format of the .fits value given
by the ``observation_date`` key.

For ``fiber_state``, the NRES and HARPS store this in a single string in 'OBJECTS' and 'ESO DPR TYPE', respectively.
For NRES the value of the header looks like ``thar&thar&none`` for a frame with Thorium-Argon (ThAr) lit on fibers 0,1 and
fiber 2 unlit. For HARPS, the same configuration (but no third fiber since it does not exist) would be
``WAVE,WAVE,THAR2``. We will convert ``WAVE,WAVE,THAR2`` to ``thar&thar&none`` with the type_keys next.

type_keys
~~~~~~~~~

``type_keys`` is by far the most confusing part of configuring an instrument. This may change in a future release.
``type_keys`` is a dictionary which takes the value of any .fits header value and converts it in place. E.g. if my .fits
header for my raw data for the ``fiber_state`` was: ``{'ESO DPR TYPE': 'WAVE,WAVE,THAR2'}`` and I set
``type_keys = {..., 'WAVE,WAVE,THAR2': 'thar&thar&none'}``, then any time ``xwavecal`` reads the ``fiber_state`` item
it will read 'thar&thar&none'

fiber_state
~~~~~~~~~~~
A note on ``fiber_state``: One must convert whatever ``fiber_state`` value in your .fits file to be
of the string format interpretable by ``wavecal``. This format is always ``fiber0lamp&fiber1lamp&fiber2lamp``.
Where ``fiberxlamp`` is the type of light coming through that lamp. E.g. if I had a fictional instrument with two
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
TODO: rename this to ``cal`` so as not to cause confusion.

Important note
~~~~~~~~~~~~~~
None of these translations will ever be saved onto the fits header of your output data product. The fits
header of your data will *not* have ``read_noise`` etc appended as extra headers. Setting header_keys and type_keys
builds a translator which understands how to interpret your fits header, ``xwavecal`` does not modify existing header keys.


time_format
~~~~~~~~~~~

In [reduction], ``time_format`` is the time format of the observation date from
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
- ``template_trace_id ``
- ``overlap_min_peak_snr ``
- ``max_red_overlap ``
- ``max_blue_overlap ``
- ``global_scale_range ``
- ``min_peak_snr ``
- ``approx_detector_range_angstroms ``
- ``approx_num_orders ``
- ``principle_order_number``
- ``m0_range ``
- ``flux_tol``

There are several other parameters you will most likely not need to change.
Let us go through the pertinant ones in the list above one-by-one:

- ``main_spectrum_name`` : this is the name of the .fits extension that contains
  the BinTableHDU of the spectrum that ``xwavecal`` will calibrate.
- ``blaze_corrected_spectrum_name`` : this is the name of the .fits extension that contains
  the BinTableHDU of the blaze corrected spectrum that ``xwavecal`` will use to aid its
  calibration of ``main_spectrum_name``. If you do not have a blaze corrected spectrum, set
  this to some string (that is not in the data) such as ``'None'``.
- ``template_trace_id `` : this is the trace id (id column in the input spectrum) for the
  diffraction order that you want to save as a template. This template will be used to identify this same
  diffraction order in all subsequent spectra you reduce. It will have a ref_id associated with it
  such that the diffraction order number understood by ``xwavecal`` is ``ref_id + m0`` where
  ``m0`` is the principle order number. I recommend setting the id to some middle order on the detector.
- ``ref_id`` : this is the reference id you wish to assign the template spectrum such that the
  diffraction order number understood by ``xwavecal`` for the template spectrum is ``ref_id + m0`` where
  ``m0`` is the principle order number.
- ``overlap_min_peak_snr `` : the minimum signal to noise for an emission peak to be considered in the overlap algorithm.
  see Brandt et al. 2019 for a discussion of the overlap algorithm. I recommend this be set to something low like 5. In
  general, overlap fitting works better if more peaks are detected. For NRES we use 5 and detect ~4000 peaks.
- ``flux_tol`` : If two emission peaks from neighboring orders have flux f1 and f2, ``flux_tol`` is
  the maximum allowed value of abs(f1 - f2)/(mean(f1, f2)) for two peaks to be considered
  a matched pair in the overlap algorithm.
- ``min_peak_snr `` : the minimum signal to noise for an emission peak to be used to constrain the wavelength
  solution after overlap detection. This should be something reasonable like 10 or 20 so
  as to detect between 1000 and 2000 emission lines. Weak lines are often contamination from trace elements
  (which are not in reference line lists and so would throw off our algorithm).
- ``max_red_overlap `` : The maximum allowed pixel coordinate for a peak to be considered for our overlap algorithm.
- ``max_blue_overlap `` :
  * The overlap algorithm will try to match peaks from
    (0, max_red_overlap) to (max_pixel, max_pixel - max_blue_overlap). Where max_pixel is the width of
    your detector (in x, i.e. number of columns, e.g. 4096).
- ``approx_detector_range_angstroms ``: If the spectrograph covers the spectral range 3000A to 9000A, then
  ``approx_detector_range_angstroms = 5000``. Note this value does not need to be precise.
- ``approx_num_orders `` : approximate number of distinct diffraction orders in the spectrum. E.g. 67 for NRES.
  Note this is not the number of traces (visible light streaks on the echelle detector) but the number of diffraction orders.
  I.e. num_of_traces/num_of_lit_fibers. This does not need to be precise either.
- ``global_scale_range ``: See Brandt et al. 2019 for a discussion of the global scale.
  This is the range about the initial guess where ``xwavecal`` will search for the global scale.
  * For example: if the guess generated by ``xwavecal`` is ``K`` and if ``global_scale_range = (0.8, 1.2)``
    then ``xwavecal`` will search for the global scale between ``0.8K`` and ``1.2K``.
- ``principle_order_number``: This needs to exactly correct. This is the true diffraction order
  number of the diffraction order with ref_id = 0. If you do not know this, insert the m0 identification stage
  (I will cover how to do this in a following section), and set ``m0_range`` to a reasonable range of values.
- ``m0_range ``: the range of possible ``m0`` (principle order number) values. This is only used if you
  are searching for ``m0`` (i.e. if you have included 'xwavecal.wavelength.IdentifyPrincipleOrderNumber' in
  the set of stages for wavecal frames). I will discuss this more later.


Reducing a directory of data
----------------------------
To reduce a batch of example data containing lampflats and wavelength calibrations (hereafter wavecal),
we would run (if in the root directory of this repo):

.. code-block:: bash

    xwavecal_reduce_dir --input-dir xwavecal/tests/data/
     --output-dir ~/Downloads --config-file xwavecal/data/nres_config.ini

This will output the reduced data files and intermediate data products (e.g. Trace files) into
~/Downloads. A .db file will be created in the place specified in :code:`nres_config.ini`. If you
re-reduce the same data, the entries in the .db will be updated appropriately.

When reducing wavecals, ``xwavecal`` will automatically select the trace files created
from lampflats which have the nearest observation date.

If you want to fpack (.fz) the output files. You must first install :code:`libcfitsio`.
E.g. via :code:`sudo apt install libcfitsio-bin` on linux.
Then run the xwavecal reduction command with the added flag: :code:`--fpack`. The files
are fpacked with a quantization of 10^6 by default. This gives an error of roughly 10^(-7) on a frame
consisting of gaussian noise only.

Reducing select files
---------------------
To reduce files by specifying paths, specify the data paths separated by spaces:

.. code-block:: bash

    xwavecal_reduce --data-paths
     xwavecal/tests/data/nres_test_data/cptnrs03-fa13-20190405-0004-w00.fits.fz
      xwavecal/tests/data/nres_test_data/cptnrs03-fa13-20190405-0014-a00.fits.fz
       --output-dir ~/Downloads --config-file xwavecal/data/nres_config.ini

For clarity, w00 is a lampflat and a00 is a ThAr exposure. Again, ``xwavecal`` will automatically reduce lampflats and
generate trace files first.
Note that if the lampflat specified is further from the wavecal in observation date than another lampflat
you already reduced which is in the database, ``xwavecal`` will find the closest lampflat
in the data base and use that instead. You would want to specify a different (blank) database in order
to force using a lampflat which is very far away. Again, files can be compressed with fpack (after installing
:code:`libcfitsio`) by adding :code:`--fpack` to the command line call.


Configuring for full data reduction (experimental)
==================================================

One can use ``xwavecal`` to fully reduce their data by adding stages to the [stages] section, and
by adding options to the [reduction] section of the config.ini file. The pipeline is
automatic, however you have to change roughly twice the number of options in the config.ini file and so
errors are more likely to occur. Example configuration files for IRD (Subaru), HARPS, and NRES spectrographs
are in the ``xwavecal/example_config/``. Those configuration files are meant to be examples only: they were made
on a limited set of IRD and HARPS data. The pipeline may not function well on all data from those instruments
using my example configuration files. The value of each configuration parameter will in those example files will
change often as I tweak the files.


Configuring a new instrument
----------------------------


Indicating header keywords
--------------------------
We need to tell ``xwavecal`` where the read_noise, etc... lies in the fits headers
of the input raw data files.

We first copy one of the example config.ini files inside of :code:`xwavecal/data/`. Next
we uncomment out the stage :code:`MakeFiberTemplate` in the section [stages].

In the section [data] of the config file, specify in header_keys which header keys
in the fits file correspond to which observables (e.g. read_noise for harps is RON).

In the type_keys, specify which outputs of the :code:`type` header key correspond to
a lampflat or a wavecal. E.g. for nres, wavecal frames have the value :code:`DOUBLE` under the header key :code:`OBSTYPE`. Therefore in type_keys, I would
have an entry :code:`{'DOUBLE': 'wavecal'}`, and in header_keys, I would have an entry
:code:`{'type': 'OBSTYPE'}`. One can insert tuples into header_keys. I.e. if you need information
from more than one field. E.g. for HARPS, I made my unique identifier (mjd-obs, chip id) because
each raw harps frame has both the blue and the red parts of the spectra as different chips.

Orientating the frames
----------------------
In section [stages] are all the reduction stages. For the finished HARPS config file,
you will notice some of the first stages are Rot90 and FlipHoriz, which rotate the frame
90 degrees counter-clockwise and flip it about the vertical axis. We do this so that the dispersion
of the frame agrees with the NRES (the nres_config.ini file does not have these flips accordingly).
Prior to tracing, but after overscan trimming, every frame must be orientated so that:
The wavelength of any given diffraction order increases from left to right in pixel (x=0 to x=Nx), and:
The diffraction orders become overall bluer as one heads up the detector (bottom to top, y=0 to y=Ny).

Making the template prior to first reduction
--------------------------------------------
In section [reduction], :code:`template_trace_id` gives the trace id (:code:`id` in the trace.fits files created)
for the diffraction order
that ``xwavecal`` will use to make a template from on the first wavecal frame you reduce. For HARPS,
I set :code:`template_trace_id = 10` arbitrarily. I recommend you don't select diffraction orders
that are known to be problematic (e.g. are near the edge). Specify the paths in the config.ini file
so that they are where you want them. Namely, you need to specify the line list path and the .db database path.

Next, reduce a lampflat and wavecal via :code:`xwavecal_reduce_dir`, or with  :code:`xwavecal_reduce`. The lampflat
reduction will make a trace file, a blaze file, and a processed lampflat file.

Reducing any wavecal will produce a template. The template is a . For all wavecal files which resemble
those you just processed, for all of time (provided you don't delete the database or the fibers.fits file)
you will never need to make another template. This template is just the 1d spectrum of the
order specified by :code:`template_trace_id`. Echelle looks for an order with a matching spectrum, and labels
it with the reference id (:code:`ref_id`) given in [reduction] of the config.ini. This template, along with
any processed files (e.g. the trace files etc) will be saved in the database .db file at the path
specified in the config.ini file.

Reduction
=========

Lampflats must always be reduced before wavelength calibration frames (e.g. Thorium-Argon (ThAr) exposures).
This is because lampflats are used to determine where the light falls, which is in turn
used to extract data. This ordering is handled for you if you supply at least one lampflat in the
data to reduce.

Lampflats
---------

Wavelength calibration files
----------------------------

Notes on reduction
------------------

The ``xwavecal`` database handles instruments independently. You can safely reduce data from
separate instruments simulataneously, provided the .fits keywords provided in :code:`config.ini` are enough
to specify each input .fits file to the appropriate instrument. By default, ``xwavecal`` uses the instrument
name (nres03 for instance) and the site name (cpt for instance). One sets in the :code:`config.ini` where
to find these specifiers in a .fits header and under what keywords.

Accessing Data Products
=======================
In this section I cover how to access the various output data products.

Traces
------
Traces are the y positions, as a function of x, of the center of flux for a given diffraction order. E.g. the ladder-rungs
on an echelle spectrograph. If your input lampflats have 67 visible orders, and are 4096 pixels wide, then the output
trace files that ``xwavecal`` generates are tables with 67 rows and 4096 + 1 columns. The additional column contains
the trace id. The column headers are :code:`id` for the trace id, and :code:`centers` for the y positions of the trace.

Trace files by default have :code:`_trace` appended onto the end of the filename (but before the filetype extension).

Assume the output trace file is named :code:`X_trace.fits.fz`. You can access the table of traces by doing the following.

.. code-block:: python

    from astropy.io import fits
    from astropy.table import table

    trace = Table(fits.open('X_trace.fits.fz')['TRACE'].data)

You could do the following to plot the trace centers atop the raw data.

.. code-block:: python

    import matplotlib.pyplot as plt

    trace = Table(fits.open('X_trace.fits.fz')['TRACE'].data)
    im = fits.open('lampflat.fits.fz')[1].data

    plt.imshow(im)

    for tr in trace['centers']:
        plt.plot(tr)

The output will look like:

Blaze
-----

Wavecal (wavelength calibration) data products
----------------------------------------------
Here we address how see the extracted spectra and other products from a wavecal lamp file,
including the spectrum's wavelength solution, and
the fluxes and associated standard 1-sigma uncertainties. The data products associated with
a calibration file are

.. code-block:: python

    import matplotlib.pyplot as plt

    im = fits.open('wavecal.fits.fz')
    im.info()

TODO ...
The wavelength model and wavelength coefficients are saved in the fits header
for each spectrum extension that has wavelengths. The model and coefficients have
keywords MODEL and MCOEFFS, respectively, in the header.
ID keywords: IDTRACE, IDBLAZE, IDLIST, IDTEMPL

What is :code:`ref_id`


The reference line list
-----------------------
We include the original ThAr (Thorium-Argon) atlas from the European Southern Observatory (ESO). This was retrieved
from http://www.eso.org/sci/facilities/paranal/instruments/uves/tools/tharatlas.html in late
2019. This line list was designed for spectrographs with a resolving power (R) of 100,000, and thus
it may not be suited for your instrument if it has a lower or larger R. Moreover, the wavelengths are air wavelengths.
It is up to you to download a line list suitable for your instrument (if the ThAr atlas is not suitable)
and correct the line list for the index of refraction of air if necessary.

Contributions
=============
We encourage and welcome contributions to ``xwavecal``. The master branch is protected
so the workflow for contributing is first to open a branch and then make a pull request.
One approving review from an administrator is required before the branch can be merged.

License
=======
MIT license, see LICENSE for more details.
