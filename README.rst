echelle
=======
.. image:: https://coveralls.io/repos/github/gmbrandt/echelle/badge.svg?branch=master
    :target: https://coveralls.io/github/gmbrandt/echelle?branch=master

.. image:: https://travis-ci.org/gmbrandt/echelle.svg?branch=master
    :target: https://travis-ci.org/gmbrandt/echelle


A library of routines for wavelength calibrating echelle 
spectrographs for high precision radial velocity work. Included is
a limited data reduction pipeline which does: overscan subtraction and trimming, gain
normalization, tracing, extraction and wavelength calibration. All output products are
saved in a user-set sqlite3 database.

At best, using ``echelle`` only requires editing a config.ini file for your data.
I cover how to do that in this readme.

Installation
============
:code:`echelle` is installed via pip by running

.. code-block:: bash

    pip install .

While in the root directory of this repository. It can also be installed by running

.. code-block:: bash

    pip install git+https://github.com/gmbrandt/echelle

Wavelength Calibrating Spectrum
===============================
This section covers how to wavelength calibrate data which already have a spectrum, and a blaze
corrected spectrum. Using ``echelle`` with spectra is preferred.

If you have raw data only and extracting a spectrum is difficult, you may try the experimental data
reduction pipeline included with :code:`echelle`, see section "Reducing Data". However, I highly recommend extracting
the spectrum first and running ``echelle`` in the preferred way.

Information about the input data products are to
be set by users via a :code:`config.ini` file. See the files
:code:`echelle/tests/data/test_config.ini` or :code:`echelle/example_config/nres_config.ini`
for the configuration to reduce NRES data. I now cover how to wavelength calibrate data, using the
Network of Robotic Echelle Spectrographs (NRES) from Las Cumbres Observatory
as an example. Configuring this wavelength solution to work for your instrument only involves (except
for rare instances) making a new config.ini file. In this example, I will cover how to set
every option in the "nres_config_wcs_only.ini" file.

Pointing to the database
------------------------
Before reducing, copy the nres_config.ini file to a new location (or edit it in place) and
change :code:`database_path` under the :code:`reduction` section to the path where you
want to the database to exist. The parent folder for the database must already exist. E.g. for myself,
this is :code:`"/home/gmbrandt/Downloads/pipeline.db"` . The surrounding :code:`" "` quotes must be there for
the config file to process properly.

The database will keep track of all your processed files. All processed calibration files are saved under the
table :code:`caldata` in the .db file specified.


Data settings
---------------
Here we tell :code:`echelle` via the config file where various information lies in the header of
your data.

In section [data] we will need to edit:

- ``primary_data_extension``
- ``files_contain``
- ``header_keys``
- ``type_keys``

data_class is also editable, but most likely will not need to be changed. data_class is the
Python object used to load in your data. The default ``echelle.images.Image`` should be fine for your data.

``primary_data_extension``

In [reduction], we will need to change ``time_format``. This is the time format of the observation date from
the fits header. This must be a string contained in double quotes ``" "`` and understood by
``datetime.datetime.strftime``. Then replace single ``%`` with ``%%`` (to fix a quirk of using a config file).

Wavelength calibration settings
-------------------------------
To wavelength calibrate your data, the following settings in config.ini may need to be changed:

- ``line_list_path``
- ``data_base_path``
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

The ``principle_order_number`` (which we call ``m0``) is the true diffraction order index of the diffraction order
with reference id 0 (:code:`ref_id`) in the extracted sectrum. If you know it for your instrument,
great. If not: go to [stages] inside of the config.ini file and uncomment the stage
:code:`IdentifyPrincipleOrderNumber`. Then set ``m0_range`` to a suitably wide range
which encompasses your guess for where ``m0`` likely lies. If you have no idea, set ``m0_range = (5, 200)``.
Most echelle spectrographs have ``m0`` between 10 and 100.

``flux_tol`` is the tolerance (float between 0 and 1) to which two emission
peaks must agree to be considered a true match in the overlap algorithm.
Thus, if your blaze correction is poor (or non-existent) you should change ``flux_tol`` to 0.5.


Reducing a directory of data
----------------------------
To reduce a batch of example data containing lampflats and wavelength calibrations (hereafter wavecal),
we would run (if in the root directory of this repo):

.. code-block:: bash

    echelle_reduce_dir --input-dir echelle/tests/data/
     --output-dir ~/Downloads --config-file echelle/data/nres_config.ini

This will output the reduced data files and intermediate data products (e.g. Trace files) into
~/Downloads. A .db file will be created in the place specified in :code:`nres_config.ini`. If you
re-reduce the same data, the entries in the .db will be updated appropriately.

When reducing wavecals, :code:`echelle` will automatically select the trace files created
from lampflats which have the nearest observation date.

If you want to fpack (.fz) the output files. You must first install :code:`libcfitsio`.
E.g. via :code:`sudo apt install libcfitsio-bin` on linux.
Then run the echelle reduction command with the added flag: :code:`--fpack`. The files
are fpacked with a quantization of 10^6 by default. This gives an error of roughly 10^(-7) on a frame
consisting of gaussian noise only.

Reducing select files
---------------------
To reduce files by specifying paths, specify the data paths separated by spaces:

.. code-block:: bash

    echelle_reduce --data-paths
     echelle/tests/data/nres_test_data/cptnrs03-fa13-20190405-0004-w00.fits.fz
      echelle/tests/data/nres_test_data/cptnrs03-fa13-20190405-0014-a00.fits.fz
       --output-dir ~/Downloads --config-file echelle/data/nres_config.ini

For clarity, w00 is a lampflat and a00 is a ThAr exposure. Again, :code:`echelle` will automatically reduce lampflats and
generate trace files first.
Note that if the lampflat specified is further from the wavecal in observation date than another lampflat
you already reduced which is in the database, :code:`echelle` will find the closest lampflat
in the data base and use that instead. You would want to specify a different (blank) database in order
to force using a lampflat which is very far away. Again, files can be compressed with fpack (after installing
:code:`libcfitsio`) by adding :code:`--fpack` to the command line call.


Reducing Raw Data (experimental)
================================

One can use :code:`echelle` to fully reduce their data by adding stages to the [stages] section, and
by adding options to the [reduction] section of the config.ini file. The pipeline is
automatic, however you have to change roughly twice the number of options in the config.ini file and so
errors are more likely to occur. Example configuration files for IRD (Subaru), HARPS, and NRES spectrographs
are in the ``echelle/example_config/``. Those configuration files are meant to be examples only: they were made
on a limited set of IRD and HARPS data. The pipeline may not function well on all data from those instruments
using my example configuration files. The value of each configuration parameter will in those example files will
change often as I tweak the files.


Configuring a new instrument
----------------------------


Indicating header keywords
--------------------------
We need to tell :code:`echelle` where the read_noise, etc... lies in the fits headers
of the input raw data files.

We first copy one of the example config.ini files inside of :code:`echelle/data/`. Next
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
that :code:`echelle` will use to make a template from on the first wavecal frame you reduce. For HARPS,
I set :code:`template_trace_id = 10` arbitrarily. I recommend you don't select diffraction orders
that are known to be problematic (e.g. are near the edge). Specify the paths in the config.ini file
so that they are where you want them. Namely, you need to specify the line list path and the .db database path.

Next, reduce a lampflat and wavecal via :code:`echelle_reduce_dir`, or with  :code:`echelle_reduce`. The lampflat
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

The :code:`echelle` database handles instruments independently. You can safely reduce data from
separate instruments simulataneously, provided the .fits keywords provided in :code:`config.ini` are enough
to specify each input .fits file to the appropriate instrument. By default, :code:`echelle` uses the instrument
name (nres03 for instance) and the site name (cpt for instance). One sets in the :code:`config.ini` where
to find these specifiers in a .fits header and under what keywords.

Accessing Data Products
=======================
In this section I cover how to access the various output data products.

Traces
------
Traces are the y positions, as a function of x, of the center of flux for a given diffraction order. E.g. the ladder-rungs
on an echelle spectrograph. If your input lampflats have 67 visible orders, and are 4096 pixels wide, then the output
trace files that :code:`echelle` generates are tables with 67 rows and 4096 + 1 columns. The additional column contains
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
We encourage and welcome contributions to :code:`echelle`. The master branch is protected
so the workflow for contributing is first to open a branch and then make a pull request.
One approving review from an administrator is required before the branch can be merged.

License
=======
MIT license, see LICENSE for more details.
