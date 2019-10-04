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

Installation
============
:code:`echelle` is installed via pip by running

.. code-block:: bash

    pip install .

While in the root directory of this repository.

Reducing Data
=============
Information about the input data products are to
be set by users via a :code:`config.ini` file. See the files
:code:`echelle/tests/data/test_config.ini` or :code:`echelle/data/nres_config.ini`
for the configuration to reduce NRES data. I now cover how to reduce data, using the
Network of Robotic Echelle Spectrographs (NRES) from Las Cumbres Observatory
as an example.

Lampflats must always be reduced before wavelength calibration frames (e.g. Thorium-Argon (ThAr) exposures).
This is because lampflats are used to determine where the light falls, which is in turn
used to extract data. This ordering handled for you if you supply at least one lampflat in the
data to reduce.

Pointing to the database
------------------------
Before reducing, copy the nres_config.ini file to a new location (or edit it in place) and
change :code:`database_path` under the :code:`reduction` section to the path where you
want to the database to exist. The parent folder for the database must already exist. E.g. for myself,
this is :code:`"/home/gmbrandt/Downloads/pipeline.db"` . The surrounding :code:`" "` quotes must be there for
the config file to process properly.

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

Arc (wavecal) data products
-----------------
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

Configuring a new instrument
============================
In this section I cover how to reduce data from a new instrument. If overscan subtraction and overscan
trimming is not required, all that is needed (ideally) is a new config.ini file. I show this
procedure by example on HARPS. The final config file we will make is See :code:`echelle/data/HARPS_chipX_config.ini`.
where X is either 1 or 2.

Indicating header keywords
--------------------------
We need to tell :code:`echelle` where the read_noise, etc... lies in the fits headers
of the input raw data files.

We first copy one of the example config.ini files inside of :code:`echelle/data/`. Next
we uncomment out the stage :code:`MakeFiberTemplate` in the section [stages].

In the section [data] of the config file, specify in header_keys which header keys
in the fits file correspond to which observables (e.g. read_noise for harps is RON).

In the type_keys, specify which outputs of the :code:`type` header key correspond to
a lampflat or a wavecal. E.g. for nres, wavecal frames have the value
 :code:`DOUBLE` under the header key :code:`OBSTYPE`. Therefore in type_keys, I would
have an entry :code:`{'DOUBLE': 'wavecal'}`, and in header_keys, I would have an entry
:code:`{'type': 'OBSTYPE'}`. One can insert tuples into header_keys. I.e. if you need information
from more than one field. E.g. for HARPS, I made my unique identifier (mjd-obs, chip id) because
each raw harps frame has both the blue and the red parts of the spectra as different chips.

Orientating the frames
----------------------
In section [stages] are all the reduction stages. For the finished HARPS config file,
you will notice some of the first stages are Rot90 and FlipHoriz, which rotate the frame
90 degrees counter-clockwise and flip it about the vertical axis. We do this so that the dispersion
of the frame agrees with nres (the nres_config.ini file does not have these flips accordingly).
Prior to tracing, but after overscan trimming, every frame must be orientated so that:
The wavelength of any given diffraction order increases from left to right in pixel, and:
The diffraction orders become overall bluer as one heads up the detector (bottom to top).

Making the template prior to first reduction
--------------------------------------------
In section [reduction], :code:`template_trace_id`
gives the trace id (:code:`id` in the trace.fits files created) for the diffraction order
that :code:`echelle` will use to make a template from on the first arc frame you reduce. For HARPS,
I set :code:`template_trace_id`=10 arbitrarily. I recommend you don't select diffraction orders
that are known to be problematic (e.g. are near the edge).

Next, reduce a lampflat and wavecal via :code:`echelle_reduce_dir`, or with  :code:`echelle_reduce`. You will notice
that :code:`echelle` will abort the wavelength solution on the first arc you input. That is because no
template exists in the database. However, you just created one. This template is just the 1d spectrum of the
order specified by :code:`template_trace_id`. Echelle looks for an order with a matching spectrum, and labels
it with the reference id (:code:`ref_id`) given in [reduction] of the config.ini. Re-run the reduction
and the frame will be wavelength calibrated (although incorrectly because you have not changed the required settings in
the config.ini file!).

Configuring settings to wavelength calibrate your instrument
------------------------------------------------------------

 Each raw
HARPS frame has two chips (chip 1 is bluer and chip 2 is redder). The config file
for the two chips are identical except for two items: the number of diffraction orders on
the detector and the principle order number. We cover those first.

The principle order number is the true diffraction order index of the diffraction order
with reference id 0 (:code:`ref_id`) in the extracted sectrum. If you know it for your instrument,
great. If not: go to [stages] inside of the config.ini file and uncomment the stage
:code:`IdentifyPrincipleOrderNumber`.



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
