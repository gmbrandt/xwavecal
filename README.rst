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

Lampflats must always be reduced before wavelength calibration frames (e.g. ThAr exposures).
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
To reduce a batch of data containing lampflats and arcs, (using the test data supplied with this
repo) we would run (if in the root directory of this repo):

.. code-block:: bash

    echelle_reduce_dir --input-dir echelle/tests/data/
     --output-dir ~/Downloads --config-file echelle/data/nres_config.ini

This will output the reduced data files and intermediate data products (e.g. Trace files) into
~/Downloads. A .db file will be created in the place specified in :code:`nres_config.ini`. If you
re-reduce the same data, the entries in the .db will be updated appropriately.

When reducing arc lamps, :code:`echelle` will automatically select the trace files created
from lampflats which have the nearest observation date.

Reducing select files
---------------------
To reduce a file(s) by specifying paths, we simply specify the data paths separated by spaces.

.. code-block:: bash

    echelle_reduce --data-paths echelle/tests/data/lampflat.fits.fz echelle/tests/data/arc.fits.fz
     --output-dir ~/Downloads --config-file echelle/data/nres_config.ini

Again, :code:`echelle` will automatically handle the ordering. Note that if the lampflat specified is
further away in observation date from the arc specified, :code:`echelle` will find the closest lampflat
in the data base and use that instead. You would want to specify a different (blank) database in order
to force using a lampflat which is very far away.

The :code:`echelle` database handles instruments independently. You can safely reduce data from
separate instruments simulataneously, provided the .fits keywords provided in :code:`config.ini` are enough
to specify each input .fits file to the appropriate instrument. By default, :code:`echelle` uses the instrument
name (nres03 for instance) and the site name (cpt for instance). One sets in the :code:`config.ini` where
to find these specifiers in a .fits header and under what keywords.

Accessing Data Products
=======================
TODO

Configuring a new instrument
============================
TODOD
- tracing, how to measure the different params required for the config.ini file.

Contributions
=============
Contributions to the code are encouraged. The master branch is protected
so the workflow for contributing is first to open a branch, and then make a pull request.
One approving review from an administrator is required.

License
=======
MIT license, see LICENSE for more details.
