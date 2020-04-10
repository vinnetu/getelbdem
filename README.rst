getelbdem
=========

A small script to download and marge raster images from the tile server of the Estonian Land Board
(`<https://www.maaamet.ee/>`_). Those tiles are normally available the geoportal of Estonian Land Board
(`<https://geoportaal.maaamet.ee/>`_).

The script takes input coordinates for the two opposing corners of the required area and downloads
required tiles and merges them together into required shape. WARNING: too big range of coordinates
attempts to download huge amount of raster data and thus is not allowed resulting in your IP be banned.

Install
-------
To install download it to your work directory

.. code-block:: bash

    $ git clone https://github.com/vinnetu/getelbdem


Usage
-----

Example for downloading digital elevation model (DEM) of Tartu with 5 m resolution. For defining extents
all coordinates need to be given in L-EST97 (EPSG:3301) coordinate system (`<https://epsg.io/3301/>`_).

.. code-block:: python

    import getelbdem

    # Download area specified by coordinates

    # Location of shapefile of required grid structure
    gridshp = "epk10T.shp"

    # Set the filename

    fn = "tartu.tif"
    # 657556.3,6473499.7 : 661012.6,6475770.9
    x1 = 6473500
    y1 = 657500
    x2 = 6475700
    y2 = 661000

    # Download the file to wanted location defined by L-EST97 coordinate system
    # Resolution can be either 1, 5, or 10
    getelbdem.dlarea(x1, y1, x2, y2, fn, gridshp, resolution = 10)


Example for downloading digital elevation model (DEM) of Tartu with 5 m resolution.

.. code-block:: python

    import pandas as pd
    import getelbdem

    data = pd.read_csv(conf.sitesfn, sep=',')
    data = data.dropna(subset=['X', 'Y'])

    z = data['X'].astype(float)
    data['X'] = data['Y'].astype(float)
    data['Y'] = z
    conf = {'maxrange': 5000,
            'resolution' : 5,
            'gridshp' : '',
            'rasterdir' : ''
            }
    #
    getelbdem.getptsmaps(data, conf)  # Put the a call to the main function in the file.

LICENSE
-------





