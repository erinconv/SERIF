=====
SERIF
=====

**SERIF** is a Python toolbox for performing automatic river segmentation.
This repository also contains different scripts that will help you in the process
of segmentation.

Requirements
------------

To use this toolbox, ArcGIS pro is required. A complete step-by-step to install ArcGIS pro can be found  `here <https://pro.arcgis.com/en/pro-app/latest/get-started/install-and-sign-in-to-arcgis-pro.htm>`_


SERIF toolbox Installation
--------------------------

#. Create a new conda environment.

    .. code-block:: bash

        conda create -n SERIF python=3.9.21 -y

#. Add conda-forge chanel

    .. code-block:: bash

        conda config --add channels conda-forge

#. Now Activate the conda environment

    .. code-block:: bash

        conda activate SERIF

#. Now, install the necessary dependencies

    .. code-block:: bash

        conda install -c conda-forge geopandas rasterio matplotlib gdal=3.10.0 python=3.9.21 -y

Detailed documentation
----------------------

Detailed information and tutorials with the step-by-step can be found in the following link: