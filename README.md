## EnKF-C ##

**EnKF-C** provides a compact generic framework for off-line data assimilation (DA) into large-scale layered geophysical models with the ensemble Kalman filter (EnKF).
Following are its other main features:

- coded in C for GNU/Linux platform;

- model-agnostic;

- can conduct DA in EnKF, ensemble optimal interpolation (EnOI), or hybrid EnKF/EnOI modes;

- permits multiple model grids;

- can handle rectangular, curvilinear, or unstructured horizontal grids, z, sigma or hybrid vertical grids.

**EnKF-C** is coded for simplicity, scalability and robustness. To handle as large systems as possible it uses shared memory capabilities of MPI-3. Here is a snapshot of ensemble spread of sea surface temperature from the 96-member EnKF ocean forecasting system with MOM5 based OFAM3 model (51 x 1500 x 3600 grid), assimilating about 14M super-observations at each 3-day cycle.

![](sst-spread.png)

For more information see [README](https://github.com/sakov/enkf-c/blob/master/enkf/README) and [user guide](https://github.com/sakov/enkf-c/blob/master/enkf/doc/enkf-userguide.pdf). (An older version of the user guide is also available from [arXiv](http://arxiv.org/abs/1410.1233).) Have a feel for how the code works by running the included example.

Checkout **EnKF-C** by running `git clone https://github.com/sakov/enkf-c`
or `svn checkout https://github.com/sakov/enkf-c`.
