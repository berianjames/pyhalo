******
pyhalo
******

:Date: March 10, 2012
:Version: 0.1
:Authors: Berian James
:Web site: https://github.com/berianjames/pyhalo
:Copyright: This document has been placed in the public domain.
:License: This code is released under the MIT license.

=====================================
PyHalo - halo modelling for cosmology
=====================================

This is a Python implementation of the halokit_ library. This code implements the power spectrum backend for the halo modelling and the release version does not include code to perform halo density profile calculations. Much of the release code is also (and better) implemented in Roban Kramer's cosmolopy_ library.

It is plausible that further work on this code will be performed in the future; a fork by Henrik Brink (brinkar) that extends the code considerably already exists and may be posted separately.

.. _halokit: https://github.com/berianjames/matlab-science-functions

.. _cosmolopy: http://roban.github.com/CosmoloPy/

* pyhalo.py contains the fundamental methods for power spectrum calculation
* ips_test.py is a driver script to test computation of the power spectrum
* ps_lin_0.3.dat is a reference power spectrum for comparison
* xi_lin_0.3.dat is the correlation function pair of the reference power spectrum
* nfw.py provides a naive implementation of the NFW halo density profile.
