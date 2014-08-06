# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function, division
from astropy.tests.helper import pytest
from .. import catalog

try:
    from scipy.ndimage import convolve
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


@pytest.mark.skipif('not HAS_SCIPY')
def test_catalog_image():
    # TODO: implement me
    pass


@pytest.mark.skipif('not HAS_SCIPY')
def test_catalog_table():
    # TODO: implement me
    pass