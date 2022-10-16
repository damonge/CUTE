import pytest
import numpy as np
import cutepy as cute


def test_bins_smoke():
    b = cute.CuteBin(10, 10.)
    assert (b._bin.nbins == 10)
    assert (np.fabs(b._bin.r_max-10.) < 1E-5)
    assert (np.fabs(b._bin.n_logint) < 1E-5)
    assert (np.fabs(b._bin.i_r_max - 0.1) < 1E-5)
    assert (np.fabs(b._bin.log_r_max - 1.) < 1E-5)

    b = cute.CuteBin(10, 10., is_deg=True)
    assert (b._bin.nbins == 10)
    assert (np.fabs(b._bin.r_max-np.radians(10.)) < 1E-5)
    assert (np.fabs(b._bin.n_logint) < 1E-5)
    assert (np.fabs(b._bin.i_r_max - 1/np.radians(10.)) < 1E-5)
    assert (np.fabs(b._bin.log_r_max - np.log10(np.radians(10.))) < 1E-5)

    # Log binning but no r_min
    with pytest.raises(ValueError):
        b = cute.CuteBin(10, 10., is_log=True)

    # r_min > r_max
    with pytest.raises(ValueError):
        b = cute.CuteBin(10, 10., rmin_log=11., is_log=True)

    b = cute.CuteBin(10, 10., rmin_log=1., is_log=True)
    assert (b._bin.nbins == 10)
    assert (np.fabs(b._bin.r_max-10.) < 1E-5)
    assert (np.fabs(b._bin.n_logint-10.) < 1E-5)
    assert (np.fabs(b._bin.i_r_max - 0.1) < 1E-5)
    assert (np.fabs(b._bin.log_r_max - 1.) < 1E-5)


def test_bins2d_smoke():
    # Log binning but no r_min
    with pytest.raises(ValueError):
        b = cute.CuteBin2D(10, 10., 10, is_mu2=True, is_log1=True)
    # r_min > r_max
    with pytest.raises(ValueError):
        b = cute.CuteBin2D(10, 10., 10, is_mu2=True, is_log1=True,
                           rmin_log1=11)
    # Not clear what to do about second dimension
    with pytest.raises(ValueError):
        b = cute.CuteBin2D(10, 10., 10)

    b = cute.CuteBin2D(10, 10., 5, is_mu2=True)
    assert (b._bin.nbins == 10)
    assert (np.fabs(b._bin.r_max-10.) < 1E-5)
    assert (np.fabs(b._bin.n_logint) < 1E-5)
    assert (np.fabs(b._bin.i_r_max - 0.1) < 1E-5)
    assert (np.fabs(b._bin.log_r_max - 1.) < 1E-5)
    assert (np.fabs(b._bin.nbins2 == 5))
    assert (np.fabs(b._bin.r_max2-1.) < 1E-5)
    assert (b._bin.is_mu2 == 1)
    assert (b._bin.mu2_from_zero == 0)
