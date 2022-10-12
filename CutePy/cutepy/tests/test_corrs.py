import pytest
import numpy as np
import cutepy as cute
import healpy as hp
from scipy.special import eval_legendre


def test_angcorr():
    nside = 128
    npix = hp.nside2npix(nside)
    ls = np.arange(3*nside)
    # 1-degree smoothing
    sig_beam = np.radians(1.)
    sig_pix = np.radians(17.71/nside)
    cl = np.exp(-ls*(ls+1)*sig_beam**2)/(ls+10.)
    cl_pixwin = cl*np.exp(-ls*(ls+1)*sig_pix**2)
    ra, dec = hp.pix2ang(nside, np.arange(npix), lonlat=True)

    b = cute.CuteBin(20, 10., is_deg=True)
    xis_A = []
    xis_X = []
    for i in range(10):
        mp = hp.synfast(cl, nside)
        # Auto-correlation
        corr = cute.CuteAngularCorrelator(b, [ra, dec], weights=mp)
        corr.run()
        xi = corr.get_xi()[1]
        xis_A.append(xi)
        # Cross-correlation
        corr = cute.CuteAngularCorrelator(b, [ra, dec], weights=mp,
                                          coords2=[ra, dec], weights2=mp)
        corr.run()
        xi = corr.get_xi()[1]
        xis_X.append(xi)
    xis_A = np.array(xis_A)
    xiA_mean = np.mean(xis_A, axis=0)
    xiA_std = np.std(xis_A, axis=0)
    xis_X = np.array(xis_X)
    xiX_mean = np.mean(xis_X, axis=0)
    xiX_std = np.std(xis_X, axis=0)

    th = b.get_r_values()
    cth = np.cos(np.radians(th))
    xi_th = np.array([np.sum(cl_pixwin*(2*ls+1)*eval_legendre(ls.astype(int), c))
                     for c in cth])/(4*np.pi)

    # Uncomment if you'd like to plot the results
    import matplotlib.pyplot as plt
    plt.figure()
    plt.errorbar(th, xiA_mean, yerr=xiA_std, fmt='k.')
    plt.plot(th, xi_th, 'r-')
    for xi in xis_A:
        plt.plot(th, xi, 'b-', alpha=0.1)
    plt.figure()
    plt.errorbar(th, xiX_mean, yerr=xiX_std, fmt='k.')
    plt.plot(th, xi_th, 'r-')
    for xi in xis_X:
        plt.plot(th, xi, 'b-', alpha=0.1)
    plt.show()
    assert np.all(np.fabs(xi_th-xiA_mean) < 2*xiA_std)
    assert np.all(np.fabs(xi_th-xiX_mean) < 2*xiX_std)
