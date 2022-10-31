import numpy as np
import cutepy as cute
import healpy as hp
from scipy.special import eval_legendre
from .utils import get_comm_world


def test_angcorr_mpi():
    comm = get_comm_world()

    nside = 128
    npix = hp.nside2npix(nside)
    ls = np.arange(3*nside).astype(int)
    # 1-degree smoothing
    sig_beam = np.radians(1.)
    cl = np.exp(-ls*(ls+1)*sig_beam**2)/(ls+10.)
    ra, dec = hp.pix2ang(nside, np.arange(npix), lonlat=True)

    b = cute.CuteBin(20, 10., is_deg=True)
    mp = hp.synfast(cl, nside)
    # Auto-correlation (no MPI)
    corr1 = cute.CuteAngularCorrelator(b, [ra, dec], weights=mp)
    corr1.run()
    # Auto-correlation (with MPI)
    corr2 = cute.CuteAngularCorrelator(b, [ra, dec], weights=mp)
    corr2.run(comm=comm)

    assert np.allclose(corr1.ww, corr2.ww, atol=1E-8, rtol=1E-5)
    assert np.allclose(corr1.ct, corr2.ct, atol=1E-8, rtol=1E-5)


def test_angcorr_tracer():
    nside = 128
    npix = hp.nside2npix(nside)
    ls = np.arange(3*nside).astype(int)
    # 1-degree smoothing
    sig_beam = np.radians(1.)
    sig_pix = np.radians(17.71/nside)
    cl = 0.001*np.exp(-ls*(ls+1)*sig_beam**2)/(ls+10.)
    cl_pixwin = cl*np.exp(-ls*(ls+1)*sig_pix**2)
    ra, dec = hp.pix2ang(nside, np.arange(npix), lonlat=True)

    b = cute.CuteBin(20, 10., is_deg=True)

    # Random catalog
    ra_r = ra
    dec_r = dec
    rr = cute.CuteAngularCorrelator(b, [ra_r, dec_r])
    rr.run()

    xis = []
    for i in range(10):
        # Data catalog
        mp = hp.synfast(cl, nside)
        nmap = np.random.poisson(1+mp)
        coords = [np.repeat(ra, nmap), np.repeat(dec, nmap)]
        # DD
        dd = cute.CuteAngularCorrelator(b, coords)
        dd.run()
        # DR
        dr = cute.CuteAngularCorrelator(b, coords, coords2=[ra_r, dec_r])
        dr.run()
        # L&S
        xi = dd.get_xi(dr=dr, rr=rr, kind='tracer')[1]
        xis.append(xi)
    xis = np.array(xis)
    xi_mean = np.mean(xis, axis=0)
    xi_std = np.std(xis, axis=0)

    th = b.get_r_values()
    cth = np.cos(np.radians(th))
    xi_th = np.array([np.sum(cl_pixwin*(2*ls+1)*eval_legendre(ls, c))
                     for c in cth])/(4*np.pi)

    # Uncomment if you'd like to plot the results
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.errorbar(th[1:], xi_mean[1:], yerr=xi_std[1:], fmt='k.')
    # plt.plot(th[1:], xi_th[1:], 'r-')
    # for xi in xis:
    #     plt.plot(th[1:], xi[1:], 'b-', alpha=0.1)
    # plt.show()
    assert np.all(np.fabs(xi_th-xi_mean)[1:] < 2*xi_std[1:])


def test_angcorr_sample():
    nside = 128
    npix = hp.nside2npix(nside)
    ls = np.arange(3*nside).astype(int)
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
    xi_th = np.array([np.sum(cl_pixwin*(2*ls+1)*eval_legendre(ls, c))
                     for c in cth])/(4*np.pi)

    # Uncomment if you'd like to plot the results
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.errorbar(th, xiA_mean, yerr=xiA_std, fmt='k.')
    # plt.plot(th, xi_th, 'r-')
    # for xi in xis_A:
    #     plt.plot(th, xi, 'b-', alpha=0.1)
    # plt.figure()
    # plt.errorbar(th, xiX_mean, yerr=xiX_std, fmt='k.')
    # plt.plot(th, xi_th, 'r-')
    # for xi in xis_X:
    #     plt.plot(th, xi, 'b-', alpha=0.1)
    # plt.show()
    assert np.all(np.fabs(xi_th-xiA_mean) < 2*xiA_std)
    assert np.all(np.fabs(xi_th-xiX_mean) < 2*xiX_std)
