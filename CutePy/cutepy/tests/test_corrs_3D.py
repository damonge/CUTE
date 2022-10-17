import numpy as np
import cutepy as cute
from scipy.integrate import quad
from scipy.special import spherical_jn


class GauSim3D(object):
    def __init__(self, lbox, ng, sigma=0.2, gamma=2.5, ncell_smooth=2):
        self.k_pivot = 2.0/lbox
        self.p0 = 1.0
        self.gamma = gamma
        self.lbox = lbox
        self.vol = self.lbox**3
        self.ng = ng
        self.dk = 2*np.pi/self.lbox
        self.dx = self.lbox/self.ng
        self.rs = self.dx*ncell_smooth
        self.kmax = self.ng*self.dk

        sig2 = quad(self._sig2_integ, np.log(self.dk),
                    np.log(self.kmax))[0]
        self.p0 = sigma**2/sig2
        sig2 = quad(self._sig2_integ, np.log(self.dk),
                    np.log(self.kmax))[0]

        kv = self.ng*self.dk*np.fft.fftfreq(self.ng)
        k3d = np.sqrt(kv[:, None, None]**2 +
                      kv[None, :, None]**2 +
                      kv[None, None, :self.ng//2+1]**2)
        self.sigk = np.sqrt(0.5*self.pk(k3d)*self.vol)*self.ng**3/self.vol

    def _sig2_integ(self, lk):
        k = np.exp(lk)
        return k**3*self.pk(k)/(2*np.pi**2)

    def pk(self, k):
        x = k/self.k_pivot
        bm = np.exp(-0.5*(k*self.rs))
        return self.p0*bm**2*x/(1+x**(self.gamma+1))

    def gen_sim(self):
        dre, dim = np.random.randn(2, self.ng, self.ng,
                                   self.ng//2+1)*self.sigk
        return np.fft.irfftn(dre + 1j*dim)

    def xyz2ijk(self, x, y, z):
        return [np.floor(x/self.dx).astype(int),
                np.floor(y/self.dx).astype(int),
                np.floor(z/self.dx).astype(int)]


def test_mono_sample():
    ng = 64
    lbox = 1.
    gs = GauSim3D(lbox=lbox, ng=ng)

    def pk2xi(lk, r):
        k = np.exp(lk)
        x = k*r
        j0 = spherical_jn(0, x)
        return k**3*gs.pk(k)*j0/(2*np.pi**2)

    def get_xi(r):
        return quad(lambda lk: pk2xi(lk, r), np.log(gs.dk), np.log(gs.kmax))[0]

    b = cute.CuteBin(20, 0.2)
    xis_A = []
    xis_X = []
    for i in range(10):
        mp = gs.gen_sim().flatten()
        x, y, z = gs.lbox*np.random.rand(3, ng**3//2)
        ix, iy, iz = gs.xyz2ijk(x, y, z)
        ibox = iz+ng*(iy+ng*ix)
        mp = mp[ibox]
        # Auto-correlation
        corr = cute.CuteCorrelator(b, [x, y, z], weights=mp)
        corr.run()
        xi = corr.get_xi()[1]
        xis_A.append(xi)
        # Cross-correlation
        corr = cute.CuteCorrelator(b, [x, y, z], weights=mp,
                                   coords2=[x, y, z], weights2=mp)
        corr.run()
        xi = corr.get_xi()[1]
        xis_X.append(xi)
    xis_A = np.array(xis_A)
    xiA_mean = np.mean(xis_A, axis=0)
    xiA_std = np.std(xis_A, axis=0)
    xis_X = np.array(xis_X)
    xiX_mean = np.mean(xis_X, axis=0)
    xiX_std = np.std(xis_X, axis=0)

    r = b.get_r_values()
    xi_th = np.array([get_xi(rr) for rr in r])

    # Uncomment if you'd like to plot the results
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.errorbar(r, xiA_mean, yerr=xiA_std, fmt='k.')
    # plt.plot(r, xi_th, 'r-')
    # for xi in xis_A:
    #     plt.plot(r, xi, 'b-', alpha=0.1)
    # plt.figure()
    # plt.errorbar(r, xiX_mean, yerr=xiX_std, fmt='k.')
    # plt.plot(r, xi_th, 'r-')
    # for xi in xis_X:
    #     plt.plot(r, xi, 'b-', alpha=0.1)
    # plt.show()
    assert np.all(np.fabs(xi_th-xiA_mean) < 2*xiA_std)
    assert np.all(np.fabs(xi_th-xiX_mean) < 2*xiX_std)


def test_2D_rmu_sample():
    ng = 64
    lbox = 1.
    gs = GauSim3D(lbox=lbox, ng=ng)

    mp = gs.gen_sim().flatten()
    x, y, z = gs.lbox*np.random.rand(3, ng**3//2)
    ix, iy, iz = gs.xyz2ijk(x, y, z)
    ibox = iz+ng*(iy+ng*ix)
    mp = mp[ibox]

    b_2D = cute.CuteBin2D(20, 0.2, 10, is_mu2=True)
    corr_2D = cute.CuteCorrelator2D(b_2D, [x, y, z], weights=mp)
    corr_2D.run()
    r, _, xi_2D = corr_2D.get_xi()

    b_mono = cute.CuteBin(20, 0.2)
    corr_mono = cute.CuteCorrelator(b_mono, [x, y, z], weights=mp)
    corr_mono.run()
    xi_mono = corr_mono.get_xi()[1]

    # Check that mean over mu is similar to monopole
    xi_mean = np.mean(xi_2D, axis=1)
    xi_std = np.std(xi_2D, axis=1)
    assert np.all(np.fabs(xi_mono-xi_mean) < 3*xi_std)

    # Check that we recover monopole by summing weights and counts over mu
    wwsum = np.sum(corr_2D.ww.reshape([b_2D.nbins1, b_2D.nbins2]), axis=1)
    ctsum = np.sum(corr_2D.ct.reshape([b_2D.nbins1, b_2D.nbins2]), axis=1)
    xi_mono_B = wwsum/ctsum
    assert np.all(np.fabs(xi_mono_B/xi_mono-1) < 1E-5)


def test_2D_ps_sample():
    ng = 64
    lbox = 1.
    gs = GauSim3D(lbox=lbox, ng=ng)

    mp = gs.gen_sim().flatten()
    x, y, z = gs.lbox*np.random.rand(3, ng**3//2)
    ix, iy, iz = gs.xyz2ijk(x, y, z)
    ibox = iz+ng*(iy+ng*ix)
    mp = mp[ibox]

    b_2D = cute.CuteBin2D(20, 0.2, 20, rmax2=0.2)
    corr_2D = cute.CuteCorrelator2D(b_2D, [x, y, z], weights=mp)
    corr_2D.run()
    rt, rp, xi_2D = corr_2D.get_xi()

    b_mono = cute.CuteBin(10, 0.2)
    corr_mono = cute.CuteCorrelator(b_mono, [x, y, z], weights=mp)
    corr_mono.run()
    r, xi_mono = corr_mono.get_xi()

    r2d = np.sqrt(rt[:, None]**2+rp[None, :]**2)
    sumxi = np.histogram(r2d.flatten(), bins=10, range=[0, 0.2],
                         weights=xi_2D.flatten())[0]
    sumxi2 = np.histogram(r2d.flatten(), bins=10, range=[0, 0.2],
                          weights=xi_2D.flatten()**2)[0]
    sumct = np.histogram(r2d.flatten(), bins=10, range=[0, 0.2])[0]
    xi_mean = sumxi/sumct
    xi_std = np.sqrt(sumxi2/sumct-xi_mean**2)

    # Check that averaging over the same r is similar to monopole
    assert np.all(np.fabs(xi_mean-xi_mono) < 2*xi_std)
