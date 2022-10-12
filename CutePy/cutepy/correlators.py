import numpy as np
from cutepy import cutelib as lib


def _get_C_correlator(kind):
    if kind == 'xi_r':
        return lib.xi_r_Acorr_bf_C, lib.xi_r_Xcorr_bf_C
    if kind == 'xi_theta':
        return lib.xi_th_Acorr_bf_C, lib.xi_th_Xcorr_bf_C
    else:
        raise ValueError(f"Unknown correlation kind {kind}.")


class CuteCorrelator(object):
    kind = 'xi_r'

    def __init__(self, bins, coords, weights=None,
                 coords2=None, weights2=None):
        self._check_bins(bins)
        self.bins = bins
        self.coords = coords
        self.weights = weights
        self.coords2 = coords2
        self.weights2 = weights2
        self.ww = None
        self.ct = None
        self.ww_norm = self._get_ww_norm()
        self._cA, self._cX = _get_C_correlator(self.kind)

    def _get_ww_norm(self):
        n1 = len(self.coords[0])
        if self.weights is None:
            sum_w1 = n1
            sum_w1_2 = n1
        else:
            sum_w1 = np.sum(self.weights)
            sum_w1_2 = np.sum(self.weights**2)
        if self.coords2 is None:
            return 0.5*(sum_w1**2-sum_w1_2)
        else:
            n2 = len(self.coords2[0])
            if self.weights2 is None:
                sum_w2 = n2
            else:
                sum_w2 = np.sum(self.weights2)
            return sum_w1*sum_w2

    def _check_bins(self, bins):
        pass

    def _cute_driver(self):
        if self.weights is None:
            w = np.ones(len(self.coords[0]))
        else:
            w = self.weights

        if self.coords2 is None:  # Auto-correlation
            d = self._cA(self.bins._bin, True,
                         self.coords[0], self.coords[1], self.coords[2],
                         w, 2*self.bins.nbins)
        else:  # Cross-correlation
            if self.weights2 is None:
                w2 = np.ones(len(self.coords2[0]))
            else:
                w2 = self.weights2
            d = self._cX(self.bins._bin, True,
                         self.coords[0], self.coords[1], self.coords[2], w,
                         self.coords2[0], self.coords2[1], self.coords2[2], w2,
                         2*self.bins.nbins)
        self.ww = d[:self.bins.nbins]
        self.ct = d[self.bins.nbins:]

    def run(self):
        if self.ww is not None:
            return
        self._cute_driver()

    def get_xi(self, kind='sample', dr=None, rr=None, rd=None):
        # Check other correlators use same bins
        for c in [dr, rr, rd]:
            if (c is not None) and (not (c.bins == self.bins)):
                raise ValueError("Other correlator uses different bins")

        if kind == 'sample':  # Mean of weights
            goodbins = self.ct > 0
            xi = np.zeros(self.bins.nbins)
            xi[goodbins] = self.ww[goodbins]/self.ct[goodbins]
            if dr is not None:
                goodbins = dr.ct > 0
                xi_r = np.zeros(self.bins.nbins)
                xi_r[goodbins] = dr.ww[goodbins]/dr.ct[goodbins]
                xi = xi - xi_r
        elif kind == 'tracer':  # Landy-Szalay
            if (rr is None) or (dr is None):
                raise ValueError("Need RR and DR for L&S estimzator")
            xi = np.zeros(self.bins.nbins)
            ddh = self.ww/self.ww_norm
            rrh = rr.ww/rr.ww_norm
            drh = dr.ww/dr.ww_norm
            goodbins = rrh > 0
            if rd is not None:
                rdh = rd.ww/rd.ww_norm
            else:
                rdh = drh
            # Landy-Szalay estimator
            xi[goodbins] = (ddh-drh-rdh+rrh)[goodbins]/rrh[goodbins]
        else:
            raise ValueError(f"Unknown correlation type {kind}.")
        r = self.bins.get_r_values()
        return r, xi


class CuteAngularCorrelator(CuteCorrelator):
    kind = 'xi_theta'

    def _cute_driver(self):
        if self.weights is None:
            w = np.ones(len(self.coords[0]))
        else:
            w = self.weights

        if self.coords2 is None:  # Auto-correlation
            d = self._cA(self.bins._bin, True,
                         np.cos(np.radians(90-self.coords[1])),
                         np.radians(self.coords[0]),
                         w, 2*self.bins.nbins)
        else:  # Cross-correlation
            if self.weights2 is None:
                w2 = np.ones(len(self.coords2[0]))
            else:
                w2 = self.weights2
            d = self._cX(self.bins._bin, True,
                         np.cos(np.radians(90-self.coords[1])),
                         np.radians(self.coords[0]), w,
                         np.cos(np.radians(90-self.coords2[1])),
                         np.radians(self.coords2[0]), w2,
                         2*self.bins.nbins)
        self.ww = d[:self.bins.nbins]
        self.ct = d[self.bins.nbins:]

    def _check_bins(self, bins):
        if not bins.is_deg:
            raise ValueError("Binning must be for an angle variable.")
