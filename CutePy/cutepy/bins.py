import numpy as np
from cutepy import cutelib as lib


class CuteBin2D(object):
    def __init__(self, nbins1, rmax1, nbins2, rmax2=None, is_log1=False,
                 rmin_log1=None, is_mu2=False):
        self._bin = None

        if is_log1:
            if (rmin_log1 is None):
                raise ValueError("I need an `rmin_log1` if using log binning")
            if (rmin_log1 >= rmax1):
                raise ValueError("`rmin_log1` >= `rmax1`!")
        if rmin_log1 is None:
            rmin_log1 = 0.
        if is_mu2:
            rmax2 = 1.
        else:
            if rmax2 is None:
                raise ValueError("I need `rmax2` if not binning in mu.")

        self._bin = lib.get_bins_2d_C(nbins1, int(is_log1), rmax1,  rmin_log1,
                                      nbins2, rmax2, int(is_mu2))

        self.nbins1 = nbins1
        self.nbins2 = nbins2
        self.nbins_total = nbins1*nbins2
        self.is_log1 = is_log1
        self.is_mu2 = is_mu2
        self.r_max1 = self._bin.r_max
        self.log_r_max1 = self._bin.log_r_max
        self.n_logint1 = self._bin.n_logint
        self.r_max2 = self._bin.r_max2

    def get_r_values(self):
        i1_s = np.arange(self.nbins1)+0.5
        i2_s = np.arange(self.nbins2)+0.5

        if self.is_log1:
            lrs = (i1_s-self.nbins1)/self.n_logint1+self.log_r_max1
            r1 = 10.**lrs
        else:
            r1 = i1_s*self.r_max1/self.nbins1
        if self.is_mu2:
            x2 = i2_s/self.nbins2
        else:
            x2 = i2_s*self.r_max2/self.nbins2
        return r1, x2

    def __del__(self):
        if self._bin is not None:
            if lib.cute_bin_free is not None:
                lib.cute_bin_free(self._bin)
            self._bin = None

    def __eq__(self, other):
        if isinstance(other, CuteBin2D):
            return ((self.nbins1 == other.nbins1) and
                    (self.r_max1 == other.r_max1) and
                    (self.is_log1 == other.is_log1) and
                    (self.n_logint1 == other.n_logint1) and
                    (self.nbins2 == other.nbins2) and
                    (self.r_max2 == other.r_max2) and
                    (self.is_mu2 == other.is_mu2))
        return False


class CuteBin(object):
    def __init__(self, nbins, rmax, is_log=False, rmin_log=None, is_deg=False):
        self._bin = None

        if is_log:
            if (rmin_log is None):
                raise ValueError("I need an `rmin_log` if using log binning")
            if (rmin_log >= rmax):
                raise ValueError("`rmin_log` >= `rmax`!")
        if rmin_log is None:
            rmin_log = 0.

        self._bin = lib.get_bins_C(nbins, int(is_log), rmax,  rmin_log,
                                   int(is_deg))
        self.nbins = nbins
        self.is_log = is_log
        self.is_deg = is_deg
        self.r_max = self._bin.r_max
        self.log_r_max = self._bin.log_r_max
        self.n_logint = self._bin.n_logint

    def get_r_values(self):
        i_s = np.arange(self.nbins)+0.5

        if self.is_log:
            lrs = (i_s-self.nbins)/self.n_logint+self.log_r_max
            rs = 10.**lrs
        else:
            rs = i_s*self.r_max/self.nbins
        if self.is_deg:
            rs = np.degrees(rs)
        return rs

    def __del__(self):
        if self._bin is not None:
            if lib.cute_bin_free is not None:
                lib.cute_bin_free(self._bin)
            self._bin = None

    def __eq__(self, other):
        if isinstance(other, CuteBin):
            return ((self.nbins == other.nbins) and
                    (self.r_max == other.r_max) and
                    (self.is_log == other.is_log) and
                    (self.n_logint == other.n_logint) and
                    (self.is_deg == other.is_deg))
        return False
