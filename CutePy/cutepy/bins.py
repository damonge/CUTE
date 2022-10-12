import numpy as np
from cutepy import cutelib as lib


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
