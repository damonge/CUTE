# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.1
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    import _cutelib
else:
    import _cutelib

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


DTORAD = _cutelib.DTORAD
class CuteBin(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    nbins = property(_cutelib.CuteBin_nbins_get, _cutelib.CuteBin_nbins_set)
    is_log = property(_cutelib.CuteBin_is_log_get, _cutelib.CuteBin_is_log_set)
    n_logint = property(_cutelib.CuteBin_n_logint_get, _cutelib.CuteBin_n_logint_set)
    r_max = property(_cutelib.CuteBin_r_max_get, _cutelib.CuteBin_r_max_set)
    i_r_max = property(_cutelib.CuteBin_i_r_max_get, _cutelib.CuteBin_i_r_max_set)
    log_r_max = property(_cutelib.CuteBin_log_r_max_get, _cutelib.CuteBin_log_r_max_set)
    is_2D = property(_cutelib.CuteBin_is_2D_get, _cutelib.CuteBin_is_2D_set)
    nbins2 = property(_cutelib.CuteBin_nbins2_get, _cutelib.CuteBin_nbins2_set)
    is_mu2 = property(_cutelib.CuteBin_is_mu2_get, _cutelib.CuteBin_is_mu2_set)
    r_max2 = property(_cutelib.CuteBin_r_max2_get, _cutelib.CuteBin_r_max2_set)
    i_r_max2 = property(_cutelib.CuteBin_i_r_max2_get, _cutelib.CuteBin_i_r_max2_set)
    nbins_total = property(_cutelib.CuteBin_nbins_total_get, _cutelib.CuteBin_nbins_total_set)

    def __init__(self):
        _cutelib.CuteBin_swiginit(self, _cutelib.new_CuteBin())
    __swig_destroy__ = _cutelib.delete_CuteBin

# Register CuteBin in _cutelib:
_cutelib.CuteBin_swigregister(CuteBin)


def cute_bin_new(nr, is_log, rmax, rmin_log, unit):
    return _cutelib.cute_bin_new(nr, is_log, rmax, rmin_log, unit)

def cute_bin_2d_new(nr, is_log, rmax, rmin_log, nr2, is_mu2, rmax2):
    return _cutelib.cute_bin_2d_new(nr, is_log, rmax, rmin_log, nr2, is_mu2, rmax2)

def cute_bin_free(bin):
    return _cutelib.cute_bin_free(bin)

def cute_get_rmax(bin):
    return _cutelib.cute_get_rmax(bin)

def cute_angular_corr_bf(bin, get_counts, DD, counts, n1, cth1, phi1, w1, n2, cth2, phi2, w2, NodeThis, NNodes):
    return _cutelib.cute_angular_corr_bf(bin, get_counts, DD, counts, n1, cth1, phi1, w1, n2, cth2, phi2, w2, NodeThis, NNodes)

def cute_xi_r_corr_bf(bin, get_counts, DD, counts, n1, x1, y1, z1, w1, n2, x2, y2, z2, w2, NodeThis, NNodes):
    return _cutelib.cute_xi_r_corr_bf(bin, get_counts, DD, counts, n1, x1, y1, z1, w1, n2, x2, y2, z2, w2, NodeThis, NNodes)

def get_bins_C(nbins, islog, rmax, rmin_log, isdeg):
    return _cutelib.get_bins_C(nbins, islog, rmax, rmin_log, isdeg)

def get_bins_2d_C(nbins1, islog1, rmax1, rmin_log1, nbins2, rmax2, ismu2):
    return _cutelib.get_bins_2d_C(nbins1, islog1, rmax1, rmin_log1, nbins2, rmax2, ismu2)

def xi_th_Xcorr_bf_C(bin, get_counts, n1c, n1p, n1w, n2c, n2p, n2w, NodeThis, NNodes, dout):
    return _cutelib.xi_th_Xcorr_bf_C(bin, get_counts, n1c, n1p, n1w, n2c, n2p, n2w, NodeThis, NNodes, dout)

def xi_th_Acorr_bf_C(bin, get_counts, n1c, n1p, n1w, NodeThis, NNodes, dout):
    return _cutelib.xi_th_Acorr_bf_C(bin, get_counts, n1c, n1p, n1w, NodeThis, NNodes, dout)

def xi_r_Xcorr_bf_C(bin, get_counts, n1x, n1y, n1z, n1w, n2x, n2y, n2z, n2w, NodeThis, NNodes, dout):
    return _cutelib.xi_r_Xcorr_bf_C(bin, get_counts, n1x, n1y, n1z, n1w, n2x, n2y, n2z, n2w, NodeThis, NNodes, dout)

def xi_r_Acorr_bf_C(bin, get_counts, n1x, n1y, n1z, n1w, NodeThis, NNodes, dout):
    return _cutelib.xi_r_Acorr_bf_C(bin, get_counts, n1x, n1y, n1z, n1w, NodeThis, NNodes, dout)


