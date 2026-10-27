"""
writers/centre.py
Basis Set Converter — the ppm centre a basis is written for

Formats that store no ppm reference are read at a fixed one: 4.65 ppm for 1H, LCModel's
PPMSEP, FSL-MRS's nucleus default (for its .json and NIfTI-MRS), the centre jMRUI is read
at and the ppm axis the FID-A and MRSCloud writers store. A basis function stating another
(centerFreq, e.g. 3.0 ppm in Osprey's own basis sets) is shifted by the difference before
it is written, so its peaks stay at their ppm.
Used by write_lcmodel, write_niftimrs, write_fsLmrs, write_jmrui, write_fida and write_mrscloud,
and, through at_centre, by the writers of prior knowledge (MIDAS, GAVA, VeSPA, PyAMARES, OXSA).
"""

import numpy as np

# the ppm at 0 Hz of a format that states none
CENTRE = 4.65


def fid_at_centre(core, centre=CENTRE):
    """
    The FID referenced to `centre`: shifted by the difference from the centre the basis
    function states (centerFreq); one stating none is taken as it is.
    """
    fid    = np.asarray(core['fid']).ravel()
    stated = core.get('centerFreq')
    if stated is None or float(stated) == centre:
        return fid
    t = np.arange(len(fid)) / float(core['sw'])
    return fid * np.exp(2j * np.pi * (float(stated) - centre) * float(core['sf']) * t)


def at_centre(basis_list, centre=CENTRE):
    """
    The basis set referenced to `centre`, for writers that read peaks off its spectrum: each
    function stating another centre shifted (fid_at_centre) and then stating none. A set
    stating none, or only a ppmCalib, is left as it is.
    """
    return [core if core.get('centerFreq') is None else
            dict(core, fid=fid_at_centre(core, centre), centerFreq=None) for core in basis_list]
