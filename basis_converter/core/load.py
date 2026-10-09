"""
core/load.py
Basis Set Converter — reading a basis set

Reads basis functions with the reader for their format, as core.detect names it.
Used by the GUI (screen4_convert.py) and by scripts, without a GUI.
"""

import os


def load_basis(files, fmt, sw=None, sf=None):
    """
    Load basis functions using the appropriate reader.

    Parameters
    ----------
    files : list of str, the basis files
    fmt   : str, the format as core.detect.detect_format names it
    sw    : float or None, spectral width in Hz (SpinWizard files do not hold it)
    sf    : float or None, spectrometer frequency in MHz (likewise)

    Returns
    -------
    list of core struct dicts, one per metabolite
    """
    if 'MARSS .mat' in fmt:
        from readers.read_marss import read_marss_file, read_marss_folder
        if len(files) == 1:
            return [read_marss_file(files[0])]
        return [read_marss_file(f) for f in files
                if f.lower().endswith('.mat')]

    if 'LCModel .BASIS' in fmt:
        from readers.read_lcmodel import read_lcmodel_basis
        return read_lcmodel_basis(files[0])

    if 'LCModel .raw' in fmt or 'MARSS .raw' in fmt:
        from readers.read_lcmodel import read_lcmodel_raw
        return [read_lcmodel_raw(f) for f in files
                if f.lower().endswith('.raw')]

    if 'Osprey' in fmt:
        from readers.read_osprey import read_osprey
        return read_osprey(files[0])

    if 'FSL-MRS' in fmt:
        from readers.read_fsLmrs import read_fsLmrs_file
        return [read_fsLmrs_file(f) for f in files
                if f.lower().endswith('.json')]

    if 'INSPECTOR' in fmt:
        from readers.read_inspector import read_inspector
        return read_inspector(files[0])

    if 'jMRUI' in fmt:
        from readers.read_jmrui import read_jmrui_file
        return [read_jmrui_file(f) for f in files
                if f.lower().endswith('.txt')]

    if 'NIfTI-MRS' in fmt:
        from readers.read_niftimrs import read_niftimrs
        results = []
        for f in files:
            if f.lower().endswith('.nii.gz'):
                results.extend(read_niftimrs(f))
        return results

    if 'ProFit' in fmt:
        from readers.read_profit import read_profit
        return read_profit(files[0])

    if 'MRSCloud' in fmt:
        from readers.read_mrscloud import read_mrscloud_file
        return [read_mrscloud_file(f) for f in files
                if f.lower().endswith('.mat')]

    if 'FID-A' in fmt:
        from readers.read_fida import read_fida_file
        return [read_fida_file(f) for f in files
                if f.lower().endswith('.mat')]

    if 'PyAMARES' in fmt:
        from readers.read_pyamares import read_pyamares
        # PyAMARES CSV doesn't contain FIDs — can only convert to OXSA
        raise RuntimeError(
            "PyAMARES CSV contains prior knowledge only, not FIDs. "
            "Cannot convert to FID-based formats."
        )

    if 'SpinWizard' in fmt:
        from readers.read_spinwizard import read_spinwizard_file, read_spinwizard_folder
        if len(files) == 1:
            return [read_spinwizard_file(files[0], sw=sw, sf=sf)]
        return read_spinwizard_folder(os.path.dirname(files[0]), sw=sw, sf=sf)

    if 'MIDAS' in fmt:
        from readers.read_midas import read_midas
        return read_midas(files[0])

    if 'GAVA' in fmt:
        from readers.read_gava import read_gava
        return read_gava(files[0])

    raise RuntimeError(f"No reader available for format: {fmt}")
