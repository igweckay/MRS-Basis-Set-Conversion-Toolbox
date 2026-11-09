"""
core/write.py
Basis Set Converter — writing a basis set

Writes basis functions with the writer for the tool chosen, as core.formats lists the
tools and their output formats.
Used by the GUI (screen4_convert.py) and by scripts, without a GUI.
"""

import os


def write_basis(basis_list, tool_id, out_fmt, outdir, te=None, seq=None):
    """
    Write basis functions using the appropriate writer.

    Parameters
    ----------
    basis_list : list of core struct dicts
    tool_id    : str, a tool's 'id' in core.formats.ALL_TOOLS
    out_fmt    : str, one of that tool's format ids, or 'default'
    outdir     : str, the output folder
    te         : float or None, echo time in ms (30 if None; 26 for OXSA)
    seq        : str or None, sequence name ('unedited' if None)

    Returns
    -------
    list of str, the paths written
    """
    stated_te = te
    te  = te or 30.0
    seq = seq or 'unedited'

    if tool_id == 'lcmodel':
        from writers.write_lcmodel import write_lcmodel
        result = write_lcmodel(basis_list, outdir,
                               mode=out_fmt if out_fmt in ('raw','basis','both') else 'both',
                               te=te, seq=seq)
        paths = result.get('raw_files', []) + \
                ([result['basis_file']] if 'basis_file' in result else [])
        return paths

    if tool_id in ('spant', 'abfit'):
        # SPANT uses LCModel .raw format
        from writers.write_lcmodel import write_lcmodel
        result = write_lcmodel(basis_list, outdir,
                               mode='raw' if out_fmt != 'nifti' else 'raw',
                               te=te, seq=seq)
        if out_fmt == 'nifti':
            from writers.write_niftimrs import write_niftimrs
            return write_niftimrs(basis_list, outdir)
        return result.get('raw_files', [])

    if tool_id == 'tarquin':
        from writers.write_lcmodel import write_lcmodel
        result = write_lcmodel(basis_list, outdir, mode='basis', te=te, seq=seq)
        return [result.get('basis_file', '')]

    if tool_id == 'osprey':
        from writers.write_osprey import write_osprey
        out_file = os.path.join(outdir, 'basis.mat')
        if out_fmt == 'nifti':
            from writers.write_niftimrs import write_niftimrs
            return write_niftimrs(basis_list, outdir)
        if out_fmt == 'basis':
            from writers.write_lcmodel import write_lcmodel
            result = write_lcmodel(basis_list, outdir, mode='basis', te=te, seq=seq)
            return [result.get('basis_file', '')]
        write_osprey(basis_list, out_file, te=float(te or 30))
        return [out_file]

    if tool_id == 'fsLmrs':
        from writers.write_fsLmrs import write_fsLmrs_folder
        return write_fsLmrs_folder(basis_list, outdir)

    if tool_id == 'inspector':
        from writers.write_inspector import write_inspector
        out_file = os.path.join(outdir, 'basis.mat')
        write_inspector(basis_list, out_file)
        return [out_file]

    if tool_id in ('jmrui', 'quest', 'nmbscope', 'aqses'):
        from writers.write_jmrui import write_jmrui_folder
        return write_jmrui_folder(basis_list, outdir)

    if tool_id == 'niftimrs':
        from writers.write_niftimrs import write_niftimrs
        return write_niftimrs(basis_list, outdir)

    if tool_id == 'marss':
        from writers.write_marss import (write_marss_folder,
                                          write_marss_combined)
        from writers.write_lcmodel import write_lcmodel_raw_folder
        if out_fmt == 'mat_combined':
            out_file = os.path.join(outdir, 'basis_combined.mat')
            write_marss_combined(basis_list, out_file)
            return [out_file]
        if out_fmt == 'raw':
            return write_lcmodel_raw_folder(basis_list, outdir)
        return write_marss_folder(basis_list, outdir)

    if tool_id == 'mrscloud':
        from writers.write_mrscloud import write_mrscloud_folder
        return write_mrscloud_folder(basis_list, outdir)

    if tool_id in ('fida', 'fida_fit', 'spinach'):
        from writers.write_fida import write_fida_folder
        return write_fida_folder(basis_list, outdir)

    if tool_id == 'pyamares':
        from writers.write_pyamares import write_pyamares
        out_file = os.path.join(outdir, 'pyamares_priors.csv')
        write_pyamares(basis_list, out_file)
        return [out_file]

    if tool_id == 'oxsa':
        from writers.write_oxsa import write_oxsa
        out_file = os.path.join(outdir, 'run_oxsa_amares.m')
        write_oxsa(basis_list, out_file, seqte=float(stated_te or 26.0))
        return [out_file]

    if tool_id == 'midas':
        from writers.write_midas import write_midas
        out_file = os.path.join(outdir, 'midas_priors.xml')
        write_midas(basis_list, out_file)
        return [out_file]

    if tool_id == 'gava':
        from writers.write_gava import write_gava
        out_file = os.path.join(outdir, 'gava_priors.txt')
        write_gava(basis_list, out_file)
        return [out_file]

    if tool_id in ('vespa_analysis', 'vespa_gen2'):
        from writers.write_vespa import write_vespa
        out_file = os.path.join(outdir, 'vespa_priors.xml')
        write_vespa(basis_list, out_file)
        return [out_file]

    if tool_id in ('jet', 'spinwizard'):
        from writers.write_spinwizard import write_spinwizard
        return write_spinwizard(basis_list, outdir)

    if tool_id == 'profit':
        from writers.write_profit import write_profit
        out_file = os.path.join(outdir, 'profit_struct.mat')
        write_profit(basis_list, out_file)
        return [out_file]

    raise RuntimeError(f"No writer available for tool: {tool_id}")
