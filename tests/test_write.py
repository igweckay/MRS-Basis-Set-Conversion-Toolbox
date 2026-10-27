"""
tests/test_write.py
Converting a basis set without the GUI: core.write and core.load.fill_missing, as
screen4_convert runs them.

Run from the repository root: python -m pytest tests
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAMPLES = os.path.join(ROOT, 'sample_basis_sets')
sys.path.insert(0, os.path.join(ROOT, 'basis_converter'))

from core.detect import detect_format                   # noqa: E402
from core.load import fill_missing, load_basis          # noqa: E402
from core.write import write_basis                      # noqa: E402


def _jmrui():
    return load_basis([os.path.join(SAMPLES, 'JMRUI', name)
                       for name in ('NAA.txt', 'Cr391.txt')], 'jMRUI .txt')


def test_a_basis_is_written_for_the_tool_chosen_without_the_gui(tmp_path):
    written = write_basis(_jmrui(), 'fsLmrs', 'default', str(tmp_path))

    assert sorted(os.path.basename(p) for p in written) == ['Cr391.json', 'NAA.json']
    assert detect_format(written)['format'] == 'FSL-MRS .json'
    assert 'tkinter' not in sys.modules


def test_parameters_the_files_lack_are_filled_in_and_kept_where_stated(tmp_path):
    files = [os.path.join(SAMPLES, 'FSL-MRS', name) for name in ('NAA.json', 'Cr391.json')]
    basis = load_basis(files, detect_format(files)['format'])
    sw = basis[0]['sw']

    fill_missing(basis, {'sf': 123.2, 'sw': 1.0, 'te': 30})            # sf, sw stated: kept
    written = write_basis(basis, 'lcmodel', 'basis', str(tmp_path))

    assert [b['sf'] for b in basis] == [123.261803, 123.261803] and basis[0]['sw'] == sw
    assert [b['te'] for b in basis] == [30, 30]
    assert detect_format(written)['data']['sf'] == 123.261803


def test_raw_files_without_their_parameters_take_those_given(tmp_path):
    folder = os.path.join(SAMPLES, 'LCModel', 'raw_files')
    files = sorted(os.path.join(folder, f) for f in os.listdir(folder) if f.endswith('.raw'))
    basis = load_basis(files, detect_format(files)['format'], sw=4000.0, sf=123.2)

    written = write_basis(basis, 'lcmodel', 'basis', str(tmp_path))

    assert detect_format(written)['data']['sw'] == 4000


def test_fsl_mrs_fields_are_the_frequency_and_the_linewidth(tmp_path):
    """FSL-MRS's fsl_io reads basis_centre as the spectrometer frequency in MHz and
    basis_width as the linewidth; the ppm reference is its nucleus' default, not stored."""
    import json
    from readers.read_fsLmrs import read_fsLmrs_file
    from writers.write_fsLmrs import write_fsLmrs_file

    asc = read_fsLmrs_file(os.path.join(SAMPLES, 'FSL-MRS', 'Asc.json'))
    assert asc['sf'] == 123.261803 and asc['centerFreq'] is None and asc['linewidth'] is None

    out = str(tmp_path / 'Asc.json')
    write_fsLmrs_file(dict(asc, linewidth=2.5, centerFreq=3.0), out)
    written = json.load(open(out))['basis']
    assert written['basis_centre'] == 123.261803 and written['basis_width'] == 2.5
