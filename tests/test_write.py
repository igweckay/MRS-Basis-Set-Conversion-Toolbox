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


def test_mrscloud_files_are_read_as_mrscloud_saves_them():
    """MRSCloud saves each metabolite as the one struct outA, not as top-level fields."""
    folder = os.path.join(SAMPLES, 'MRSCLOUD')
    files = sorted(os.path.join(folder, f) for f in os.listdir(folder)
                   if f.startswith('Siemens_') and f.endswith('.mat'))
    found = detect_format(files)
    basis = load_basis(files, found['format'])

    assert found['format'] == 'MRSCloud .mat' and detect_format(files[:1])['format'] == 'MRSCloud .mat'
    assert len(basis) == len(files) == 25
    assert (basis[0]['name'], basis[0]['sw'], basis[0]['sf']) == ('Ace', 4000.0, 123.04753)


def test_a_basis_stating_no_linewidth_is_written_with_the_default(tmp_path):
    """FSL-MRS states the linewidth as null when it has none; every writer takes its default."""
    files = [os.path.join(SAMPLES, 'FSL-MRS', name) for name in ('NAA.json', 'Cr391.json')]
    basis = load_basis(files, detect_format(files)['format'])
    assert [b['linewidth'] for b in basis] == [None, None]

    for tool in ('inspector', 'fida', 'mrscloud'):
        outdir = tmp_path / tool
        outdir.mkdir()
        written = write_basis(basis, tool, 'default', str(outdir))
        found = detect_format(written)

        assert found['status'] == 'full', tool
        assert len(load_basis(written, found['format'])) == 2, tool


def test_every_format_read_at_4_65_ppm_keeps_a_stated_centre(tmp_path):
    """FSL-MRS, jMRUI, FID-A and MRSCloud state no ppm centre and are read at 4.65, as LCModel is:
    a basis put at 3.0 ppm (Osprey's own sets) is shifted, so NAA stays at 2.01 ppm."""
    import numpy as np
    naa = load_basis([os.path.join(SAMPLES, 'JMRUI', 'NAA.txt')], 'jMRUI .txt')[0]
    t = np.arange(len(naa['fid'])) / naa['sw']
    moved = dict(naa, fid=naa['fid'] * np.exp(2j * np.pi * 1.65 * naa['sf'] * t), centerFreq=3.0)

    for tool in ('fsLmrs', 'jmrui', 'fida', 'mrscloud', 'lcmodel'):
        outdir = tmp_path / tool
        outdir.mkdir()
        written = write_basis([moved], tool, 'raw' if tool == 'lcmodel' else 'default', str(outdir))
        back = load_basis(written, detect_format(written)['format'], sw=naa['sw'], sf=naa['sf'])[0]

        assert np.allclose(back['fid'], naa['fid'], atol=1e-4 * np.abs(naa['fid']).max()), tool


def test_osprey_creatines_hold_both_their_parts(tmp_path):
    """MARSS splits Cr and PCr into the 3.9 ppm CH2 and the rest; Osprey's own Cr and PCr hold
    both, and its CrCH2 is a negative CH2 correction, so the parts are summed into Cr and PCr."""
    import numpy as np
    parts = ('Cr391', 'CrNo391', 'PCr393', 'PCrNo393', 'NAA')
    basis = load_basis([os.path.join(SAMPLES, 'JMRUI', f'{p}.txt') for p in parts], 'jMRUI .txt')
    fid = {b['name']: b['fid'] for b in basis}

    written = write_basis(basis, 'osprey', 'mat', str(tmp_path))
    back = {b['name']: b['fid'] for b in load_basis(written, detect_format(written)['format'])}

    assert 'CrCH2' not in back
    for name, summed in (('Cr', fid['Cr391'] + fid['CrNo391']), ('PCr', fid['PCr393'] + fid['PCrNo393'])):
        got = back[name][:len(summed)]
        scale = np.vdot(summed, got) / np.vdot(summed, summed)       # Osprey normalises the set
        assert np.linalg.norm(got - scale * summed) < 1e-6 * np.linalg.norm(got), name


def test_osprey_states_the_centre_a_basis_stating_none_is_at(tmp_path):
    """A basis read from a format with no ppm centre is at 4.65 ppm (writers/centre.py); Osprey
    reads its peaks against the centerFreq it states."""
    basis = load_basis([os.path.join(SAMPLES, 'JMRUI', 'NAA.txt')], 'jMRUI .txt')
    written = write_basis(basis, 'osprey', 'mat', str(tmp_path))

    assert load_basis(written, detect_format(written)['format'])[0]['centerFreq'] == 4.65


def test_nifti_mrs_keeps_a_stated_centre(tmp_path):
    """NIfTI-MRS states no ppm centre either (FSL-MRS reads 1H at 4.65); its .raw step, before
    spec2nii, shifts a basis put at 3.0 ppm as the other writers do."""
    import numpy as np
    from readers.read_lcmodel import read_lcmodel_raw
    from writers.write_niftimrs import _write_raw
    naa = load_basis([os.path.join(SAMPLES, 'JMRUI', 'NAA.txt')], 'jMRUI .txt')[0]
    t = np.arange(len(naa['fid'])) / naa['sw']
    moved = dict(naa, fid=naa['fid'] * np.exp(2j * np.pi * 1.65 * naa['sf'] * t), centerFreq=3.0)

    _write_raw(moved, str(tmp_path / 'NAA.raw'))
    back = read_lcmodel_raw(str(tmp_path / 'NAA.raw'), sw=naa['sw'], sf=naa['sf'])

    assert np.allclose(back['fid'], naa['fid'], atol=1e-4 * np.abs(naa['fid']).max())


def test_mrscloud_converted_to_basis_is_mrscloud_own_basis(tmp_path):
    """MRSCloud's basis sets are centred near 3.0 ppm, as the ppm axis each file stores says
    (2.97 here); converted to .BASIS they are MRSCloud's own .BASIS export of the same set,
    up to one scale for the set and the 5 digits a .BASIS keeps."""
    import re
    import numpy as np
    folder = os.path.join(SAMPLES, 'MRSCLOUD')
    files = sorted(os.path.join(folder, f) for f in os.listdir(folder)
                   if f.startswith('Siemens_') and f.endswith('.mat'))
    written = write_basis(load_basis(files, 'MRSCloud .mat'), 'lcmodel', 'basis', str(tmp_path))

    def spectra(path):
        out = {}
        for block in re.split(r'\$BASIS(?![0-9])', open(path).read())[1:]:
            name = re.search(r"METABO\s*=\s*'([^']*)'", block).group(1).strip()
            body = re.split(r'(?m)^\s*\$END\s*$', block, maxsplit=1)[1].split('$')[0]
            v = np.array(body.split(), dtype=float)
            out[name] = v[0::2] + 1j * v[1::2]
        return out

    own = spectra(os.path.join(SAMPLES, 'LCModel', 'BASIS_file', 'LCModel_Siemens_UnEdited_sLASER_105_.BASIS'))
    ours = spectra(written[0])
    scales = []
    for name, spec in own.items():
        a = np.fft.ifft(spec)                              # MRSCloud's is not zero-filled
        b = np.fft.ifft(ours[name])[:len(a)]
        k = np.vdot(a, b) / np.vdot(a, a)
        scales.append(k)
        assert np.linalg.norm(b - k * a) < 1e-4 * np.linalg.norm(b), name
    assert np.allclose(scales, scales[0], rtol=1e-4)


def test_spant_and_abfit_get_the_format_asked_for(tmp_path):
    """spant and ABfit read a .BASIS (BasiC's recommended output for them) and .raw files."""
    basis = _jmrui()
    for tool in ('spant', 'abfit'):
        for out_fmt, wanted in (('basis', ['basis_set.BASIS']), ('raw', ['Cr391.raw', 'NAA.raw'])):
            outdir = tmp_path / tool / out_fmt
            written = write_basis(basis, tool, out_fmt, str(outdir))

            assert sorted(os.path.basename(p) for p in written) == wanted, (tool, out_fmt)
