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


def test_raw_files_without_seqpar_are_lcmodel_raw_and_ask_for_what_they_lack():
    """A .raw holding only $NMID is plain LCModel .raw, whoever wrote it (MARSS does too):
    the label says LCModel, the spectral width, frequency and echo time are asked for."""
    folder = os.path.join(SAMPLES, 'LCModel', 'raw_files')
    files = sorted(os.path.join(folder, f) for f in os.listdir(folder) if f.endswith('.raw'))
    result = detect_format(files)

    assert result['status'] == 'full'
    assert result['format'] == 'LCModel .raw'
    assert result['compatible_tools'] == ['LCModel', 'SPANT', 'TARQUIN']
    assert sorted(result['missing']) == ['sf', 'sw', 'te']
    assert 'MARSS' in result['message']
    assert len(load_basis(files, result['format'], sw=4000.0, sf=123.2)) == 23


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


def _folder(*parts):
    folder = os.path.join(SAMPLES, *parts)
    return [os.path.join(folder, f) for f in os.listdir(folder) if not f.startswith('.')]


def test_a_folder_holding_two_lcmodel_basis_files_names_both_and_each_loads():
    found = detect_format(_folder('LCModel', 'BASIS_file'))

    assert found['status'] == 'several' and found['compatible_tools'] == []
    sets = {s['name']: s for s in found['sets']}
    assert list(sets) == ['LCModel_Siemens_UnEdited_sLASER_105_.BASIS', 'SLASER_Siemens26.basis']
    for name, n, te in (('LCModel_Siemens_UnEdited_sLASER_105_.BASIS', 25, 105.0),
                        ('SLASER_Siemens26.basis', 22, 26.0)):
        picked = sets[name]
        assert [os.path.basename(f) for f in picked['files']] == [name]
        assert picked['format'] == 'LCModel .BASIS' and 'LCModel' in picked['compatible_tools']
        basis = load_basis(picked['files'], picked['format'])
        assert len(basis) == len(picked['data']['names']) == n and basis[0]['te'] == te
    assert 'LCModel_Siemens_UnEdited_sLASER_105_.BASIS (LCModel .BASIS, 25 metabolites)' in found['message']


def test_mrscloud_singles_beside_their_combined_basis_are_two_sets_however_picked():
    """Browse folder lists the files in the folder's order, Browse files sorted: both are two sets."""
    files = _folder('MRSCLOUD')
    for picked_as in (files, sorted(files)):
        found = detect_format(picked_as)
        assert found['status'] == 'several'
        sets = {s['name']: s for s in found['sets']}
        assert list(sets) == ['BASIS_Siemens_UnEdited_sLASER_105_wMM.mat', 'Siemens_UnEdited_sLASER_TE105_*.mat']

        combined, singles = sets.values()
        assert combined['format'] == 'Osprey .mat' and singles['format'] == 'MRSCloud .mat'
        assert len(singles['files']) == 25 and all(os.path.basename(f).startswith('Siemens_')
                                                   for f in singles['files'])
        assert len(load_basis(singles['files'], singles['format'])) == 25
        assert len(load_basis(combined['files'], combined['format'])) == len(combined['data']['names'])


def test_a_folder_holding_one_set_is_detected_as_before():
    singles = [f for f in _folder('MRSCLOUD') if os.path.basename(f).startswith('Siemens_')]
    for files, fmt in ((_folder('MARSS', 'basis_functions_metab_mm_mat'), 'MARSS .mat'),
                       (_folder('FIDA'), 'FID-A .mat'), (_folder('OSPREY'), 'Osprey .mat'),
                       (singles, 'MRSCloud .mat'), (_folder('LCModel', 'BASIS_file')[:1], 'LCModel .BASIS')):
        found = detect_format(files)
        assert (found['status'], found['format'], 'sets' in found) == ('full', fmt, False)


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

    names, fids, _, _ = _osprey(basis, tmp_path)                    # as Osprey reads the .mat
    back = dict(zip(names, fids.T))

    assert 'CrCH2' not in back
    for name, summed in (('Cr', fid['Cr391'] + fid['CrNo391']), ('PCr', fid['PCr393'] + fid['PCrNo393'])):
        # from the second point: Osprey's DC correction (op_dccorr) takes a constant off the spectrum
        got, summed = back[name][1:len(summed)], summed[1:]
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
    """spant and ABfit read a .BASIS (BasiC's recommended output for them, and its default) and
    .raw files."""
    basis = _jmrui()
    for tool in ('spant', 'abfit'):
        for out_fmt, wanted in (('basis', ['basis_set.BASIS']), ('raw', ['Cr391.raw', 'NAA.raw']),
                                ('default', ['basis_set.BASIS'])):
            outdir = tmp_path / tool / out_fmt
            written = write_basis(basis, tool, out_fmt, str(outdir))

            assert sorted(os.path.basename(p) for p in written) == wanted, (tool, out_fmt)


def test_osprey_names_are_a_cell_of_strings(tmp_path):
    """Osprey finds its metabolites by strcmp on BASIS.name, a cell of strings in its own sets;
    a char matrix pads the shorter names with blanks, 'NAA  ', which no strcmp finds."""
    import numpy as np
    import scipy.io as sio
    written = write_basis(_jmrui(), 'osprey', 'mat', str(tmp_path))
    names = sio.loadmat(written[0])['BASIS']['name'][0, 0]

    assert names.dtype == object
    assert {'NAA', 'Cr'} <= {str(np.ravel(n)[0]) for n in names.ravel()}


def test_vespa_priors_are_read_back(tmp_path):
    """core.detect names a VeSPA .xml, and readers/read_vespa reads it; load_basis takes it too."""
    basis = _jmrui()
    written = write_basis(basis, 'vespa_analysis', 'default', str(tmp_path))
    found = detect_format(written)

    assert found['format'] == 'VeSPA .xml'
    assert sorted(b['name'] for b in load_basis(written, found['format'])) == ['Cr391', 'NAA']


_VESPA = {}


def _vespa_jmrui():
    """The 23 jMRUI sample functions, and the VeSPA prior written for them (once: it takes a while)."""
    import tempfile
    if not _VESPA:
        folder = os.path.join(SAMPLES, 'JMRUI')
        basis = load_basis([os.path.join(folder, f) for f in sorted(os.listdir(folder))], 'jMRUI .txt')
        _VESPA['basis'] = basis
        _VESPA['path'] = write_basis(basis, 'vespa_analysis', 'default', tempfile.mkdtemp())[0]
    return _VESPA['basis'], _VESPA['path']


def _vespa_lines(path):
    """Each prior metabolite's lines as (ppm, area, phase) arrays, its spins, and the prior's comment."""
    import numpy as np
    import xml.etree.ElementTree as ET
    prior = ET.parse(path).getroot().find('prior')
    lines = {m.findtext('name'): np.array([[float(l.findtext(k)) for k in ('ppm', 'area', 'phase')]
                                           for l in m.findall('line')])
             for m in prior.findall('prior_metabolite')}
    spins = {m.findtext('name'): int(m.findtext('spins')) for m in prior.findall('prior_metabolite')}
    return lines, spins, prior.findtext('comment')


def test_vespa_priors_list_each_line_once():
    """VeSPA adds up every line of a prior metabolite: no line is written twice."""
    _, path = _vespa_jmrui()
    lines, _, comment = _vespa_lines(path)

    assert len(lines) == 21
    for name, found in lines.items():
        ppms = list(found[:, 0])
        assert len(ppms) == len(set(ppms)), (name, ppms)
    naa = lines['NAA']
    assert 2.00 <= naa[naa[:, 1].argmax(), 0] <= 2.01
    assert 'PPM from literature' not in comment


def _vespa_fid(lines, n, sw, sf, shape=lambda t: 1.0, centre=4.65):
    """VeSPA's create_fid (common/util/generic_spectral.py) for prior lines (ppm, area, phase), at a
    dataset of n points, sw Hz, sf MHz and resppm `centre`, times a lineshape shape(t)."""
    import numpy as np
    t = np.arange(n)[:, None] / sw
    hz = sw / 2 - sf * (lines[:, 0] - centre)                            # mrs_dataset.ppm2hz
    return (lines[:, 1] * np.exp(1j * np.radians(lines[:, 2]))
            * np.exp(1j * 2 * np.pi * hz * t)).sum(axis=1) * shape(t[:, 0])


def _as_vespa_holds(fid):
    """A basis FID as VeSPA holds the same spectrum: conjugated (its ppm axis runs the other way)
    and chopped (ppm2hz's sw / 2)."""
    import numpy as np
    return np.conj(fid) * (-1.0) ** np.arange(len(fid))


def test_vespa_lines_are_the_lines_the_fid_holds(tmp_path):
    """A set whose lines share a Voigt lineshape (2 Hz Lorentzian, 3 Hz Gaussian), as glycine's
    singlet shows it: a function of known lines (a singlet, a doublet, a line with a phase) is
    written line for line, undamped, at its ppm, its amplitude at t = 0 and its phase in VeSPA's
    sign (VeSPA's FID is BasiC's conjugate). A copy with its first point halved, as FID-A's sets
    hold it, gives the same lines at its own scale. A 30 Hz hump is no lines of that lineshape:
    it is left out, and the prior says so and states the lineshape. The echo time given is the
    prior's."""
    import numpy as np
    sw, sf, n = 4000.0, 123.2618, 2048
    t = np.arange(n) / sw
    voigt = np.exp(-np.pi * 2.0 * t - (np.pi * 3.0 * t) ** 2 / (4 * np.log(2)))

    def fid(lines):                                                      # ppm, area, phase
        return voigt * (lines[:, 1] * np.exp(1j * np.radians(lines[:, 2]))
                        * np.exp(2j * np.pi * np.outer(t, (lines[:, 0] - 4.65) * sf))).sum(axis=1)
    known = np.array([[3.027, 3.0, 0.0],
                      [1.31 + 3.45 / sf, 1.5, 0.0],
                      [1.31 - 3.45 / sf, 1.5, 0.0],
                      [3.913, 0.7, 40.0]])
    halved = 2 * fid(known)
    halved[0] /= 2
    hump = 5.0 * np.exp(-np.pi * 30.0 * t + 2j * np.pi * (0.9 - 4.65) * sf * t)
    basis = [{'name': name, 'fid': f, 'sw': sw, 'sf': sf, 'n': n} for name, f in
             (('Gly', fid(np.array([[3.548, 2.0, 0.0]]))), ('Known', fid(known)),
              ('Halved', halved), ('Hump', hump))]
    path = write_basis(basis, 'vespa_analysis', 'default', str(tmp_path), te=35.0)[0]
    lines, spins, comment = _vespa_lines(path)

    assert set(lines) == {'Gly', 'Known', 'Halved'}
    assert 'Voigt 2.000 Hz Lorentzian, 3.000 Hz Gaussian' in comment and 'Hump' in comment
    for name, scale in (('Known', 1.0), ('Halved', 2.0)):
        got = lines[name][np.argsort(lines[name][:, 0])]
        want = known[np.argsort(known[:, 0])]
        assert len(got) == 4 and spins[name] == 4, (name, got)
        assert np.allclose(got[:, 0], want[:, 0], atol=0.001), name
        assert np.allclose(got[:, 1], scale * want[:, 1], rtol=0.001), name
        assert np.allclose(got[:, 2], -want[:, 2], atol=0.5), name
    import xml.etree.ElementTree as ET
    assert float(ET.parse(path).getroot().find('prior').findtext('seqte')) == 35.0


def test_an_undamped_set_keeps_its_multiplets(tmp_path):
    """A set simulated undamped (the SpinWizard mouse set) measures a 0 Hz lineshape on its singlet.
    A multiplet's components then come out a hair wide (0.0004 Hz): it is lines, not a hump, while
    a 30 Hz hump still is one."""
    import numpy as np
    sw, sf, n = 4000.0, 123.2618, 2048
    t = np.arange(n) / sw

    def fid(lines):                                                      # ppm, area
        return (lines[:, 1] * np.exp(2j * np.pi * np.outer(t, (lines[:, 0] - 4.65) * sf))).sum(axis=1)
    doublet = np.array([[1.31 + 3.45 / sf, 1.5], [1.31 - 3.45 / sf, 1.5]])
    hump = 5.0 * np.exp(-np.pi * 30.0 * t + 2j * np.pi * (0.9 - 4.65) * sf * t)
    basis = [{'name': name, 'fid': f, 'sw': sw, 'sf': sf, 'n': n} for name, f in
             (('Gly', fid(np.array([[3.548, 2.0]]))), ('Lac', fid(doublet) * np.exp(-np.pi * 0.0004 * t)),
              ('MM', hump))]
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('error')                 # no 0 / 0 on a 0 Hz lineshape
        path = write_basis(basis, 'vespa_analysis', 'default', str(tmp_path))[0]
    lines, _, comment = _vespa_lines(path)

    assert set(lines) == {'Gly', 'Lac'} and 'wrong): MM.' in comment
    assert 'Lorentzian 0.000 Hz' in comment            # not a Gaussian of 1e-6 Hz
    largest = lines['Lac'][np.argsort(lines['Lac'][:, 1])[-2:]]          # with lines at 1e-4 of them
    assert np.allclose(np.sort(largest[:, 0]), np.sort(doublet[:, 0]), atol=0.001)  # for the 0.0004 Hz
    assert np.allclose(largest[:, 1], 1.5, rtol=0.001)


def test_a_zero_filled_set_gives_the_lines_of_the_set_it_was_filled_from():
    """An FSL-MRS file converted from a .BASIS keeps its 2n points, the second n near 0 (1e-6), as the
    FSL-MRS sample does. HLSVD fits no damped lines to that edge, so the zero-filling is left off:
    the lines and the lineshape are those of the FIDs before they were filled."""
    import numpy as np
    from writers.write_vespa import vespa_lines
    folder = os.path.join(SAMPLES, 'JMRUI')
    basis = load_basis([os.path.join(folder, f'{n}.txt') for n in ('Cr391', 'Lac')], 'jMRUI .txt')
    rng = np.random.default_rng(0)
    filled = [dict(b, fid=np.r_[b['fid'], 1e-6 * np.abs(b['fid']).max() * rng.standard_normal(len(b['fid']))],
                   n=2 * len(b['fid'])) for b in basis]
    want, got = vespa_lines(basis), vespa_lines(filled)

    assert got['lineshape'] == want['lineshape'] and not got['short'] and not got['humps']
    for (name, lines, _), (_, filled_lines, _) in zip(want['metabolites'], got['metabolites']):
        assert np.allclose(np.array(filled_lines)[:, :2], np.array(lines)[:, :2], rtol=1e-3, atol=1e-3), name


def test_the_priors_writers_refuse_two_functions_of_one_name(tmp_path):
    """The priors name each function's lines after it, and VeSPA's importer, PyAMARES and OXSA keep
    them by name: with two functions named NAA, the second one's lines took the first one's place.
    Every priors writer now says so."""
    import pytest
    folder = os.path.join(SAMPLES, 'JMRUI')
    basis = load_basis([os.path.join(folder, f'{n}.txt') for n in ('NAA', 'Cr391')], 'jMRUI .txt')
    twice = [basis[0], dict(basis[1], name='NAA')]
    for tool in ('vespa_analysis', 'gava', 'midas', 'pyamares', 'oxsa'):
        with pytest.raises(RuntimeError, match='more than one function is named NAA'):
            write_basis(twice, tool, 'default', str(tmp_path / tool))


def test_an_oxsa_hump_component_narrower_than_the_reference_has_no_negative_offset(tmp_path):
    """OXSA holds a hump's components as rows whose linewidth is the reference line's plus an offset.
    A component narrower than the reference (a sharp component too small to count against the hump)
    gave a negative offset, which OXSA's initializePriorKnowledge cannot bound: the offset is now
    never below 0."""
    import numpy as np
    sw, sf, n = 4000.0, 123.2618, 2048
    t = np.arange(n) / sw
    lorentz = lambda ppm, width, amp: amp * np.exp(-np.pi * width * t + 2j * np.pi * (ppm - 4.65) * sf * t)
    basis = [{'name': name, 'fid': f, 'sw': sw, 'sf': sf, 'n': n} for name, f in
             (('Gly', lorentz(3.548, 2.0, 2.0)),
              ('MM', lorentz(0.9, 30.0, 5.0) + lorentz(1.2, 1.0, 0.02)))]      # a sharp 0.4 % beside the hump
    path = write_basis(basis, 'oxsa', 'default', str(tmp_path))[0]
    _, comment, fields, peaks = _oxsa_pk(path)

    rows = [p for p in peaks if p['IV']['peakName'].startswith('MM')]
    assert len(rows) >= 2 and 'damped component' in ' '.join(comment)
    assert all(p['PK']['base_linewidth'] is None or p['PK']['base_linewidth'] >= 0 for p in rows)
    assert any(p['PK']['base_linewidth'] is None and p['PK']['G_linewidth'] is not None for p in rows[1:])


def test_vespa_lines_rebuild_each_jmrui_function():
    """The jMRUI sample's lines share a 1 Hz Lorentzian (its singlets show it). VeSPA's create_fid
    formula, given the written lines and times that lineshape, as VeSPA's fit applies one
    lineshape to all of a metabolite's lines, gives back each function as VeSPA holds it, to the
    writer's criterion (with 10 % for the ppm, areas and phases written rounded, which move the
    residual by about 2e-6). The lines do not cancel each other: their areas add up to at most
    twice the FID at t = 0. MM, a broad hump, and MM09, a 10 Hz singlet, hold no line of the
    set's lineshape: they are left out, and the prior says so and points to VeSPA's own
    macromolecule model. NAA's lines add up to NAA's FID at t = 0 as creatine's add up to
    creatine's: one scale for the whole set."""
    import numpy as np
    from writers.centre import at_centre
    from writers.lines import TOL
    basis, path = _vespa_jmrui()
    lines, spins, comment = _vespa_lines(path)

    assert set(lines) == {b['name'] for b in basis} - {'MM', 'MM09'}
    assert 'Lorentzian 1.000 Hz' in comment
    assert comment.rstrip().endswith('analysis/constants.py:548).') and ': MM, MM09.' in comment
    residual, cancel = {}, {}
    for core in at_centre(basis):
        name = core['name']
        if name in ('MM', 'MM09'):
            continue
        fid, sw, sf = core['fid'], core['sw'], core['sf']
        rebuilt = _vespa_fid(lines[name], len(fid), sw, sf, lambda t: np.exp(-np.pi * 1.0 * t))
        held = _as_vespa_holds(fid)
        residual[name] = np.linalg.norm((rebuilt - held)[1:]) / np.linalg.norm(held[1:])
        cancel[name] = lines[name][:, 1].sum() / abs(fid[0])
    assert max(residual.values()) <= 1.1 * TOL, residual
    assert max(cancel.values()) <= 2, cancel

    def at_zero(name):
        found = lines[name]
        return np.sum(found[:, 1] * np.exp(-1j * np.radians(found[:, 2])))
    fids = {b['name']: b['fid'] for b in basis}
    assert abs(at_zero('NAA') / at_zero('CrNo391') / (fids['NAA'][0] / fids['CrNo391'][0]) - 1) < 1e-3
    assert spins['NAA'] == 5 and spins['Cho'] == 5


def _long_fid(n, sw, components, sf=123.2618):
    """A FID of n points at sw Hz: (ppm, width in Hz, amplitude) Lorentzian components."""
    import numpy as np
    t = np.arange(n) / sw
    return sum(amp * np.exp(-np.pi * width * t + 2j * np.pi * (ppm - 4.65) * sf * t)
               for ppm, width, amp in components)


def test_a_long_vespa_set_is_written_in_seconds(tmp_path):
    """MRSCloud's sets are 2 s, 8192 points. The order search of the lines' HLSVD computes the
    Hankel matrix's singular vectors for the orders it tries, not for the 300 it may try: a
    singlet's took 5 s on its own, a 32-function set's three minutes. Three such functions are
    written in well under the time one took (the bound is generous: the machine may be busy),
    each to its lines."""
    import time
    sw, n = 4000.0, 8192
    basis = [{'name': name, 'fid': _long_fid(n, sw, parts), 'sw': sw, 'sf': 123.2618, 'n': n}
             for name, parts in (('Cr', [(3.027, 1.0, 3.0)]), ('sI', [(3.34, 1.0, 1.0)]),
                                 ('Lac', [(1.31 + 3.45 / 123.2618, 1.0, 1.5),
                                          (1.31 - 3.45 / 123.2618, 1.0, 1.5)]))]
    start = time.time()
    path = write_basis(basis, 'vespa_analysis', 'default', str(tmp_path))[0]
    took = time.time() - start
    lines, _, _ = _vespa_lines(path)

    assert {name: len(found) for name, found in lines.items()} == {'Cr': 1, 'sI': 1, 'Lac': 2}
    assert took < 60.0, took                      # 5-7 s here unloaded; minutes before


def test_the_line_refinement_takes_p_plus_one_rows_not_the_fid(monkeypatch):
    """The refinement of the lines' frequencies is scipy's trust-region step, which depends on the
    Jacobian J only through J^H J and J^H r. It is handed those as p + 1 rows for p lines, not
    as the residual's 2n rows (16384 for MRSCloud: an SVD of 16384 x p each step), and the lines
    still meet the criterion. Components 1.0 to 1.4 Hz wide take many lines of a 1 Hz
    lineshape, and the refinement moves them."""
    import numpy as np
    from writers import lines as fit
    sw, n = 4000.0, 4096
    fid = _long_fid(n, sw, [(2.30, 1.4, 1.0), (2.32, 1.2, 0.7), (2.35, 1.0, 0.5)])
    rows = []
    least_squares = fit.least_squares

    def spy(fun, x0, jac, **options):
        rows.append((len(x0), jac(x0).shape, fun(x0).shape))
        return least_squares(fun, x0, jac=jac, **options)
    monkeypatch.setattr(fit, 'least_squares', spy)
    found = fit.fit_lines(fid, sw, {'lorentz': 1.0, 'gauss': 0.0})

    assert rows and all(jac == (p + 1, p) and res == (p + 1,) for p, jac, res in rows), rows[:3]
    assert found['residual'] <= fit.TOL and len(found['freq']) > 10


def test_a_lorentzian_s_line_sums_are_the_products_of_its_lines():
    """For a Lorentzian set the Gram matrix of the lines and the sums the refinement's normal
    equations take (t**l L**2 exp(2 pi i d t) summed over the FID, l = 0, 1, 2) are in closed form,
    geometric sums and their derivatives: they are the products of the lines' matrix, for a
    0.5 s and a 2 s FID, lines 0.3 Hz apart and the same line twice included. The lines' matrix,
    built in blocks, is the lines'."""
    import numpy as np
    from writers.lines import _Lines
    rng = np.random.default_rng(1)
    for n, sw in ((2047, 4000.0), (8191, 4096.0)):
        model = _Lines(n, sw, {'lorentz': 1.0, 'gauss': 0.0})
        f = np.r_[rng.uniform(-sw / 2, sw / 2, 30), 10.0, 10.3, 10.3]
        A = model.L[:, None] * np.exp(2j * np.pi * np.outer(model.t, f))
        cols = np.arange(len(f)) % 2 == 0
        closed = model.sums(f, None, 3, cols)
        products = [A.conj().T @ A] + [A.conj().T @ (model.t[:, None] ** l * A[:, cols]) for l in (1, 2)]

        assert model.closed
        assert np.abs(model.atoms(f) - A).max() < 1e-11
        for got, want in zip(closed, products):
            assert np.abs(got - want).max() < 1e-11 * np.abs(want).max()


def test_a_vespa_export_reports_each_metabolite_as_it_is_fitted(tmp_path):
    """Screen 4 lists each metabolite as its lines are fitted: write_basis hands the VeSPA writer a
    progress callback (screen 4 gives its step list's), called before the functions are decomposed,
    with the set's lineshape, and once a metabolite, in order, a hump said to be left out. The GUI's list is not tested here."""
    folder = os.path.join(SAMPLES, 'JMRUI')
    basis = load_basis([os.path.join(folder, f) for f in ('NAA.txt', 'Cr391.txt', 'MM09.txt')],
                       'jMRUI .txt')
    steps = []
    write_basis(basis, 'vespa_analysis', 'default', str(tmp_path),
                progress=lambda text, sub=None: steps.append((text, sub)))

    assert [text for text, _ in steps] == ['Decomposing each function', "Lineshape of the set's lines",
                                           'Lines 1/3: NAA', 'Lines 2/3: Cr391', 'Lines 3/3: MM09']
    assert steps[0][1].startswith('3 functions') and steps[1][1].startswith('Lorentzian 1.000 Hz')
    import re
    assert re.match(r'\d+ lines, residual', steps[2][1]) and 'hump' in steps[4][1]


def test_gava_and_midas_hold_the_vespa_prior_s_lines_as_vespa_simulation_writes_them(tmp_path):
    """VeSPA-Simulation writes its simulation's same ppms, areas and phases to its Analysis prior,
    to GAVA and to MIDAS (simulation/dialog_mixed_metabolite_output.py). So the GAVA rows and the
    MIDAS fitt_PriorLine values of the jMRUI sample are the VeSPA prior's lines, number for
    number: absolute ppm, areas on the set's one scale (factor 1), phases in VeSPA's sign, the
    humps MM and MM09 left out and named. The layout is VeSPA's: GAVA's ';' comment lines ending
    with its results heading, rows abbr, dims 0.0, the line's index from 0, ppm, area, phase,
    tab-separated, no newline after the last; MIDAS's comments line0000 on and params numbered
    over the file, indented with tabs by ElementTree, as VeSPA's util_xml.indent and
    tree.write leave them."""
    import numpy as np
    import xml.etree.ElementTree as ET
    basis, path = _vespa_jmrui()
    prior, _, _ = _vespa_lines(path)
    (tmp_path / 'gava').mkdir()
    (tmp_path / 'midas').mkdir()
    gava = open(write_basis(basis, 'gava', 'default', str(tmp_path / 'gava'))[0]).read()
    midas_path = write_basis(basis, 'midas', 'default', str(tmp_path / 'midas'))[0]

    header = [line for line in gava.split('\n') if line.startswith(';')]
    rows = [line.split('\t') for line in gava.split('\n') if not line.startswith(';')]
    assert gava.split('\n')[:len(header)] == header and not gava.endswith('\n')
    assert header[-3:] == [';', ';Simulation Spectral Results', ';' + '-' * 75]
    root = ET.parse(midas_path).getroot()
    comments = root.find('VESPA_SIMULATION_MIDAS_EXPORT').findall('comment')
    params = root.find('PRIOR_METABOLITE_INFORMATION').findall('param')
    assert [c.get('line') for c in comments] == [f'line{i:04d}' for i in range(len(comments))]
    assert [p.get('name') for p in params] == [f'fitt_PriorLine{i:05d}' for i in range(1, len(params) + 1)]
    midas = [p.get('value').split('++') for p in params]

    for fields in (rows, midas):
        assert all(len(f) == 8 and f[1:4] == ['0.0', '0.0', '0.0'] for f in fields)
        names = list(dict.fromkeys(f[0] for f in fields))
        assert names == list(prior)
        for name in names:
            mine = [f for f in fields if f[0] == name]
            assert [int(f[4]) for f in mine] == list(range(len(mine)))
            assert np.array_equal(np.array([[float(x) for x in f[5:]] for f in mine]), prior[name]), name
    for text in ('\n'.join(header), '\n'.join(c.get('value') for c in comments)):
        assert 'Lorentzian 1.000 Hz' in text and 'factor 1' in text and 'MM, MM09' in text
    xml = open(midas_path).read()
    assert xml.startswith('<FITT_Generic_XML Creation_date="') and '\n\t<VESPA_SIMULATION' in xml
    assert '\n\t\t<param name="fitt_PriorLine00001" value="' in xml and xml.endswith('</FITT_Generic_XML>')


def test_gava_and_midas_state_the_echo_time_and_report_each_metabolite(tmp_path):
    """GAVA's and MIDAS's headers state the echo time given (they always said 26 ms), and screen 4
    lists each metabolite as its lines are fitted, as for VeSPA."""
    import xml.etree.ElementTree as ET
    for tool in ('gava', 'midas'):
        steps = []
        (tmp_path / tool).mkdir()
        path = write_basis(_jmrui(), tool, 'default', str(tmp_path / tool), te=35.0,
                           progress=lambda text, sub=None: steps.append(text))[0]
        if tool == 'gava':
            header = [line[1:] for line in open(path).read().split('\n') if line.startswith(';')]
        else:
            header = [c.get('value') for c in ET.parse(path).getroot().iter('comment')]

        assert 'TE: 35.0 ms' in header, tool
        assert steps == ['Decomposing each function', "Lineshape of the set's lines",
                         'Lines 1/2: NAA', 'Lines 2/2: Cr391'], tool


def test_gava_and_midas_name_the_source_file_without_its_folders(tmp_path):
    """The GAVA and MIDAS headers name the file the set came from, not its full path, which would
    show the user's folders in a file they may share."""
    import xml.etree.ElementTree as ET
    for tool in ('gava', 'midas'):
        (tmp_path / tool).mkdir()
        path = write_basis(_jmrui(), tool, 'default', str(tmp_path / tool))[0]
        if tool == 'gava':
            header = [line[1:] for line in open(path).read().split('\n') if line.startswith(';')]
        else:
            header = [c.get('value') for c in ET.parse(path).getroot().iter('comment')]

        assert 'Source: NAA.txt' in header, tool
        assert SAMPLES not in open(path).read() and 'sample_basis_sets' not in open(path).read(), tool


def test_jmrui_files_name_themselves_without_the_output_folder(tmp_path):
    """Each jMRUI .txt states its own file name, as FSL-MRS's writejMRUItxt does, not the full
    output path, which would show the user's folders in a file they may share. FSL-MRS names a
    function without SignalNames after that name, so it is the metabolite's; BasiC's reader
    reads the files back as before."""
    import numpy as np
    outdir = tmp_path / 'my_folders' / 'jmrui_out'
    basis = _jmrui()
    written = write_basis(basis, 'jmrui', 'default', str(outdir))

    for path in written:
        text = open(path).read()
        assert f'\nFilename: {os.path.basename(path)}\n' in text
        assert str(tmp_path) not in text and 'my_folders' not in text and 'jmrui_out' not in text
    back = load_basis(written, detect_format(written)['format'])
    assert [b['name'] for b in back] == [b['name'] for b in basis]
    for a, b in zip(basis, back):
        assert np.allclose(a['fid'], b['fid'], atol=1e-12 * np.abs(a['fid']).max())


def test_the_tool_lists_say_what_oxsa_gets_and_reads():
    """OXSA gets a PK_<set>.m prior knowledge function (write_oxsa), not a script, and reads that,
    not a PyAMARES CSV: the PyAMARES peaks handed to write_oxsa are gone, from the code and from
    what it says."""
    from core.formats import ALL_TOOLS, get_compatible_tools
    assert next(t for t in ALL_TOOLS if t['id'] == 'oxsa')['ext'] == '.m priors'
    assert get_compatible_tools('PyAMARES .csv') == ['PyAMARES']
    texts = {name: open(os.path.join(ROOT, name), encoding='utf-8').read()
             for name in ('README.md', 'basis_converter/readers/read_pyamares.py',
                          'basis_converter/core/detect.py', 'basis_converter/core/load.py',
                          'basis_converter/core/formats.py', 'basis_converter/screen2_output.py')}
    for name, text in texts.items():
        for stale in ('write_oxsa()', 'for write_oxsa', 'convert to OXSA', '.m script',
                      'run_oxsa_amares', 'literature chemical shifts', 'relative areas'):
            assert stale not in text, (name, stale)


def _pyamares_peaks(path, move=None):
    """The peaks of a PyAMARES prior as PyAMARES evaluates them (kernel/PriorKnowledge.py 0.3.28),
    in its internal units: ak, freq (Hz), dk (pi * Hz), phi (radians) and g. A number cell is in
    the CSV's units (ppm times the MHz the comment states, Hz, degrees, unitconverter); an
    expression names a peak to its left and is in the internal units: <peak>*<ratio>,
    <peak>+-<offset>Hz, <peak> or <peak>+-<offset in dk>, <peak>+-<radians>, <peak>. Numbers in
    expressions are positional (refindall reads one number and its unit; 1e-05 would be two).
    A sub-peak (a name with digits) is tied to its base (the name without them) in every row.
    `move`: offsets added to the first column's freq, dk and phi before the ties are followed,
    as a fit moving the free parameters would. Returns the peaks in order, the comment lines,
    the bounds by peak and each peak's tied cells as {row: the peak named}."""
    import csv
    import re
    import numpy as np
    text = open(path).read().splitlines()
    comment = [line for line in text if line.startswith('#')]
    assert text[:len(comment)] == comment                  # find_header_row: comments first
    rows = list(csv.reader(text[len(comment):]))
    assert [r[0] for r in rows] == ['Index', 'Initial Values', 'amplitude', 'chemicalshift',
                                    'linewidth', 'phase', 'g', 'Bounds', 'amplitude',
                                    'chemicalshift', 'linewidth', 'phase', 'g']
    mhz = float(re.search(r'MHz=([0-9.]+)', comment[0]).group(1))
    number = r'([+-][0-9]+(?:\.[0-9]+)?)'
    keys = ('ak', 'freq', 'dk', 'phi', 'g')
    forms = (r'\*([0-9]+(?:\.[0-9]+)?)', number + 'Hz', number + '?', number, '')
    units = (1.0, mhz, np.pi, np.pi / 180, 1.0)
    peaks, bounds, ties = {}, {}, {}
    for i, name in enumerate(rows[0][1:], 1):
        bounds[name] = [rows[k][i] for k in range(8, 13)]
        peak, ties[name] = {}, {}
        for key, form, unit, cell in zip(keys, forms, units, (rows[k][i] for k in range(2, 7))):
            try:
                peak[key] = float(cell) * unit
                continue
            except ValueError:
                pass
            other = re.match(r'[A-Za-z_]+', cell).group(0)
            tied = re.fullmatch(re.escape(other) + form, cell)
            assert other in peaks and tied, (name, key, cell)          # KeyError: not to its left
            x = float(tied.group(1)) if tied.groups() and tied.group(1) else (1.0 if key == 'ak' else 0.0)
            peak[key] = peaks[other][key] * x if key == 'ak' else peaks[other][key] + x
            ties[name][key] = other
        base = ''.join(c for c in name if not c.isdigit())
        if name != base:
            assert ties[name] == dict.fromkeys(keys, base), (name, ties[name])
        if move and not peaks:
            peak = dict(peak, **{k: peak[k] + move.get(k, 0.0) for k in ('freq', 'dk', 'phi')})
        peaks[name] = peak
    return peaks, comment, bounds, ties


def _pyamares_fid(peaks, n, sw):
    """PyAMARES' model (kernel/fid.py equation6) of peaks, at deadtime 0."""
    import numpy as np
    t = np.arange(n) / sw
    return sum(p['ak'] * np.exp(1j * p['phi']) * np.exp(-p['dk'] * (1 - p['g'] + p['g'] * t) * t)
               * np.exp(2j * np.pi * p['freq'] * t) for p in peaks)


def _voigt_set():
    """Functions whose lines share a Voigt lineshape (2 Hz Lorentzian, 3 Hz Gaussian): singlets Gly,
    Cr391 and Cr (the largest), a function of four known lines (two phased), a two-component
    hump named MM09 and the jMRUI sample's MM. Returns the set and the known lines."""
    import numpy as np
    sw, sf, n = 4000.0, 123.2618, 2048
    t = np.arange(n) / sw
    voigt = np.exp(-np.pi * 2.0 * t - (np.pi * 3.0 * t) ** 2 / (4 * np.log(2)))

    def fid(lines):                                                      # ppm, area, phase
        return voigt * (lines[:, 1] * np.exp(1j * np.radians(lines[:, 2]))
                        * np.exp(2j * np.pi * np.outer(t, (lines[:, 0] - 4.65) * sf))).sum(axis=1)
    known = np.array([[3.027, 3.0, 0.0],
                      [1.31 + 5.0 / sf, 1.5, 10.0],
                      [1.31 - 5.0 / sf, 1.4, -20.0],
                      [3.913, 0.7, 40.0]])
    hump = (5.0 * np.exp(-np.pi * 30.0 * t + 2j * np.pi * (0.9 - 4.65) * sf * t)
            + 2.0j * np.exp(-np.pi * 45.0 * t + 2j * np.pi * (1.3 - 4.65) * sf * t))
    mm = load_basis([os.path.join(SAMPLES, 'JMRUI', 'MM.txt')], 'jMRUI .txt')[0]
    basis = [{'name': name, 'fid': f, 'sw': sw, 'sf': sf, 'n': n} for name, f in
             (('Gly', fid(np.array([[3.548, 2.0, 0.0]]))), ('Known', fid(known)),
              ('Cr391', fid(np.array([[3.913, 1.0, 5.0]]))), ('Cr', fid(np.array([[3.027, 3.0, 0.0]]))),
              ('MM09', hump))] + [mm]
    return basis, known


def test_pyamares_gets_each_metabolite_s_lines_as_tied_sub_peaks(tmp_path):
    """PyAMARES fits each peak as ak exp(i phi) exp(-dk (1 - g + g t) t) exp(2 pi i f t), with the
    data shifted by the carrier (kernel/fid.py, PriorKnowledge.initialize_FID): with
    carrier = 4.65 and deadtime = 0 a basis' line c exp(2 pi i f t) is a peak at 4.65 + f / sf ppm
    with phase +angle(c), the basis' own sign, not VeSPA's. A set whose lines share a Voigt
    lineshape (2 Hz Lorentzian, 3 Hz Gaussian): each metabolite is its largest line as a base peak
    with a digit-free name and its other lines as <base>2, <base>3, ... tied to it (amplitude ratio,
    shift offset in Hz, phase offset in radians, linewidth and g the base's), and the Voigt is
    AMARES' damping exactly (dk = pi * 2 + (pi * 3)**2 / (4 ln 2), g its Gaussian share).
    PyAMARES' model of the written peaks gives back each function from its second point to the
    criterion. A hump of two damped components is written as its own metabolite, linewidths tied
    as offsets; the jMRUI sample's MM, which HLSVD components do not rebuild, is left out and
    named. MM09 and MM, Cr391 and Cr would be one row in result_sum: MM09 is MM_a, Cr391 is
    Cr_a. The first line states MHz, the carrier, deadtime=0 and each metabolite's total."""
    import re
    import numpy as np
    from writers.lines import TOL
    basis, known = _voigt_set()
    sw, sf, n = 4000.0, 123.2618, 2048
    t = np.arange(n) / sw
    steps = []
    path = write_basis(basis, 'pyamares', 'default', str(tmp_path), te=35.0,
                       progress=lambda text, sub=None: steps.append((text, sub)))[0]
    peaks, comment, bounds, ties = _pyamares_peaks(path)

    groups = {}
    for name in peaks:
        groups.setdefault(''.join(c for c in name if not c.isdigit()), []).append(name)
    assert list(groups) == ['Cr', 'Gly', 'Known', 'Cr_a', 'MM_a']        # the reference first
    assert groups['Known'] == ['Known', 'Known2', 'Known3', 'Known4']
    assert [p for p in groups['MM_a']] == ['MM_a', 'MM_a2']
    source = {'Cr_a': 'Cr391', 'MM_a': 'MM09'}
    for name in groups:
        core = next(b for b in basis if b['name'] == source.get(name, name))
        want = core['fid'] * np.exp(2j * np.pi * 4.65 * sf * t)              # the carrier's shift
        got = _pyamares_fid([peaks[p] for p in groups[name]], n, sw)
        assert np.linalg.norm((got - want)[1:]) <= 1.1 * TOL * np.linalg.norm(want[1:]), name

    k = sorted(groups['Known'], key=lambda p: peaks[p]['freq'])        # 1.27, 1.35, 3.027, 3.913
    assert np.allclose([peaks[p]['freq'] / sf for p in k], np.sort(known[:, 0]), atol=0.001)
    assert np.allclose([peaks[p]['ak'] for p in k], [1.4, 1.5, 3.0, 0.7], rtol=0.001)
    assert np.allclose(np.degrees([peaks[p]['phi'] for p in k]), [-20.0, 10.0, 0.0, 40.0], atol=0.5)
    dk = np.pi * 2.0 + (np.pi * 3.0) ** 2 / (4 * np.log(2))
    assert np.allclose([peaks[p]['dk'] for p in k], dk, rtol=1e-3) and \
        np.allclose([peaks[p]['g'] for p in k], (dk - np.pi * 2.0) / dk, rtol=1e-3)
    assert np.allclose(sorted(peaks[p]['dk'] / np.pi for p in groups['MM_a']), [30.0, 45.0], atol=1e-3)
    assert bounds['Known'][0] == '(0,' and bounds['Known2'] == ['', '', '(0,', '', '']

    assert 'MHz=123.2618' in comment[0] and 'carrier=4.65' in comment[0] and 'deadtime=0' in comment[0]
    assert 'TE 35 ms' in comment[0]
    total = dict(re.findall(r'([A-Za-z_]+) ([0-9.]+)(?:;|\.$)', comment[0].split('function): ')[1]))
    for name, members in groups.items():
        assert np.isclose(float(total[name]), sum(peaks[p]['ak'] for p in members), rtol=1e-8), name
    text = ' '.join(comment)
    assert 'Left out: MM:' in text and 'g_global=False' in text
    assert 'Cr_a = Cr391' in text and 'MM_a = MM09' in text
    assert [s for s, _ in steps][-1] == 'Lines 6/6: MM' and 'left out' in steps[-1][1]
    assert steps[-2][1].startswith('2 damped components')


def test_pyamares_prior_of_the_jmrui_sample_keeps_each_metabolite_apart(tmp_path):
    """The jMRUI sample as PyAMARES reads it: 22 metabolites (MM left out, MM09 as its 10 Hz peak),
    each one row of result_sum (names distinct once digits are removed), each rebuilt from its
    second point by PyAMARES' model to the criterion, NAA where it is (2.008 ppm, its largest
    line 29.62 with phase -3.2 degrees, as HSVD gives it), at the set's 1 Hz Lorentzian."""
    import re
    import numpy as np
    from writers.lines import TOL
    folder = os.path.join(SAMPLES, 'JMRUI')
    basis = load_basis([os.path.join(folder, f) for f in sorted(os.listdir(folder))], 'jMRUI .txt')
    path = write_basis(basis, 'pyamares', 'default', str(tmp_path))[0]
    peaks, comment, bounds, ties = _pyamares_peaks(path)

    groups = {}
    for name in peaks:
        groups.setdefault(''.join(c for c in name if not c.isdigit()), []).append(name)
    source = {b['name']: b for b in basis}
    renamed = dict(re.findall(r'([A-Za-z_]+) = ([A-Za-z0-9]+)', next(c for c in comment if 'Names' in c)))
    assert len(groups) == 22 and 'MM' not in groups and renamed['MM_a'] == 'MM09'
    assert set(renamed.get(g, g) for g in groups) == set(source) - {'MM'}
    for name, members in groups.items():
        core = source[renamed.get(name, name)]
        t = np.arange(len(core['fid'])) / core['sw']
        want = core['fid'] * np.exp(2j * np.pi * 4.65 * core['sf'] * t)
        got = _pyamares_fid([peaks[p] for p in members], len(t), core['sw'])
        assert np.linalg.norm((got - want)[1:]) <= 1.1 * TOL * np.linalg.norm(want[1:]), name
    naa = peaks['NAA']
    assert abs(naa['freq'] / 123.2618 - 2.008) < 1e-3 and abs(naa['ak'] - 29.62) < 0.01
    assert abs(np.degrees(naa['phi']) + 3.22) < 0.01 and np.isclose(naa['dk'], np.pi, rtol=1e-6)
    assert np.isclose(peaks['MM_a']['dk'], 10 * np.pi, rtol=1e-6) and len(groups['MM_a']) == 1


def test_pyamares_fits_one_phase_one_lineshape_and_one_shift_for_the_set(tmp_path):
    """The data have one zero-order phase, and the basis says where the metabolites sit and how
    broad they are relative to one another. The largest singlet (Cr) is the first column and holds
    the set's free shift (+-0.1 ppm), linewidth (0 to its own + 80 Hz) and phase (+-180 degrees);
    every other metabolite's first peak is tied to it as the basis holds them (shift offset in Hz,
    phase offset in radians, linewidth and g equal; a hump's linewidth an offset in dk and its g
    fixed at 0), its amplitude free. So moving the reference by a shift, a broadening and a phase
    moves every peak by the same: each function comes back times that shift, phase and (for
    the lines) the reference's lineshape. Free: one amplitude per metabolite and the three."""
    import numpy as np
    basis, known = _voigt_set()
    sw, sf, n = 4000.0, 123.2618, 2048
    t = np.arange(n) / sw
    path = write_basis(basis, 'pyamares', 'default', str(tmp_path))[0]
    peaks, comment, bounds, ties = _pyamares_peaks(path)
    bases = [p for p in peaks if not any(c.isdigit() for c in p)]
    assert bases == ['Cr', 'Gly', 'Known', 'Cr_a', 'MM_a']
    ref = peaks['Cr']
    lw = ref['dk'] / np.pi
    assert ties['Cr'] == {} and bounds['Cr'][0] == '(0,' and bounds['Cr'][3:] == ['(-180,180)', '(0,1)']
    lo, hi = (float(x) for x in bounds['Cr'][1].strip('()').split(','))
    assert np.isclose(lo, 2.927, atol=1e-6) and np.isclose(hi, 3.127, atol=1e-6)
    lo, hi = (float(x) for x in bounds['Cr'][2].strip('()').split(','))
    assert lo == 0 and np.isclose(hi, lw + 80, atol=1e-6)
    for b in ('Gly', 'Known', 'Cr_a'):
        assert ties[b] == {'freq': 'Cr', 'dk': 'Cr', 'phi': 'Cr', 'g': 'Cr'}, b
        assert bounds[b] == ['(0,', '', '(0,', '', ''], b
    assert ties['MM_a'] == {'freq': 'Cr', 'dk': 'Cr', 'phi': 'Cr'} and bounds['MM_a'][4] == '(0,0)'
    assert peaks['MM_a']['g'] == 0 and np.isclose(peaks['MM_a']['dk'], 30 * np.pi, atol=1e-6)
    fixed = lambda cell: cell.startswith('(') and len(set(cell.strip('()').split(','))) == 1
    free = [(p, k) for p in peaks for k, key in enumerate(('ak', 'freq', 'dk', 'phi', 'g'))
            if key not in ties[p] and key != 'g' and not fixed(bounds[p][k])]
    assert free == [(b, 0) for b in bases[:1]] + [('Cr', 1), ('Cr', 2), ('Cr', 3)] + \
        [(b, 0) for b in bases[1:]]
    assert f'{len(bases) + 3} free parameters' in ' '.join(comment)
    assert 'the first column, Cr' in ' '.join(comment)

    shift, broad, turn = 0.02 * sf, np.pi * 3.0, np.radians(30.0)
    moved = _pyamares_peaks(path, move={'freq': shift, 'dk': broad, 'phi': turn})[0]
    common = np.exp(2j * np.pi * shift * t + 1j * turn)
    for b in bases:
        members = [p for p in peaks if ''.join(c for c in p if not c.isdigit()) == b]
        before = _pyamares_fid([peaks[p] for p in members], n, sw)
        after = _pyamares_fid([moved[p] for p in members], n, sw)
        if b == 'MM_a':                                       # its own Lorentzian components
            want = before * common * np.exp(-broad * t)
        else:                                                 # the reference's lineshape
            want = before * common
            assert all(np.isclose(moved[p]['dk'], moved['Cr']['dk']) and
                       moved[p]['g'] == moved['Cr']['g'] for p in members), b
            reshape = (lambda q: np.exp(-q['dk'] * (1 - q['g'] + q['g'] * t) * t))
            want = want * reshape(moved['Cr']) / reshape(peaks['Cr'])
        assert np.allclose(after, want, atol=1e-9 * np.abs(before).max()), b


def test_basic_reads_back_its_own_pyamares_prior(tmp_path):
    """BasiC's reader skips the '#' lines and evaluates the tied cells as PyAMARES does, so the
    peaks it lists are the CSV's: NAA where it is and its sub-peaks where their offsets put them,
    in ppm, Hz and degrees, the expressions kept."""
    import numpy as np
    from readers.read_pyamares import read_pyamares
    path = write_basis(_jmrui(), 'pyamares', 'default', str(tmp_path))[0]
    peaks = _pyamares_peaks(path)[0]
    read = read_pyamares(path)

    assert read['names'] == list(peaks) and read['mhz'] == 123.2618
    assert np.allclose(read['amplitude'], [p['ak'] for p in peaks.values()], rtol=1e-12)
    assert np.allclose(read['chemicalshift'], [p['freq'] / 123.2618 for p in peaks.values()], atol=1e-9)
    assert np.allclose(read['linewidth'], [p['dk'] / np.pi for p in peaks.values()], atol=1e-9)
    assert np.allclose(read['phase'], [np.degrees(p['phi']) for p in peaks.values()], atol=1e-7)
    i = read['names'].index('NAA2')
    assert read['expressions']['chemicalshift'][i].startswith('NAA+')
    assert read['expressions']['chemicalshift'][0] is None                  # the reference
    assert abs(read['chemicalshift'][read['names'].index('NAA')] - 2.008) < 1e-3


def test_a_pyamares_prior_basic_wrote_is_detected_as_pyamares(tmp_path):
    """Screen 1 detects the CSV BasiC writes for PyAMARES as a PyAMARES prior, with its peaks."""
    path = write_basis(_jmrui(), 'pyamares', 'default', str(tmp_path))[0]
    found = detect_format([path])

    assert found['format'] == 'PyAMARES .csv' and found['status'] == 'full'
    assert found['data']['names'] == list(_pyamares_peaks(path)[0])


def _oxsa_pk(path):
    """An OXSA prior knowledge function as preparePriorKnowledge assembles it: the function's name,
    its header comment, the three field lists and each peak's bounds, initial values and prior
    knowledge as dicts by field. A cell is a name ('...'), [] (None), [lo,hi] (a list) or a number.
    Checks the file's form: one function named as its file, returning
    AMARES.priorKnowledge.preparePriorKnowledge(fields, values)."""
    import re
    lines = open(path).read().splitlines()
    head = re.fullmatch(r'function pk = ([A-Za-z]\w*)\(\)', lines[0])
    assert head and head.group(1) == os.path.splitext(os.path.basename(path))[0], lines[0]
    comment = [l[2:] if l.startswith('% ') else '' for l in lines[1:]
               if l.startswith('%') and not l.startswith('%%')]

    def cell(text):
        text = text.strip()
        if text.startswith("'"):
            return text[1:-1].replace("''", "'")
        if text == '[]':
            return None
        if text.startswith('['):
            return [float(x) for x in text[1:-1].split(',')]
        return float(text)

    def split(row):                                      # commas outside [] and ''
        cells, depth, quoted, start = [], 0, False, 0
        for i, ch in enumerate(row):
            quoted ^= ch == "'"
            depth += (ch == '[') - (ch == ']') if not quoted else 0
            if ch == ',' and not depth and not quoted:
                cells.append(row[start:i])
                start = i + 1
        return cells + [row[start:]]

    fields, values, i = {}, {}, 0
    while i < len(lines):
        f = re.fullmatch(r'fields\.(\w+) = \{(.*)\};', lines[i])
        v = re.fullmatch(r'values\.(\w+) = \{\.\.\.', lines[i])
        if f:
            fields[f.group(1)] = [cell(c) for c in split(f.group(2))]
        elif v:
            rows = []
            i += 1
            while lines[i] != '};':
                assert lines[i].endswith(';'), lines[i]
                rows.append([cell(c) for c in split(lines[i][:-1])])
                i += 1
            values[v.group(1)] = rows
        i += 1
    assert lines[-2:] == ['pk = AMARES.priorKnowledge.preparePriorKnowledge(fields, values);', 'end']
    arrays = {'Bounds': 'boundsCellArray', 'IV': 'IVCellArray', 'PK': 'PKCellArray'}
    assert all(len(r) == len(fields[k]) for k, a in arrays.items() for r in values[a])
    peaks = [{'B': dict(zip(fields['Bounds'], b)), 'IV': dict(zip(fields['IV'], iv)),
              'PK': dict(zip(fields['PK'], pk))}
             for b, iv, pk in zip(*(values[a] for a in arrays.values()))]
    return head.group(1), comment, fields, peaks


def _oxsa_check_constraints(peaks):
    """checkConstraints.m's rules (5-59), and initializeOptimization's (41-46: a group's peaks with
    bounds stand together). Returns the free parameters as initializeOptimization counts them."""
    import re
    names = [p['B']['peakName'] for p in peaks]
    assert all(p['IV']['peakName'] == p['PK']['peakName'] == p['B']['peakName'] for p in peaks)
    assert len(set(names)) == len(names) and all(re.fullmatch(r'[A-Za-z]\w{0,62}', n) for n in names)
    free = []
    for field in peaks[0]['B']:
        if field == 'peakName':
            continue
        group = 0
        for k, p in enumerate(peaks, 1):
            bound = p['B'][field]
            if bound is None:
                continue
            assert len(bound) == 2 and bound[0] < bound[1], (p['B']['peakName'], field, bound)
            value = p['IV'][field] if field in p['IV'] else p['PK'][field]
            assert bound[0] <= value <= bound[1], (p['B']['peakName'], field, value, bound)
            g = p['PK'].get('G_' + field)
            if g is not None and g == group:
                continue
            group = g if g is not None else group
            free.append((k, field))
    for field in ('G_linewidth', 'G_amplitude', 'G_phase', 'G_chemShiftDelta', 'G_sigma'):
        numbers = [p['PK'].get(field) for p in peaks]
        for g in set(n for n in numbers if n is not None):
            assert numbers.index(g) + 1 == g, (field, g)
            members = [k for k, n in enumerate(numbers, 1) if n == g]
            with_bounds = [k for k in members if peaks[k - 1]['B'].get(field[2:]) is not None]
            between = [k for k in range(members[0], members[-1] + 1) if k not in members and
                       peaks[k - 1]['B'].get(field[2:]) is not None]
            assert not (with_bounds and between), (field, g, between)
    return free


def _oxsa_model(peaks, move=None):
    """Each peak's parameters as createModelConstraints.m and applyModelConstraints.m make them
    from the free ones (here their initial values, plus `move` for row 1's), for the fields the
    writer uses: a free parameter is a peak's own where its bounds are given; a group (G_*) takes
    its first peak's; a peak in a multiplet with a ratio has its amplitude times amplitudeRatio
    (108-128); RelPhase adds degrees to the phase (139-142), base_linewidth Hz to the linewidth
    (99-101); empty bounds and no group fix a parameter at its initial value (46-64)."""
    move = move or {}
    out = []
    for k, p in enumerate(peaks, 1):
        q = {}
        for field in ('chemShift', 'linewidth', 'amplitude', 'phase', 'sigma'):
            if field not in p['IV']:
                q[field] = 0.0
                continue
            g = p['PK'].get('G_' + field) if field != 'chemShift' else None
            lead = peaks[int(g) - 1] if g is not None else p
            if field == 'sigma' and g is not None and p['B']['sigma'] is None:
                q[field] = 0.0                  # a sigma group's peak needs bounds of its own (157)
            elif lead['B'].get(field) is not None:
                q[field] = lead['IV'][field] + ((g or k) == 1) * move.get(field, 0.0)
            else:
                q[field] = p['IV'][field]
        if p['PK']['multiplet'] is not None and p['PK']['amplitudeRatio'] is not None:
            q['amplitude'] *= p['PK']['amplitudeRatio']
        q['phase'] += p['PK']['RelPhase'] or 0.0
        q['linewidth'] += p['PK'].get('base_linewidth') or 0.0
        q['name'], q['group'] = p['PK']['peakName'], p['PK']['G_amplitude']
        out.append(q)
    return out


def _oxsa_fid(params, n, sw, mhz):
    """makeModelFid.m: each peak amplitude exp(i phase pi / 180) exp(-|t| linewidth pi)
    exp(-2 pi^2 t^2 sigma^2) exp(i 2 pi t chemShift imagingFrequency), summed; t = k / sw."""
    import numpy as np
    t = np.arange(n) / sw
    return sum(q['amplitude'] * np.exp(1j * np.radians(q['phase']))
               * np.exp(-np.pi * q['linewidth'] * t) * np.exp(-2 * np.pi ** 2 * t ** 2 * q['sigma'] ** 2)
               * np.exp(2j * np.pi * t * q['chemShift'] * mhz) for q in params)


def _oxsa_groups(params):
    groups = {}
    for q in params:
        groups.setdefault(params[int(q['group']) - 1]['name'], []).append(q)
    return groups


def test_oxsa_gets_a_prior_knowledge_function_in_oxsa_s_own_form(tmp_path):
    """OXSA fits each peak as amplitude exp(i phase) exp(-|t| linewidth pi) exp(-2 pi^2 t^2 sigma^2)
    exp(i 2 pi t chemShift MHz) (makeModelFid.m): a basis' line c exp(2 pi i f t) is a peak at
    chemShift = ppm - 4.65, the carrier, with phase +angle(c) in degrees, the basis' own sign, and
    the set's Voigt (2 Hz Lorentzian, 3 Hz Gaussian) is linewidth 2 Hz and sigma 3 / 2.3548. The
    prior is PK_<set>.m returning AMARES.priorKnowledge.preparePriorKnowledge(fields, values) with
    OXSA's fields, and checkConstraints' rules hold. Each metabolite's lines are rows tied to its
    largest: amplitude by amplitudeRatio in one-component multiplets of its G_amplitude group, phase
    by G_phase 1 and RelPhase, linewidth by G_linewidth 1 (a hump's own widths by base_linewidth),
    sigma by G_sigma; chemical shifts fixed. OXSA's model of the rows gives back each function from
    its second point to the criterion; with VeSPA's phase sign it would not. Moving row 1's
    linewidth and phase moves every metabolite alike. The jMRUI MM is left out and named."""
    import re
    import numpy as np
    from writers.lines import TOL
    basis, known = _voigt_set()
    sw, sf, n = 4000.0, 123.2618, 2048
    t = np.arange(n) / sw
    steps = []
    path = write_basis(basis, 'oxsa', 'default', str(tmp_path), te=35.0,
                       progress=lambda text, sub=None: steps.append((text, sub)))[0]
    assert os.path.basename(path) == 'PK_basis.m'
    name, comment, fields, peaks = _oxsa_pk(path)
    assert fields == {
        'Bounds': ['peakName', 'chemShift', 'linewidth', 'amplitude', 'phase', 'chemShiftDelta',
                   'amplitudeRatio', 'sigma'],
        'IV': ['peakName', 'chemShift', 'linewidth', 'amplitude', 'phase', 'sigma'],
        'PK': ['peakName', 'multiplet', 'chemShiftDelta', 'amplitudeRatio', 'G_linewidth',
               'G_amplitude', 'G_phase', 'RelPhase', 'G_chemShiftDelta', 'refPeak', 'base_linewidth',
               'G_sigma']}
    free = _oxsa_check_constraints(peaks)
    params = _oxsa_model(peaks)
    groups = _oxsa_groups(params)
    assert list(groups) == ['Cr', 'Gly', 'Known', 'Cr391', 'MM09']           # the reference first
    assert [q['name'] for q in groups['Known']] == ['Known', 'Known_2', 'Known_3', 'Known_4']
    assert free == [(1, 'linewidth')] + [(k, 'amplitude') for k in (1, 2, 3, 7, 8)] + \
        [(1, 'phase'), (1, 'sigma')]
    assert all(p['B']['chemShift'] is None for p in peaks)                   # shifts fixed
    assert all(p['PK']['G_phase'] == 1 and p['PK']['G_linewidth'] == 1 for p in peaks)
    assert peaks[0]['B']['phase'] == [0.0, 360.0] and peaks[0]['B']['amplitude'] == [0.0, np.inf]
    assert peaks[0]['PK']['refPeak'] == 1 and all(p['PK']['refPeak'] == 0 for p in peaks[1:])
    assert all(p['PK']['multiplet'] == 1 and p['B']['amplitude'] is None and
               p['PK']['G_amplitude'] == 3 for p in peaks[3:6])              # Known_2-4

    source = {b['name']: b for b in basis}
    for group, members in groups.items():
        want = source[group]['fid']
        got = _oxsa_fid(members, n, sw, sf)
        assert np.linalg.norm((got - want)[1:]) <= 1.1 * TOL * np.linalg.norm(want[1:]), group
    vespa = [dict(q, phase=-q['phase']) for q in groups['Known']]          # VeSPA's sign
    assert np.linalg.norm((_oxsa_fid(vespa, n, sw, sf) - source['Known']['fid'])[1:]) > \
        0.1 * np.linalg.norm(source['Known']['fid'][1:])

    k = sorted(groups['Known'], key=lambda q: q['chemShift'])               # 1.27, 1.35, 3.027, 3.913
    assert np.allclose([q['chemShift'] + 4.65 for q in k], np.sort(known[:, 0]), atol=0.001)
    assert np.allclose([q['amplitude'] for q in k], [1.4, 1.5, 3.0, 0.7], rtol=0.001)
    assert np.allclose([(q['phase'] + 180) % 360 - 180 for q in k], [-20.0, 10.0, 0.0, 40.0], atol=0.5)
    assert all(np.isclose(q['linewidth'], 2.0, rtol=1e-3) and
               np.isclose(q['sigma'], 3.0 / (2 * np.sqrt(2 * np.log(2))), rtol=1e-3) for q in k)
    assert np.allclose(sorted(q['linewidth'] for q in groups['MM09']), [30.0, 45.0], atol=1e-3)
    assert all(q['sigma'] == 0 for q in groups['MM09'])

    moved = _oxsa_groups(_oxsa_model(peaks, move={'linewidth': 3.0, 'phase': 30.0}))
    for group, members in groups.items():
        before, after = _oxsa_fid(members, n, sw, sf), _oxsa_fid(moved[group], n, sw, sf)
        want = before * np.exp(1j * np.radians(30.0) - np.pi * 3.0 * t)
        assert np.allclose(after, want, atol=1e-9 * np.abs(before).max()), group

    text = ' '.join(' '.join(comment).split())
    assert name == 'PK_basis' and 'pk = PK_basis();' in text and 'AMARES.amaresFit(' in text
    assert "'imagingFrequency', 123.2618" in text and 'relative to the carrier, 4.65 ppm' in text
    assert 'TE 35 ms' in text and 'Free: 8 parameters' in text
    assert 'Left out (approximated by nothing): MM:' in text
    listed = text.split('the amount of the basis function): ')[1].split('. ')[0]
    totals = {g: float(v) for g, v in re.findall(r'(\w+) \d+(?:-\d+)? ([0-9.]*[0-9])', listed)}
    for group, members in groups.items():
        assert np.isclose(totals[group], sum(q['amplitude'] for q in members), rtol=1e-8), group
    assert [s for s, _ in steps][-1] == 'Lines 6/6: MM' and 'left out' in steps[-1][1]


def test_oxsa_prior_of_the_jmrui_sample_rebuilds_each_metabolite(tmp_path):
    """The jMRUI sample as OXSA reads it: PK_JMRUI.m, 22 metabolites (MM left out, MM09 as its 10 Hz
    peak), each rebuilt from its second point by OXSA's model of its rows to the criterion; sI's
    largest line is row 1; NAA's first row is its 2.008 ppm line, -2.642 from the carrier, 29.62 at
    -3.22 degrees; 24 free parameters: 22 amplitudes, one linewidth and one phase."""
    import numpy as np
    from writers.lines import TOL
    folder = os.path.join(SAMPLES, 'JMRUI')
    basis = load_basis([os.path.join(folder, f) for f in sorted(os.listdir(folder))], 'jMRUI .txt')
    path = write_basis(basis, 'oxsa', 'default', str(tmp_path))[0]
    assert os.path.basename(path) == 'PK_JMRUI.m'
    name, comment, fields, peaks = _oxsa_pk(path)
    free = _oxsa_check_constraints(peaks)
    groups = _oxsa_groups(_oxsa_model(peaks))

    source = {b['name']: b for b in basis}
    assert len(groups) == 22 and set(groups) == set(source) - {'MM'} and list(groups)[0] == 'sI'
    for group, members in groups.items():
        core = source[group]
        got = _oxsa_fid(members, len(core['fid']), core['sw'], core['sf'])
        want = core['fid']
        assert np.linalg.norm((got - want)[1:]) <= 1.1 * TOL * np.linalg.norm(want[1:]), group
    naa = groups['NAA'][0]
    assert abs(naa['chemShift'] + 2.642) < 1e-3 and abs(naa['amplitude'] - 29.62) < 0.01
    assert abs((naa['phase'] + 180) % 360 - 180 + 3.22) < 0.01 and np.isclose(naa['linewidth'], 1.0)
    assert len(groups['MM09']) == 1 and np.isclose(groups['MM09'][0]['linewidth'], 10.0, rtol=1e-6)
    assert len(free) == 24 and 'Free: 24 parameters' in ' '.join(' '.join(comment).split())


def test_no_prior_writer_reads_peaks_off_the_spectrum_at_literature_shifts():
    """The priors are the basis set's own lines (write_vespa.vespa_lines): the FFT peak picking at
    literature shifts that OXSA last used is gone from the writers."""
    import writers
    import writers.write_vespa as write_vespa
    folder = os.path.dirname(writers.__file__)
    for helper in ('_peaks_at_ppms', '_fft_peaks', 'PPM_CALIB'):
        assert not hasattr(write_vespa, helper), helper
        for f in os.listdir(folder):
            if f.endswith('.py'):
                assert helper not in open(os.path.join(folder, f)).read(), (f, helper)


def test_prior_writers_keep_a_stated_centre(tmp_path):
    """The priors read peaks off the spectrum at 4.65 ppm: a basis put at 3.0 ppm (Osprey's and
    MRSCloud's sets) is shifted first, so NAA is listed at 2.01 ppm."""
    import numpy as np
    naa = load_basis([os.path.join(SAMPLES, 'JMRUI', 'NAA.txt')], 'jMRUI .txt')[0]
    t = np.arange(len(naa['fid'])) / naa['sw']
    moved = dict(naa, fid=naa['fid'] * np.exp(2j * np.pi * 1.65 * naa['sf'] * t), centerFreq=3.0)

    for tool in ('pyamares', 'oxsa'):
        (tmp_path / tool).mkdir()
        for core, label in ((naa, 'stated none'), (moved, 'at 3.0')):
            path = write_basis([core], tool, 'default', str(tmp_path / tool))[0]
            if tool == 'pyamares':
                ppm = _pyamares_peaks(path)[0]['NAA']['freq'] / naa['sf']
            else:
                ppm = _oxsa_groups(_oxsa_model(_oxsa_pk(path)[3]))['NAA'][0]['chemShift'] + 4.65
            assert abs(ppm - 2.008) < 0.01, (tool, label, ppm)


def test_raw_files_keep_the_echo_time(tmp_path):
    """LCModel's .raw states the echo time in $SEQPAR (ECHOT); it is written and read back."""
    written = write_basis(_jmrui(), 'lcmodel', 'raw', str(tmp_path), te=26.0)

    assert [b['te'] for b in load_basis(written, detect_format(written)['format'])] == [26.0, 26.0]


def test_spinwizard_files_keep_the_frequency_they_state(tmp_path):
    """A frequency passed to load_basis is for files that state none; SpinWizard's state theirs."""
    folder = os.path.join(SAMPLES, 'SPINWIZARD')
    files = sorted(os.path.join(folder, f) for f in os.listdir(folder) if '.' not in f)
    basis = load_basis(files, 'SpinWizard', sw=2000.0, sf=63.9)

    assert {(b['sw'], b['sf']) for b in basis} == {(4000.0, 123.2618)}


def test_a_basis_keeps_the_linewidth_lcmodel_reads_in_ppm(tmp_path):
    """LCModel's FWHMBA is the basis' linewidth in ppm, a basis' linewidth is in Hz (FSL-MRS's
    basis_width, FID-A's linewidth); a .BASIS read and written again states the FWHMBA it did."""
    from readers.read_lcmodel import read_lcmodel_basis
    from writers.write_lcmodel import write_lcmodel_basis
    cows = os.path.join(SAMPLES, 'LCModel', 'BASIS_file', 'SLASER_Siemens26.basis')
    basis = read_lcmodel_basis(cows)                           # FWHMBA= 9.99999978E-03
    sf = basis[0]['sf']
    assert abs(basis[0]['linewidth'] - 0.01 * sf) < 1e-6

    def fwhmba(path):
        return float(next(l for l in open(path) if 'FWHMBA' in l).split('=')[1].strip(' ,\n'))

    write_lcmodel_basis([dict(b, linewidth=2.0) for b in basis], str(tmp_path / 'hz.BASIS'))
    write_lcmodel_basis([dict(b, linewidth=None, fwhmba=-1.0) for b in basis], str(tmp_path / 'stated.BASIS'))
    assert abs(fwhmba(tmp_path / 'hz.BASIS') - 2.0 / sf) < 1e-9
    assert fwhmba(tmp_path / 'stated.BASIS') == -1.0


def test_formats_without_a_centre_field_keep_a_stated_centre(tmp_path):
    """INSPECTOR, ProFit, MARSS and SpinWizard store no ppm centre either and are read at 4.65:
    a basis put at 3.0 ppm is shifted as for the others, so NAA stays at 2.01 ppm."""
    import numpy as np
    basis = _jmrui()
    t = np.arange(len(basis[0]['fid'])) / basis[0]['sw']
    moved = [dict(b, fid=b['fid'] * np.exp(2j * np.pi * 1.65 * b['sf'] * t), centerFreq=3.0) for b in basis]

    for tool in ('inspector', 'profit', 'marss', 'spinwizard'):
        outdir = tmp_path / tool
        outdir.mkdir()
        written = write_basis(moved, tool, 'default', str(outdir))
        back = load_basis(written, detect_format(written)['format'], sw=basis[0]['sw'], sf=basis[0]['sf'])

        for b, c in zip(basis, sorted(back, key=lambda c: [x['name'] for x in basis].index(c['name']))):
            assert np.allclose(c['fid'], b['fid'], atol=1e-4 * np.abs(b['fid']).max()), (tool, b['name'])


def test_nifti_mrs_keeps_the_echo_time(tmp_path):
    """NIfTI-MRS states the echo time (EchoTime, in s); spec2nii's raw route leaves it out, the
    writer adds it, and it is read back in ms."""
    written = write_basis(_jmrui(), 'niftimrs', 'default', str(tmp_path), te=26.0)

    assert [b['te'] for b in load_basis(written, detect_format(written)['format'])] == [26.0, 26.0]


def test_mrscloud_and_osprey_keep_the_linewidth(tmp_path):
    """MRSCloud's outA.linewidth and Osprey's BASIS.linewidth (Hz) are written and read back."""
    basis = [dict(b, linewidth=2.5) for b in _jmrui()]
    for tool in ('mrscloud', 'osprey'):
        outdir = tmp_path / tool
        outdir.mkdir()
        written = write_basis(basis, tool, 'default', str(outdir), te=26.0)
        back = load_basis(written, detect_format(written)['format'])

        assert {b['linewidth'] for b in back if b['name'] in ('NAA', 'Cr391', 'Cr')} == {2.5}, tool


def test_nifti_mrs_holds_the_fids_as_fsl_mrs_reads_them(tmp_path):
    """NIfTI-MRS FIDs are the conjugate of a basis' own (FSL-MRS reads NAA at 2.01 ppm only
    then): they are written conjugated, in the file the writer returns, and read back."""
    import nibabel as nib
    import numpy as np
    basis = _jmrui()
    written = write_basis(basis, 'niftimrs', 'default', str(tmp_path))

    assert sorted(os.listdir(tmp_path)) == ['Cr391.json', 'Cr391.nii.gz', 'NAA.json', 'NAA.nii.gz']
    for b, path in zip(basis, written):
        stored = np.asarray(nib.load(path).dataobj).squeeze()
        assert np.allclose(stored, np.conj(b['fid']), atol=1e-4 * np.abs(b['fid']).max())
    for b, c in zip(basis, load_basis(written, detect_format(written)['format'])):
        assert np.allclose(c['fid'], b['fid'], atol=1e-4 * np.abs(b['fid']).max())


def test_the_nifti_mrs_samples_are_the_simulation_the_others_hold():
    """The NIfTI-MRS sample set is the conjugated files of the conj step: read, they are the
    jMRUI sample's FIDs."""
    import numpy as np
    files = [os.path.join(SAMPLES, 'NIFTI-MRS', name) for name in ('NAA.nii.gz', 'Cr391.nii.gz')]
    for b, c in zip(_jmrui(), load_basis(files, 'NIfTI-MRS')):
        assert np.allclose(c['fid'], b['fid'], atol=1e-4 * np.abs(b['fid']).max()), b['name']


def test_a_ppmcalib_is_not_taken_as_the_centre(tmp_path):
    """MARSS files state a ppmCalib of 4.675, but their functions are at 4.65 (read there, NAA is at
    2.008 ppm): Osprey and SpinWizard state 4.65 for them, and their priors are the jMRUI sample's,
    which holds the same FIDs and states no ppmCalib."""
    folder = os.path.join(SAMPLES, 'MARSS', 'basis_functions_metab_mm_mat')
    marss = load_basis([os.path.join(folder, f'{name}.mat') for name in ('NAA', 'Cr391')], 'MARSS .mat')
    assert {b['ppmCalib'] for b in marss} == {4.675}

    (tmp_path / 'osprey').mkdir()
    written = write_basis(marss, 'osprey', 'default', str(tmp_path / 'osprey'), te=26.0)
    assert {b['centerFreq'] for b in load_basis(written, detect_format(written)['format'])} == {4.65}
    (tmp_path / 'spinwizard').mkdir()
    write_basis(marss, 'spinwizard', 'default', str(tmp_path / 'spinwizard'))
    assert float(open(tmp_path / 'spinwizard' / 'BasisRFOffset.txt').read()) == 4.65

    import re
    import numpy as np
    number = r'-?\d+\.\d+(?:[eE][-+]?\d+)?'

    def priors(basis, tool):
        """The text of the priors without paths and dates, and their numbers."""
        outdir = tmp_path / f'{tool}_{basis[0]["source"].split(".")[-1]}'
        outdir.mkdir()
        path = write_basis(basis, tool, 'default', str(outdir))[0]
        lines = [l for l in open(path).read().splitlines() if os.sep not in l and 'PK_' not in l and
                 not any(w in l.lower() for w in ('source', 'date', 'time', 'created', 'generated'))]
        return [re.sub(number, '#', l) for l in lines], [float(x) for l in lines for x in re.findall(number, l)]
    for tool in ('vespa_analysis', 'midas', 'gava', 'oxsa'):
        (text, values), (text_j, values_j) = priors(marss, tool), priors(_jmrui(), tool)
        # OXSA's RelPhase is the difference of two phases each written to 6 decimals
        assert text == text_j and np.allclose(values, values_j, rtol=1e-6,
                                              atol=1e-5 if tool == 'oxsa' else 1e-9), tool


def _osprey(basis, tmp_path):
    """The BASIS struct the Osprey writer saves: its names, FIDs, spectra and ppm axis."""
    import numpy as np
    import scipy.io as sio
    path = write_basis(basis, 'osprey', 'default', str(tmp_path), te=26.0)[0]
    b = sio.loadmat(path, squeeze_me=True, struct_as_record=False)['BASIS']
    names = [str(n) for n in np.ravel(b.name)]
    return names, np.asarray(b.fids), np.asarray(b.specs), np.asarray(b.ppm).ravel()


def test_osprey_macromolecules_are_osprey_own(tmp_path):
    """The Osprey writer builds Osprey's MM and lipids as its fit_makeBasis does: in proton
    equivalents of a third of creatine's 3.027 ppm singlet, Gaussians of decay 1 / (pi lw)."""
    import numpy as np
    folder = os.path.join(SAMPLES, 'JMRUI')
    basis = load_basis([os.path.join(folder, f'{p}.txt') for p in ('Cr391', 'CrNo391', 'NAA')], 'jMRUI .txt')
    names, fids, specs, ppm = _osprey(basis, tmp_path)
    cr, mm12 = np.real(specs[:, names.index('Cr')]), np.real(specs[:, names.index('MM12')])

    # three protons for creatine's singlet, two for MM12 (Osprey's Lorentzian area, near the sum)
    singlet = cr[(ppm > 3.027 - 0.1) & (ppm < 3.027 + 0.1)].sum()
    assert abs(mm12.sum() / singlet / (2 / 3) - 1) < 0.15
    n = 16 * len(mm12)                                                       # its width, finely sampled
    fine = np.real(np.fft.fftshift(np.fft.fft(fids[:, names.index('MM12')], n)))
    half = np.flatnonzero(fine > fine.max() / 2)
    width = (half[-1] - half[0]) * (ppm[1] - ppm[0]) * len(mm12) / n
    assert abs(width - 2 * np.sqrt(2 * np.log(2)) / 2 * 0.15) < 0.005     # decay 1 / (pi lw), lw 0.15 ppm


def test_a_macromolecule_the_set_brings_takes_osprey_own_place(tmp_path):
    """A set's MM09 is a macromolecule: it is written as MM09, in place of Osprey's Gaussian, not as
    a lipid; Osprey's own Lip09 is added beside it."""
    import numpy as np
    folder = os.path.join(SAMPLES, 'JMRUI')
    basis = load_basis([os.path.join(folder, f'{p}.txt') for p in ('Cr391', 'CrNo391', 'MM09')], 'jMRUI .txt')
    names, fids, _, _ = _osprey(basis, tmp_path)
    given = next(b for b in basis if b['name'] == 'MM09')['fid'][1:]

    assert names.count('MM09') == 1 and 'Lip09' in names
    written = fids[1:len(given) + 1, names.index('MM09')]                    # from the second point:
    scale = np.vdot(given, written) / np.vdot(given, given)                  # its DC correction
    assert np.linalg.norm(written - scale * given) < 1e-6 * np.linalg.norm(written)


def test_creatine_parts_named_as_lcmodel_keeps_them_are_summed_for_osprey(tmp_path):
    """The FSL-MRS sample names creatine's parts as LCModel keeps names, six characters (CrNo39,
    PCrNo3): they are still the parts Osprey's Cr and PCr are summed from, CH3 and CH2 in one."""
    import numpy as np
    folder = os.path.join(SAMPLES, 'FSL-MRS')
    basis = load_basis([os.path.join(folder, f'{p}.json') for p in ('Cr391', 'CrNo39', 'PCr393', 'PCrNo3')],
                       'FSL-MRS .json')
    names, _, specs, ppm = _osprey(basis, tmp_path)

    assert 'CrNo39' not in names and 'PCrNo3' not in names
    for name in ('Cr', 'PCr'):
        spec = np.real(specs[:, names.index(name)])
        ch3, ch2 = (spec[(ppm > c - 0.1) & (ppm < c + 0.1)].sum() for c in (3.027, 3.913))
        assert abs(ch3 / ch2 / (3 / 2) - 1) < 0.15, name


def test_osprey_macromolecules_can_be_left_out_as_fit_makebasis_leaves_them(tmp_path):
    """add_mm=False is fit_makeBasis's addMMFlag 0: the set's own MM and MM09 are kept as they are, none
    of Osprey's MM and lipids are added; the H2O Osprey's water fit takes is still there."""
    import numpy as np
    import scipy.io as sio
    folder = os.path.join(SAMPLES, 'JMRUI')
    basis = load_basis([os.path.join(folder, f'{p}.txt') for p in ('Cr391', 'CrNo391', 'NAA', 'MM', 'MM09')],
                       'jMRUI .txt')
    path = write_basis(basis, 'osprey', 'mat', str(tmp_path), te=26.0, add_mm=False)[0]
    b = sio.loadmat(path, squeeze_me=True, struct_as_record=False)['BASIS']
    names = [str(n) for n in np.ravel(b.name)]

    assert sorted(names) == ['Cr', 'H2O', 'MM', 'MM09', 'NAA'] and b.flags.addedMM == 0
    given = {x['name']: np.asarray(x['fid'])[1:] for x in basis}
    scale = np.vdot(given['NAA'], b.fids[1:, names.index('NAA')]) / np.vdot(given['NAA'], given['NAA'])
    for name in ('MM', 'MM09'):                                          # unscaled: as NAA is
        written = b.fids[1:, names.index(name)]
        assert np.linalg.norm(written - scale * given[name]) < 1e-6 * np.linalg.norm(written), name


def test_matlab_structs_state_their_size_as_a_double(tmp_path):
    """FID-A, MRSCloud and Osprey files hold one struct, stating sz as MATLAB's size() gives it, a
    double: Osprey's fit_makeBasis loads one struct per file, and an integer sz rounds its arithmetic."""
    import numpy as np
    import scipy.io as sio
    for tool in ('fida', 'mrscloud', 'osprey'):
        outdir = tmp_path / tool
        outdir.mkdir()
        path = write_basis(_jmrui(), tool, 'default', str(outdir), te=26.0)[0]
        saved = {k: v for k, v in sio.loadmat(path, squeeze_me=True, struct_as_record=False).items()
                 if not k.startswith('_')}
        assert len(saved) == 1, tool
        assert np.asarray(next(iter(saved.values())).sz).dtype == float, tool


def test_a_profit_file_with_one_metabolite_is_profit(tmp_path):
    """A 1x1 cell loads as the array it holds: ProFit is told by its fields, not by a cell array."""
    import numpy as np
    naa = load_basis([os.path.join(SAMPLES, 'JMRUI', 'NAA.txt')], 'jMRUI .txt')
    written = write_basis(naa, 'profit', 'default', str(tmp_path))
    found = detect_format(written)

    assert found['format'] == 'ProFit .mat'
    assert np.allclose(load_basis(written, found['format'])[0]['fid'], naa[0]['fid'])


def test_a_writer_that_writes_nothing_says_so(tmp_path):
    """The writers that write a file per function carry on past one that fails, but not past all:
    into a folder that cannot be written, none is written, and they say so."""
    import pytest
    for tool in ('fida', 'marss', 'mrscloud', 'fsLmrs', 'lcmodel'):
        outdir = tmp_path / tool
        outdir.mkdir()
        outdir.chmod(0o500)
        try:
            with pytest.raises(RuntimeError, match='Nothing written'):
                write_basis(_jmrui(), tool, 'raw' if tool == 'lcmodel' else 'default', str(outdir))
        finally:
            outdir.chmod(0o700)


def test_the_osprey_sample_states_the_centre_it_is_at(tmp_path):
    """The OSPREY sample is the COWS MARSS set, at 4.65 ppm: it states 4.65, and read there NAA is
    at 2.008 ppm."""
    import numpy as np
    naa = next(b for b in load_basis([os.path.join(SAMPLES, 'OSPREY', 'ospreybasis_marss_to_cows_26.mat')],
                                     'Osprey .mat') if b['name'] == 'NAA')
    fid, sw, sf = naa['fid'], naa['sw'], naa['sf']
    n = 16 * len(fid)
    spec = np.real(np.fft.fftshift(np.fft.ifft(fid * np.exp(-np.pi * np.arange(len(fid)) / sw), n)))
    ppm = naa['centerFreq'] - np.fft.fftshift(np.fft.fftfreq(n, 1 / sw)) / sf
    window = (ppm > 1.9) & (ppm < 2.1)

    assert naa['centerFreq'] == 4.65
    assert abs(ppm[window][np.argmax(spec[window])] - 2.008) < 0.003


def test_the_osprey_sample_holds_creatine_and_water_as_osprey_does(tmp_path):
    """The OSPREY sample is the MARSS set as the Osprey writer writes it: its Cr and PCr hold the CH3
    and the CH2 (three protons to two), and its H2O is below its metabolites."""
    import numpy as np
    basis = {b['name']: b for b in load_basis([os.path.join(SAMPLES, 'OSPREY', 'ospreybasis_marss_to_cows_26.mat')],
                                              'Osprey .mat')}
    def spectrum(b):
        fid, sw, sf = b['fid'], b['sw'], b['sf']
        ppm = b['centerFreq'] - np.fft.fftshift(np.fft.fftfreq(len(fid), 1 / sw)) / sf
        return np.real(np.fft.fftshift(np.fft.ifft(fid * np.exp(-np.pi * np.arange(len(fid)) / sw)))), ppm

    for name in ('Cr', 'PCr'):
        spec, ppm = spectrum(basis[name])
        ch3, ch2 = (spec[(ppm > c - 0.1) & (ppm < c + 0.1)].sum() for c in (3.027, 3.913))
        assert abs(ch3 / ch2 / (3 / 2) - 1) < 0.15, name
    peak = lambda name: spectrum(basis[name])[0].max()
    assert peak('H2O') < max(peak(n) for n in ('NAA', 'Cr', 'PCh'))


def test_text_formats_keep_a_basis_in_double_precision(tmp_path):
    """.BASIS, .raw, NIfTI-MRS (written by spec2nii from a .raw), jMRUI and SpinWizard files are
    text, or made from text: written with 17 digits, they give back the FIDs they were given."""
    import numpy as np
    folder = os.path.join(SAMPLES, 'MARSS', 'basis_functions_metab_mm_mat')
    basis = load_basis([os.path.join(folder, f'{name}.mat') for name in ('NAA', 'GPC')], 'MARSS .mat')
    for tool, fmt in [('lcmodel', 'basis'), ('lcmodel', 'raw'), ('niftimrs', 'default'),
                      ('jmrui', 'default'), ('spinwizard', 'default')]:
        written = write_basis(basis, tool, fmt, str(tmp_path / f'{tool}_{fmt}'), te=26.0)
        back = {b['name']: b for b in load_basis(written, detect_format(written)['format'],
                                                  sw=basis[0]['sw'], sf=basis[0]['sf'])}
        for b in basis:
            a, c = np.asarray(b['fid']), np.asarray(back[b['name']]['fid'])[:len(b['fid'])]
            assert np.linalg.norm(c - a) < 1e-12 * np.linalg.norm(a), (tool, fmt, b['name'])


def test_osprey_water_is_left_out_of_other_formats_as_osprey_exports_them(tmp_path):
    """Osprey's H2O is for its water-reference fit: kept in an Osprey .mat, and left out of the
    formats other tools fit with, as Osprey's io_writelcmBASIS and io_writeTarquinBASIS leave it out."""
    basis = load_basis([os.path.join(SAMPLES, 'OSPREY', 'ospreybasis_marss_to_cows_26.mat')], 'Osprey .mat')
    assert 'H2O' in [b['name'] for b in basis]

    lcm = write_basis(basis, 'lcmodel', 'basis', str(tmp_path / 'lcm'))
    fsl = write_basis(basis, 'fsLmrs', 'default', str(tmp_path / 'fsl'))
    osp_nifti = write_basis(basis, 'osprey', 'nifti', str(tmp_path / 'osp_nifti'))
    osp = write_basis(basis, 'osprey', 'mat', str(tmp_path / 'osp'), te=26.0)

    for written in (lcm, fsl, osp_nifti):
        names = [b['name'] for b in load_basis(written, detect_format(written)['format'])]
        assert 'H2O' not in names and 'NAA' in names
    assert 'H2O' in [b['name'] for b in load_basis(osp, 'Osprey .mat')]


def test_a_set_written_for_osprey_reads_back_as_it_was_given(tmp_path):
    """Osprey's model changes a set (its names, creatine's parts summed, the first point, the scale)
    and adds to it (MM, lipids, water); the .mat holds the set as given, which is read back."""
    import numpy as np
    folder = os.path.join(SAMPLES, 'MARSS', 'basis_functions_metab_mm_mat')
    files = [os.path.join(folder, f) for f in sorted(os.listdir(folder)) if f.endswith('.mat')]
    basis = load_basis(files, 'MARSS .mat')
    for add_mm in (True, False):
        written = write_basis(basis, 'osprey', 'mat', str(tmp_path / str(add_mm)), te=26.0, add_mm=add_mm)
        back = load_basis(written, 'Osprey .mat')

        assert [b['name'] for b in back] == [b['name'] for b in basis]
        for a, b in zip(basis, back):
            assert np.array_equal(np.asarray(a['fid']), b['fid']) and not b['water_ref']


def test_an_osprey_set_written_for_osprey_again_keeps_its_water_out_of_other_formats(tmp_path):
    """An Osprey set's H2O, read back from the .mat BasiC wrote of it, is still Osprey's water."""
    basis = load_basis([os.path.join(SAMPLES, 'OSPREY', 'ospreybasis_marss_to_cows_26.mat')], 'Osprey .mat')
    mat = write_basis(basis, 'osprey', 'mat', str(tmp_path / 'osp'), te=26.0)
    back = load_basis(mat, 'Osprey .mat')
    assert [b['water_ref'] for b in back] == [b['water_ref'] for b in basis] and any(b['water_ref'] for b in back)

    lcm = write_basis(back, 'lcmodel', 'basis', str(tmp_path / 'lcm'))
    assert 'H2O' not in [b['name'] for b in load_basis(lcm, 'LCModel .BASIS')]


def test_an_osprey_set_changed_since_basic_wrote_it_is_read_as_osprey_holds_it(tmp_path):
    import numpy as np
    import scipy.io as sio
    written = write_basis(_jmrui(), 'osprey', 'mat', str(tmp_path), te=26.0)[0]
    mat = sio.loadmat(written)
    mat['BASIS'][0, 0]['fids'][:, 0] *= 2
    sio.savemat(written, {'BASIS': mat['BASIS']})

    back = {b['name']: b for b in load_basis([written], 'Osprey .mat')}
    assert {'NAA', 'Cr', 'H2O'} <= set(back) and back['H2O']['water_ref']


def test_fida_and_mrscloud_files_state_the_centre_they_are_at(tmp_path):
    """Osprey's fit_makeBasis puts a file's ppm axis at the centerFreq it states, and at 3.0 ppm
    when it states none: BasiC's FID-A and MRSCloud files are at 4.65 and say so, in centerFreq
    and in their ppm axis. A file stating a centerFreq is read at it, as Osprey reads it."""
    import numpy as np
    import scipy.io as sio
    naa = load_basis([os.path.join(SAMPLES, 'JMRUI', 'NAA.txt')], 'jMRUI .txt')
    for tool, struct in (('fida', 'out'), ('mrscloud', 'outA')):
        written = write_basis(naa, tool, 'default', str(tmp_path / tool))
        out = sio.loadmat(written[0], simplify_cells=True)[struct]
        assert out['centerFreq'] == 4.65, tool
        assert abs((out['ppm'][0] + out['ppm'][-1]) / 2 - 4.65) < 1e-3, tool
        assert load_basis(written, detect_format(written)['format'])[0]['centerFreq'] == 4.65, tool

        out['centerFreq'] = 3.0
        sio.savemat(written[0], {struct: out})
        assert load_basis(written, detect_format(written)['format'])[0]['centerFreq'] == 3.0, tool


def test_fsl_mrs_files_hold_the_conjugate_as_fsl_mrs_stores_it(tmp_path):
    """FSL-MRS saves the conjugate of its FID and conjugates it back on reading (fsl_io's
    write_fsl_basis_file and readFSLBasisFiles); unconjugated, FSL-MRS finds NAA at 7.3 ppm. The
    FSL-MRS sample, converted from SLASER_Siemens26.basis, reads as that basis, and BasiC's
    files hold the conjugate of the FID it has."""
    import json
    import numpy as np
    sample = load_basis([os.path.join(SAMPLES, 'FSL-MRS', 'NAA.json')], 'FSL-MRS .json')[0]
    lcmodel = [b for b in load_basis([os.path.join(SAMPLES, 'LCModel', 'BASIS_file', 'SLASER_Siemens26.basis')],
                                     'LCModel .BASIS') if b['name'] == 'NAA'][0]
    fid, got = lcmodel['fid'], sample['fid'][:len(lcmodel['fid'])]          # the sample is zero-filled
    scale = np.vdot(fid, got) / np.vdot(fid, fid)
    assert np.linalg.norm(got - scale * fid) < 1e-9 * np.linalg.norm(got)

    naa = load_basis([os.path.join(SAMPLES, 'JMRUI', 'NAA.txt')], 'jMRUI .txt')
    written = write_basis(naa, 'fsLmrs', 'default', str(tmp_path))
    with open(written[0]) as f:
        stored = json.load(f)['basis']
    assert np.allclose(np.array(stored['basis_re']) + 1j * np.array(stored['basis_im']), np.conj(naa[0]['fid']))
    assert np.allclose(load_basis(written, 'FSL-MRS .json')[0]['fid'], naa[0]['fid'])


def test_each_writer_runs_from_the_command_line(tmp_path):
    """Each writer's docstring gives its own command line: run as a script from another folder,
    it finds the writers package and says how it is used, and OXSA's writes a function it has
    no literature lines for."""
    import glob
    import shutil
    import subprocess
    for path in sorted(glob.glob(os.path.join(ROOT, 'basis_converter', 'writers', 'write_*.py'))):
        run = subprocess.run([sys.executable, path], cwd=str(tmp_path), capture_output=True, text=True)
        assert run.returncode == 1 and run.stdout.startswith('Usage:'), (path, run.stderr[-300:])

    unknown = tmp_path / 'unknown'
    unknown.mkdir()
    shutil.copy(os.path.join(SAMPLES, 'MARSS', 'basis_functions_metab_mm_mat', 'NAA.mat'), unknown / 'Xyz.mat')
    oxsa = os.path.join(ROOT, 'basis_converter', 'writers', 'write_oxsa.py')
    run = subprocess.run([sys.executable, oxsa, str(unknown), str(tmp_path)], cwd=str(tmp_path),
                         capture_output=True, text=True)
    assert run.returncode == 0 and list(tmp_path.glob('*.m')), run.stderr[-300:]


def test_the_echo_time_written_is_the_one_given_else_the_one_stated(tmp_path):
    """Every file written for a set holds the echo time given, else the one the set states; a .raw
    states none where the set states none (FID-A and MRSCloud 0), and a .BASIS or Osprey's .mat,
    which need one, 30. Without one given, each .raw keeps its own function's."""
    import scipy.io as sio
    naa = load_basis([os.path.join(SAMPLES, 'JMRUI', 'NAA.txt')], 'jMRUI .txt')
    cases = (('stated', [dict(b, te=26.0) for b in naa], None, 26.0, 26.0),
             ('given', [dict(b, te=26.0) for b in naa], 35.0, 35.0, 35.0),
             ('none', [dict(b, te=None) for b in naa], None, None, 30.0))
    for case, basis, te, raw_te, set_te in cases:
        out = tmp_path / case
        written = write_basis(basis, 'lcmodel', 'both', str(out / 'lcmodel'), te=te)
        raw = [p for p in written if p.endswith('.raw')]
        assert [b['te'] for b in load_basis(raw, 'LCModel .raw')] == [raw_te], case
        assert load_basis([out / 'lcmodel' / 'basis_set.BASIS'], 'LCModel .BASIS')[0]['te'] == set_te, case

        written = write_basis(basis, 'osprey', 'mat', str(out / 'osprey'), te=te)
        assert sio.loadmat(written[0], simplify_cells=True)['BASIS']['te'] == set_te, case
        for tool in ('fida', 'mrscloud'):
            written = write_basis(basis, tool, 'default', str(out / tool), te=te)
            back = load_basis(written, detect_format(written)['format'])[0]['te']
            assert (back or None) == raw_te, (case, tool)      # they hold te 0 where none is stated

    mixed = [dict(naa[0], te=26.0), dict(naa[0], name='Cr', te=35.0)]
    written = write_basis(mixed, 'lcmodel', 'raw', str(tmp_path / 'mixed'))
    assert sorted(b['te'] for b in load_basis(written, 'LCModel .raw')) == [26.0, 35.0]


def test_an_osprey_set_stating_two_centres_is_written_at_one(tmp_path):
    """MRSCloud's v2 sets hold their metabolites with 4.68 ppm at 0 Hz and their MM at centerFreq.
    Osprey's .mat states one centre, the first function's, so each function is moved to it:
    NAA and MM09 stay where they are on the .mat's ppm axis, and in the set given, read back. A
    function stating no centre is at 4.65 and is moved too."""
    import numpy as np
    import scipy.io as sio
    sw, n, Bo = 4000.0, 2048, 2.89
    sf, t = Bo * 42.577, np.arange(n) / sw
    line = lambda ppm, centre: np.exp(2j * np.pi * (ppm - centre) * sf * t - t / 0.05)
    path = str(tmp_path / 'v2.mat')
    sio.savemat(path, {'BASIS': {
        'fids': np.stack([line(2.008, 4.68), line(0.91, 3.0)], axis=1),
        'name': np.array(['NAA', 'MM09'], dtype=object), 'nMets': 1, 'spectralwidth': sw, 'n': n,
        'Bo': Bo, 'te': 30.0, 'seq': 'unedited', 'centerFreq': 3.0, 'mrsCloudVersion': 'v2'}})
    written = write_basis(load_basis([path], 'Osprey .mat'), 'osprey', 'mat', str(tmp_path / 'out'), add_mm=False)

    b = sio.loadmat(written[0], squeeze_me=True, struct_as_record=False)['BASIS']
    names, specs, ppm = [str(x) for x in np.ravel(b.name)], np.asarray(b.specs), np.asarray(b.ppm).ravel()
    for name, want in (('NAA', 2.008), ('MM09', 0.91)):
        assert abs(ppm[np.argmax(np.abs(specs[:, names.index(name)]))] - want) < 0.02, name
    for core in load_basis(written, 'Osprey .mat'):                     # the set given, read back
        freq = np.fft.fftshift(np.fft.fftfreq(len(core['fid']), 1 / core['sw']))
        peak = core['centerFreq'] + freq[np.argmax(np.abs(np.fft.fftshift(np.fft.fft(core['fid']))))] / core['sf']
        assert abs(peak - {'NAA': 2.008, 'MM09': 0.91}[core['name']]) < 0.02, core['name']

    # a function stating no centre is at 4.65, beside one stating 3.0
    mixed = [{'name': 'NAA', 'fid': line(2.008, 3.0), 'sw': sw, 'sf': sf, 'n': n, 'centerFreq': 3.0},
             {'name': 'Cr', 'fid': line(3.027, 4.65), 'sw': sw, 'sf': sf, 'n': n}]
    written = write_basis(mixed, 'osprey', 'mat', str(tmp_path / 'mixed'), add_mm=False)
    b = sio.loadmat(written[0], squeeze_me=True, struct_as_record=False)['BASIS']
    names, specs, ppm = [str(x) for x in np.ravel(b.name)], np.asarray(b.specs), np.asarray(b.ppm).ravel()
    assert abs(ppm[np.argmax(np.abs(specs[:, names.index('Cr')]))] - 3.027) < 0.02


def test_mrscloud_files_hold_the_fields_fit_makebasis_reads(tmp_path):
    """Osprey's fit_makeBasis reads FID-A's te and seq of each file, which MRSCloud's own files
    leave out: written, an MRSCloud file holds them, and its echo time is read back."""
    import scipy.io as sio
    basis = [dict(b, te=26.0) for b in _jmrui()]
    written = write_basis(basis, 'mrscloud', 'default', str(tmp_path))

    out = sio.loadmat(written[0], simplify_cells=True)['outA']
    assert out['seq'] == 'unedited' and out['te'] == 26.0
    assert [b['te'] for b in load_basis(written, 'MRSCloud .mat')] == [b['te'] for b in basis]


def test_a_basis_states_the_sequence_lcmodel_truncated_in_full(tmp_path):
    """The sample .BASIS states SEQ='sLASE', cut short as LCModel cuts names; it is read as
    sLASER, so screen 3 shows it and the FID-A and MRSCloud files fit_makeBasis reads hold it."""
    import scipy.io as sio
    cows = [os.path.join(SAMPLES, 'LCModel', 'BASIS_file', 'SLASER_Siemens26.basis')]
    assert detect_format(cows)['data']['seq'] == 'sLASER'
    basis = load_basis(cows, 'LCModel .BASIS')
    assert {b['seq'] for b in basis} == {'sLASER'}

    written = write_basis(basis[:1], 'mrscloud', 'default', str(tmp_path), seq=basis[0]['seq'])
    assert sio.loadmat(written[0], simplify_cells=True)['outA']['seq'] == 'sLASER'


def test_fsl_mrs_keeps_the_echo_time_and_nifti_mrs_the_linewidth(tmp_path):
    """Neither format has a field for them: FSL-MRS's meta, which FSL-MRS does not read, holds the
    echo time (EchoTime, s, as NIfTI-MRS states it), and a user-defined NIfTI-MRS key the linewidth."""
    basis = [dict(b, te=26.0, linewidth=2.5) for b in _jmrui()]
    for tool, field in (('fsLmrs', 'te'), ('niftimrs', 'linewidth')):
        written = write_basis(basis, tool, 'default', str(tmp_path / tool), te=26.0)
        back = load_basis(written, detect_format(written)['format'])
        assert [b[field] for b in back] == [b[field] for b in basis], tool


def test_osprey_water_is_two_protons_beside_the_metabolites(tmp_path):
    """Osprey scales the data to its metabolites' peak, and its own sets peak in a metabolite: the
    synthetic H2O is two protons by area, below them, so the set is not normalised to the water."""
    import numpy as np
    folder = os.path.join(SAMPLES, 'JMRUI')
    basis = load_basis([os.path.join(folder, f'{p}.txt') for p in ('Cr391', 'CrNo391', 'NAA', 'GPC')], 'jMRUI .txt')
    names, _, specs, ppm = _osprey(basis, tmp_path)
    peaks = np.max(np.real(specs), axis=0)
    water, cr = np.real(specs[:, names.index('H2O')]), np.real(specs[:, names.index('Cr')])

    assert max(peaks[names.index(n)] for n in ('Cr', 'NAA', 'GPC')) == peaks.max() == 1.0
    singlet = cr[(ppm > 3.027 - 0.1) & (ppm < 3.027 + 0.1)].sum()            # three protons
    assert abs(water.sum() / singlet / (2 / 3) - 1) < 0.15


def test_every_priors_writer_writes_a_set_of_broad_humps_alone(tmp_path):
    """jMRUI's MM alone, or beside Glu: no singlet, and lines of a lineshape so broad that their
    Gram matrix rounds short of positive definite while their frequencies are refined. Each
    priors writer still writes its file."""
    folder = os.path.join(SAMPLES, 'JMRUI')
    for names in (('MM',), ('Glu', 'MM')):
        basis = load_basis([os.path.join(folder, f'{n}.txt') for n in names], 'jMRUI .txt')
        for tool in ('vespa_analysis', 'gava', 'midas', 'pyamares', 'oxsa'):
            written = write_basis(basis, tool, 'default', str(tmp_path / f'{tool}-{len(names)}'))
            assert written and all(os.path.getsize(p) for p in written), (names, tool)


def test_a_set_without_singlets_takes_the_lineshape_of_its_narrowest_lines(tmp_path):
    """Lac and MM09 hold no singlet: the lineshape is Lac's 1 Hz, not a median between Lac and
    the 10 Hz MM09, so MM09 is left out as a hump and Lac is rebuilt from a few lines."""
    from writers.lines import set_lineshape
    from writers.write_vespa import vespa_lines
    folder = os.path.join(SAMPLES, 'JMRUI')
    basis = load_basis([os.path.join(folder, f'{n}.txt') for n in ('Lac', 'MM09')], 'jMRUI .txt')
    shape = set_lineshape(basis)
    assert shape['source'] == 'lines' and abs(shape['lorentz'] - 1.0) < 0.01

    found = vespa_lines(basis)
    assert found['humps'] == ['MM09'] and not found['short']
    assert [name for name, _, _ in found['metabolites']] == ['Lac'] and len(found['metabolites'][0][1]) < 30


def test_the_pyamares_prior_names_its_reference_a_singlet_only_where_it_is_one(tmp_path):
    """The comment says which peak carries the set's free shift, linewidth and phase: sI, a
    singlet, in the jMRUI sample; Glu alone has none, and its reference is called a metabolite."""
    folder = os.path.join(SAMPLES, 'JMRUI')

    def comment(names):
        basis = load_basis([os.path.join(folder, f'{n}.txt') for n in names], 'jMRUI .txt')
        written = write_basis(basis, 'pyamares', 'default', str(tmp_path / '-'.join(names)))
        return open(written[0]).read()

    assert '(the largest line of the largest singlet)' in comment(('NAA', 'sI'))
    assert '(the largest line of the largest metabolite)' in comment(('Glu',))


def test_a_function_its_lines_rebuild_only_in_part_is_named_with_their_number(tmp_path):
    """jMRUI's MM alone takes its own 950 Hz lineshape and is rebuilt only in part, by far fewer
    lines than the most allowed: each priors file names it with the lines it has, not the most."""
    from writers.write_vespa import vespa_lines
    basis = load_basis([os.path.join(SAMPLES, 'JMRUI', 'MM.txt')], 'jMRUI .txt')
    found = vespa_lines(basis)
    (name, lines, residual), = found['metabolites']
    assert found['short'] == [f'MM ({len(lines)} lines, {residual:.0e})'] and len(lines) < 100

    for tool in ('vespa_analysis', 'gava', 'pyamares', 'oxsa'):
        written = write_basis(basis, tool, 'default', str(tmp_path / tool))
        text = open(written[0]).read()
        assert f'Rebuilt only to: {found["short"][0]}.' in text and '1000 lines' not in text, tool


def test_a_set_of_functions_without_signal_says_so_for_each_priors_writer(tmp_path):
    """A function that is all zeros has no lines; beside others it is skipped, but a set of
    only such functions has nothing to write, and the writer says why rather than numpy."""
    import numpy as np
    import pytest
    basis = load_basis([os.path.join(SAMPLES, 'JMRUI', 'Lac.txt')], 'jMRUI .txt')
    basis[0]['fid'] = np.zeros_like(basis[0]['fid'])
    for tool in ('vespa_analysis', 'gava', 'midas', 'pyamares', 'oxsa'):
        with pytest.raises(RuntimeError, match='No basis function holds any signal'):
            write_basis(basis, tool, 'default', str(tmp_path / tool))


def test_raw_files_that_cannot_be_read_are_called_unreadable(tmp_path):
    """Each .raw label reads with the same reader: when it fails, screen 1 says why rather than
    offering a choice of names that would fail the same way on screen 4."""
    for name in ('a.raw', 'b.raw'):
        (tmp_path / name).write_text('not a basis function\n')
    result = detect_format(sorted(str(p) for p in tmp_path.iterdir()))

    assert result['status'] == 'failed' and result['format'] == 'Unreadable'
    assert result['message'].startswith('Could not read .raw file')


def test_a_singlet_without_signal_is_skipped_as_any_function_without_signal(tmp_path):
    """Cr391 all zeros has no line to measure the lineshape on: beside Lac it is skipped and the
    lineshape is Lac's; alone, the writer says there is no signal."""
    import numpy as np
    import pytest
    folder = os.path.join(SAMPLES, 'JMRUI')
    basis = load_basis([os.path.join(folder, f'{n}.txt') for n in ('Cr391', 'Lac')], 'jMRUI .txt')
    basis[0]['fid'] = np.zeros_like(basis[0]['fid'])
    for tool in ('vespa_analysis', 'gava', 'midas', 'pyamares', 'oxsa'):
        written = write_basis(basis, tool, 'default', str(tmp_path / tool))
        assert 'Lac' in open(written[0]).read(), tool
    with pytest.raises(RuntimeError, match='No basis function holds any signal'):
        write_basis(basis[:1], 'vespa_analysis', 'default', str(tmp_path / 'alone'))


def test_each_priors_file_names_the_functions_it_found_no_lines_for(tmp_path):
    """Cr391 all zeros beside Lac is left out of every priors file, and each one says so."""
    import numpy as np
    folder = os.path.join(SAMPLES, 'JMRUI')
    basis = load_basis([os.path.join(folder, f'{n}.txt') for n in ('Cr391', 'Lac')], 'jMRUI .txt')
    basis[0]['fid'] = np.zeros_like(basis[0]['fid'])
    for tool in ('vespa_analysis', 'gava', 'midas', 'pyamares', 'oxsa'):
        written = write_basis(basis, tool, 'default', str(tmp_path / tool))
        text = open(written[0]).read()
        assert 'Cr391' in text and ('No lines found for: Cr391' in text or
                                    'Cr391: no lines found' in text), tool
