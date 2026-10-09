"""
tests/test_lcmodel_basis.py
LCModel .BASIS read/write against LCModel's own files.

A .BASIS stores each metabolite as a spectrum, fft(conj(FID)); a .raw stores conj(FID)
(FSL-MRS reads both so). The jMRUI and LCModel .raw samples are the same simulation, so
a jMRUI basis written to .BASIS and read back must be that .raw's FID.

Run from the repository root: python -m pytest tests
"""

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAMPLES = os.path.join(ROOT, 'sample_basis_sets')
sys.path.insert(0, os.path.join(ROOT, 'basis_converter'))

from readers.read_jmrui import read_jmrui_file          # noqa: E402
from readers.read_lcmodel import read_lcmodel_basis     # noqa: E402
from writers.write_lcmodel import write_lcmodel_basis   # noqa: E402


def _raw_fid(path):
    """The FID in an LCModel .raw: the pairs after its last $END, conjugated."""
    with open(path) as f:
        values = np.array(f.read().split('$END')[-1].split(), dtype=float)
    return np.conj(values[0::2] + 1j * values[1::2])


def test_a_basis_written_and_read_back_is_lcmodels_own(tmp_path):
    naa = read_jmrui_file(os.path.join(SAMPLES, 'JMRUI', 'NAA.txt'))
    out = str(tmp_path / 'naa.BASIS')
    write_lcmodel_basis([naa], out)

    fid = read_lcmodel_basis(out)[0]['fid']
    lcmodel = _raw_fid(os.path.join(SAMPLES, 'LCModel', 'raw_files', 'NAA.raw'))[:len(fid)]

    scale = np.vdot(lcmodel, fid) / np.vdot(lcmodel, lcmodel)
    assert np.linalg.norm(fid - scale * lcmodel) < 1e-4 * np.linalg.norm(fid)


def test_a_basis_survives_a_round_trip(tmp_path):
    first, second = str(tmp_path / 'a.BASIS'), str(tmp_path / 'b.BASIS')
    write_lcmodel_basis([read_jmrui_file(os.path.join(SAMPLES, 'JMRUI', name))
                         for name in ('NAA.txt', 'Cr391.txt')], first)
    write_lcmodel_basis(read_lcmodel_basis(first), second)

    # equal to print precision (6E13.5): the zero-filled half's rounding moves last digits
    for a, b in zip(read_lcmodel_basis(first), read_lcmodel_basis(second)):
        assert np.linalg.norm(a['fid'] - b['fid']) < 1e-5 * np.linalg.norm(a['fid'])


def test_a_basis_is_read_by_its_detected_format_without_the_gui():
    from core.detect import detect_format
    from core.load import load_basis

    files = [os.path.join(SAMPLES, 'JMRUI', name) for name in ('NAA.txt', 'Cr391.txt')]
    found = detect_format(files)
    basis = load_basis(files, found['format'])

    assert 'jMRUI' in found['format'] and [m['name'] for m in basis] == ['NAA', 'Cr391']
    assert 'tkinter' not in sys.modules


def _peak_ppm(core, centre):
    fid = np.asarray(core['fid']).ravel()
    freq = np.fft.fftshift(np.fft.fftfreq(len(fid), 1 / core['sw']))
    return centre + freq[np.argmax(np.abs(np.fft.fftshift(np.fft.fft(fid))))] / core['sf']


def test_a_basis_stating_another_centre_keeps_its_peaks_where_they_are(tmp_path):
    """Osprey's own basis sets put 3.0 ppm at 0 Hz; LCModel's is written for 4.65."""
    naa = read_jmrui_file(os.path.join(SAMPLES, 'JMRUI', 'NAA.txt'))      # NAA at 2.01
    moved = dict(naa, fid=naa['fid'] * np.exp(2j * np.pi * 1.65 * naa['sf']
                                                 * np.arange(len(naa['fid'])) / naa['sw']),
                 centerFreq=3.0)                                       # the same, at 3.0
    out = str(tmp_path / 'naa.BASIS')
    write_lcmodel_basis([moved, dict(naa, name='same')], out)

    back = read_lcmodel_basis(out)

    assert abs(_peak_ppm(moved, 3.0) - _peak_ppm(naa, 4.65)) < 0.01
    assert [round(_peak_ppm(b, 4.65), 2) for b in back] == [round(_peak_ppm(naa, 4.65), 2)] * 2
    from readers.read_lcmodel import read_lcmodel_raw
    from writers.write_lcmodel import write_lcmodel_raw
    write_lcmodel_raw(moved, str(tmp_path / 'naa.raw'))
    raw = read_lcmodel_raw(str(tmp_path / 'naa.raw'))
    assert round(_peak_ppm(raw, 4.65), 2) == round(_peak_ppm(naa, 4.65), 2)
