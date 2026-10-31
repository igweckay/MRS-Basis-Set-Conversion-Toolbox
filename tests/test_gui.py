"""
tests/test_gui.py
The screens, without anyone at the window.

The step bar: BaseScreen.__init__ runs with the Tk frame and the widget builders
stubbed out, so no display is needed. The parameters summary: screen 3 is built on a
hidden Tk root and typed into; the test is skipped when Tk cannot open a window.
Osprey's macromolecules: screen 2 is built on a hidden Tk root and its tools picked,
skipped likewise; screen 4 writes the .mat with its widgets stubbed out, as for the
step bar. Several basis sets in one folder: screen 1 is built on a hidden Tk root, given
each folder and one of its sets picked, skipped likewise.

The screens are imported in a separate Python, so that tkinter stays out of this
process and the other tests can still check that core runs without it.

Run from the repository root: python -m pytest tests
"""

import importlib.util
import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

STEPS_SEEN = """
import tkinter as tk
from base_screen import BaseScreen
from screen1_upload import Screen1Upload
from screen2_output import Screen2Output
from screen3_params import Screen3Params
from screen4_convert import Screen4Convert

seen = []
tk.Frame.__init__ = lambda self, *args, **kwargs: None
BaseScreen._build_header = lambda self: None
BaseScreen._build_step_bar = lambda self: seen.append(self.current_step)
for screen in (Screen1Upload, Screen2Output, Screen3Params, Screen4Convert):
    screen.build_body = lambda self: None
    screen(None, {}, None)
print(getattr(BaseScreen, 'current_step', None), *seen)
"""


@pytest.mark.skipif(importlib.util.find_spec('tkinter') is None, reason='no tkinter')
def test_the_step_bar_highlights_each_screens_own_step():
    run = subprocess.run([sys.executable, '-c', STEPS_SEEN], capture_output=True, text=True,
                         cwd=os.path.join(ROOT, 'basis_converter'))

    assert run.returncode == 0, run.stderr
    assert run.stdout.split() == ['1', '1', '2', '3', '4']


SUMMARY_SEEN = """
import json
import tkinter as tk
from tkinter import ttk

try:
    root = tk.Tk()
except tk.TclError as error:
    print(error)
    raise SystemExit(3)
root.withdraw()

from screen3_params import Screen3Params


def widgets(widget):
    yield widget
    for child in widget.winfo_children():
        yield from widgets(child)


def summary(name):
    # The summary card is built last; each of its cells holds the name, then the value.
    cells = [w for w in widgets(screen) if isinstance(w, tk.Frame)
             and [c.cget('text') for c in w.winfo_children()
                  if isinstance(c, tk.Label)][:1] == [name]]
    return cells[-1].winfo_children()[1].cget('text')


def seen_now():
    return [summary('Echo time (TE)'), summary('Sequence')]


# Built as the App builds it, then shown with what detect_format finds in
# sample_basis_sets/LCModel/raw_files: no SW, Larmor frequency or TE, no sequence.
state = {}
screen = Screen3Params(root, state, root)
state.update(core={'sw': None, 'sf': None, 'Bo': None},
             missing_fields=['sw', 'sf', 'te'])
screen.on_show()

te_var = str(screen._entries['te'])
te = next(w for w in widgets(screen)
          if isinstance(w, tk.Entry) and str(w.cget('textvariable')) == te_var)
seq = next(w for w in widgets(screen) if isinstance(w, ttk.Combobox))

seen = {'shown': seen_now()}
te.delete(0, 'end')
seen['cleared'] = seen_now()
te.insert('end', '2')
seen['typed 2'] = seen_now()
te.insert('end', '6')
seen['typed 26'] = seen_now()
seq.current([str(v) for v in seq.cget('values')].index('sLASER'))
seen['picked sLASER'] = seen_now()
print(json.dumps(seen))
"""


@pytest.mark.skipif(importlib.util.find_spec('tkinter') is None, reason='no tkinter')
def test_the_parameters_summary_follows_what_the_user_types_and_picks():
    run = subprocess.run([sys.executable, '-c', SUMMARY_SEEN], capture_output=True, text=True,
                         cwd=os.path.join(ROOT, 'basis_converter'))
    if run.returncode == 3:
        pytest.skip(f'Tk cannot open a window: {run.stdout.strip()}')

    assert run.returncode == 0, run.stderr
    assert json.loads(run.stdout) == {
        'shown': ['30', 'unedited'],
        'cleared': ['—', 'unedited'],
        'typed 2': ['2', 'unedited'],
        'typed 26': ['26', 'unedited'],
        'picked sLASER': ['26', 'sLASER'],
    }


ADD_MM_SHOWN = """
import json
import tkinter as tk

try:
    root = tk.Tk()
except tk.TclError as error:
    print(error)
    raise SystemExit(3)
root.withdraw()

from core.formats import ALL_TOOLS
from screen2_output import Screen2Output


def widgets(widget):
    yield widget
    for child in widget.winfo_children():
        yield from widgets(child)


def shown(widget):
    # packed, and so is every frame it sits in, up to the scrolled frame
    while widget is not screen.scroll_frame:
        if not widget.winfo_manager():
            return False
        widget = widget.nametowidget(widget.winfo_parent())
    return True


def checkbox():
    box, = [w for w in widgets(screen) if isinstance(w, tk.Checkbutton)]
    return box


def pick(tool_id, fmt=None):
    screen._on_tool_selected(next(t for t in ALL_TOOLS if t['id'] == tool_id))
    if fmt:
        next(w for w in widgets(screen) if isinstance(w, tk.Radiobutton)
             and w.cget('value') == fmt).invoke()


def seen_now():
    box = checkbox()
    if not shown(box):
        return 'hidden'
    return {'checked': root.getboolean(box.getvar(box.cget('variable'))), 'text': box.cget('text')}


class App:                               # what the screen asks of main.App
    def quit(self):
        pass

    def show_screen(self, name):
        pass


state = {}
screen = Screen2Output(root, state, App())
seen = {'none picked': seen_now()}
pick('osprey')
seen['osprey'] = seen_now()
checkbox().invoke()
pick('osprey')                           # its card clicked again: no tool picked
seen['osprey dropped'] = seen_now()
pick('osprey', 'mat')
seen['osprey again'] = seen_now()
checkbox().invoke()
seen['unchecked'] = seen_now()
screen._on_next()
seen['next, unchecked'] = state.get('add_mm')
pick('lcmodel', 'basis')
seen['lcmodel'] = seen_now()
screen._on_next()
seen['next, lcmodel'] = state.get('add_mm', 'not set')
pick('fsLmrs')
seen['fsl-mrs'] = seen_now()
print(json.dumps(seen))
"""


@pytest.mark.skipif(importlib.util.find_spec('tkinter') is None, reason='no tkinter')
def test_osprey_alone_offers_its_macromolecules_checked_and_screen_2_keeps_the_choice():
    run = subprocess.run([sys.executable, '-c', ADD_MM_SHOWN], capture_output=True, text=True,
                         cwd=os.path.join(ROOT, 'basis_converter'))
    if run.returncode == 3:
        pytest.skip(f'Tk cannot open a window: {run.stdout.strip()}')

    assert run.returncode == 0, run.stderr
    text = "Add Osprey's macromolecules and lipids (fit_makeBasis addMMFlag)"
    assert json.loads(run.stdout) == {
        'none picked': 'hidden',
        'osprey': {'checked': True, 'text': text},
        'osprey dropped': 'hidden',
        'osprey again': {'checked': True, 'text': text},
        'unchecked': {'checked': False, 'text': text},
        'next, unchecked': False,
        'lcmodel': 'hidden',
        'next, lcmodel': 'not set',
        'fsl-mrs': 'hidden',
    }


PYAMARES_TEXT = """
import json
import tkinter as tk

try:
    root = tk.Tk()
except tk.TclError as error:
    print(error)
    raise SystemExit(3)
root.withdraw()

from core.formats import ALL_TOOLS
from screen2_output import Screen2Output


class App:
    def quit(self):
        pass

    def show_screen(self, name):
        pass


screen = Screen2Output(root, {}, App())
screen._on_tool_selected(next(t for t in ALL_TOOLS if t['id'] == 'pyamares'))
print(json.dumps([screen.info_panel_title.cget('text'), screen.info_panel_body.cget('text')]))
"""


@pytest.mark.skipif(importlib.util.find_spec('tkinter') is None, reason='no tkinter')
def test_screen_2_says_what_the_pyamares_prior_holds():
    run = subprocess.run([sys.executable, '-c', PYAMARES_TEXT], capture_output=True, text=True,
                         cwd=os.path.join(ROOT, 'basis_converter'))
    if run.returncode == 3:
        pytest.skip(f'Tk cannot open a window: {run.stdout.strip()}')

    assert run.returncode == 0, run.stderr
    title, body = json.loads(run.stdout)
    assert title == 'PyAMARES — CSV prior output'
    for said in ("each metabolite's lines", 'tied', 'one phase, one lineshape and one shift',
                 'largest singlet', 'carrier', 'truncate_initial_points=1', 'uncertainty'):
        assert said in body, said
    assert 'dominant peak' not in body


WRITTEN_FOR_OSPREY = """
import os
import sys
import tkinter as tk
from base_screen import BaseScreen
from core.load import load_basis
from screen4_convert import Screen4Convert

tk.Frame.__init__ = lambda self, *args, **kwargs: None
BaseScreen._build_header = lambda self: None
BaseScreen._build_step_bar = lambda self: None
Screen4Convert.build_body = lambda self: None

folder = os.path.join('..', 'sample_basis_sets', 'JMRUI')
basis = load_basis([os.path.join(folder, f'{p}.txt') for p in ('Cr391', 'CrNo391', 'NAA', 'MM', 'MM09')],
                   'jMRUI .txt')
for outdir, state in zip(sys.argv[1:], ({'core': {'te': 26.0}, 'add_mm': False},
                                        {'core': {'te': 26.0}})):
    Screen4Convert(None, state, None)._write(basis, 'osprey', 'mat', outdir)
"""


@pytest.mark.skipif(importlib.util.find_spec('tkinter') is None, reason='no tkinter')
def test_screen_4_leaves_osprey_macromolecules_out_when_screen_2_says_so(tmp_path):
    """state['add_mm'] False is fit_makeBasis's addMMFlag 0: none of Osprey's MM and lipids are
    added and the set's own MM and MM09 are kept; the H2O BasiC adds for Osprey's water fit is
    there either way. Without add_mm in the state, Osprey's are added, as before."""
    import numpy as np
    import scipy.io as sio
    left_out, added = tmp_path / 'left_out', tmp_path / 'added'
    run = subprocess.run([sys.executable, '-c', WRITTEN_FOR_OSPREY, str(left_out), str(added)],
                         capture_output=True, text=True, cwd=os.path.join(ROOT, 'basis_converter'))
    assert run.returncode == 0, run.stderr

    def written(outdir):
        b = sio.loadmat(str(outdir / 'basis.mat'), squeeze_me=True, struct_as_record=False)['BASIS']
        return sorted(str(n) for n in np.ravel(b.name)), b.flags.addedMM

    ospreys = {'MM12', 'MM14', 'MM17', 'MM20', 'Lip09', 'Lip13', 'Lip20'}
    assert written(left_out) == (['Cr', 'H2O', 'MM', 'MM09', 'NAA'], 0)
    names, flag = written(added)
    assert ospreys <= set(names) and flag == 1


SETS_PICKED = """
import json
import os
import tkinter as tk
from tkinter import ttk

try:
    root = tk.Tk()
except tk.TclError as error:
    print(error)
    raise SystemExit(3)
root.withdraw()

from core.load import load_basis
from screen1_upload import Screen1Upload


def widgets(widget):
    yield widget
    for child in widget.winfo_children():
        yield from widgets(child)


def seen_now():
    labels = [w.cget('text') for w in widgets(screen.result_frame) if isinstance(w, tk.Label)]
    boxes = [w for w in widgets(screen.result_frame) if isinstance(w, ttk.Combobox)]
    return {
        'titles': [t for t in labels if t in ('Basis set', 'Detected format')],
        'badges': [t[1:].strip() for t in labels if t[:1] in '✓⚠✗'],
        'dropdowns': [[str(v) for v in box.cget('values')] for box in boxes],
        'shown': [box.get() for box in boxes],
        'next': str(screen.next_btn.cget('state')),
    }


def pick(name):
    box = next(w for w in widgets(screen.result_frame) if isinstance(w, ttk.Combobox)
               and name in [str(v) for v in w.cget('values')])
    box.set(name)
    box.event_generate('<<ComboboxSelected>>')
    seen = seen_now()
    seen['files'] = len(state['files'])
    seen['format'] = state['detected_format']
    seen['loaded'] = len(load_basis(state['files'], state['detected_format']))
    return seen


class App:                               # what the screen asks of main.App
    def quit(self):
        pass

    def show_screen(self, name):
        pass


samples = os.path.join('..', 'sample_basis_sets')
seen = {}
for folder, picks in ((('LCModel', 'BASIS_file'), ['SLASER_Siemens26.basis',
                                                   'LCModel_Siemens_UnEdited_sLASER_105_.BASIS']),
                      (('MRSCLOUD',), ['Siemens_UnEdited_sLASER_TE105_*.mat',
                                       'BASIS_Siemens_UnEdited_sLASER_105_wMM.mat']),
                      (('JMRUI',), [])):
    state = {}
    screen = Screen1Upload(root, state, App())
    path = os.path.join(samples, *folder)
    screen._load_files([os.path.join(path, f) for f in os.listdir(path)])
    seen[folder[0]] = [seen_now()] + [pick(name) for name in picks]
print(json.dumps(seen))
"""


@pytest.mark.skipif(importlib.util.find_spec('tkinter') is None, reason='no tkinter')
def test_screen_1_lets_the_user_pick_one_of_several_basis_sets_and_converts_only_that_one():
    run = subprocess.run([sys.executable, '-c', SETS_PICKED], capture_output=True, text=True,
                         cwd=os.path.join(ROOT, 'basis_converter'))
    if run.returncode == 3:
        pytest.skip(f'Tk cannot open a window: {run.stdout.strip()}')

    assert run.returncode == 0, run.stderr
    seen = json.loads(run.stdout.splitlines()[-1])
    lcmodel = ['LCModel_Siemens_UnEdited_sLASER_105_.BASIS', 'SLASER_Siemens26.basis']
    mrscloud = ['BASIS_Siemens_UnEdited_sLASER_105_wMM.mat', 'Siemens_UnEdited_sLASER_TE105_*.mat']
    unpicked = {'titles': ['Basis set'], 'badges': ['Pick a set'], 'shown': [''], 'next': 'disabled'}
    picked = {'titles': ['Basis set', 'Detected format'], 'badges': ['Picked', 'Recognized'],
              'next': 'normal'}
    assert seen['LCModel'] == [
        dict(unpicked, dropdowns=[lcmodel]),
        dict(picked, dropdowns=[lcmodel], shown=['SLASER_Siemens26.basis'],
             files=1, format='LCModel .BASIS', loaded=22),
        dict(picked, dropdowns=[lcmodel], shown=['LCModel_Siemens_UnEdited_sLASER_105_.BASIS'],
             files=1, format='LCModel .BASIS', loaded=25),
    ]
    assert seen['MRSCLOUD'][0] == dict(unpicked, dropdowns=[mrscloud])
    assert seen['MRSCLOUD'][1] == dict(picked, dropdowns=[mrscloud], shown=['Siemens_UnEdited_sLASER_TE105_*.mat'],
                                       files=25, format='MRSCloud .mat', loaded=25)
    assert seen['MRSCLOUD'][2] == dict(picked, dropdowns=[mrscloud], shown=['BASIS_Siemens_UnEdited_sLASER_105_wMM.mat'],
                                       files=1, format='Osprey .mat', loaded=33)
    assert seen['JMRUI'] == [{'titles': ['Detected format'], 'badges': ['Recognized'], 'dropdowns': [],
                              'shown': [], 'next': 'normal'}]
