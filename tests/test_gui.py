"""
tests/test_gui.py
The screens, without anyone at the window.

The step bar: BaseScreen.__init__ runs with the Tk frame and the widget builders
stubbed out, so no display is needed. The parameters summary: screen 3 is built on a
hidden Tk root and typed into; the test is skipped when Tk cannot open a window.
Osprey's macromolecules: screen 2 is built on a hidden Tk root and its tools picked,
skipped likewise; screen 4 writes the .mat with its widgets stubbed out, as for the
step bar. Several basis sets in one folder: screen 1 is built on a Tk root shown fully
transparent, as Aqua's Tk 8.6 sends no events to a window never shown, given each folder and
one of its sets picked, skipped likewise. New conversion: main.App is built
hidden, converts, and is reset, skipped likewise. A failed conversion and a VeSPA export:
main.App is built hidden and screen 4 converts, skipped likewise; these three also check that
the conversion thread never calls Tk.

The screens are imported in a separate Python, so that tkinter stays out of this
process and the other tests can still check that core runs without it.

Run from the repository root: python -m pytest tests
"""

import importlib.util
import json
import os
import re
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
                 'largest singlet', 'Broad humps', 'left out', 'carrier',
                 'truncate_initial_points=1', 'uncertainty'):
        assert said in body, said
    assert 'dominant peak' not in body


@pytest.mark.skipif(importlib.util.find_spec('tkinter') is None, reason='no tkinter')
def test_screen_2_says_what_the_oxsa_prior_holds():
    run = subprocess.run([sys.executable, '-c', PYAMARES_TEXT.replace("== 'pyamares'", "== 'oxsa'")],
                         capture_output=True, text=True, cwd=os.path.join(ROOT, 'basis_converter'))
    if run.returncode == 3:
        pytest.skip(f'Tk cannot open a window: {run.stdout.strip()}')

    assert run.returncode == 0, run.stderr
    title, body = json.loads(run.stdout)
    assert title == 'OXSA — prior knowledge output'
    for said in ('PK_<set>.m', "OXSA's own form", "each metabolite's lines", 'tied',
                 'one phase and one linewidth', 'Broad humps', 'left out', 'chemical shifts',
                 'fixed', 'AMARES.amaresFit'):
        assert said in body, said
    assert 'run_oxsa_amares.m' not in body and 'script' not in body


@pytest.mark.skipif(importlib.util.find_spec('tkinter') is None, reason='no tkinter')
def test_screen_2_says_what_the_vespa_gava_and_midas_priors_hold():
    """VeSPA's prior, and GAVA and MIDAS as VeSPA-Simulation exports them, hold each metabolite's
    lines fitted from the basis set on the source's scale (write_vespa.vespa_lines), with the
    broad humps left out: not literature shifts, not relative areas."""
    for tool, title, own in (
            ('vespa_analysis', 'VeSPA Analysis — prior XML output',
             ("'Basis Set from File...'", "VeSPA's own macromolecule model")),
            ('gava', 'GAVA — text prior output', ('VeSPA-Simulation', "the VeSPA prior's lines")),
            ('midas', 'MIDAS — prior XML output', ('VeSPA-Simulation', "the VeSPA prior's lines"))):
        run = subprocess.run([sys.executable, '-c', PYAMARES_TEXT.replace("== 'pyamares'", f"== '{tool}'")],
                             capture_output=True, text=True, cwd=os.path.join(ROOT, 'basis_converter'))
        if run.returncode == 3:
            pytest.skip(f'Tk cannot open a window: {run.stdout.strip()}')

        assert run.returncode == 0, run.stderr
        shown, body = json.loads(run.stdout)
        assert shown == title
        for said in ("each metabolite's lines fitted from the basis set", "source's scale",
                     'lineshape', 'Broad humps', 'left out') + own:
            assert said in body, (tool, said)
        for stale in ('literature', 'relative areas', 'simulated FID', 'echo time'):
            assert stale not in body, (tool, stale)


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
root.attributes('-alpha', 0.0)       # shown, but not seen: Aqua's Tk 8.6 sends no events to a
                                     # window never shown

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
    root.update()
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
    screen.pack()
    screen._load_files([os.path.join(path, f) for f in os.listdir(path)])
    root.update()
    seen[folder[0]] = [seen_now()] + [pick(name) for name in picks]
    screen.destroy()
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




# Screen 4 converts in a thread. Tk may only be called from the main thread: a tkinter that takes
# Tcl for unthreaded (3.11 and 3.12 with Tcl 9) runs a call from the thread as it comes, which
# corrupts Tcl's memory or hangs. So the scripts below build the App on a tkapp that notes every
# call made from any other thread, and run screen 4 as it runs for the user, with no stand-ins.
TK_WATCHED = """
import threading
import tkinter
import _tkinter

off_main = []
tk_create = _tkinter.create


class Watched:
    def __init__(self, tkapp):
        self._tkapp = tkapp

    def __getattr__(self, name):
        attr = getattr(self._tkapp, name)
        if not callable(attr):
            return attr

        def call(*args, **kwargs):
            if threading.current_thread() is not threading.main_thread():
                off_main.append(name)
            return attr(*args, **kwargs)
        return call


tkinter._tkinter.create = lambda *args: Watched(tk_create(*args))
"""


def run_hidden(script, *args, timeout=180):
    """The script in a separate Python, as the tests above run theirs; hanging counts as failing."""
    try:
        return subprocess.run([sys.executable, '-c', script, *args], capture_output=True, text=True,
                              cwd=os.path.join(ROOT, 'basis_converter'), timeout=timeout)
    except subprocess.TimeoutExpired as hung:
        stderr = hung.stderr.decode(errors='replace') if isinstance(hung.stderr, bytes) else hung.stderr
        pytest.fail(f'still running after {timeout} s; stderr: {stderr}')


NEW_CONVERSION = TK_WATCHED + """
import json
import os
import sys
import time
import tkinter as tk

import screen4_convert
from core.formats import ALL_TOOLS
from main import App

try:
    app = App()
except tk.TclError as error:
    print(error)
    raise SystemExit(3)
app.withdraw()

asked = []
screen4_convert.filedialog.askdirectory = lambda **kwargs: asked.append(1) or sys.argv[1]
screen1, screen2, screen4 = (app.screens[name] for name in ('screen1', 'screen2', 'screen4'))


def widgets(widget):
    yield widget
    for child in widget.winfo_children():
        yield from widgets(child)


def convert(seen, key):
    # The conversion thread reaches the widgets through the main loop, so this runs in it.
    app.show_screen('screen4')               # starts the conversion
    back = next(w for w in widgets(screen4) if isinstance(w, tk.Button) and w.cget('text') == '← Back')
    started = [str(screen4.new_btn.cget('state')), str(back.cget('state'))]
    end = time.time() + 60
    while str(screen4.new_btn.cget('state')) != 'normal' and time.time() < end:
        yield
    labels = [w.cget('text') for w in widgets(screen4.scroll_frame) if isinstance(w, tk.Label)]
    seen[key] = {'title': screen4.progress_title.cget('text'), 'error': labels[-1], 'asked': len(asked),
                 'new conversion, back while converting': started,
                 'back after': str(back.cget('state'))}


def seen_now():
    return {
        'card': len(screen1.result_frame.winfo_children()),
        'next 1': str(screen1.next_btn.cget('state')),
        'tool': screen2.selected_tool and screen2.selected_tool['id'],
        'next 2': str(screen2.next_btn.cget('state')),
        'files': len(app.state.get('files', [])),
    }


def clicks(seen):
    folder = os.path.join('..', 'sample_basis_sets', 'JMRUI')
    screen1._load_files([os.path.join(folder, f) for f in ('NAA.txt', 'Cr391.txt')])
    app.show_screen('screen2')
    screen2._on_tool_selected(next(t for t in ALL_TOOLS if t['id'] == 'lcmodel'))
    next(w for w in widgets(screen2) if isinstance(w, tk.Radiobutton) and w.cget('value') == 'basis').invoke()
    seen['before'] = seen_now()
    screen2._on_next()
    app.state['params'] = {'te': 26.0}
    yield from convert(seen, 'converted')
    screen4.new_btn.invoke()
    yield
    seen['after'] = seen_now()
    del asked[:]
    yield from convert(seen, 'empty')


def tick():
    try:
        next(steps)
        app.after(20, tick)
    except StopIteration:
        app.quit()


seen = {}
steps = clicks(seen)
app.after(0, tick)
app.mainloop()
seen['tk off the main thread'] = sorted(set(off_main))
print(json.dumps(seen))
"""


@pytest.mark.skipif(importlib.util.find_spec('tkinter') is None, reason='no tkinter')
def test_new_conversion_starts_from_an_empty_screen_1_and_screen_4_wants_files(tmp_path):
    """After a conversion, New conversion leaves no card on screen 1, its Next disabled, and no
    tool picked on screen 2, as at start-up. Screen 4 shown with no files says so without asking
    for an output folder. New conversion and Back are disabled while screen 4 converts, so that
    Back and Next cannot start a second conversion beside the first, and enabled at its end."""
    run = run_hidden(NEW_CONVERSION, str(tmp_path))
    if run.returncode == 3:
        pytest.skip(f'Tk cannot open a window: {run.stdout.strip()}')

    assert run.returncode == 0, run.stderr
    seen = json.loads(run.stdout.splitlines()[-1])
    assert seen['before'] == {'card': 1, 'next 1': 'normal', 'tool': 'lcmodel', 'next 2': 'normal', 'files': 2}
    assert seen['converted']['title'] == 'Conversion complete'
    assert seen['converted']['new conversion, back while converting'] == ['disabled', 'disabled']
    assert seen['converted']['back after'] == 'normal'
    assert [f for f in os.listdir(tmp_path) if f.endswith('.BASIS')]
    assert seen['after'] == {'card': 0, 'next 1': 'disabled', 'tool': None, 'next 2': 'disabled', 'files': 0}
    assert seen['empty'] == {'title': 'Conversion failed', 'asked': 0,
                             'new conversion, back while converting': ['disabled', 'disabled'],
                             'back after': 'normal',
                             'error': 'No input files to convert. Go back to screen 1 and upload a basis set.'}
    assert seen['tk off the main thread'] == []


CONVERTED = TK_WATCHED + """
import json
import sys
import time
import tkinter as tk

import screen4_convert
from main import App

try:
    app = App()
except tk.TclError as error:
    print(error)
    raise SystemExit(3)
app.withdraw()

asked = []                               # for each folder asked: on the main thread?
screen4_convert.filedialog.askdirectory = lambda **kwargs: (
    asked.append(threading.current_thread() is threading.main_thread()) or sys.argv[1])
screen4 = app.screens['screen4']


def widgets(widget):
    yield widget
    for child in widget.winfo_children():
        yield from widgets(child)


def labels(widget):
    return [w.cget('text') for w in widgets(widget) if isinstance(w, tk.Label)]


end = time.time() + 120


def wait():
    if str(screen4.new_btn.cget('state')) == 'normal' or time.time() > end:
        app.quit()
    else:
        app.after(20, wait)


app.state.update(json.loads(sys.argv[2]))
app.show_screen('screen4')
app.after(20, wait)
app.mainloop()
print(json.dumps({'title': screen4.progress_title.cget('text'),
                  'steps': [labels(row)[1:] for row in screen4.step_list.winfo_children()],
                  'window': labels(screen4.scroll_frame), 'asked on the main thread': asked,
                  'tk off the main thread': sorted(set(off_main))}))
"""


@pytest.mark.skipif(importlib.util.find_spec('tkinter') is None, reason='no tkinter')
def test_a_failed_conversion_prints_its_traceback_and_shows_its_message(tmp_path):
    """Read as single MRSCloud files, MRSCloud's combined BASIS_*.mat fails in the reader: the
    window shows the message alone, as before, and stderr holds the traceback to it."""
    combined = os.path.join('..', 'sample_basis_sets', 'MRSCLOUD', 'BASIS_Siemens_UnEdited_sLASER_105_wMM.mat')
    run = run_hidden(CONVERTED, str(tmp_path), json.dumps({
        'files': [combined], 'detected_format': 'MRSCloud .mat',
        'output_tool': 'LCModel', 'output_tool_id': 'lcmodel', 'output_format': 'basis'}))
    if run.returncode == 3:
        pytest.skip(f'Tk cannot open a window: {run.stdout.strip()}')

    assert run.returncode == 0, run.stderr
    seen = json.loads(run.stdout.splitlines()[-1])
    message = "Missing field 'txfrq' in BASIS_Siemens_UnEdited_sLASER_105_wMM.mat"
    assert seen['title'] == 'Conversion failed' and seen['window'][-1] == message
    assert seen['steps'] == [['Loading input files', 'MRSCloud .mat']]
    assert seen['asked on the main thread'] == [True] and seen['tk off the main thread'] == []
    assert 'Traceback (most recent call last)' in run.stderr
    assert 'in _convert' in run.stderr and 'read_mrscloud.py' in run.stderr
    assert run.stderr.rstrip().endswith(f'RuntimeError: {message}')


@pytest.mark.skipif(importlib.util.find_spec('tkinter') is None, reason='no tkinter')
def test_screen_4_lists_each_vespa_metabolite_as_its_lines_are_fitted(tmp_path):
    """The VeSPA writer's progress reaches screen 4's step list from the conversion thread
    through the main loop, a row as the functions are decomposed, one for the lineshape and one a
    metabolite, in order."""
    folder = os.path.join('..', 'sample_basis_sets', 'JMRUI')
    run = run_hidden(CONVERTED, str(tmp_path), json.dumps({
        'files': [os.path.join(folder, f) for f in ('NAA.txt', 'Cr391.txt')], 'detected_format': 'jMRUI .txt',
        'output_tool': 'VeSPA Analysis', 'output_tool_id': 'vespa_analysis', 'output_format': 'default'}))
    if run.returncode == 3:
        pytest.skip(f'Tk cannot open a window: {run.stdout.strip()}')

    assert run.returncode == 0, run.stderr
    seen = json.loads(run.stdout.splitlines()[-1])
    assert seen['title'] == 'Conversion complete'
    steps = seen['steps']
    assert [step[0] for step in steps] == ['Loading input files', 'Loaded 2 metabolites',
                                           'Writing VeSPA Analysis format', 'Decomposing each function',
                                           "Lineshape of the set's lines", 'Lines 1/2: NAA',
                                           'Lines 2/2: Cr391', 'Done']
    assert steps[3][1].startswith('2 functions') and steps[4][1].startswith('Lorentzian 1.000 Hz')
    assert re.match(r'\d+ lines, residual', steps[5][1])
    assert steps[7] == ['Done', f'Output: {tmp_path}']
    assert [f for f in os.listdir(tmp_path) if f.endswith('.xml')]
    assert seen['asked on the main thread'] == [True] and seen['tk off the main thread'] == []
