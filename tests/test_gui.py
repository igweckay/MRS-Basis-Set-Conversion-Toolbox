"""
tests/test_gui.py
The screens, without anyone at the window.

The step bar: BaseScreen.__init__ runs with the Tk frame and the widget builders
stubbed out, so no display is needed. The parameters summary: screen 3 is built on a
hidden Tk root and typed into; the test is skipped when Tk cannot open a window.

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
