"""
tests/test_gui.py
The screens' step bar, without a display: BaseScreen.__init__ runs with the Tk frame
and the widget builders stubbed out.

The screens are imported in a separate Python, so that tkinter stays out of this
process and the other tests can still check that core runs without it.

Run from the repository root: python -m pytest tests
"""

import importlib.util
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
