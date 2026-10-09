"""
screen4_convert.py
Basis Set Converter — Screen 4: Conversion and output

Runs the conversion in a background thread so the GUI stays responsive.
Shows a step-by-step progress list and the output files when done.

State read by this screen:
    state['files']          : list of input file paths
    state['detected_format']: input format string
    state['detection_status']: 'full' | 'partial' | 'ambiguous'
    state['core']           : dict with FIDs and parameters
    state['output_tool']    : e.g. 'LCModel'
    state['output_tool_id'] : e.g. 'lcmodel'
    state['output_format']  : e.g. 'both' | 'raw' | 'default'
    state['params']         : user-supplied parameter overrides

State written by this screen:
    state['output_path']    : path to output folder or file
"""

import os
import tkinter as tk
from tkinter import ttk, filedialog
import threading

from base_screen import BaseScreen


class Screen4Convert(BaseScreen):

    current_step = 4

    def build_body(self):
        self._status     = 'idle'   # idle | running | done | error
        self._steps_done = []
        self._output_paths = []

        # Scrollable body
        outer = tk.Frame(self, bg="white")
        outer.pack(fill="both", expand=True)

        self._canvas = tk.Canvas(outer, bg="white", highlightthickness=0)
        scrollbar    = ttk.Scrollbar(outer, orient="vertical",
                                     command=self._canvas.yview)
        self._canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        self._canvas.pack(side="left", fill="both", expand=True)

        self.scroll_frame = tk.Frame(
            self._canvas, bg="white", padx=24, pady=12
        )
        self._cw = self._canvas.create_window(
            (0, 0), window=self.scroll_frame, anchor="nw"
        )
        self.scroll_frame.bind("<Configure>",
            lambda e: self._canvas.configure(
                scrollregion=self._canvas.bbox("all")))
        self._canvas.bind("<Configure>",
            lambda e: self._canvas.itemconfig(self._cw, width=e.width))

        self._build_footer()

    def on_show(self):
        """Start conversion when this screen is shown."""
        for w in self.scroll_frame.winfo_children():
            w.destroy()
        self._status       = 'idle'
        self._steps_done   = []
        self._output_paths = []
        self._build_context_bar()
        self._build_progress_area()
        # Start conversion in background thread
        thread = threading.Thread(target=self._run_conversion, daemon=True)
        thread.start()

    #################### Context bar ####################

    def _build_context_bar(self):
        in_fmt   = self.state.get('detected_format', '—')
        out_tool = self.state.get('output_tool', '—')
        core     = self.state.get('core', {})
        n_mets   = len(core.get('names', []))
        Bo       = core.get('Bo')
        te       = core.get('te')

        ctx = tk.Frame(
            self.scroll_frame, bg="#f3f4f6",
            padx=12, pady=8
        )
        ctx.pack(fill="x", pady=(0, 12))

        info = f"{in_fmt}  →  {out_tool}"
        if n_mets:
            info += f"   ·   {n_mets} metabolites"
        if Bo:
            info += f"   ·   {Bo:.3f} T"
        if te:
            info += f"   ·   TE {te:.0f} ms"

        tk.Label(
            ctx,
            text=info,
            font=("Helvetica", 10),
            bg="#f3f4f6", fg="#374151"
        ).pack(side="left")

    #################### Progress area #################### 

    def _build_progress_area(self):
        self.progress_card = tk.Frame(
            self.scroll_frame,
            bg="white",
            highlightbackground="#e5e7eb",
            highlightthickness=1
        )
        self.progress_card.pack(fill="x")

        # Header
        self.progress_header = tk.Frame(self.progress_card, bg="white",
                                        padx=14, pady=10)
        self.progress_header.pack(fill="x")

        self.progress_title = tk.Label(
            self.progress_header,
            text="Preparing…",
            font=("Helvetica", 11, "bold"),
            bg="white", fg="#111827"
        )
        self.progress_title.pack(side="left")

        self.progress_count = tk.Label(
            self.progress_header,
            text="",
            font=("Helvetica", 9),
            bg="white", fg="#9ca3af"
        )
        self.progress_count.pack(side="right")

        # Progress bar
        self.progress_bar = ttk.Progressbar(
            self.progress_card,
            orient="horizontal",
            mode="indeterminate",
            length=400
        )
        self.progress_bar.pack(fill="x", padx=0, pady=0)
        self.progress_bar.start(15)

        # Step list frame
        self.step_list = tk.Frame(self.progress_card, bg="white",
                                  padx=14, pady=8)
        self.step_list.pack(fill="x")

        # Output area (shown when done)
        self.output_area = tk.Frame(self.scroll_frame, bg="white")
        # packed later when done

    def _add_step(self, text, state='done', sub=None):
        """Add a step row to the progress list (called from main thread)."""
        row = tk.Frame(self.step_list, bg="white")
        row.pack(fill="x", pady=2)

        # Icon
        if state == 'done':
            icon_text, icon_bg, icon_fg = "✓", "#d1fae5", "#065f46"
        elif state == 'error':
            icon_text, icon_bg, icon_fg = "✗", "#fecaca", "#991b1b"
        else:
            icon_text, icon_bg, icon_fg = "·", "#e5e7eb", "#9ca3af"

        icon = tk.Label(
            row,
            text=icon_text,
            font=("Helvetica", 9, "bold"),
            bg=icon_bg, fg=icon_fg,
            width=2, padx=3
        )
        icon.pack(side="left")

        # Text
        col_frame = tk.Frame(row, bg="white")
        col_frame.pack(side="left", padx=(6, 0))

        tk.Label(
            col_frame,
            text=text,
            font=("Helvetica", 10),
            bg="white", fg="#111827"
        ).pack(anchor="w")

        if sub:
            tk.Label(
                col_frame,
                text=sub,
                font=("Helvetica", 9),
                bg="white", fg="#9ca3af"
            ).pack(anchor="w")

        self.step_list.update_idletasks()

    def _set_done(self, output_paths):
        """Called when conversion succeeds."""
        self.progress_bar.stop()
        self.progress_bar.configure(mode="determinate", value=100)
        self.progress_title.config(text="Conversion complete")
        self.progress_count.config(text="All steps done")
        self.progress_card.config(highlightbackground="#bbf7d0")
        self._build_output_panel(output_paths)
        self._enable_footer_buttons()

    def _set_error(self, message):
        """Called when conversion fails."""
        self.progress_bar.stop()
        self.progress_bar.configure(mode="determinate", value=0)
        self.progress_title.config(text="Conversion failed", fg="#dc2626")
        self.progress_card.config(highlightbackground="#fecaca")

        err_card = tk.Frame(
            self.scroll_frame,
            bg="#fef2f2",
            highlightbackground="#fecaca",
            highlightthickness=1,
            padx=14, pady=10
        )
        err_card.pack(fill="x", pady=(8, 0))

        tk.Label(
            err_card,
            text="✗  Conversion failed",
            font=("Helvetica", 11, "bold"),
            bg="#fef2f2", fg="#dc2626"
        ).pack(anchor="w")

        tk.Label(
            err_card,
            text=str(message),
            font=("Helvetica", 9),
            bg="#fef2f2", fg="#991b1b",
            wraplength=560, justify="left"
        ).pack(anchor="w", pady=(6, 0))

        self._enable_footer_buttons()

    def _build_output_panel(self, output_paths):
        """Show output files when done."""
        out_card = tk.Frame(
            self.scroll_frame,
            bg="#f0fdf4",
            highlightbackground="#bbf7d0",
            highlightthickness=1
        )
        out_card.pack(fill="x", pady=(10, 0))

        hdr = tk.Frame(out_card, bg="#f0fdf4", padx=14, pady=10)
        hdr.pack(fill="x")

        tk.Label(
            hdr,
            text="✓  Output ready",
            font=("Helvetica", 11, "bold"),
            bg="#f0fdf4", fg="#065f46"
        ).pack(side="left")

        n_files = len(output_paths) if isinstance(output_paths, list) \
                  else 1
        tk.Label(
            hdr,
            text=f"{n_files} file(s) written",
            font=("Helvetica", 9),
            bg="#f0fdf4", fg="#065f46"
        ).pack(side="right")

        ttk.Separator(out_card, orient="horizontal").pack(fill="x")

        # File list
        files = output_paths if isinstance(output_paths, list) \
                else [output_paths]
        for fpath in files[:8]:
            row = tk.Frame(out_card, bg="#f0fdf4", padx=14, pady=4)
            row.pack(fill="x")
            tk.Label(
                row,
                text=f"📄  {os.path.basename(fpath)}",
                font=("Helvetica", 10),
                bg="#f0fdf4", fg="#111827"
            ).pack(side="left")

            size = self._file_size_str(fpath)
            tk.Label(
                row,
                text=size,
                font=("Helvetica", 9),
                bg="#f0fdf4", fg="#9ca3af"
            ).pack(side="right")

        if len(files) > 8:
            tk.Label(
                out_card,
                text=f"  + {len(files)-8} more files",
                font=("Helvetica", 9),
                bg="#f0fdf4", fg="#9ca3af",
                padx=14
            ).pack(anchor="w")

        ttk.Separator(out_card, orient="horizontal").pack(fill="x")

        # Actions
        actions = tk.Frame(out_card, bg="#f0fdf4", padx=14, pady=8)
        actions.pack(fill="x")

        outdir = os.path.dirname(files[0]) if files else ''

        tk.Button(
            actions,
            text="📁  Open output folder",
            command=lambda: os.system(f"open '{outdir}'"),
            relief="flat", bg="#dbeafe", fg="#1d4ed8",
            padx=12, pady=5, cursor="hand2"
        ).pack(side="left", padx=(0, 8))

        tk.Button(
            actions,
            text="Copy path",
            command=lambda: self._copy_to_clipboard(outdir),
            relief="flat", bg="#f3f4f6", fg="#374151",
            padx=12, pady=5, cursor="hand2"
        ).pack(side="left")

        self._output_paths = files
        self.state['output_path'] = outdir

    #################### Conversion logic #################### 

    def _run_conversion(self):
        """
        Run the full conversion pipeline in a background thread.
        Uses self.after() to update the UI safely from the thread.
        """
        try:
            self._convert()
        except Exception as ex:
            self.after(0, self._set_error, str(ex))

    def _step(self, text, sub=None):
        """Add a step to the UI from the background thread."""
        self.after(0, self._add_step, text, 'done', sub)

    def _convert(self):
        """Main conversion pipeline."""
        import sys
        sys.path.insert(0, os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))

        files    = self.state.get('files', [])
        fmt      = self.state.get('detected_format', '')
        core     = self.state.get('core', {})
        tool_id  = self.state.get('output_tool_id', '')
        out_fmt  = self.state.get('output_format', 'default')

        #################### Ask user for output folder #################### 
        outdir = self.after(0, self._ask_output_folder)
        # Since filedialog must run on main thread, use a synchronized approach
        outdir_holder = [None]
        done_event    = threading.Event()

        def ask():
            d = filedialog.askdirectory(title="Select output folder")
            outdir_holder[0] = d or os.path.expanduser("~/Desktop/basis_output")
            done_event.set()

        self.after(0, ask)
        done_event.wait()
        outdir = outdir_holder[0]
        os.makedirs(outdir, exist_ok=True)

        self.after(0, self.progress_title.config, {'text': 'Converting…'})

        ##################### Step 1: Load all basis functions #################### 
        self._step("Loading input files", sub=fmt)
        basis_list = self._load_basis(files, fmt)

        if not basis_list:
            raise RuntimeError("No basis functions could be loaded.")

        self._step(f"Loaded {len(basis_list)} metabolites")

        # ── Step 2: Apply user params #################### 
        params = self.state.get('params', {})
        if params:
            self._step("Applying user parameters")
            for core_item in basis_list:
                for key, val in params.items():
                    # Overwrite if missing OR if currently None
                    if val and (key not in core_item or core_item[key] is None):
                        core_item[key] = val

        #################### Step 3: Run writer #################### 
        self._step(f"Writing {self.state.get('output_tool', '')} format",
                   sub=f"output format: {out_fmt}")

        output_paths = self._write(basis_list, tool_id, out_fmt, outdir)

        self._step("Done", sub=f"Output: {outdir}")

        ####################  Show output #################### 
        self.after(0, self._set_done, output_paths)

    def _load_basis(self, files, fmt):
        """Load basis functions using the appropriate reader."""
        from core.load import load_basis
        sw = self.state.get('core', {}).get('sw') or \
             self.state.get('params', {}).get('sw')
        sf = self.state.get('core', {}).get('sf') or \
             self.state.get('params', {}).get('sf')
        return load_basis(files, fmt, sw=sw, sf=sf)

    def _write(self, basis_list, tool_id, out_fmt, outdir):
        """Run the appropriate writer."""
        te  = self.state.get('core', {}).get('te') or \
              self.state.get('params', {}).get('te') or 30.0
        seq = self.state.get('core', {}).get('seq') or \
              self.state.get('params', {}).get('seq') or 'unedited'

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
            te = float(self.state.get('core', {}).get('te') or
                      self.state.get('params', {}).get('te') or 26.0)
            write_oxsa(basis_list, out_file, seqte=te)
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

    def _ask_output_folder(self):
        pass   # placeholder — actual call is in _convert via threading

    #################### Footer ####################

    def _build_footer(self):
        footer = tk.Frame(self, bg="white", pady=12, padx=20)
        footer.pack(fill="x", side="bottom")
        ttk.Separator(self, orient="horizontal").pack(
            fill="x", side="bottom", before=footer
        )

        tk.Button(
            footer, text="← Back",
            command=lambda: self.app.show_screen("screen3"),
            relief="flat", bg="#f3f4f6", fg="#374151",
            padx=14, pady=6, cursor="hand2"
        ).pack(side="left")

        self.new_btn = tk.Button(
            footer, text="New conversion →",
            command=self._new_conversion,
            relief="flat", bg="#e5e7eb", fg="#9ca3af",
            padx=14, pady=6, state="disabled"
        )
        self.new_btn.pack(side="right")

    def _enable_footer_buttons(self):
        self.new_btn.config(
            bg="#dbeafe", fg="#1d4ed8",
            state="normal", cursor="hand2"
        )

    def _new_conversion(self):
        """Reset state and go back to screen 1."""
        for key in ['files', 'detected_format', 'detection_status',
                    'core', 'missing_fields', 'ambiguous_options',
                    'compatible_tools', 'format_confirmed',
                    'output_tool', 'output_tool_id', 'output_format',
                    'output_path', 'params']:
            self.state.pop(key, None)
        self.app.show_screen("screen1")

    #################### Utilities ##############################

    def _copy_to_clipboard(self, text):
        self.clipboard_clear()
        self.clipboard_append(text)

    def _file_size_str(self, path):
        try:
            size = os.path.getsize(path)
            if size < 1024:
                return f"{size} B"
            elif size < 1024**2:
                return f"{size/1024:.1f} KB"
            else:
                return f"{size/1024**2:.1f} MB"
        except Exception:
            return ""