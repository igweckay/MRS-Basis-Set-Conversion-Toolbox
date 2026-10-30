"""
screen1_upload.py
Basis Set Converter — Screen 1: Upload

The user drops or browses for their basis set files.
This screen detects the format automatically and shows the result.

What this screen does:
    1. Lets the user select files or a folder via Browse buttons
    2. Calls detect_format() to identify what was uploaded
    3. Shows the detected format, file list, and extracted metadata
    4. If format is ambiguous, shows a dropdown to let the user confirm
       (likewise, if the files hold several basis sets, to pick the one to convert)
    5. If partial data is found, warns the user but allows proceeding
    6. If nothing usable, shows an error and blocks Next

State written by this screen:
    state['files']            : list of file paths
    state['detected_format']  : string name of detected format
    state['detection_status'] : 'full' | 'partial' | 'ambiguous' | 'several' | 'failed'
    state['core']             : dict with whatever fields could be extracted
    state['missing_fields']   : list of field names not found
    state['uploaded_files']   : every file uploaded ('files' holds the picked set's only)
    state['basis_sets']       : the sets found, when the upload holds more than one
    state['set_picked']       : index of the set picked among them, None until then
"""

import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from base_screen import BaseScreen
from core.detect import detect_format


class Screen1Upload(BaseScreen):

    current_step = 1

    def build_body(self):
        # Main scrollable content area
        self.body = tk.Frame(self, bg="white", padx=24, pady=16)
        self.body.pack(fill="both", expand=True)

        self._build_drop_zone()
        self._build_result_area()

        # Footer — Next disabled until files are loaded
        self.next_btn = None
        self._build_footer()

    def on_show(self):
        """Refresh when navigating back to this screen."""
        # If files were already loaded, re-show the result
        if self.state.get('files'):
            self._show_result()

    ############# Drop zone #############

    def _build_drop_zone(self):
        self.drop_frame = tk.Frame(
            self.body,
            bg="#f9fafb",
            relief="flat",
            bd=1,
            padx=20,
            pady=30
        )
        self.drop_frame.pack(fill="x", pady=(0, 12))

        tk.Label(
            self.drop_frame,
            text="⬆  Drop your basis set here",
            font=("Helvetica", 13, "bold"),
            bg="#f9fafb",
            fg="#374151"
        ).pack()

        tk.Label(
            self.drop_frame,
            text="Drop a single file, multiple files, or an entire folder.\n"
                 "The format will be detected automatically.",
            font=("Helvetica", 10),
            bg="#f9fafb",
            fg="#6b7280",
            justify="center"
        ).pack(pady=(4, 12))

        # Browse buttons row
        btn_row = tk.Frame(self.drop_frame, bg="#f9fafb")
        btn_row.pack()

        tk.Button(
            btn_row,
            text="Browse files",
            command=self._browse_files,
            relief="flat",
            bg="white",
            fg="#374151",
            padx=14, pady=6,
            cursor="hand2",
            bd=1
        ).pack(side="left", padx=4)

        tk.Button(
            btn_row,
            text="Browse folder",
            command=self._browse_folder,
            relief="flat",
            bg="white",
            fg="#374151",
            padx=14, pady=6,
            cursor="hand2",
            bd=1
        ).pack(side="left", padx=4)

        # Accepted formats pills
        tk.Label(
            self.drop_frame,
            text="Accepted: .mat  .BASIS  .raw  .json  .txt  .nii.gz  .csv  .xml  folder",
            font=("Helvetica", 9),
            bg="#f9fafb",
            fg="#9ca3af"
        ).pack(pady=(12, 0))

    ############# Result area #############

    def _build_result_area(self):
        """Container where detection result card will be drawn."""
        self.result_frame = tk.Frame(self.body, bg="white")
        self.result_frame.pack(fill="both", expand=True)

    def _clear_result(self):
        for widget in self.result_frame.winfo_children():
            widget.destroy()

    ############# Footer #############

    def _build_footer(self):
        footer = tk.Frame(self, bg="white", pady=12, padx=20)
        footer.pack(fill="x", side="bottom")

        ttk.Separator(self, orient="horizontal").pack(
            fill="x", side="bottom", before=footer
        )

        # Clear button (hidden until files loaded)
        self.clear_btn = tk.Button(
            footer,
            text="Clear",
            command=self._clear_all,
            relief="flat",
            bg="#f3f4f6",
            fg="#6b7280",
            padx=14, pady=6,
            cursor="hand2"
        )
        self.clear_btn.pack(side="left")
        self.clear_btn.pack_forget()   # hidden initially

        # Next button — disabled until detection succeeds
        self.next_btn = tk.Button(
            footer,
            text="Next →",
            command=lambda: self.app.show_screen("screen2"),
            relief="flat",
            bg="#e5e7eb",
            fg="#9ca3af",
            padx=14, pady=6,
            state="disabled"
        )
        self.next_btn.pack(side="right")

    def _enable_next(self):
        self.next_btn.config(
            bg="#dbeafe", fg="#1d4ed8",
            state="normal", cursor="hand2"
        )
        self.clear_btn.pack(side="left")

    def _disable_next(self):
        self.next_btn.config(
            bg="#e5e7eb", fg="#9ca3af",
            state="disabled", cursor=""
        )
        self.clear_btn.pack_forget()

    ############## Browse handlers #############

    def _browse_files(self):
        paths = filedialog.askopenfilenames(
            title="Select basis set files",
            filetypes=[
                ("All supported", "*.mat *.BASIS *.basis *.raw *.json *.txt *.gz *.csv *.xml"),
                ("MATLAB files", "*.mat"),
                ("LCModel basis", "*.BASIS *.basis"),
                ("LCModel raw", "*.raw"),
                ("JSON files", "*.json"),
                ("Text files", "*.txt"),
                ("NIfTI", "*.gz"),
                ("PyAMARES priors", "*.csv"),
                ("VeSPA priors", "*.xml"),
                ("All files", "*.*"),
            ]
        )
        if paths:
            self._load_files(list(paths))

    def _browse_folder(self):
        folder = filedialog.askdirectory(title="Select basis set folder")
        if folder:
            IGNORED_NO_EXT = {'Icon'}
            paths = [
                os.path.join(folder, f)
                for f in os.listdir(folder)
                if os.path.isfile(os.path.join(folder, f))
                and not f.startswith('.')
                and f not in ('Thumbs.db', 'desktop.ini')
                and not (os.path.splitext(f)[1] == ''
                         and f.strip() in IGNORED_NO_EXT)
            ]
            if paths:
                self._load_files(paths)
            else:
                messagebox.showwarning(
                    "Empty folder",
                    "The selected folder contains no files."
                )

    ############## Core: load and detect #############

    def _load_files(self, paths):
        """
        Called after the user selects files or a folder.
        Runs detection and updates the shared state.
        """
        self.state['files'] = paths

        # Run format detection
        result = detect_format(paths)

        # Several basis sets: each is a detect_format result of its own, applied once picked
        self.state['uploaded_files'] = paths
        self.state['basis_sets']     = result.get('sets', [])
        self.state['set_picked']     = None

        self._apply_detection(result)
        self._show_result()

    def _apply_detection(self, result):
        """Write a detect_format result to the shared state."""
        self.state['detected_format']    = result.get('format', 'Unknown')
        self.state['detection_status']   = result.get('status', 'failed')
        self.state['core']               = result.get('data', {})
        self.state['missing_fields']     = result.get('missing', [])
        self.state['ambiguous_options']  = result.get('ambiguous_options', [])
        self.state['compatible_tools']   = result.get('compatible_tools', [])
        self.state['format_confirmed']   = False   # reset on new load

    def _show_result(self):
        """Draw the detection result card based on current state."""
        self._clear_result()

        if self.state.get('basis_sets'):
            self._show_sets_card()
            if self.state.get('set_picked') is None:
                self._disable_next()
                return

        status    = self.state.get('detection_status', 'failed')
        confirmed = self.state.get('format_confirmed', False)

        if status == 'failed':
            self._show_error_card()
            self._disable_next()
        elif status == 'ambiguous' and not confirmed:
            # Still ambiguous — show dropdown for user to confirm
            self._show_ambiguous_card()
            self._enable_next()
        elif status == 'partial' and not confirmed:
            self._show_partial_card()
            self._enable_next()
        else:
            # 'full', or user confirmed ambiguous/partial — show green card
            self._show_success_card()
            self._enable_next()

    ############# Result cards #############

    def _card_frame(self, border_color="#e5e7eb", bg_color="white"):
        """Helper — returns a framed card widget."""
        outer = tk.Frame(
            self.result_frame,
            bg=border_color,
            padx=1, pady=1
        )
        outer.pack(fill="x", pady=(0, 8))

        inner = tk.Frame(outer, bg=bg_color, padx=14, pady=12)
        inner.pack(fill="x")

        return inner

    def _card_header(self, card, format_text, badge_text, badge_color,
                     badge_fg, show_dropdown=False, dropdown_options=None,
                     title="Detected format", on_select=None, value=None,
                     show_tools=True):
        """Draw the top row of a result card. The dropdown confirms the format,
        or, given on_select, calls it when an option is picked."""
        header = tk.Frame(card, bg=card["bg"])
        header.pack(fill="x", pady=(0, 4))

        tk.Label(
            header,
            text=title,
            font=("Helvetica", 9),
            bg=card["bg"],
            fg="#6b7280"
        ).pack(side="left")

        if show_dropdown and dropdown_options:
            # Ambiguous — show a dropdown
            var = tk.StringVar(value=dropdown_options[0] if value is None else value)
            dropdown = ttk.Combobox(
                header,
                textvariable=var,
                values=dropdown_options,
                state="readonly",
                width=max(22, max(len(o) for o in dropdown_options)),
                font=("Helvetica", 10)
            )
            dropdown.pack(side="left", padx=(8, 0))
            if on_select is None:
                self.format_var = var
                dropdown.bind("<<ComboboxSelected>>", self._on_format_confirmed)
            else:
                self.set_var = var
                dropdown.bind("<<ComboboxSelected>>", on_select)
        else:
            tk.Label(
                header,
                text=format_text,
                font=("Helvetica", 11, "bold"),
                bg=card["bg"],
                fg="#111"
            ).pack(side="left", padx=(8, 0))

        # Badge
        tk.Label(
            header,
            text=badge_text,
            font=("Helvetica", 9),
            bg=badge_color,
            fg=badge_fg,
            padx=8, pady=2
        ).pack(side="right")

        # Compatible tools row
        tools = self.state.get('compatible_tools', []) if show_tools else []
        if tools:
            tools_row = tk.Frame(card, bg=card["bg"])
            tools_row.pack(fill="x", pady=(0, 6))
            tk.Label(
                tools_row,
                text="Compatible with",
                font=("Helvetica", 9),
                bg=card["bg"],
                fg="#6b7280"
            ).pack(side="left")
            tk.Label(
                tools_row,
                text="  ·  ".join(tools),
                font=("Helvetica", 9, "bold"),
                bg=card["bg"],
                fg="#374151"
            ).pack(side="left", padx=(8, 0))

    def _card_files(self, card, files):
        """
        Draw a scrollable file list inside a card.
        All files shown — scrollable if more than 6 are present.
        """
        ttk.Separator(card, orient="horizontal").pack(fill="x", pady=(0, 6))

        row_height   = 20
        visible_rows = 6
        canvas_height = row_height * min(len(files), visible_rows)

        canvas_frame = tk.Frame(card, bg=card["bg"])
        canvas_frame.pack(fill="x")

        canvas = tk.Canvas(
            canvas_frame,
            bg=card["bg"],
            height=canvas_height,
            highlightthickness=0
        )
        scrollbar = ttk.Scrollbar(
            canvas_frame,
            orient="vertical",
            command=canvas.yview
        )
        canvas.configure(yscrollcommand=scrollbar.set)

        if len(files) > visible_rows:
            scrollbar.pack(side="right", fill="y")

        canvas.pack(side="left", fill="x", expand=True)

        inner = tk.Frame(canvas, bg=card["bg"])
        canvas_window = canvas.create_window((0, 0), window=inner, anchor="nw")

        for path in files:
            fname = os.path.basename(path)
            size  = self._file_size_str(path)
            row   = tk.Frame(inner, bg=card["bg"])
            row.pack(fill="x", pady=1)
            tk.Label(
                row,
                text=f"📄  {fname}",
                font=("Helvetica", 10),
                bg=card["bg"],
                fg="#374151"
            ).pack(side="left")
            tk.Label(
                row,
                text=size,
                font=("Helvetica", 9),
                bg=card["bg"],
                fg="#9ca3af"
            ).pack(side="right")

        inner.update_idletasks()
        canvas.configure(scrollregion=canvas.bbox("all"))

        def _on_canvas_resize(event):
            canvas.itemconfig(canvas_window, width=event.width)
        canvas.bind("<Configure>", _on_canvas_resize)

        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        canvas.bind("<MouseWheel>", _on_mousewheel)
        canvas.bind("<Button-4>", lambda e: canvas.yview_scroll(-1, "units"))
        canvas.bind("<Button-5>", lambda e: canvas.yview_scroll( 1, "units"))

        if len(files) > visible_rows:
            tk.Label(
                card,
                text=f"{len(files)} files total — scroll to see all",
                font=("Helvetica", 9),
                bg=card["bg"],
                fg="#9ca3af"
            ).pack(anchor="w", pady=(2, 0))

        ttk.Separator(card, orient="horizontal").pack(fill="x", pady=(6, 0))

    def _card_metadata(self, card):
        """Draw the 2x2 metadata grid inside a card."""
        core    = self.state.get('core', {})
        missing = self.state.get('missing_fields', [])
        fmt     = self.state.get('detected_format', '')

        # PyAMARES CSV shows peaks instead of spectral parameters
        if 'PyAMARES' in fmt or fmt.endswith('.csv'):
            n_peaks = len(core.get('names', []))
            meta = {
                "Peaks":          str(n_peaks) if n_peaks else "—",
                "Field strength": "—",
                "Bandwidth":      "—",
                "Points":         "—",
            }
        else:
            # Only show a value if present and not None/zero
            n_metabs = len(core.get('names', []))
            Bo       = core.get('Bo')
            sw       = core.get('sw')
            n        = core.get('n')
            meta = {
                "Metabolites":    str(n_metabs) if n_metabs else "—",
                "Field strength": f"{Bo:.3f} T"  if Bo      else "—",
                "Bandwidth":      f"{sw:.0f} Hz" if sw      else "—",
                "Points":         str(n)          if n       else "—",
            }

        grid = tk.Frame(card, bg=card["bg"])
        grid.pack(fill="x", pady=(8, 0))

        for i, (key, val) in enumerate(meta.items()):
            col = i % 2
            row = i // 2
            cell = tk.Frame(grid, bg=card["bg"])
            cell.grid(row=row, column=col, sticky="w", padx=(0, 24), pady=2)
            tk.Label(
                cell,
                text=key,
                font=("Helvetica", 9),
                bg=card["bg"],
                fg="#9ca3af"
            ).pack(anchor="w")
            color = "#f59e0b" if key.lower().replace(" ", "_") in missing else "#111"
            tk.Label(
                cell,
                text=val,
                font=("Helvetica", 10, "bold"),
                bg=card["bg"],
                fg=color
            ).pack(anchor="w")

        ############## Scrollable metabolite name list #############
        names = core.get('names', [])
        if names and 'PyAMARES' not in fmt:
            tk.Label(
                card,
                text="Metabolites",
                font=("Helvetica", 9),
                bg=card["bg"],
                fg="#9ca3af"
            ).pack(anchor="w", pady=(10, 2))

            row_height    = 18
            visible_rows  = 6
            canvas_height = row_height * min(len(names), visible_rows)

            canvas_frame = tk.Frame(card, bg=card["bg"])
            canvas_frame.pack(fill="x")

            canvas = tk.Canvas(
                canvas_frame,
                bg=card["bg"],
                height=canvas_height,
                highlightthickness=0
            )
            scrollbar = ttk.Scrollbar(
                canvas_frame,
                orient="vertical",
                command=canvas.yview
            )
            canvas.configure(yscrollcommand=scrollbar.set)

            if len(names) > visible_rows:
                scrollbar.pack(side="right", fill="y")

            canvas.pack(side="left", fill="x", expand=True)

            inner = tk.Frame(canvas, bg=card["bg"])
            canvas_window = canvas.create_window((0, 0), window=inner, anchor="nw")

            for name in names:
                tk.Label(
                    inner,
                    text=f"  {name}",
                    font=("Helvetica", 9),
                    bg=card["bg"],
                    fg="#374151",
                    anchor="w"
                ).pack(fill="x", pady=1)

            inner.update_idletasks()
            canvas.configure(scrollregion=canvas.bbox("all"))

            def _on_resize(event, cw=canvas_window, cv=canvas):
                cv.itemconfig(cw, width=event.width)
            canvas.bind("<Configure>", _on_resize)

            def _on_wheel(event, cv=canvas):
                cv.yview_scroll(int(-1 * (event.delta / 120)), "units")
            canvas.bind("<MouseWheel>", _on_wheel)
            canvas.bind("<Button-4>", lambda e, cv=canvas: cv.yview_scroll(-1, "units"))
            canvas.bind("<Button-5>", lambda e, cv=canvas: cv.yview_scroll( 1, "units"))

            if len(names) > visible_rows:
                tk.Label(
                    card,
                    text=f"{len(names)} metabolites — scroll to see all",
                    font=("Helvetica", 9),
                    bg=card["bg"],
                    fg="#9ca3af"
                ).pack(anchor="w", pady=(2, 0))

    def _card_peaks(self, card):
        """
        Draw a scrollable peak list for PyAMARES CSV inside a card.
        All peaks are shown — user can scroll to see them all.
        """
        core   = self.state.get('core', {})
        params = core.get('peak_params', {})
        names  = params.get('names', [])
        cs     = params.get('chemicalshift', [])
        amp    = params.get('amplitude', [])
        lw     = params.get('linewidth', [])

        if not names:
            return

        ttk.Separator(card, orient="horizontal").pack(fill="x", pady=(0, 6))

        ############## Column headers #############
        header_row = tk.Frame(card, bg=card["bg"])
        header_row.pack(fill="x", pady=(0, 2))
        for txt, w in [("Name", 12), ("Chem shift", 14),
                       ("Amplitude", 12), ("Linewidth", 10)]:
            tk.Label(
                header_row, text=txt,
                font=("Helvetica", 9), bg=card["bg"], fg="#9ca3af",
                width=w, anchor="w"
            ).pack(side="left")

        ttk.Separator(card, orient="horizontal").pack(fill="x", pady=(0, 3))

        ############## Scrollable canvas containing all peaks #############
        row_height = 18
        visible_rows = 8
        canvas_height = row_height * visible_rows

        canvas_frame = tk.Frame(card, bg=card["bg"])
        canvas_frame.pack(fill="x")

        canvas = tk.Canvas(
            canvas_frame,
            bg=card["bg"],
            height=canvas_height,
            highlightthickness=0
        )
        scrollbar = ttk.Scrollbar(
            canvas_frame,
            orient="vertical",
            command=canvas.yview
        )
        canvas.configure(yscrollcommand=scrollbar.set)

        # Only show scrollbar if there are more peaks than visible rows
        if len(names) > visible_rows:
            scrollbar.pack(side="right", fill="y")

        canvas.pack(side="left", fill="x", expand=True)

        # Inner frame that holds all peak rows
        inner = tk.Frame(canvas, bg=card["bg"])
        canvas_window = canvas.create_window((0, 0), window=inner, anchor="nw")

        ############## Draw every peak #############
        for name, c, a, l in zip(names, cs, amp, lw):
            row = tk.Frame(inner, bg=card["bg"])
            row.pack(fill="x", pady=1)
            tk.Label(row, text=name,           font=("Helvetica", 9),
                     bg=card["bg"], fg="#111",     width=12, anchor="w").pack(side="left")
            tk.Label(row, text=f"{c:.4f} ppm", font=("Helvetica", 9),
                     bg=card["bg"], fg="#374151",  width=14, anchor="w").pack(side="left")
            tk.Label(row, text=f"{a:.4f}",     font=("Helvetica", 9),
                     bg=card["bg"], fg="#374151",  width=12, anchor="w").pack(side="left")
            tk.Label(row, text=f"{l:.1f} Hz",  font=("Helvetica", 9),
                     bg=card["bg"], fg="#374151",  width=10, anchor="w").pack(side="left")

        # Update scroll region once inner frame is fully built
        inner.update_idletasks()
        canvas.configure(scrollregion=canvas.bbox("all"))

        # Make canvas width follow the card width
        def _on_canvas_resize(event):
            canvas.itemconfig(canvas_window, width=event.width)
        canvas.bind("<Configure>", _on_canvas_resize)

        # Mouse wheel scrolling
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        canvas.bind("<MouseWheel>", _on_mousewheel)   # Windows / macOS
        canvas.bind("<Button-4>",                      # Linux scroll up
            lambda e: canvas.yview_scroll(-1, "units"))
        canvas.bind("<Button-5>",                      # Linux scroll down
            lambda e: canvas.yview_scroll(1, "units"))

        # Peak count label below the list
        tk.Label(
            card,
            text=f"{len(names)} peaks total — scroll to see all",
            font=("Helvetica", 9),
            bg=card["bg"],
            fg="#9ca3af"
        ).pack(anchor="w", pady=(4, 0))

    def _show_success_card(self):
        fmt  = self.state.get('detected_format', '')
        card = self._card_frame(border_color="#bbf7d0", bg_color="#f0fdf4")
        self._card_header(
            card,
            format_text=fmt,
            badge_text="✓  Recognized",
            badge_color="#bbf7d0",
            badge_fg="#065f46"
        )
        self._card_files(card, self.state['files'])
        self._card_metadata(card)

        # For PyAMARES, also show the peak table
        if 'PyAMARES' in fmt or fmt.endswith('.csv'):
            self._card_peaks(card)

    def _show_ambiguous_card(self):
        # Determine which formats are plausible
        options = self.state.get('ambiguous_options',
                                  ["MARSS .mat", "FID-A .mat", "MRSCloud .mat"])
        # If user already confirmed a format, show green card
        if self.state.get('format_confirmed'):
            card = self._card_frame(border_color="#bbf7d0", bg_color="#f0fdf4")
            badge_text  = "✓  Confirmed"
            badge_color = "#bbf7d0"
            badge_fg    = "#065f46"
        else:
            card = self._card_frame(border_color="#fde68a", bg_color="#fffbeb")
            badge_text  = "⚠  Confirm format"
            badge_color = "#fde68a"
            badge_fg    = "#92400e"
        self._card_header(
            card,
            format_text="",
            badge_text=badge_text,
            badge_color=badge_color,
            badge_fg=badge_fg,
            show_dropdown=True,
            dropdown_options=options
        )
        self._card_files(card, self.state['files'])
        self._card_metadata(card)

    def _show_sets_card(self):
        """The files hold several basis sets: a dropdown picks the one to convert.
        Once picked, the card stays above the set's own card, to pick another."""
        sets   = self.state['basis_sets']
        picked = self.state.get('set_picked')
        if picked is None:
            card = self._card_frame(border_color="#fde68a", bg_color="#fffbeb")
            badge_text, badge_color, badge_fg = "⚠  Pick a set", "#fde68a", "#92400e"
            note = (f"These files hold {len(sets)} basis sets. Pick the one to "
                    f"convert; only its files are converted.")
        else:
            card = self._card_frame(border_color="#bbf7d0", bg_color="#f0fdf4")
            badge_text, badge_color, badge_fg = "✓  Picked", "#bbf7d0", "#065f46"
            note = (f"Set {picked + 1} of the {len(sets)} these files hold; "
                    f"only its files are converted.")
        self._card_header(
            card,
            format_text="",
            badge_text=badge_text,
            badge_color=badge_color,
            badge_fg=badge_fg,
            show_dropdown=True,
            dropdown_options=[s['name'] for s in sets],
            title="Basis set",
            on_select=self._on_set_picked,
            value="" if picked is None else sets[picked]['name'],
            show_tools=False
        )
        if picked is None:
            self._card_files(card, self.state.get('uploaded_files', []))
            for s in sets:
                tk.Label(
                    card,
                    text=f"{s['name']}:  {s['format']}, "
                         f"{len(s['data'].get('names') or [])} metabolites",
                    font=("Helvetica", 9),
                    bg=card["bg"],
                    fg="#374151"
                ).pack(anchor="w", pady=(4, 0))
        tk.Label(
            card,
            text=note,
            font=("Helvetica", 9),
            bg=card["bg"],
            fg=badge_fg,
            wraplength=560,
            justify="left"
        ).pack(anchor="w", pady=(8, 0))

    def _show_partial_card(self):
        card = self._card_frame(border_color="#fde68a", bg_color="#fffbeb")
        self._card_header(
            card,
            format_text=self.state['detected_format'],
            badge_text="⚠  Partial data",
            badge_color="#fde68a",
            badge_fg="#92400e"
        )
        self._card_files(card, self.state['files'])
        self._card_metadata(card)

        # Warning message
        missing = self.state.get('missing_fields', [])
        msg = (
            f"Format not recognized, but core data was found. "
            f"Missing fields ({', '.join(missing)}) will be asked "
            f"on the next screen."
        )
        tk.Label(
            card,
            text=msg,
            font=("Helvetica", 9),
            bg="#fffbeb",
            fg="#92400e",
            wraplength=560,
            justify="left"
        ).pack(anchor="w", pady=(8, 0))

    def _show_error_card(self):
        card = self._card_frame(border_color="#fecaca", bg_color="#fef2f2")
        self._card_header(
            card,
            format_text="Unreadable",
            badge_text="✗  Cannot proceed",
            badge_color="#fecaca",
            badge_fg="#991b1b"
        )
        self._card_files(card, self.state['files'])

        tk.Label(
            card,
            text="No usable data could be extracted. "
                 "Please check that your files are valid basis set files.",
            font=("Helvetica", 9),
            bg="#fef2f2",
            fg="#991b1b",
            wraplength=560,
            justify="left"
        ).pack(anchor="w", pady=(8, 0))

    ############## Dropdown callbacks #############

    def _on_set_picked(self, event):
        """User picked one of several basis sets: from now on only its files
        and its detected format are in the state."""
        sets  = self.state.get('basis_sets', [])
        names = [s['name'] for s in sets]
        if self.set_var.get() not in names:
            return
        picked = names.index(self.set_var.get())
        self.state['set_picked'] = picked
        self.state['files']      = sets[picked]['files']
        self._apply_detection(sets[picked])
        self._show_result()

    def _on_format_confirmed(self, event):
        """
        User confirmed format from the ambiguous dropdown.
        Re-run detection with the chosen format so metadata is correct,
        then redraw the card green.
        """
        chosen = self.format_var.get()
        self.state['detected_format']  = chosen
        self.state['format_confirmed'] = True

        # Re-read the first file using the chosen format's reader
        # so metadata (sw, sf, Bo, n, names) reflects the correct format
        files = self.state.get('files', [])
        if files:
            result = self._detect_with_format(chosen, files)
            if result:
                self.state['detection_status'] = 'full'
                self.state['core']             = result.get('data', {})
                self.state['missing_fields']   = result.get('missing', [])

        self._show_result()

    def _detect_with_format(self, fmt, files):
        """
        Run the appropriate reader for the user-chosen format.
        Returns a detect-style result dict, or None on failure.
        """
        import os
        try:
            if 'MARSS' in fmt:
                from readers.read_marss import read_marss_file
                r     = read_marss_file(files[0])
                names = [os.path.splitext(os.path.basename(p))[0]
                         for p in files]
                return {
                    'data': {
                        'sw':    r['sw'],
                        'sf':    r['sf'],
                        'n':     r['n'],
                        'Bo':    round(r['sf'] / 42.577, 3),
                        'names': names,
                    },
                    'missing': [],
                }

            if 'FID-A' in fmt:
                from readers.read_fida import read_fida_file
                r     = read_fida_file(files[0])
                names = [os.path.splitext(os.path.basename(p))[0]
                         for p in files]
                return {
                    'data': {
                        'sw':    r['sw'],
                        'sf':    r['sf'],
                        'n':     r['n'],
                        'Bo':    r.get('Bo'),
                        'names': names,
                    },
                    'missing': [] if r.get('te') else ['te'],
                }

            if 'MRSCloud' in fmt:
                from readers.read_mrscloud import read_mrscloud_file
                r     = read_mrscloud_file(files[0])
                names = [os.path.splitext(os.path.basename(p))[0]
                         for p in files]
                return {
                    'data': {
                        'sw':    r['sw'],
                        'sf':    r['sf'],
                        'n':     r['n'],
                        'Bo':    r.get('Bo'),
                        'names': names,
                    },
                    'missing': [],
                }

            if 'Osprey' in fmt:
                from readers.read_osprey import read_osprey
                results = read_osprey(files[0])
                return {
                    'data': {
                        'sw':    results[0]['sw'],
                        'sf':    results[0]['sf'],
                        'n':     results[0]['n'],
                        'Bo':    results[0].get('Bo'),
                        'names': [x['name'] for x in results],
                    },
                    'missing': [],
                }

            if 'INSPECTOR' in fmt:
                from readers.read_inspector import read_inspector
                results = read_inspector(files[0])
                return {
                    'data': {
                        'sw':    results[0]['sw'],
                        'sf':    results[0]['sf'],
                        'n':     results[0]['n'],
                        'Bo':    round(results[0]['sf'] / 42.577, 3),
                        'names': [x['name'] for x in results],
                    },
                    'missing': [],
                }

            if 'LCModel .raw' in fmt:
                from readers.read_lcmodel import read_lcmodel_raw
                r     = read_lcmodel_raw(files[0])
                names = [os.path.splitext(os.path.basename(p))[0]
                         for p in files]
                return {
                    'data': {
                        'sw':    r['sw'],
                        'sf':    r['sf'],
                        'n':     r['n'],
                        'Bo':    round(r['sf'] / 42.577, 3) if r['sf'] else None,
                        'names': names,
                    },
                    'missing': ['te'],
                }

            if 'LCModel .BASIS' in fmt:
                from readers.read_lcmodel import read_lcmodel_basis
                results = read_lcmodel_basis(files[0])
                return {
                    'data': {
                        'sw':    results[0]['sw'],
                        'sf':    results[0]['sf'],
                        'n':     results[0]['n'],
                        'Bo':    round(results[0]['sf'] / 42.577, 3),
                        'names': [x['name'] for x in results],
                    },
                    'missing': [] if results[0].get('te') else ['te'],
                }

        except Exception as ex:
            print(f"  Warning: could not re-read with format '{fmt}': {ex}")
        return None

    ############## Clear #############

    def _clear_all(self):
        """Reset screen and shared state."""
        for key in ['files', 'detected_format', 'detection_status',
                    'core', 'missing_fields', 'ambiguous_options',
                    'compatible_tools', 'format_confirmed',
                    'uploaded_files', 'basis_sets', 'set_picked']:
            self.state.pop(key, None)
        self._clear_result()
        self._disable_next()

    ############## Utilities #############

    def _file_size_str(self, path):
        """Return human-readable file size string."""
        try:
            size = os.path.getsize(path)
            if size < 1024:
                return f"{size} B"
            elif size < 1024 ** 2:
                return f"{size / 1024:.1f} KB"
            else:
                return f"{size / 1024 ** 2:.1f} MB"
        except OSError:
            return ""