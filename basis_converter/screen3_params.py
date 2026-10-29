"""
screen3_params.py
Basis Set Converter — Screen 3: Missing parameters

This screen only appears when required fields could not be found
in the uploaded basis set. It asks the user to supply them.

Required fields that may be missing:
    sw   - spectral width in Hz
    sf   - Larmor frequency in MHz
    te   - echo time in ms
    Bo   - field strength in Tesla (derived from sf if present)

Fields are split into three tiers:
    Required  - must be filled before converting
    Optional  - defaults available, shown in amber
    Derived   - calculated automatically, shown in green

State written by this screen:
    state['params'] : dict of all parameter values (from file + user input)
"""

import tkinter as tk
from tkinter import ttk

from base_screen import BaseScreen


############## Parameter definitions #############

# Each param: (label, unit, default, description, tier)
# tier: 'required' | 'optional' | 'derived'

PARAMS = [
    {
        'key'    : 'sw',
        'label'  : 'Spectral width (SW)',
        'unit'   : 'Hz',
        'default': '4000',
        'desc'   : 'Bandwidth of the acquisition.',
        'tier'   : 'required',
    },
    {
        'key'    : 'sf',
        'label'  : 'Larmor frequency',
        'unit'   : 'MHz',
        'default': '123.26',
        'desc'   : 'Proton Larmor frequency at your field strength.',
        'tier'   : 'required',
    },
    {
        'key'    : 'te',
        'label'  : 'Echo time (TE)',
        'unit'   : 'ms',
        'default': '30',
        'desc'   : 'Echo time used during simulation or acquisition.',
        'tier'   : 'optional',
    },
    {
        'key'    : 'seq',
        'label'  : 'Sequence',
        'unit'   : '',
        'default': 'unedited',
        'desc'   : 'Pulse sequence name.',
        'tier'   : 'optional',
        'choices': ['unedited', 'sLASER', 'PRESS', 'STEAM',
                    'LASER', 'MEGA-PRESS', 'MEGA-sLASER'],
    },
]


class Screen3Params(BaseScreen):

    current_step = 3

    def build_body(self):
        self._entries   = {}   # key → tk.StringVar
        self._tier_data = []   # list of (key, tier, var)
        self._summary   = {}   # key → its value label in the summary card

        # Scrollable body
        outer = tk.Frame(self, bg="white")
        outer.pack(fill="both", expand=True)

        canvas    = tk.Canvas(outer, bg="white", highlightthickness=0)
        scrollbar = ttk.Scrollbar(outer, orient="vertical",
                                  command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        self.scroll_frame = tk.Frame(canvas, bg="white", padx=24, pady=12)
        self._cw = canvas.create_window(
            (0, 0), window=self.scroll_frame, anchor="nw"
        )

        self.scroll_frame.bind("<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>",
            lambda e: canvas.itemconfig(self._cw, width=e.width))
        canvas.bind("<MouseWheel>",
            lambda e: canvas.yview_scroll(int(-1*(e.delta/120)), "units"))

        self._build_content()
        self._build_footer()

    def on_show(self):
        """Called every time this screen is shown — rebuild with current state."""
        for w in self.scroll_frame.winfo_children():
            w.destroy()
        self._entries   = {}
        self._tier_data = []
        self._summary   = {}
        self._build_content()

    ############## Content #############

    def _build_content(self):
        missing  = self.state.get('missing_fields', [])
        core     = self.state.get('core', {})
        in_fmt   = self.state.get('detected_format', '')
        out_tool = self.state.get('output_tool', '')

        # Context bar
        ctx = tk.Frame(self.scroll_frame, bg="#f3f4f6", padx=12, pady=8)
        ctx.pack(fill="x", pady=(0, 12))
        tk.Label(
            ctx,
            text=f"{in_fmt}  →  {out_tool}",
            font=("Helvetica", 10),
            bg="#f3f4f6",
            fg="#374151"
        ).pack(side="left")

        # Determine which params to show
        required_params = []
        optional_params = []
        derived_params  = []

        for p in PARAMS:
            key = p['key']
            # Already in core and not in missing → derived / already known
            if key in core and core[key] and key not in missing:
                derived_params.append(p)
            elif key in missing or core.get(key) is None:
                if p['tier'] == 'required':
                    required_params.append(p)
                else:
                    optional_params.append(p)
            else:
                optional_params.append(p)

        # If nothing is missing — show a confirmation summary
        if not required_params and not optional_params:
            self._build_all_found(core, derived_params)
        else:
            # Intro message
            msg = (
                "Some parameters could not be found in your file. "
                "Required fields must be filled in before converting."
            )
            msg_frame = tk.Frame(
                self.scroll_frame, bg="#eff6ff",
                highlightbackground="#bfdbfe", highlightthickness=1,
                padx=12, pady=8
            )
            msg_frame.pack(fill="x", pady=(0, 12))
            tk.Label(
                msg_frame,
                text=msg,
                font=("Helvetica", 10),
                bg="#eff6ff", fg="#1d4ed8",
                wraplength=560, justify="left"
            ).pack(anchor="w")

            if required_params:
                self._section_label("Required", "#dc2626")
                for p in required_params:
                    self._build_field(p, tier='required')

            if optional_params:
                self._section_label("Optional — defaults available", "#92400e")
                for p in optional_params:
                    self._build_field(p, tier='optional')

            if derived_params:
                self._section_label("Found in file — no action needed", "#065f46")
                for p in derived_params:
                    self._build_derived(p, core)

        # Summary card
        self._build_summary(core, required_params, optional_params,
                            derived_params)

    def _section_label(self, text, color):
        tk.Label(
            self.scroll_frame,
            text=text,
            font=("Helvetica", 10, "bold"),
            bg="white",
            fg=color
        ).pack(anchor="w", pady=(12, 4))

    def _build_field(self, p, tier):
        """Build an input field for a missing parameter."""
        key = p['key']

        # Card
        border_color = "#fecaca" if tier == 'required' else "#fde68a"
        bg_color     = "#fef2f2" if tier == 'required' else "#fffbeb"

        card = tk.Frame(
            self.scroll_frame,
            bg=bg_color,
            highlightbackground=border_color,
            highlightthickness=1,
            padx=14, pady=10
        )
        card.pack(fill="x", pady=3)

        # Header row
        hdr = tk.Frame(card, bg=bg_color)
        hdr.pack(fill="x")

        tk.Label(
            hdr,
            text=p['label'],
            font=("Helvetica", 11, "bold"),
            bg=bg_color, fg="#111827"
        ).pack(side="left")

        badge_text  = "Required" if tier == 'required' else "Optional"
        badge_bg    = "#fecaca" if tier == 'required' else "#fde68a"
        badge_fg    = "#991b1b" if tier == 'required' else "#92400e"
        tk.Label(
            hdr,
            text=badge_text,
            font=("Helvetica", 9),
            bg=badge_bg, fg=badge_fg,
            padx=8, pady=2
        ).pack(side="right")

        tk.Label(
            card,
            text=p['desc'],
            font=("Helvetica", 9),
            bg=bg_color, fg="#6b7280"
        ).pack(anchor="w", pady=(2, 6))

        # Input row
        input_row = tk.Frame(card, bg=bg_color)
        input_row.pack(fill="x")

        var = tk.StringVar(value=self.state.get('params', {}).get(
            key, p['default']))
        self._entries[key] = var

        if 'choices' in p:
            # Dropdown
            combo = ttk.Combobox(
                input_row,
                textvariable=var,
                values=p['choices'],
                state="readonly",
                font=("Helvetica", 10),
                width=20
            )
            combo.pack(side="left")
        else:
            # Text entry
            entry = tk.Entry(
                input_row,
                textvariable=var,
                font=("Helvetica", 11),
                width=14,
                relief="flat",
                highlightbackground="#d1d5db",
                highlightthickness=1
            )
            entry.pack(side="left")

        if p['unit']:
            tk.Label(
                input_row,
                text=p['unit'],
                font=("Helvetica", 10),
                bg=bg_color, fg="#6b7280"
            ).pack(side="left", padx=(6, 0))

        if tier == 'optional' and p.get('default'):
            tk.Label(
                card,
                text=f"Default: {p['default']} {p['unit']}".strip(),
                font=("Helvetica", 9),
                bg=bg_color, fg="#9ca3af"
            ).pack(anchor="w", pady=(4, 0))

    def _build_derived(self, p, core):
        """Build a read-only card for a field already found in the file."""
        key   = p['key']
        value = core.get(key)
        if value is None:
            return

        card = tk.Frame(
            self.scroll_frame,
            bg="#f0fdf4",
            highlightbackground="#bbf7d0",
            highlightthickness=1,
            padx=14, pady=8
        )
        card.pack(fill="x", pady=3)

        hdr = tk.Frame(card, bg="#f0fdf4")
        hdr.pack(fill="x")

        tk.Label(
            hdr,
            text=p['label'],
            font=("Helvetica", 11, "bold"),
            bg="#f0fdf4", fg="#111827"
        ).pack(side="left")

        tk.Label(
            hdr,
            text="Found in file",
            font=("Helvetica", 9),
            bg="#bbf7d0", fg="#065f46",
            padx=8, pady=2
        ).pack(side="right")

        val_str = f"{value:.4f} {p['unit']}".strip() \
                  if isinstance(value, float) else f"{value} {p['unit']}".strip()
        tk.Label(
            card,
            text=val_str,
            font=("Helvetica", 12, "bold"),
            bg="#f0fdf4", fg="#065f46"
        ).pack(anchor="w", pady=(4, 0))

    def _build_all_found(self, core, derived):
        """Show a green confirmation when nothing is missing."""
        card = tk.Frame(
            self.scroll_frame,
            bg="#f0fdf4",
            highlightbackground="#bbf7d0",
            highlightthickness=1,
            padx=14, pady=12
        )
        card.pack(fill="x", pady=(0, 12))

        tk.Label(
            card,
            text="✓  All parameters found",
            font=("Helvetica", 12, "bold"),
            bg="#f0fdf4", fg="#065f46"
        ).pack(anchor="w")

        tk.Label(
            card,
            text="Nothing extra is needed — you can proceed directly to conversion.",
            font=("Helvetica", 10),
            bg="#f0fdf4", fg="#065f46"
        ).pack(anchor="w", pady=(4, 0))

    def _build_summary(self, core, required, optional, derived):
        """Parameter summary card at the bottom."""
        ttk.Separator(self.scroll_frame, orient="horizontal").pack(
            fill="x", pady=(16, 8)
        )

        tk.Label(
            self.scroll_frame,
            text="Parameters summary",
            font=("Helvetica", 10, "bold"),
            bg="white", fg="#374151"
        ).pack(anchor="w", pady=(0, 6))

        grid = tk.Frame(self.scroll_frame, bg="white")
        grid.pack(fill="x")

        all_params = required + optional + derived
        for i, p in enumerate(all_params):
            key = p['key']
            row = i // 2
            col = i % 2

            cell = tk.Frame(grid, bg="white")
            cell.grid(row=row, column=col, sticky="w",
                      padx=(0, 24), pady=3)

            tk.Label(
                cell,
                text=p['label'],
                font=("Helvetica", 9),
                bg="white", fg="#9ca3af"
            ).pack(anchor="w")

            self._summary[key] = tk.Label(
                cell,
                font=("Helvetica", 10, "bold"),
                bg="white"
            )
            self._summary[key].pack(anchor="w")
            self._show_summary_value(p, core)

            # Follow what the user types or picks
            if key in self._entries:
                self._entries[key].trace_add(
                    'write', lambda *_, p=p: self._show_summary_value(p, core))

        for c in range(2):
            grid.columnconfigure(c, weight=1)

    def _show_summary_value(self, p, core):
        """Write a parameter's current value into its label in the summary card."""
        key = p['key']
        if key in self._entries:
            val_str = self._entries[key].get() or "—"
            color   = "#f59e0b"   # amber for user-supplied
        elif key in core and core[key]:
            v = core[key]
            val_str = f"{v:.4f} {p['unit']}".strip() \
                      if isinstance(v, float) else f"{v} {p['unit']}".strip()
            color   = "#065f46"   # green for found in file
        else:
            val_str = "—"
            color   = "#9ca3af"

        self._summary[key].config(text=val_str, fg=color)

    ############## Footer #############

    def _build_footer(self):
        footer = tk.Frame(self, bg="white", pady=12, padx=20)
        footer.pack(fill="x", side="bottom")
        ttk.Separator(self, orient="horizontal").pack(
            fill="x", side="bottom", before=footer
        )

        tk.Button(
            footer, text="← Back",
            command=lambda: self.app.show_screen("screen2"),
            relief="flat", bg="#f3f4f6", fg="#374151",
            padx=14, pady=6, cursor="hand2"
        ).pack(side="left")

        tk.Button(
            footer, text="Convert →",
            command=self._on_convert,
            relief="flat", bg="#dbeafe", fg="#1d4ed8",
            padx=14, pady=6, cursor="hand2"
        ).pack(side="right")

    def _on_convert(self):
        """Collect all parameter values and proceed to screen 4."""
        core   = self.state.get('core', {}).copy()
        params = self.state.get('params', {}).copy()

        # Merge user-supplied values into params
        for key, var in self._entries.items():
            val = var.get().strip()
            if val:
                try:
                    params[key] = float(val)
                except ValueError:
                    params[key] = val   # keep as string for seq etc.

        # Merge params into core (params override)
        for key, val in params.items():
            if val:
                core[key] = val

        # Derive Bo from sf if not present
        if core.get('sf') and not core.get('Bo'):
            core['Bo'] = float(core['sf']) / 42.577

        # Derive sw from dt if not present
        if core.get('dt') and not core.get('sw'):
            core['sw'] = 1.0 / float(core['dt'])

        self.state['core']   = core
        self.state['params'] = params

        self.app.show_screen("screen4")
