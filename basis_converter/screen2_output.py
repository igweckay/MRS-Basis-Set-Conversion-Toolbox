"""
screen2_output.py
Basis Set Converter — Screen 2: Output tool selection

The user selects which software tool they want to convert to.
Tools are grouped into three sections:
    1. Fitting and quantification tools
    2. Basis set generators
    3. Peak fitting tools (PyAMARES, OXSA)

Tools with multiple output format options show a sub-panel when selected.

State written by this screen:
    state['output_tool']    : str, e.g. 'LCModel'
    state['output_tool_id'] : str, e.g. 'lcmodel'
    state['output_format']  : str, e.g. 'both' / 'mat' / 'default'
"""

import tkinter as tk
from tkinter import ttk
from base_screen import BaseScreen
from core.formats import FITTING_TOOLS, GENERATOR_TOOLS, PEAK_FITTING_TOOLS


class Screen2Output(BaseScreen):

    current_step = 2

    def build_body(self):
        self.selected_tool   = None
        self.selected_format = None
        self._tool_buttons   = {}

        # Scrollable area
        outer = tk.Frame(self, bg="white")
        outer.pack(fill="both", expand=True)

        self._canvas = tk.Canvas(outer, bg="white", highlightthickness=0)
        scrollbar = ttk.Scrollbar(outer, orient="vertical",
                                  command=self._canvas.yview)
        self._canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        self._canvas.pack(side="left", fill="both", expand=True)

        self.scroll_frame = tk.Frame(self._canvas, bg="white", padx=24, pady=12)
        self._cw = self._canvas.create_window(
            (0, 0), window=self.scroll_frame, anchor="nw"
        )

        self.scroll_frame.bind("<Configure>", self._on_frame_configure)
        self._canvas.bind("<Configure>", self._on_canvas_configure)
        self._canvas.bind("<MouseWheel>",
            lambda e: self._canvas.yview_scroll(int(-1*(e.delta/120)), "units"))
        self._canvas.bind("<Button-4>",
            lambda e: self._canvas.yview_scroll(-1, "units"))
        self._canvas.bind("<Button-5>",
            lambda e: self._canvas.yview_scroll(1, "units"))

        self._build_sections()
        self._build_format_panel()
        self._build_info_panel()
        self._build_footer()

    def _on_frame_configure(self, event):
        self._canvas.configure(scrollregion=self._canvas.bbox("all"))

    def _on_canvas_configure(self, event):
        self._canvas.itemconfig(self._cw, width=event.width)

    def on_show(self):
        pass

    ############## Sections #############

    def _build_sections(self):
        self._build_section("Fitting and quantification tools", FITTING_TOOLS)
        ttk.Separator(self.scroll_frame, orient="horizontal").pack(fill="x", pady=10)
        self._build_section("Basis set generators", GENERATOR_TOOLS)
        ttk.Separator(self.scroll_frame, orient="horizontal").pack(fill="x", pady=10)
        self._build_section(
            "Peak fitting tools  —  produce prior knowledge, not basis sets",
            PEAK_FITTING_TOOLS
        )

    def _build_section(self, label, tools):
        tk.Label(
            self.scroll_frame,
            text=label,
            font=("Helvetica", 10, "bold"),
            bg="white",
            fg="#6b7280"
        ).pack(anchor="w", pady=(4, 6))

        grid = tk.Frame(self.scroll_frame, bg="white")
        grid.pack(fill="x")

        col_count = 4
        for i, tool in enumerate(tools):
            self._build_tool_card(grid, tool, i // col_count, i % col_count)
        for c in range(col_count):
            grid.columnconfigure(c, weight=1)

    def _build_tool_card(self, parent, tool, row, col):
        is_dimmed = tool['dimmed']
        bg  = "#f3f4f6" if is_dimmed else "white"
        fg  = "#9ca3af" if is_dimmed else "#111827"

        frame = tk.Frame(
            parent, bg=bg, relief="flat", bd=0,
            highlightbackground="#e5e7eb",
            highlightthickness=1,
            padx=10, pady=8
        )
        frame.grid(row=row, column=col, padx=4, pady=4, sticky="nsew")

        name_lbl = tk.Label(frame, text=tool['name'],
                            font=("Helvetica", 11, "bold"), bg=bg, fg=fg)
        name_lbl.pack(anchor="w")

        ext_lbl = tk.Label(frame, text=tool['ext'],
                           font=("Helvetica", 9), bg=bg, fg="#9ca3af")
        ext_lbl.pack(anchor="w")

        if not is_dimmed:
            for widget in [frame, name_lbl, ext_lbl]:
                widget.configure(cursor="hand2")
                widget.bind("<Button-1>",
                    lambda e, t=tool: self._on_tool_selected(t))

        self._tool_buttons[tool['id']] = frame

    ############## Format sub-panel #############

    def _build_format_panel(self):
        self.format_panel = tk.Frame(
            self.scroll_frame, bg="#eff6ff",
            highlightbackground="#bfdbfe", highlightthickness=1,
            padx=14, pady=10
        )

        self.format_panel_title = tk.Label(
            self.format_panel, text="Output format",
            font=("Helvetica", 10, "bold"), bg="#eff6ff", fg="#1d4ed8"
        )
        self.format_panel_title.pack(anchor="w", pady=(0, 6))

        self.format_options_frame = tk.Frame(self.format_panel, bg="#eff6ff")
        self.format_options_frame.pack(fill="x")

        self._format_var = tk.StringVar()

    def _build_info_panel(self):
        self.info_panel = tk.Frame(
            self.scroll_frame, bg="#fefce8",
            highlightbackground="#fde047", highlightthickness=1,
            padx=14, pady=10
        )
        self.info_panel_title = tk.Label(
            self.info_panel, text="",
            font=("Helvetica", 10, "bold"),
            bg="#fefce8", fg="#854d0e"
        )
        self.info_panel_title.pack(anchor="w")
        self.info_panel_body = tk.Label(
            self.info_panel, text="",
            font=("Helvetica", 9),
            bg="#fefce8", fg="#713f12",
            wraplength=520, justify="left"
        )
        self.info_panel_body.pack(anchor="w", pady=(4, 0))

    def _show_info_panel(self, title, body):
        self._hide_format_panel()
        self.info_panel_title.config(text=title)
        self.info_panel_body.config(text=body)
        self.info_panel.pack(fill="x", pady=(8, 0))

    def _hide_info_panel(self):
        self.info_panel.pack_forget()

    def _show_format_panel(self, tool):
        for w in self.format_options_frame.winfo_children():
            w.destroy()

        self.format_panel_title.config(
            text=f"Output format for {tool['name']}"
        )
        self._format_var.set('')

        for fmt in tool['formats']:
            row = tk.Frame(self.format_options_frame, bg="#eff6ff")
            row.pack(fill="x", pady=2)

            tk.Radiobutton(
                row, text=fmt['label'],
                variable=self._format_var, value=fmt['id'],
                font=("Helvetica", 10), bg="#eff6ff", fg="#111827",
                activebackground="#eff6ff", cursor="hand2",
                command=lambda fid=fmt['id']: self._on_format_selected(fid)
            ).pack(side="left")

            tk.Label(
                row, text=f"  —  {fmt['sub']}",
                font=("Helvetica", 9), bg="#eff6ff", fg="#6b7280"
            ).pack(side="left")

        self.format_panel.pack(fill="x", pady=(8, 0))

    def _hide_format_panel(self):
        self.format_panel.pack_forget()
        if hasattr(self, 'info_panel'):
            self._hide_info_panel()

    ############## Selection handlers #############

    def _on_tool_selected(self, tool):
        # Deselect previous
        if self.selected_tool:
            prev = self._tool_buttons.get(self.selected_tool['id'])
            if prev:
                prev.config(bg="white", highlightbackground="#e5e7eb")
                for c in prev.winfo_children():
                    c.config(bg="white")

        # Toggle off
        if self.selected_tool and self.selected_tool['id'] == tool['id']:
            self.selected_tool   = None
            self.selected_format = None
            self._hide_format_panel()
            self._disable_next()
            return

        # Select
        self.selected_tool   = tool
        self.selected_format = None

        card = self._tool_buttons[tool['id']]
        card.config(bg="#dbeafe", highlightbackground="#3b82f6")
        for c in card.winfo_children():
            c.config(bg="#dbeafe")

        # Tools that generate scripts/priors rather than basis set files
        INFO_TOOLS = {
            'oxsa': (
                "OXSA — MATLAB script output",
                "Selecting OXSA will produce a pre-configured run_oxsa_amares.m "
                "script with your metabolite names, PPM positions, linewidth bounds "
                "and relative amplitudes already filled in from your basis set. "
                "Open the script in MATLAB, set your OXSA path and subject file, and run."
            ),
            'vespa_analysis': (
                "VeSPA Analysis — prior XML output",
                "Selecting VeSPA Analysis will produce a VIFF XML prior file "
                "loadable via 'Basis Set from File' in VeSPA Analysis. "
                "Peak positions are taken from literature chemical shifts "
                "and areas from the simulated FID at your echo time."
            ),
            'pyamares': (
                "PyAMARES — CSV prior output",
                "Selecting PyAMARES will produce a CSV prior knowledge file "
                "with metabolite names and dominant peak PPM positions "
                "estimated from the basis set."
            ),
            'midas': (
                "MIDAS — prior XML output",
                "Selecting MIDAS will produce a FITT_Generic_XML prior file "
                "compatible with the MIDAS FITT2 application. "
                "Peak positions and relative areas are derived from the basis set."
            ),
            'gava': (
                "GAVA — text prior output",
                "Selecting GAVA will produce a tab-separated text prior file "
                "with metabolite peak positions and areas from the basis set."
            ),
        }

        if tool['id'] in INFO_TOOLS:
            title, body = INFO_TOOLS[tool['id']]
            self._show_info_panel(title, body)
            self._enable_next()
        elif tool['formats']:
            self._hide_info_panel()
            self._show_format_panel(tool)
            self._disable_next()
        else:
            self._hide_format_panel()
            self.selected_format = 'default'
            self._enable_next()

    def _on_format_selected(self, fmt_id):
        self.selected_format = fmt_id
        self._enable_next()

    ############## Footer #############

    def _build_footer(self):
        footer = tk.Frame(self, bg="white", pady=12, padx=20)
        footer.pack(fill="x", side="bottom")
        ttk.Separator(self, orient="horizontal").pack(
            fill="x", side="bottom", before=footer
        )

        tk.Button(
            footer, text="← Back",
            command=lambda: self.app.show_screen("screen1"),
            relief="flat", bg="#f3f4f6", fg="#374151",
            padx=14, pady=6, cursor="hand2"
        ).pack(side="left")

        self.next_btn = tk.Button(
            footer, text="Next →",
            command=self._on_next,
            relief="flat", bg="#e5e7eb", fg="#9ca3af",
            padx=14, pady=6, state="disabled"
        )
        self.next_btn.pack(side="right")

    def _enable_next(self):
        self.next_btn.config(bg="#dbeafe", fg="#1d4ed8",
                             state="normal", cursor="hand2")

    def _disable_next(self):
        self.next_btn.config(bg="#e5e7eb", fg="#9ca3af",
                             state="disabled", cursor="")

    def _on_next(self):
        if not self.selected_tool:
            return
        self.state['output_tool']    = self.selected_tool['name']
        self.state['output_tool_id'] = self.selected_tool['id']
        self.state['output_format']  = self.selected_format or 'default'
        self.app.show_screen("screen3")