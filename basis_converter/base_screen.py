"""
base_screen.py
Basis Set Converter — base class for all screens

All four screens inherit from BaseScreen.
Provides the shared header, step indicator, and footer helpers.

Kept in its own file so that main.py and screen files can both
import from it without creating a circular import.
"""

import tkinter as tk
from tkinter import ttk


class BaseScreen(tk.Frame):
    """
    Parent class for all screens.
    Subclasses must set current_step (1–4) and implement build_body().
    """

    STEPS = ["Upload", "Output", "Parameters", "Convert"]
    current_step = 1          # override in each subclass

    def __init__(self, parent, state, app):
        super().__init__(parent, bg="white")
        self.state = state        # shared dict
        self.app   = app          # reference to App for navigation

        self._build_header()
        self._build_step_bar()
        self.build_body()

    ################ Header #####################

    def _build_header(self):
        header = tk.Frame(self, bg="white", pady=12, padx=20)
        header.pack(fill="x")

        tk.Label(
            header,
            text="Basis Set Converter",
            font=("Helvetica", 16, "bold"),
            bg="white",
            fg="#111"
        ).pack(side="left")

        tk.Label(
            header,
            text="Convert between MRS basis set formats",
            font=("Helvetica", 10),
            bg="white",
            fg="#888"
        ).pack(side="left", padx=(10, 0))

        # Exit button — top right of header
        tk.Button(
            header,
            text="✕  Exit",
            command=self.app.quit,
            relief="flat",
            bg="#fef2f2",
            fg="#991b1b",
            padx=10, pady=4,
            cursor="hand2"
        ).pack(side="right")

        ttk.Separator(self, orient="horizontal").pack(fill="x")

    ##################### Step indicator #####################

    def _build_step_bar(self):
        bar = tk.Frame(self, bg="white", pady=10, padx=20)
        bar.pack(fill="x")

        for i, step in enumerate(self.STEPS, start=1):
            if i < self.current_step:
                circle_bg, circle_fg, label_fg = "#d1fae5", "#065f46", "#666"
                circle_text = "✓"
            elif i == self.current_step:
                circle_bg, circle_fg, label_fg = "#dbeafe", "#1d4ed8", "#111"
                circle_text = str(i)
            else:
                circle_bg, circle_fg, label_fg = "#f3f4f6", "#9ca3af", "#9ca3af"
                circle_text = str(i)

            step_frame = tk.Frame(bar, bg="white")
            step_frame.pack(side="left")

            tk.Label(
                step_frame,
                text=circle_text,
                font=("Helvetica", 9, "bold"),
                bg=circle_bg,
                fg=circle_fg,
                width=2,
                relief="flat"
            ).pack(side="left")

            tk.Label(
                step_frame,
                text=f" {step}",
                font=("Helvetica", 10,
                      "bold" if i == self.current_step else "normal"),
                bg="white",
                fg=label_fg
            ).pack(side="left")

            if i < len(self.STEPS):
                tk.Label(
                    bar,
                    text="  ——  ",
                    font=("Helvetica", 9),
                    bg="white",
                    fg="#d1d5db"
                ).pack(side="left")

        ttk.Separator(self, orient="horizontal").pack(fill="x")

    ##################### Footer helper #####################

    def make_footer(self, parent, back_screen=None, next_screen=None,
                    next_text="Next →", back_text="← Back",
                    next_command=None):
        """
        Standard footer with Back and Next buttons.
        Pass back_screen / next_screen as screen name strings,
        or pass next_command for custom behavior (e.g. start conversion).
        """
        footer = tk.Frame(parent, bg="white", pady=12, padx=20)
        footer.pack(fill="x", side="bottom")

        ttk.Separator(parent, orient="horizontal").pack(
            fill="x", side="bottom", before=footer
        )

        if back_screen:
            tk.Button(
                footer,
                text=back_text,
                command=lambda: self.app.show_screen(back_screen),
                relief="flat",
                bg="#f3f4f6",
                fg="#374151",
                padx=14, pady=6,
                cursor="hand2"
            ).pack(side="left")

        if next_command:
            tk.Button(
                footer,
                text=next_text,
                command=next_command,
                relief="flat",
                bg="#dbeafe",
                fg="#1d4ed8",
                padx=14, pady=6,
                cursor="hand2"
            ).pack(side="right")
        elif next_screen:
            tk.Button(
                footer,
                text=next_text,
                command=lambda: self.app.show_screen(next_screen),
                relief="flat",
                bg="#dbeafe",
                fg="#1d4ed8",
                padx=14, pady=6,
                cursor="hand2"
            ).pack(side="right")

    ##################### Subclass interface #####################

    def build_body(self):
        """Override in each screen to add screen-specific content."""
        raise NotImplementedError

    def on_show(self):
        """
        Called every time this screen is brought to the front.
        Override to refresh content when navigating back/forward.
        """
        pass
