"""Full-width tab bar in the HCI style: navy band, white labels, accent
underline on the active tab. Same small API as CTkTabview (add/get/set plus a
click ``command``) so app.py's call sites stay unchanged.

Every tab body stays gridded in one cell and switching just tkraise()s the
chosen one -- unmapping/mapping a tab's whole subtree is a real per-widget
Cocoa operation on macOS (50-211ms stalls per switch), which raising avoids.
"""
from __future__ import annotations

from collections.abc import Callable

import customtkinter as ctk

from ui import theme


class HciTabview(ctk.CTkFrame):
    def __init__(self, master, title: str = "", logo=None, command: Callable[[], None] | None = None):
        super().__init__(master, fg_color="transparent", corner_radius=0, border_width=0)
        self._command = command
        self._frames: dict[str, ctk.CTkFrame] = {}
        self._buttons: dict[str, ctk.CTkButton] = {}
        self._underlines: dict[str, ctk.CTkFrame] = {}
        self._current: str | None = None

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)

        self._header = ctk.CTkFrame(
            self, fg_color=theme.pair("chrome"), corner_radius=0, border_width=0,
        )
        self._header.grid(row=0, column=0, sticky="ew")
        self._header.grid_columnconfigure(1, weight=1)

        if logo is not None:
            ctk.CTkLabel(self._header, image=logo, text="").grid(row=0, column=0, padx=(16, 6), pady=6)
        elif title:
            ctk.CTkLabel(
                self._header, text=title, text_color=theme.pair("chrome_fg"),
                font=ctk.CTkFont(size=15, weight="bold"),
            ).grid(row=0, column=0, padx=(18, 20), sticky="w")

        self._tab_strip = ctk.CTkFrame(
            self._header, fg_color="transparent", corner_radius=0, border_width=0,
        )
        self._tab_strip.grid(row=0, column=1, padx=(10, 0), sticky="w")

        # Right-hand slot for header controls (e.g. the light/dark toggle)
        self.header_right = ctk.CTkFrame(
            self._header, fg_color="transparent", corner_radius=0, border_width=0,
        )
        self.header_right.grid(row=0, column=2, padx=(0, 14), sticky="e")

        self._body = ctk.CTkFrame(
            self, fg_color=theme.pair("bg"), corner_radius=0, border_width=0,
        )
        self._body.grid(row=1, column=0, sticky="nsew")
        self._body.grid_columnconfigure(0, weight=1)
        self._body.grid_rowconfigure(0, weight=1)

    def add(self, name: str) -> ctk.CTkFrame:
        frame = ctk.CTkFrame(
            self._body, fg_color=theme.pair("bg"), corner_radius=0, border_width=0,
        )
        frame.grid(row=0, column=0, sticky="nsew")
        self._frames[name] = frame

        cell = ctk.CTkFrame(
            self._tab_strip, fg_color="transparent", corner_radius=0, border_width=0, width=0, height=0,
        )
        cell.pack(side="left")
        btn = ctk.CTkButton(
            cell, text=name, corner_radius=0, height=40, width=len(name) * 8 + 28,
            fg_color="transparent", hover_color=theme.pair("chrome_hover"),
            text_color=theme.pair("chrome_fg"), font=ctk.CTkFont(size=13),
            border_width=0, command=lambda n=name: self._on_click(n),
        )
        btn.pack(side="top", padx=2)
        underline = ctk.CTkFrame(
            cell, width=0, height=3, corner_radius=0, border_width=0, fg_color="transparent",
        )
        underline.pack(side="top", fill="x", padx=10)
        self._buttons[name] = btn
        self._underlines[name] = underline

        if self._current is None:
            self.set(name)
        else:
            # A newly created frame stacks above the others; keep the current tab on top
            self._frames[self._current].tkraise()
        return frame

    def get(self) -> str:
        return self._current or ""

    def set(self, name: str) -> None:
        if name not in self._frames:
            raise ValueError(f"HciTabview has no tab named '{name}'")
        self._current = name
        self._frames[name].tkraise()
        for tab_name, underline in self._underlines.items():
            active = tab_name == name
            underline.configure(fg_color=theme.pair("tab_accent") if active else "transparent")
            self._buttons[tab_name].configure(
                font=ctk.CTkFont(size=13, weight="bold" if active else "normal"),
            )

    def _on_click(self, name: str) -> None:
        self.set(name)
        if self._command is not None:
            self._command()
