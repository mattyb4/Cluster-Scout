"""Huntsman Cancer Institute color theme (light + dark) for the whole UI.

CustomTkinter widgets take ``(light, dark)`` color tuples and re-theme
themselves on ``set_appearance_mode``. Everything CustomTkinter doesn't manage
(ttk Treeview, plain Tk frames/canvases, matplotlib figures, the HTML help
pane) is registered with :func:`on_change` and re-themed from a callback when
the mode flips.

Brand colors (huntsmancancer.org): primary #285087 and white; bright
secondaries #ecd451 #f2ad49 #e65032 #96d2ff #009f7e; neutral secondaries
#f0f0f0 #707271 #505050 #5e4232 #890000. A few text-only variants below are
darkened/lightened derivatives of those, because the brights don't reach
readable contrast as small text on white (yellow, green) or on charcoal (red).
"""
from __future__ import annotations

import json
import traceback
from collections.abc import Callable
from pathlib import Path

import customtkinter as ctk

LIGHT = "light"
DARK = "dark"

NAVY = "#285087"
SKY = "#96d2ff"

# token -> (light, dark)
TOKENS: dict[str, tuple[str, str]] = {
    # Surfaces
    "bg":        ("#ffffff", "#262626"),   # window
    "surface":   ("#f0f0f0", "#303030"),   # tab content / cards
    "surface2":  ("#ffffff", "#3a3a3a"),   # frames nested on a surface
    "line":      ("#d2d3d3", "#4a4a4a"),
    # Text
    "fg":        ("#000000", "#e6e6e6"),
    "muted":     ("#4a4a4a", "#a3a5a4"),
    # Buttons
    "primary":        (NAVY, NAVY),
    "primary_hover":  ("#1d3d68", "#3a68a8"),
    # Outlined secondary button (fill + border + body-colored text), as in the HCI mockup
    "secondary":       ("#ffffff", "#3a3a3a"),
    "secondary_hover": ("#e4e5e5", "#454545"),
    # Segmented-control track / unselected segments
    "seg":       ("#707271", "#4d4d4d"),
    "seg_hover": ("#505050", "#5e5e5e"),
    # Tab bar
    "chrome":       (NAVY, "#1d1d1d"),
    "chrome_fg":    ("#ffffff", "#e6e6e6"),
    "chrome_hover": ("#1d3d68", "#333333"),
    "tab_accent":   ("#ecd451", "#f2ad49"),
    # Help "?" badge
    "badge":    ("#d2d3d3", "#4d4d4d"),
    "danger":        ("#890000", "#c93f25"),
    "danger_hover":  ("#6b0000", "#e65032"),
    "accent":        (NAVY, SKY),          # links / progress / highlights
    # Status text
    "ok":    ("#00775f", "#26c4a1"),
    "warn":  ("#8a6500", "#ecd451"),
    "err":   ("#890000", "#e65032"),
    # Tables
    "row_odd":  ("#ffffff", "#262626"),
    "row_even": ("#f6f6f6", "#2d2d2d"),
    "select_bg": (SKY, SKY),
    "select_fg": ("#000000", "#000000"),
    "head_bg":   (NAVY, NAVY),
    "head_fg":   ("#ffffff", "#ffffff"),
    "head_hover": ("#1d3d68", "#3a68a8"),
    # Plots
    "plot_bg":      ("#ffffff", "#262626"),
    "plot_axis":    ("#a3a5a4", "#555555"),
    "plot_fg":      ("#000000", "#dcdcdc"),
    "plot_strong":  ("#000000", "#ffffff"),
    "plot_legend_bg":   ("#f0f0f0", "#3a3a3a"),
    "plot_legend_edge": ("#d2d3d3", "#555555"),
    "neutral":      ("#707271", "#9a9c9b"),   # unscored needles / backbone
    "ptm_marker":   (NAVY, SKY),
    # Tooltip
    "tip_bg": ("#505050", "#3a3a3a"),
    "tip_fg": ("#ffffff", "#e6e6e6"),
}

# Static data colors (identical in both modes; drawn with a dark outline)
PP_RED = "#e65032"
PP_YELLOW = "#ecd451"
PP_GREEN = "#009f7e"

_mode: str = DARK
_listeners: list[Callable[[str], None]] = []
PREFS_FILE: Path | None = None


def mode() -> str:
    return _mode


def pair(token: str) -> tuple[str, str]:
    """(light, dark) tuple for CustomTkinter ``*_color`` arguments."""
    return TOKENS[token]


def c(token: str) -> str:
    """Hex color for *token* in the current mode (plain Tk / matplotlib)."""
    light, dark = TOKENS[token]
    return light if _mode == LIGHT else dark


def resolve(color):
    """Collapse a ``(light, dark)`` tuple to the current mode; pass strings through."""
    if isinstance(color, (tuple, list)):
        return color[0] if _mode == LIGHT else color[1]
    return color


def on_change(callback: Callable[[str], None]) -> None:
    """Call *callback(mode)* whenever :func:`set_mode` changes the mode."""
    _listeners.append(callback)


def init(prefs_file: Path) -> None:
    """Load the saved mode (default dark), install the HCI CustomTkinter
    theme, and set the appearance mode. Call once at startup, before any
    widgets are created.
    """
    global PREFS_FILE, _mode
    PREFS_FILE = prefs_file
    try:
        saved = json.loads(prefs_file.read_text()).get("appearance_mode")
        if saved in (LIGHT, DARK):
            _mode = saved
    except Exception:
        pass
    ctk.set_default_color_theme("blue")
    _install_ctk_theme()
    _install_style_patches()
    ctk.set_appearance_mode(_mode)


def set_mode(new_mode: str, persist: bool = True) -> None:
    global _mode
    if new_mode not in (LIGHT, DARK) or new_mode == _mode:
        return
    _mode = new_mode
    ctk.set_appearance_mode(new_mode)
    if persist:
        _save_prefs()
    for callback in list(_listeners):
        try:
            callback(new_mode)
        except Exception:
            traceback.print_exc()


def toggle() -> str:
    set_mode(LIGHT if _mode == DARK else DARK)
    return _mode


def _save_prefs() -> None:
    if PREFS_FILE is None:
        return
    try:
        data: dict = {}
        if PREFS_FILE.exists():
            data = json.loads(PREFS_FILE.read_text())
        data["appearance_mode"] = _mode
        PREFS_FILE.parent.mkdir(parents=True, exist_ok=True)
        PREFS_FILE.write_text(json.dumps(data, indent=2))
    except Exception:
        pass


def _install_ctk_theme() -> None:
    """Override the built-in 'blue' CustomTkinter theme with HCI colors."""
    p = pair
    white = ("#ffffff", "#ffffff")
    overrides: dict[str, dict] = {
        "CTk": {"fg_color": p("bg")},
        "CTkToplevel": {"fg_color": p("bg")},
        "CTkFrame": {
            "fg_color": p("surface"), "top_fg_color": p("surface2"),
            "border_color": p("line"),
        },
        "CTkButton": {
            "fg_color": p("primary"), "hover_color": p("primary_hover"),
            "border_color": p("line"), "text_color": white,
            "text_color_disabled": ("#d2d3d3", "#8a8c8b"),
        },
        "CTkLabel": {"text_color": p("fg")},
        "CTkEntry": {
            "fg_color": p("bg"), "border_color": p("line"), "text_color": p("fg"),
            "placeholder_text_color": p("muted"),
        },
        "CTkCheckBox": {
            "fg_color": p("accent"), "hover_color": p("accent"),
            "border_color": p("muted"), "checkmark_color": ("#ffffff", "#1d1d1d"),
            "text_color": p("fg"), "text_color_disabled": ("#a3a5a4", "#707271"),
        },
        "CTkSwitch": {
            "fg_color": ("#c4c5c5", "#4a4a4a"), "progress_color": p("accent"),
            "button_color": ("#505050", "#d5d5d5"), "button_hover_color": ("#303030", "#ffffff"),
            "text_color": p("fg"), "text_color_disabled": ("#a3a5a4", "#707271"),
        },
        "CTkRadioButton": {
            "fg_color": p("accent"), "hover_color": p("primary_hover"),
            "border_color": p("muted"), "text_color": p("fg"),
            "text_color_disabled": ("#a3a5a4", "#707271"),
        },
        "CTkProgressBar": {
            "fg_color": p("line"), "progress_color": p("accent"),
        },
        "CTkSlider": {
            "fg_color": p("line"), "progress_color": p("muted"),
            "button_color": p("accent"), "button_hover_color": p("primary_hover"),
        },
        "CTkOptionMenu": {
            "fg_color": p("primary"), "button_color": p("primary_hover"),
            "button_hover_color": ("#14294a", "#4a78b8"), "text_color": white,
            "text_color_disabled": ("#d2d3d3", "#8a8c8b"),
        },
        "CTkComboBox": {
            "fg_color": p("bg"), "border_color": p("line"),
            "button_color": p("primary"), "button_hover_color": p("primary_hover"),
            "text_color": p("fg"), "text_color_disabled": ("#a3a5a4", "#707271"),
        },
        "CTkScrollbar": {
            "button_color": ("#a3a5a4", "#555555"), "button_hover_color": ("#707271", "#707271"),
        },
        "CTkSegmentedButton": {
            "fg_color": p("seg"), "selected_color": p("primary"),
            "selected_hover_color": p("primary_hover"),
            "unselected_color": p("seg"), "unselected_hover_color": p("seg_hover"),
            "text_color": white, "text_color_disabled": ("#d2d3d3", "#8a8c8b"),
        },
        "CTkTextbox": {
            "fg_color": p("bg"), "border_color": p("line"), "text_color": p("fg"),
            "scrollbar_button_color": ("#a3a5a4", "#555555"),
            "scrollbar_button_hover_color": ("#707271", "#707271"),
        },
        "CTkScrollableFrame": {"label_fg_color": p("surface2")},
        "DropdownMenu": {
            "fg_color": p("surface2"), "hover_color": ("#d6ebfb", NAVY),
            "text_color": p("fg"),
        },
    }
    for widget, values in overrides.items():
        ctk.ThemeManager.theme.setdefault(widget, {}).update(values)


def _style_button_kwargs(kw: dict) -> None:
    """Give buttons their HCI look from fg_color alone, so the ~80 call sites
    that pass ``fg_color=_BTN2`` / ``_PRIMARY`` / ``_DANGER`` need no changes:
    secondary buttons are outlined with body-colored text; filled ones get
    white text and no border.
    """
    fg = kw.get("fg_color")
    if fg is None:
        return
    if fg == "transparent":
        kw.setdefault("text_color", pair("fg"))
    elif fg == TOKENS["secondary"]:
        kw.setdefault("text_color", pair("fg"))
        kw.setdefault("border_width", 1)
        kw.setdefault("border_color", pair("line"))
    elif fg in (TOKENS["primary"], TOKENS["danger"]):
        kw.setdefault("text_color", ("#ffffff", "#ffffff"))
        kw.setdefault("border_width", 0)


def _install_style_patches() -> None:
    """Button styling (above) plus a 1px border on every filled CTkFrame, so
    panels read as outlined cards. Transparent frames stay borderless.
    """
    orig_btn_init = ctk.CTkButton.__init__
    orig_btn_configure = ctk.CTkButton.configure
    orig_frame_init = ctk.CTkFrame.__init__

    def btn_init(self, *args, **kwargs):
        _style_button_kwargs(kwargs)
        orig_btn_init(self, *args, **kwargs)

    def btn_configure(self, require_redraw=False, **kwargs):
        _style_button_kwargs(kwargs)
        orig_btn_configure(self, require_redraw=require_redraw, **kwargs)

    def frame_init(self, *args, **kwargs):
        if kwargs.get("fg_color", None) != "transparent":
            kwargs.setdefault("border_width", 1)
        orig_frame_init(self, *args, **kwargs)

    ctk.CTkButton.__init__ = btn_init
    ctk.CTkButton.configure = btn_configure
    ctk.CTkFrame.__init__ = frame_init
