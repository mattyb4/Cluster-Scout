"""Help / Documentation tab: renders docs/help.md in-app."""
from __future__ import annotations

import customtkinter as ctk

from ui import theme
from ui.common import PROJECT_ROOT


class HelpTabMixin:
    def _build_help_tab(self, tab):
        """Load docs/help.md, convert to HTML, and display in the Help tab."""
        import markdown

        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(0, weight=1)

        help_md = PROJECT_ROOT / "docs" / "help.md"
        try:
            md_text = help_md.read_text(encoding="utf-8")
        except FileNotFoundError:
            md_text = "# Help\n\nDocumentation file not found at `docs/help.md`."

        self._help_html_body = markdown.markdown(md_text, extensions=["tables", "fenced_code"])
        self._help_frame = None
        full_html = self._help_full_html()

        try:
            import webbrowser

            from tkinterweb import HtmlFrame
            frame = HtmlFrame(
                tab, messages_enabled=False,
                # Without this, clicking a link navigates this embedded frame itself
                # away from the docs; open in the system browser instead
                on_link_click=lambda url: webbrowser.open(url),
            )
            frame.load_html(full_html)
            self._help_frame = frame
            frame.grid(row=0, column=0, sticky="nsew")
        except ImportError:
            fallback = ctk.CTkTextbox(tab, wrap="word", font=ctk.CTkFont(size=13))
            fallback.insert("1.0", md_text)
            fallback.configure(state="disabled")
            fallback.grid(row=0, column=0, sticky="nsew")

    def _help_full_html(self) -> str:
        """The help page's HTML with CSS colors for the current light/dark mode."""
        t = theme.c
        css = f"""
        <style>
            body {{ font-family: Segoe UI, Arial, sans-serif; padding: 24px;
                   background: {t("bg")}; color: {t("fg")}; line-height: 1.6; }}
            h1 {{ color: {t("accent")}; border-bottom: 2px solid {theme.NAVY}; padding-bottom: 6px; }}
            h2 {{ color: {t("accent")}; margin-top: 28px; }}
            h3 {{ color: {t("fg")}; margin-top: 20px; }}
            a {{ color: {t("accent")}; }}
            table {{ border-collapse: collapse; margin: 12px 0; width: 100%; }}
            th, td {{ border: 1px solid {t("line")}; padding: 6px 10px; text-align: left; }}
            th {{ background: {t("head_bg")}; color: {t("head_fg")}; }}
            tr:nth-child(even) {{ background: {t("row_even")}; }}
            code {{ background: {t("surface")}; padding: 2px 5px; border-radius: 3px; }}
            hr {{ border: 1px solid {t("line")}; margin: 20px 0; }}
            strong {{ color: {t("plot_strong")}; }}
            ul, ol {{ padding-left: 24px; }}
            li {{ margin: 4px 0; }}
        </style>
        """
        return f"<html><head>{css}</head><body>{self._help_html_body}</body></html>"

    def _retheme_help(self) -> None:
        if self._help_frame is not None:
            self._help_frame.load_html(self._help_full_html())
