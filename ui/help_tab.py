"""Help / Documentation tab: renders docs/help.md in-app."""
from __future__ import annotations

import customtkinter as ctk

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

        html_body = markdown.markdown(md_text, extensions=["tables", "fenced_code"])
        css = """
        <style>
            body { font-family: Segoe UI, Arial, sans-serif; padding: 24px;
                   background: #2b2b2b; color: #dcdcdc; line-height: 1.6; }
            h1 { color: #3a86ff; border-bottom: 2px solid #3a86ff; padding-bottom: 6px; }
            h2 { color: #6cb4ee; margin-top: 28px; }
            h3 { color: #a0cfff; margin-top: 20px; }
            table { border-collapse: collapse; margin: 12px 0; width: 100%; }
            th, td { border: 1px solid #555; padding: 6px 10px; text-align: left; }
            th { background: #3a3a3a; color: #dcdcdc; }
            tr:nth-child(even) { background: #333; }
            code { background: #3a3a3a; padding: 2px 5px; border-radius: 3px; }
            hr { border: 1px solid #555; margin: 20px 0; }
            strong { color: #f0f0f0; }
            ul, ol { padding-left: 24px; }
            li { margin: 4px 0; }
        </style>
        """
        full_html = f"<html><head>{css}</head><body>{html_body}</body></html>"

        try:
            import webbrowser
            from tkinterweb import HtmlFrame
            frame = HtmlFrame(
                tab, messages_enabled=False,
                # Without this, clicking a link (e.g. the COSMIC/PTMD sources)
                # navigates this embedded frame itself to that external site,
                # replacing the help content with no way back short of
                # restarting the app. Providing on_link_click replaces the
                # default in-place navigation entirely, so opening the URL in
                # the system browser instead means this frame never navigates
                # away from the docs at all.
                on_link_click=lambda url: webbrowser.open(url),
            )
            frame.load_html(full_html)
            frame.grid(row=0, column=0, sticky="nsew")
        except ImportError:
            fallback = ctk.CTkTextbox(tab, wrap="word", font=ctk.CTkFont(size=13))
            fallback.insert("1.0", md_text)
            fallback.configure(state="disabled")
            fallback.grid(row=0, column=0, sticky="nsew")
