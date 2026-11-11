"""Tests for the startup window sizing in ui/common.py (fit_window_size)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

pytest.importorskip("customtkinter")
from ui.common import fit_window_size  # noqa: E402


class TestFitWindowSize:
    def test_a_big_screen_gets_the_default_size(self):
        width, height, min_w, min_h = fit_window_size(2560, 1440, 1.0, 830)
        assert (width, height) == (1100, 820), f"a roomy screen should open at the default size, got {width}x{height}"
        assert min_w == 900 and min_h == 560, f"the usual minimums apply when the tabs fit, got {min_w}x{min_h}"

    def test_a_tab_bar_wider_than_the_usual_minimum_raises_it(self):
        _w, _h, min_w, _min_h = fit_window_size(2560, 1440, 1.0, 1030)
        assert min_w == 1030, f"the minimum width must cover the tab bar so no tab is clipped, got {min_w}"

    def test_the_window_is_never_narrower_than_the_tab_bar(self):
        width, *_ = fit_window_size(2560, 1440, 1.0, 1200)
        assert width == 1200, f"a tab bar wider than the default width should widen the window, got {width}"

    @pytest.mark.parametrize("scale", [1.0, 1.25, 1.5])
    def test_a_1080p_screen_is_never_overfilled(self, scale):
        width, height, min_w, min_h = fit_window_size(1920, 1080, scale, 1030)
        assert width * scale <= 1920 - 40 and height * scale <= 1080 - 120, (
            f"at {scale}x the window ({width}x{height}) must fit on a 1920x1080 screen with room for the taskbar"
        )
        assert min_w <= width and min_h <= height, "the minimum size can't exceed the window's own size"

    def test_a_1080p_screen_at_150_percent_shortens_the_window(self):
        _w, height, _min_w, min_h = fit_window_size(1920, 1080, 1.5, 830)
        assert 560 <= height < 820, f"only 640 units of height are free at 150%, so expected a shorter window, got {height}"
        assert min_h <= height

    def test_a_screen_too_narrow_for_the_tabs_wins_over_the_tab_bar(self):
        width, _h, min_w, _min_h = fit_window_size(1280, 720, 1.0, 1500)
        assert width <= 1240 and min_w <= width, (
            f"on a screen narrower than the tab bar the window must still fit on it, got width {width}, minimum {min_w}"
        )
