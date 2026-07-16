# macOS beta-testing handoff

Paste this to a fresh Claude Code session on the Mac where Cluster-Scout is
being tested. That session has no memory of the investigation that produced
this note — everything it needs to pick up cleanly should be here or in the
git history it points to.

## Context

Cluster-Scout is a CustomTkinter desktop app being distributed as a beta to
research collaborators via a GitHub Release (tag `v0.1.1`), downloaded as a
source zip and launched via `run.command` (which runs `uv sync` then
`uv run app.py`). It was built and has only ever been run on Windows until
this beta round — this is the first real macOS testing.

## Symptom reported

Two rounds of macOS testing so far, both on real hardware (not something the
Windows dev environment that produced these fixes could reproduce or verify):

**Round 1** — crash tracebacks on launch/interaction, including:
- `_tkinter.TclError: invalid command name` from CustomTkinter's
  `CTkTextbox._check_if_scrollbars_needed`
- `_tkinter.TclError: can't invoke "winfo" command: application has been
  destroyed` from CustomTkinter's `ScalingTracker.check_dpi_scaling`
- General "major performance issues"

**Round 2** (after the round-1 fixes below) — crashes stopped, but
performance issues persisted: "clicks to change tabs or check boxes respond
immediately, but they often do not respond" — intermittent, not a uniform
slowdown, which reads like main-thread event-loop contention rather than one
slow operation.

## What's already been fixed (verify these are actually in place first)

All in `app.py`, at module level (four `_patch_customtkinter_*` functions
called right after the `ctk.set_appearance_mode("dark")` line). Run
`git log --oneline -6` and `git show <hash>` on each of these to see the
actual diffs and reasoning in the commit messages, rather than trusting this
summary blindly:

1. **`iconbitmap()` → `iconphoto()` platform branch** — `.ico` isn't a valid
   format for `iconbitmap()` outside Windows; raises `TclError` there.
2. **Mousewheel delta scaling** (`ui/common.py`, `isolate_textbox_scroll`) —
   was dividing `event.delta` by 120 unconditionally (Windows' per-notch
   convention); macOS reports much smaller raw deltas, so scrolling in the
   log textbox was silently truncating to 0 units.
3. **`CTkTextbox._check_if_scrollbars_needed` wrapped in try/except** — the
   library's own `winfo_exists()` guard only protects the *next* scheduled
   call, not the one currently executing, so a textbox destroyed between
   scheduling and firing (e.g. a Pipeline-tab mode switch tearing down and
   rebuilding widgets) could throw before reaching that guard.
4. **`ScalingTracker.check_dpi_scaling` disabled entirely outside Windows** —
   this CustomTkinter internal polls every open window every 100ms, forever,
   to detect OS-level per-monitor DPI changes. Real feature on Windows;
   `get_window_dpi_scaling()` is hardcoded to always return `1` on macOS, so
   the loop can never detect anything there — pure wasted overhead for the
   app's entire session, confirmed by reading the actual library source at
   `.venv/lib/python*/site-packages/customtkinter/windows/widgets/scaling/scaling_tracker.py`.
5. **`AppearanceModeTracker.update` disabled entirely, all platforms** —
   near-identical loop to #4 but firing every **30ms** (3x more often).
   Detects live OS light/dark-mode toggles via `darkdetect`. This app calls
   `ctk.set_appearance_mode("dark")` once at import and never switches at
   runtime, so the loop's entire purpose never applies — verified via a real
   screenshot that disabling it has zero visual effect.

Fix #5 was the working hypothesis for round 2's "clicks often don't
respond" — a loop firing 30+ times/second regardless of what's on screen is
a plausible source of intermittent click-eating if macOS's Tk/Cocoa event
dispatch has more per-event overhead than Windows' (a commonly reported
characteristic of Tk on macOS, not something specific to this app). **This
has not been confirmed on real hardware** — that's what this testing round
is for.

## What's still unverified / not ruled out

- Whether #4 + #5 were the full explanation for the responsiveness issue, or
  just part of it.
- Whether `uv`-managed Python on macOS bundles a fully-working Tcl/Tk. This
  has been a real, documented gap in `python-build-standalone` (what `uv`
  uses) in the past. If the app is himself sluggish even doing nothing, this
  is worth checking — compare against a system/python.org Python install if
  the symptom persists.
- Whether matplotlib's `TkAgg` backend or `tkinterweb` (used for the Help tab)
  behave differently on macOS's Tk — not profiled either way.
- The app's own `_poll_queue` (`ui/pipeline_runner.py`) runs every 100ms
  forever by design (it's how background pipeline threads hand results to
  the UI) — audited and it's a cheap non-blocking `queue.get_nowait()` check
  with no Tk calls when idle, so it wasn't touched, but it hasn't been
  profiled under real load either.

## Verification steps

1. `git pull`, confirm `git log --oneline -1` shows `a0dee40` ("Eliminate
   AppearanceModeTracker's 30ms polling loop") or later, or check out tag
   `v0.1.1`.
2. Launch from a terminal directly — `uv run app.py` — rather than
   double-clicking `run.command`, so any traceback is visible live instead of
   needing a screenshot after the fact.
3. Click rapidly through all 5 tabs (Pipeline / Results / Visualization /
   Analysis Tools / Help) and toggle checkboxes (PolyPhen filter, mode radio
   buttons) for 30-60 seconds. Note whether clicks still feel dropped/delayed,
   and whether it's better, the same, or different from before.
4. Watch the terminal the whole time for any uncaught traceback — especially
   around Pipeline-tab mode switches (destroys/rebuilds widgets) and app
   close (`Cmd+Q` or the window close button).
5. **If it's still slow**, don't guess again — get real data:
   - `uv add --dev py-spy` (or `pip install py-spy` in the venv), then
     `py-spy record -o profile.svg -- uv run app.py`, reproduce the lag, quit
     the app, and open `profile.svg` — a flamegraph will show exactly where
     time is actually going during the interaction, instead of continuing to
     reason from source code alone the way this investigation has so far.
   - Alternatively `py-spy dump --pid <pid>` while it's mid-lag, a few times
     in a row, to sample what the main thread is doing at that moment.

## Reporting back

If it's fixed: good, nothing further needed. If not: bring back whichever of
the following you can get —
- Exact repro steps and whether the lag is truly continuous or intermittent.
- Any new tracebacks (full text, not paraphrased).
- A `py-spy` flamegraph or a few dump samples taken during the lag.

Real profiling data at this point would be much higher-value than another
round of source-reading and hypothesis, since two rounds of "should fix it"
based on static analysis have each turned out to be partial.
