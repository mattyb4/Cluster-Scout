@echo off
setlocal enabledelayedexpansion
title Cluster-Scout Launcher
cd /d "%~dp0"

echo ============================================
echo   Cluster-Scout - Beta Launcher
echo ============================================
echo.

where uv >nul 2>nul
if !errorlevel! neq 0 (
    echo [1/3] uv not found - installing it now, please wait...
    powershell -ExecutionPolicy ByPass -Command "irm https://astral.sh/uv/install.ps1 | iex"
    set "PATH=%USERPROFILE%\.local\bin;%PATH%"
    where uv >nul 2>nul
    if !errorlevel! neq 0 (
        echo.
        echo ERROR: Could not find or install uv automatically.
        echo Please install it manually from https://astral.sh/uv and re-run this script.
        pause
        exit /b 1
    )
    echo   uv installed successfully.
) else (
    echo [1/3] uv found.
)

echo.
echo [2/3] Installing/updating dependencies - first run may take a few minutes...
uv sync
if !errorlevel! neq 0 (
    echo.
    echo ERROR: Dependency installation failed. See the output above for details.
    pause
    exit /b 1
)

echo.
echo [3/3] Launching Cluster-Scout...
echo.
uv run app.py
if !errorlevel! neq 0 (
    echo.
    echo Cluster-Scout closed with an error ^(see above^).
    echo If you need help, send a screenshot of this window to Matt.
    pause
)
