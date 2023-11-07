@echo off
setlocal enabledelayedexpansion

:: Check for Python 3.9 or higher
echo Checking for Python version 3.9 or higher...
python --version 2>&1 | findstr /R "Python 3\.[9-9][0-9]*\..*" > nul
if errorlevel 1 (
    echo Python 3.9 or higher is not found.
    goto :end
) else (
    echo Python 3.9 or higher is installed.
)

:: Check for Git
echo Checking for Git...
git --version > nul 2>&1
if errorlevel 1 (
    echo Git is not found.
    goto :end
) else (
    echo Git is installed.
)

:: Clone the GitHub repo
echo Cloning the GitHub repo...
git clone https://github.com/Fascinat0r/corrosion-tool2.git

:: Install Jupyter Notebook
echo Installing Jupyter Notebook...
python -m pip install jupyter notebook thermopack==2.1.0

:: Launch the Jupyter Notebook
echo Launching Jupyter Notebook...
python -m jupyter notebook ./corrosion-tool2/jupiter/Main.ipynb

:end
pause
