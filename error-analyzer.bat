@echo off
setlocal enabledelayedexpansion

set "ENV_NAME=bio-work"

echo ===================================================
echo         Mutation Analysis Tool - Batch Processing
echo ===================================================
echo.

REM Activate conda environment
set "CONDA_BAT="
for /f "usebackq delims=" %%i in (`where conda.bat 2^>nul`) do (
    set "CONDA_BAT=%%i"
    goto :conda_found
)

:conda_found
if not defined CONDA_BAT (
    echo ERROR: conda.bat not found. Please ensure Conda is installed and added to PATH.
    goto :end
)

call "%CONDA_BAT%" activate %ENV_NAME% >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Failed to activate Conda environment "%ENV_NAME%".
    echo Please verify the environment exists: conda env list
    goto :end
)

REM Check if Python is available in the activated environment
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Python not found in Conda environment "%ENV_NAME%".
    goto :end
)

REM Check if analysis script exists
if not exist "%~dp0analyze_mutations.py" (
    echo ERROR: analyze_mutations.py file not found in the current directory.
    goto :end
)

echo Please select operation mode:
echo 1. Select a folder (process all mpileup CNS *.xls files in the folder)
echo 2. Process all mpileup CNS *.xls files in the current directory
echo 3. Specify a single file to process
echo.

set /p choice="Enter your choice (1/2/3): "

if "%choice%"=="1" (
    echo Opening folder picker...
    for /f "usebackq delims=" %%i in (`powershell -NoProfile -Command "Add-Type -AssemblyName System.Windows.Forms; $f = New-Object System.Windows.Forms.FolderBrowserDialog; $f.Description = 'Select folder containing mpileup CNS .xls files'; if ($f.ShowDialog() -eq 'OK') { $f.SelectedPath }"`) do set "folder_path=%%i"

    if not defined folder_path (
        echo No folder selected. Operation canceled.
        goto :end
    )

    python "%~dp0analyze_mutations.py" "!folder_path!"
    goto :end
)

if "%choice%"=="2" (
    echo Processing all mpileup CNS *.xls files in the current directory...
    
    set file_count=0
    for %%f in (*.mpileup*.cns*.xls) do (
        set /a file_count+=1
    )
    
    if !file_count! equ 0 (
        echo ERROR: No mpileup CNS *.xls files found in the current directory.
        goto :end
    )
    
    echo Found !file_count! files to process...
    
    set current=0
    for %%f in (*.mpileup*.cns*.xls) do (
        set /a current+=1
        echo.
        echo [!current!/!file_count!] Processing: %%f
        python "%~dp0analyze_mutations_single.py" "%%~dpnxf"
    )
    
    echo.
    echo Batch processing complete! Processed !file_count! files
    goto :end
)

if "%choice%"=="3" (
    echo Opening file picker...
    for /f "usebackq delims=" %%i in (`powershell -NoProfile -Command "Add-Type -AssemblyName System.Windows.Forms; $f = New-Object System.Windows.Forms.OpenFileDialog; $f.Filter = 'mpileup CNS files (*.mpileup*.cns*.xls)|*.mpileup*.cns*.xls|All files (*.*)|*.*'; $f.Multiselect = $false; if ($f.ShowDialog() -eq 'OK') { $f.FileName }"`) do set "file_path=%%i"
    
    if not defined file_path (
        echo No file selected. Operation canceled.
        goto :end
    )
    
    echo Processing file: !file_path!
    python "%~dp0analyze_mutations_single.py" "!file_path!"
    goto :end
)

echo Invalid option, please run the script again and select 1, 2, or 3.

:end
echo.
echo Press any key to exit...
pause >nul
