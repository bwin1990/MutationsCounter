@echo off
setlocal enabledelayedexpansion

echo ===================================================
echo         Mutation Analysis Tool - Batch Processing
echo ===================================================
echo.

REM Check if Python is installed
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Python not found. Please make sure Python is installed and added to PATH.
    goto :end
)

REM Check if analysis script exists
if not exist "%~dp0analyze_mutations.py" (
    echo ERROR: analyze_mutations.py file not found in the current directory.
    goto :end
)

echo Please select operation mode:
echo 1. Select a folder (process all .mpileup.cns.filter.xls files in the folder)
echo 2. Process all .mpileup.cns.filter.xls files in the current directory
echo 3. Specify a single file to process
echo.

set /p choice="Enter your choice (1/2/3): "

if "%choice%"=="1" (
    echo Starting folder selection mode...
    python "%~dp0analyze_mutations.py"
    goto :end
)

if "%choice%"=="2" (
    echo Processing all .mpileup.cns.filter.xls files in the current directory...
    
    set file_count=0
    for %%f in (*.mpileup.cns.filter.xls) do (
        set /a file_count+=1
    )
    
    if !file_count! equ 0 (
        echo ERROR: No .mpileup.cns.filter.xls files found in the current directory.
        goto :end
    )
    
    echo Found !file_count! files to process...
    
    set current=0
    for %%f in (*.mpileup.cns.filter.xls) do (
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
    echo Please drag and drop the .mpileup.cns.filter.xls file into this window, then press Enter:
    set /p file_path=
    
    if not exist "!file_path!" (
        echo ERROR: File does not exist.
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
