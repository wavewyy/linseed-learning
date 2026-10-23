@echo off
setlocal
pushd "%~dp0"
"%~dp0R\bin\Rscript.exe" "%~dp0study\01_quickstart.R"
if errorlevel 1 (
  echo ERROR: Please copy the error message for troubleshooting.
) else (
  echo SUCCESS: Open results\quickstart\proportions.png
)
popd
pause
