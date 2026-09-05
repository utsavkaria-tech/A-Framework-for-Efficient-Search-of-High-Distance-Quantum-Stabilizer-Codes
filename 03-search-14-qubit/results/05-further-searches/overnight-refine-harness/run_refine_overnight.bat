@echo off
setlocal
rem ---------------------------------------------------------------------------
rem  Overnight targeted refinement of the 5-residual incumbent (Stage-1 orbit 0006).
rem  Anchored iterative broadening: every sweep restarts from the verified incumbent
rem  and only the 3-opt breadth M grows, so the search never wanders into a worse
rem  basin.  Writes only inside this directory.
rem
rem  HOURS is the wall-clock budget.  8 hours gets through roughly M = 16 and 32.
rem ---------------------------------------------------------------------------
set EXE=%~dp0stab14_refine.exe
set LOG=%~dp0refine_overnight.log
set HOURS=9
set /a SECS=%HOURS%*3600

echo Started %DATE% %TIME%  (budget %HOURS% h) > "%LOG%"
"%EXE%" --time-limit %SECS% --top-m 16 --verbose --out "%~dp0refine_out" >> "%LOG%" 2>&1
echo Finished %DATE% %TIME% >> "%LOG%"

if exist "%~dp0refine_out\FULL_COVERAGE.txt" (
  echo *** FULL COVERAGE FOUND *** & type "%~dp0refine_out\FULL_COVERAGE.txt"
) else (
  echo Done. Final report:
  type "%~dp0refine_out\refinement_report.txt"
)
