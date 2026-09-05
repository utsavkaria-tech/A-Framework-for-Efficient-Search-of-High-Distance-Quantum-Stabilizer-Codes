@echo off
setlocal
rem ---------------------------------------------------------------------------
rem  Overnight Stage-2 attack on Stage-1 orbit 0006.
rem  Each seed runs in its own sub-directory, so no result is ever overwritten.
rem  Stops early if a full-coverage certificate (M_8 = 11186) is produced.
rem  Adjust BUDGET (seconds per seed) and the seed list to taste.
rem ---------------------------------------------------------------------------
set EXE=%~dp0stab14_2stage.exe
set LOG=%~dp0overnight_orbit0006.log
set BUDGET=1800

echo Started %DATE% %TIME%  (budget %BUDGET% s per seed) > "%LOG%"
for %%S in (1 2 3 5 8 13 21 34 55 89 144 233 377 610 987 1597) do (
  echo. >> "%LOG%"
  echo ================ seed %%S ================ >> "%LOG%"
  if not exist "%~dp0seed_%%S" mkdir "%~dp0seed_%%S"
  pushd "%~dp0seed_%%S"
  "%EXE%" --stage1-index 6 --stage2-time-limit %BUDGET% --selftest 0 --seed %%S >> "%LOG%" 2>&1
  popd
  if exist "%~dp0seed_%%S\full_coverage_solution.txt" (
    echo. >> "%LOG%"
    echo *** FULL COVERAGE FOUND at seed %%S *** >> "%LOG%"
    goto :done
  )
)
:done
echo. >> "%LOG%"
echo Finished %DATE% %TIME% >> "%LOG%"
echo Done. Summary:
findstr /C:"M_8 =" "%LOG%"
