# sitecustomize.py — enables coverage in subprocesses when COVERAGE_PROCESS_START is set
import os
if os.getenv("COVERAGE_PROCESS_START"):
    try:
        # coverage.process_startup will read COVERAGE_PROCESS_START to configure coverage
        import coverage
        coverage.process_startup()
    except Exception:
        # avoid breaking test runs if coverage isn't available
        pass
