#!/usr/bin/env python3
"""Require Faust to accept valid controls and reject unsafe table indices."""

import argparse
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--faust", required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    failures = 0
    print(subprocess.check_output([args.faust, "--version"], text=True), flush=True)
    with tempfile.TemporaryDirectory(prefix="interval-faust-rejections-") as temp:
        for backend in ["cpp", "ocpp"]:
            for precision in ["-single", "-double"]:
                for expectation in ["accept", "reject"]:
                    for program in sorted((root / expectation).glob("*.dsp")):
                        result = subprocess.run([
                            args.faust, "-lang", backend, precision,
                            str(program), "-o", str(Path(temp) / "output.cpp"),
                        ], capture_output=True, text=True)
                        label = f"{program.stem}/{backend}/{precision[1:]}"
                        # A crash is not a valid diagnostic rejection. Faust's
                        # ordinary compilation errors return status 1.
                        passed = (result.returncode == 0 if expectation == "accept"
                                  else result.returncode == 1 and bool(result.stderr.strip()))
                        if passed:
                            print(f"PASS {label}: {expectation}", flush=True)
                        else:
                            failures += 1
                            print(f"FAIL {label}: expected {expectation}, "
                                  f"exit status {result.returncode}", flush=True)
                            if result.stderr:
                                print(result.stderr.rstrip(), flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
