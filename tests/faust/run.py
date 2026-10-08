#!/usr/bin/env python3
"""Run the rain regression in both backends/precisions with an optional reference."""

import argparse
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


def same(left, right):
    return len(left) == len(right) and all(
        math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-6) for a, b in zip(left, right)
    )


def run_case(args, directory, binary, program, backend, precision, name):
    source = directory / f"{name}.cpp"
    executable = directory / name
    subprocess.run([
        binary, "-lang", backend, precision, "-I", str(args.faust_root / "libraries"),
        "-a", str(Path(__file__).with_name("runtime-architecture.cpp")),
        str(program), "-o", str(source),
    ], check=True)
    subprocess.run([
        args.cxx, "-std=c++17", "-O2", "-fwrapv", "-ffp-contract=off",
        "-I", str(args.faust_root / "architecture"), *shlex.split(args.cxxflags),
        str(source), "-o", str(executable),
    ], check=True)
    output = subprocess.run([str(executable)], check=True, capture_output=True, text=True)
    rows = [[float(value) for value in line.split()] for line in output.stdout.splitlines()]
    if len(rows) != 4096 or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise RuntimeError(f"{name}: missing or malformed output")
    if not all(math.isfinite(value) for row in rows for value in row):
        raise RuntimeError(f"{name}: non-finite output")
    if not all(any(row[channel] != 0 for row in rows) for channel in range(len(rows[0]))):
        raise RuntimeError(f"{name}: a noise channel was replaced by silence")
    return [value for row in rows for value in row]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--faust", required=True)
    parser.add_argument("--reference")
    parser.add_argument("--faust-root", required=True, type=Path)
    parser.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    parser.add_argument("--cxxflags", default="")
    args = parser.parse_args()
    # A comparator must reject a changed sample before its passing verdict matters.
    assert not same([0.0, 0.5], [0.0, 1.5])
    for binary in filter(None, [args.faust, args.reference]):
        print(subprocess.check_output([binary, "--version"], text=True), flush=True)
    programs = [Path(__file__).with_name("rain-second-noise.dsp"),
                args.faust_root / "examples/gameaudio/rain.dsp"]
    with tempfile.TemporaryDirectory(prefix="interval-faust-rain-") as temp:
        directory = Path(temp)
        for program in programs:
            for backend in ["cpp", "ocpp"]:
                for precision in ["-single", "-double"]:
                    name = f"{program.stem}-{backend}-{precision[1:]}"
                    result = run_case(args, directory, args.faust, program, backend, precision, name)
                    if args.reference:
                        reference = run_case(args, directory, args.reference, program,
                                             backend, precision, name + "-reference")
                        if not same(result, reference):
                            raise RuntimeError(f"{name}: runtime output differs from reference")
                    print(f"PASS {name}: 4096 frames, {len(result) // 4096} active channels", flush=True)


if __name__ == "__main__":
    main()
