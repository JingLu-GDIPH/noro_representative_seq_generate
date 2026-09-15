#!/usr/bin/env python3
"""Run IPHnano norovirus VP1/RdRp genotyping for final consensus FASTA files."""

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


DEFAULT_GENOTYPER_DIR = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_genotyping")
DEFAULT_PYTHON = Path("/Users/LuJ/mambaforge/envs/noro-consensus/bin/python")
DEFAULT_MAFFT = Path("/Users/LuJ/mambaforge/envs/noro-consensus/bin/mafft")


def write_mafft_wrapper(wrapper_dir, mafft_path):
    wrapper_dir.mkdir(parents=True, exist_ok=True)
    wrapper = wrapper_dir / "mafft"
    wrapper.write_text(
        "#!/bin/sh\n"
        f'exec "{mafft_path}" --quiet "$@"\n'
    )
    wrapper.chmod(0o755)
    return wrapper


def run_one(label, input_fasta, output_prefix, args, wrapper_dir):
    script = args.genotyper_dir / "genotype_norovirus.py"
    ref_dir = args.genotyper_dir / "ref_seq"
    if not script.is_file():
        raise FileNotFoundError(f"Cannot find genotyping script: {script}")
    if not ref_dir.is_dir():
        raise FileNotFoundError(f"Cannot find reference directory: {ref_dir}")
    if not input_fasta.is_file():
        raise FileNotFoundError(f"Cannot find input FASTA for {label}: {input_fasta}")

    env = os.environ.copy()
    env["PATH"] = os.pathsep.join(
        [
            str(wrapper_dir),
            str(args.python.parent),
            "/usr/local/bin",
            "/usr/bin",
            "/bin",
            "/usr/sbin",
            "/sbin",
        ]
    )

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(args.python),
        str(script),
        "--input",
        str(input_fasta),
        "--ref-dir",
        str(ref_dir),
        "--output-prefix",
        str(output_prefix),
        "--threads",
        str(args.threads),
    ]
    print(f"[{label}] running: {' '.join(command)}")
    subprocess.run(command, check=True, env=env)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run the IPHnano norovirus dual-region genotyper on final GI/GII "
            "consensus FASTA files with a MAFFT --quiet wrapper for sandbox-safe "
            "parallel execution."
        )
    )
    parser.add_argument("--gi-fasta", type=Path, default=Path("results_new/gi/07_final_results/all_final_consensus.fasta"))
    parser.add_argument("--gii-fasta", type=Path, default=Path("results_new/gii/07_final_results/all_final_consensus.fasta"))
    parser.add_argument("--outdir", type=Path, default=Path("results_new/genotyping_consensus"))
    parser.add_argument("--genotyper-dir", type=Path, default=DEFAULT_GENOTYPER_DIR)
    parser.add_argument("--python", type=Path, default=DEFAULT_PYTHON)
    parser.add_argument("--mafft", type=Path, default=DEFAULT_MAFFT)
    parser.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 8) - 2))
    parser.add_argument(
        "--group",
        choices=("gi", "gii", "both"),
        default="both",
        help="Which consensus FASTA set to genotype.",
    )
    args = parser.parse_args()

    if not args.python.is_file():
        found = shutil.which("python3")
        if not found:
            raise FileNotFoundError(f"Python not found: {args.python}")
        args.python = Path(found)
    if not args.mafft.is_file():
        found = shutil.which("mafft")
        if not found:
            raise FileNotFoundError(f"MAFFT not found: {args.mafft}")
        args.mafft = Path(found)

    with tempfile.TemporaryDirectory(prefix="noro_genotyping_mafft_") as tmpdir:
        wrapper_dir = Path(tmpdir)
        write_mafft_wrapper(wrapper_dir, args.mafft)

        if args.group in {"gi", "both"}:
            run_one(
                "GI",
                args.gi_fasta,
                args.outdir / "gi" / "gi_final_consensus",
                args,
                wrapper_dir,
            )
        if args.group in {"gii", "both"}:
            run_one(
                "GII",
                args.gii_fasta,
                args.outdir / "gii" / "gii_final_consensus",
                args,
                wrapper_dir,
            )


if __name__ == "__main__":
    main()
