#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def run(cmd, *, check=True, shell=False, **kwargs):
    """Run a command and echo it."""
    print(f"[RUN] {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    if shell:
        return subprocess.run(cmd, shell=True, check=check, **kwargs)
    return subprocess.run(cmd, shell=False, check=check, **kwargs)


def pass_deletion_positions(vcf_path):
    """Return 0-based? 1-based reference positions covered by PASS deletion alleles.

    For a VCF record POS with REF allele longer than an ALT allele, the
    deleted reference bases are POS+1 .. POS+len(REF)-len(ALT).
    """
    positions = set()
    result = subprocess.run(
        [
            "bcftools", "query",
            "-i", 'FILTER="PASS"',
            "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n",
            str(vcf_path),
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    for line in result.stdout.splitlines():
        chrom, pos_str, ref, alt = line.split("\t")
        pos = int(pos_str)
        for allele in alt.split(","):
            delta = len(ref) - len(allele)
            if delta > 0:
                # The deleted bases are the positions after POS.
                for p in range(pos + 1, pos + delta + 1):
                    positions.add((chrom, p))
    return positions


def main():
    parser = argparse.ArgumentParser(
        description="Mpox ONT CPU variant calling and consensus pipeline."
    )
    parser.add_argument(
        "--raw",
        default="hzcdcib-2.fastq.gz",
        help="Input raw FASTQ file (default: hzcdcib-2.fastq.gz).",
    )
    parser.add_argument(
        "--ref",
        default="ref.fasta",
        help="Reference FASTA file (default: ref.fasta).",
    )
    parser.add_argument(
        "--depth-filter",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable depth>=4 filtering of VCF and consensus masking (default: enabled).",
    )
    parser.add_argument(
        "--haplotype",
        "-H",
        default=None,
        help="bcftools consensus -H option for handling heterozygous genotypes "
             "(e.g. 1, 2, A, R, I). When set, -s SAMPLE is also passed so "
             "FORMAT/GT is honored. Default applies all ALT alleles.",
    )
    parser.add_argument(
        "--snp-only",
        action="store_true",
        help="Apply only SNPs to consensus (bcftools view -v snps).",
    )
    args = parser.parse_args()

    root = Path("./").resolve()
    os.chdir(root)

    raw_path = root / args.raw
    ref_path = root / args.ref

    if not raw_path.exists():
        print(f"[ERROR] Raw FASTQ not found: {raw_path}", file=sys.stderr)
        sys.exit(1)
    if not ref_path.exists():
        print(f"[ERROR] Reference FASTA not found: {ref_path}", file=sys.stderr)
        sys.exit(1)

    threads = 32
    model_name = "r1041_e82_400bps_sup_v500"
    input_dir = root / "clair3_input"
    output_dir = root / "clair3_output"
    results_dir = root / "results"

    print("[INFO] Cleaning previous outputs...")
    for d in (results_dir, input_dir, output_dir):
        if d.exists():
            shutil.rmtree(d)
    results_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== Step 1: fastplong  ===")
    run([
        "fastplong",
        "-i", str(raw_path),
        "-o", str(results_dir / "clean.fastq"),
        "-w", str(threads),
        "-j", str(results_dir / "clean.json"),
        "-h", str(results_dir / "clean.html"),
        "--length_required", "1800",
        "-m", "20",
    ])

    print("=== Step 2: sanitize reference (rename contig -> ref) ===")
    input_ref = input_dir / "ref.fa"
    with open(ref_path, "r") as fin, open(input_ref, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(">ref\n")
            else:
                fout.write(line)
    run(["samtools", "faidx", str(input_ref)])

    print("=== Step 3: minimap2 (map-ont) -> sorted BAM ===")
    minimap = subprocess.Popen(
        [
            "minimap2", "-ax", "map-ont",
            "-t", str(threads),
            "--secondary=no",
            str(input_ref),
            str(results_dir / "clean.fastq"),
        ],
        stdout=subprocess.PIPE,
    )
    sorted_bam = input_dir / "input.bam"
    with minimap.stdout as pipe_in:
        run([
            "samtools", "sort",
            "-@", str(threads),
            "-o", str(sorted_bam),
            "-",
        ], stdin=pipe_in)
    minimap.wait()
    if minimap.returncode != 0:
        sys.exit(minimap.returncode)
    run(["samtools", "index", str(sorted_bam)])

    print("=== Step 4: Clair3 SNP/Indel calling (ONT, CPU) ===")
    uid = os.getuid()
    gid = os.getgid()
    run([
        "docker", "run", "--rm",
        "--user", f"{uid}:{gid}",
        "-v", f"{input_dir}:{input_dir}",
        "-v", f"{output_dir}:{output_dir}",
        "hkubal/clair3:v2.0.1",
        "/opt/bin/run_clair3.sh",
        f"--bam_fn={input_dir}/input.bam",
        f"--ref_fn={input_dir}/ref.fa",
        f"--threads={threads}",
        "--platform=ont",
        f"--model_path=/opt/models/{model_name}",
        f"--output={output_dir}",
        "--include_all_ctgs",
        "--no_phasing_for_fa",
        "--chunk_size=50000",
    ])

    vcf = output_dir / "merge_output.vcf.gz"

    print("=== Step 5: samtools depth -> coverage masks ===")
    depth_file = results_dir / "depth.txt"
    with open(depth_file, "w") as fh:
        run(["samtools", "depth", "-aa", str(sorted_bam)], stdout=fh)

    highcov = results_dir / "highcov.bed"
    lowcov = results_dir / "lowcov.bed"

    if args.depth_filter:
        # Positions that Clair3 called as deleted should not be treated as low
        # coverage in the mask, and are kept through the depth filter.
        del_positions = pass_deletion_positions(vcf)

        with open(depth_file, "r") as fin, open(highcov, "w") as fout:
            for line in fin:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 3:
                    chrom, pos, dep = cols[0], int(cols[1]), int(cols[2])
                    if dep >= 4 or (chrom, pos) in del_positions:
                        fout.write(f"{chrom}\t{pos - 1}\t{pos}\n")

        with open(depth_file, "r") as fin, open(lowcov, "w") as fout:
            chr_name = ""
            start = 0
            end = 0
            for line in fin:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 3:
                    continue
                chrom = cols[0]
                pos = int(cols[1])
                dep = int(cols[2])
                if dep < 4 and (chrom, pos) not in del_positions:
                    if chrom == chr_name and pos == end + 1:
                        end = pos
                    else:
                        if chr_name and (end - start + 1) > 10:
                            fout.write(f"{chr_name}\t{start - 1}\t{end}\n")
                        chr_name = chrom
                        start = pos
                        end = pos
            if chr_name and (end - start + 1) > 10:
                fout.write(f"{chr_name}\t{start - 1}\t{end}\n")
    else:
        # Depth filtering disabled: create empty placeholder BEDs.
        highcov.write_text("")
        lowcov.write_text("")

    print("=== Step 6: filter Clair3 VCF (PASS" + (" + depth>=4" if args.depth_filter else "") + ") ===")
    filtered_vcf = results_dir / "filtered.vcf.gz"
    view_cmd = ["bcftools", "view", "-f", "PASS"]
    if args.snp_only:
        view_cmd.extend(["-v", "snps"])
    if args.depth_filter:
        view_cmd.extend(["-T", str(highcov)])
    view_cmd.append(str(vcf))
    view = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)
    with view.stdout as pipe_in:
        run([
            "bcftools", "norm",
            "-f", str(input_ref),
            "-Oz", "-o", str(filtered_vcf),
        ], stdin=pipe_in)
    view.wait()
    if view.returncode != 0:
        sys.exit(view.returncode)
    run(["bcftools", "index", str(filtered_vcf)])

    print("=== Step 7: bcftools consensus ===")
    consensus = results_dir / "consensus.fasta"
    consensus_cmd = [
        "bcftools", "consensus",
        "-f", str(input_ref),
    ]
    if args.depth_filter:
        consensus_cmd.extend(["-m", str(lowcov)])
    if args.haplotype:
        # Honor FORMAT/GT so -H controls het allele selection.
        consensus_cmd.extend(["-s", "SAMPLE", "-H", args.haplotype])
    consensus_cmd.append(str(filtered_vcf))
    with open(consensus, "w") as fh:
        run(consensus_cmd, stdout=fh)

    print("=== Pipeline finished ===")
    print(f"Clean reads:     {results_dir / 'clean.fastq'}")
    print(f"Alignment:       {input_dir / 'input.bam'}")
    print(f"Clair3 output:   {output_dir}/")
    print(f"Depth file:      {depth_file}")
    print(f"Low-cov mask:    {lowcov}")
    print(f"Filtered VCF:    {filtered_vcf}")
    print(f"Consensus:       {consensus}")


if __name__ == "__main__":
    main()
