"""Command-line interface: ``vcfkit stats`` and ``vcfkit filter``."""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from typing import List, Optional

from . import __version__
from .analysis import VariantFilter, VariantStats, parse_region
from .vcf import VCFFormatError, VCFReader, open_text


def _print_stats(s: dict) -> None:
    print(f"Total variants     {s['total_variants']}")
    for t, n in sorted(s["types"].items()):
        print(f"  {t:<16} {n}")
    print(f"Multi-allelic      {s['multiallelic']}")
    print(f"Ti / Tv            {s['transitions']} / {s['transversions']}  (ratio {s['ti_tv_ratio']})")
    if s["qual"]:
        q = s["qual"]
        print(f"QUAL               min {q['min']}  median {q['median']}  mean {q['mean']}  max {q['max']}")
    if s["depth"]:
        d = s["depth"]
        print(f"Depth              min {d['min']}  median {d['median']}  mean {d['mean']}  max {d['max']}")
    print("Filters            " + ", ".join(f"{k}={v}" for k, v in s["filters"].items()))
    print("Per chromosome     " + ", ".join(f"{k}={v}" for k, v in s["per_chromosome"].items()))
    for name, g in s["genotypes"].items():
        print(f"Sample {name:<11} {g}  het/hom={s['het_hom_ratio'][name]}")


def cmd_stats(args: argparse.Namespace) -> int:
    stats = VariantStats()
    with VCFReader(args.input) as reader:
        for v in reader:
            stats.add(v, reader.samples)
    summary = stats.summary()
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(summary, fh, indent=2)
    _print_stats(summary)
    return 0


def cmd_filter(args: argparse.Namespace) -> int:
    flt = VariantFilter(
        min_qual=args.min_qual,
        min_depth=args.min_depth,
        max_depth=args.max_depth,
        min_af=args.min_af,
        max_af=args.max_af,
        pass_only=args.pass_only,
        types={t.upper() for t in args.type},
        regions=[parse_region(r) for r in args.region],
        genotype=args.genotype,
    )
    rejected: Counter = Counter()
    kept = 0
    with VCFReader(args.input) as reader:
        out = open_text(args.output, "wt") if args.output != "-" else sys.stdout
        try:
            cmdline = " ".join(sys.argv[1:]) or "filter"
            for line in reader.meta:
                out.write(line + "\n")
            out.write(f"##vcfkitCommand={cmdline}\n")
            out.write(reader.header_line + "\n")
            for v in flt.apply(reader, rejected):
                out.write(v.raw + "\n")
                kept += 1
        finally:
            if out is not sys.stdout:
                out.close()
    total = kept + sum(rejected.values())
    print(f"kept {kept} / {total} variants", file=sys.stderr)
    for why, n in rejected.most_common():
        print(f"  removed ({why}): {n}", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="vcfkit", description="Summarise and filter VCF variant calls.")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("stats", help="summary statistics (types, Ti/Tv, depth, genotypes ...)")
    s.add_argument("input", help="VCF file (.vcf or .vcf.gz)")
    s.add_argument("--json", help="also write the summary as JSON")
    s.set_defaults(func=cmd_stats)

    f = sub.add_parser("filter", help="filter variants by quality, depth, type, region ...")
    f.add_argument("input")
    f.add_argument("-o", "--output", default="-", help="output VCF ('-' = stdout, .gz = gzip)")
    f.add_argument("--min-qual", type=float)
    f.add_argument("--min-depth", type=int)
    f.add_argument("--max-depth", type=int)
    f.add_argument("--min-af", type=float, help="minimum INFO/AF (first ALT)")
    f.add_argument("--max-af", type=float)
    f.add_argument("--pass-only", action="store_true", help="keep only FILTER=PASS")
    f.add_argument("--type", action="append", default=[], choices=["snp", "indel", "mnp", "other", "SNP", "INDEL", "MNP", "OTHER"],
                   help="keep only this variant type (repeatable)")
    f.add_argument("--region", action="append", default=[], help="chr[:start[-end]] (repeatable)")
    f.add_argument("--genotype", choices=["het", "hom_alt", "hom_ref"], help="genotype of the first sample")
    f.set_defaults(func=cmd_filter)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (VCFFormatError, ValueError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    except FileNotFoundError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
