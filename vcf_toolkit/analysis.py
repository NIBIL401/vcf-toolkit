"""Variant summary statistics and filtering."""
from __future__ import annotations

import re
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

from .vcf import Variant, genotype_class

_REGION_RE = re.compile(r"^([^:]+)(?::(\d[\d,]*)(?:-(\d[\d,]*))?)?$")


def parse_region(region: str) -> Tuple[str, int, float]:
    """Parse ``chr``, ``chr:start`` or ``chr:start-end`` (1-based, inclusive)."""
    m = _REGION_RE.match(region.strip())
    if not m:
        raise ValueError(f"invalid region {region!r}")
    chrom, start, end = m.groups()
    s = int(start.replace(",", "")) if start else 1
    e = float(int(end.replace(",", ""))) if end else float("inf")
    if e < s:
        raise ValueError(f"region end < start: {region!r}")
    return chrom, s, e


@dataclass
class VariantFilter:
    min_qual: Optional[float] = None
    min_depth: Optional[int] = None
    max_depth: Optional[int] = None
    min_af: Optional[float] = None
    max_af: Optional[float] = None
    pass_only: bool = False
    types: Set[str] = field(default_factory=set)  # SNP / INDEL / MNP / OTHER
    regions: List[Tuple[str, int, float]] = field(default_factory=list)
    genotype: Optional[str] = None  # het / hom_alt / hom_ref, checked on first sample

    def reason(self, v: Variant) -> Optional[str]:
        """Return why a variant fails, or None if it passes."""
        if self.pass_only and not v.is_pass:
            return "not_pass"
        if self.types and v.variant_type not in self.types:
            return "type"
        if self.min_qual is not None and (v.qual is None or v.qual < self.min_qual):
            return "low_qual"
        if self.min_depth is not None or self.max_depth is not None:
            dp = v.depth()
            if dp is None:
                return "no_depth"
            if self.min_depth is not None and dp < self.min_depth:
                return "low_depth"
            if self.max_depth is not None and dp > self.max_depth:
                return "high_depth"
        if self.min_af is not None or self.max_af is not None:
            af = v.allele_frequency()
            if af is None:
                return "no_af"
            if self.min_af is not None and af < self.min_af:
                return "low_af"
            if self.max_af is not None and af > self.max_af:
                return "high_af"
        if self.regions and not any(c == v.chrom and s <= v.pos <= e for c, s, e in self.regions):
            return "outside_region"
        if self.genotype and genotype_class(v.genotype(0)) != self.genotype:
            return "genotype"
        return None

    def apply(self, variants: Iterable[Variant], rejected: Optional[Counter] = None) -> Iterator[Variant]:
        for v in variants:
            why = self.reason(v)
            if why is None:
                yield v
            elif rejected is not None:
                rejected[why] += 1


class VariantStats:
    def __init__(self) -> None:
        self.total = 0
        self.types: Counter = Counter()
        self.per_chrom: Counter = Counter()
        self.filters: Counter = Counter()
        self.transitions = 0
        self.transversions = 0
        self.multiallelic = 0
        self.indel_lengths: Counter = Counter()
        self.quals: List[float] = []
        self.depths: List[int] = []
        self.genotypes: Dict[str, Counter] = defaultdict(Counter)
        self.substitutions: Counter = Counter()

    def add(self, v: Variant, sample_names: Optional[List[str]] = None) -> None:
        self.total += 1
        self.types[v.variant_type] += 1
        self.per_chrom[v.chrom] += 1
        for f in v.filter:
            self.filters["PASS" if f == "." else f] += 1
        if len(v.alt) > 1:
            self.multiallelic += 1
        ts = v.is_transition()
        if ts is True:
            self.transitions += 1
        elif ts is False:
            self.transversions += 1
        if ts is not None:
            self.substitutions[f"{v.ref.upper()}>{v.alt[0].upper()}"] += 1
        for length in v.indel_lengths():
            self.indel_lengths[length] += 1
        if v.qual is not None:
            self.quals.append(v.qual)
        dp = v.depth()
        if dp is not None:
            self.depths.append(dp)
        names = sample_names or [f"sample{i + 1}" for i in range(len(v.samples))]
        for name, s in zip(names, v.samples):
            self.genotypes[name][genotype_class(s.get("GT"))] += 1

    @property
    def titv(self) -> Optional[float]:
        return self.transitions / self.transversions if self.transversions else None

    def summary(self) -> dict:
        def describe(xs):
            if not xs:
                return None
            return {"min": min(xs), "median": statistics.median(xs),
                    "mean": round(statistics.fmean(xs), 2), "max": max(xs)}

        het_hom = {}
        for name, c in self.genotypes.items():
            het_hom[name] = round(c["het"] / c["hom_alt"], 3) if c["hom_alt"] else None

        return {
            "total_variants": self.total,
            "types": dict(self.types),
            "multiallelic": self.multiallelic,
            "transitions": self.transitions,
            "transversions": self.transversions,
            "ti_tv_ratio": round(self.titv, 3) if self.titv is not None else None,
            "substitutions": dict(sorted(self.substitutions.items())),
            "per_chromosome": dict(self.per_chrom),
            "filters": dict(self.filters),
            "qual": describe(self.quals),
            "depth": describe(self.depths),
            "indel_length_distribution": {str(k): v for k, v in sorted(self.indel_lengths.items())},
            "genotypes": {k: dict(v) for k, v in self.genotypes.items()},
            "het_hom_ratio": het_hom,
        }
