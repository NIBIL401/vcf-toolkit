"""Minimal, streaming VCF (v4.x) parser."""
from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from typing import IO, Dict, Iterator, List, Optional, Union

InfoValue = Union[str, bool]

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


class VCFFormatError(ValueError):
    pass


def open_text(path: str, mode: str = "rt") -> IO[str]:
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, mode)  # type: ignore[return-value]
    return open(path, mode)


def parse_info(field_str: str) -> Dict[str, InfoValue]:
    if field_str in (".", ""):
        return {}
    info: Dict[str, InfoValue] = {}
    for item in field_str.split(";"):
        if not item:
            continue
        key, sep, value = item.partition("=")
        info[key] = value if sep else True
    return info


def _is_symbolic(allele: str) -> bool:
    return allele.startswith("<") or "[" in allele or "]" in allele or allele in ("*", ".")


@dataclass
class Variant:
    chrom: str
    pos: int
    id: str
    ref: str
    alt: List[str]
    qual: Optional[float]
    filter: List[str]
    info: Dict[str, InfoValue]
    format: List[str] = field(default_factory=list)
    samples: List[Dict[str, str]] = field(default_factory=list)
    raw: str = ""

    @property
    def is_pass(self) -> bool:
        return self.filter in (["PASS"], ["."], [])

    @property
    def variant_type(self) -> str:
        alts = [a for a in self.alt if not _is_symbolic(a)]
        if not alts:
            return "OTHER"
        if all(len(a) == len(self.ref) == 1 for a in alts):
            return "SNP"
        if all(len(a) == len(self.ref) for a in alts):
            return "MNP"
        if any(len(a) != len(self.ref) for a in alts):
            return "INDEL"
        return "OTHER"

    def is_transition(self) -> Optional[bool]:
        """True/False for biallelic SNPs, None otherwise."""
        if self.variant_type != "SNP" or len(self.alt) != 1:
            return None
        return (self.ref.upper(), self.alt[0].upper()) in TRANSITIONS

    def indel_lengths(self) -> List[int]:
        """Signed lengths: positive = insertion, negative = deletion."""
        return [len(a) - len(self.ref) for a in self.alt
                if not _is_symbolic(a) and len(a) != len(self.ref)]

    def depth(self) -> Optional[int]:
        """INFO/DP, or the sum of FORMAT/DP across samples."""
        dp = self.info.get("DP")
        if isinstance(dp, str) and dp.isdigit():
            return int(dp)
        vals = [s.get("DP") for s in self.samples]
        nums = [int(v) for v in vals if v and v.isdigit()]
        return sum(nums) if nums else None

    def allele_frequency(self) -> Optional[float]:
        af = self.info.get("AF")
        if isinstance(af, str) and af not in (".", ""):
            try:
                return float(af.split(",")[0])
            except ValueError:
                return None
        return None

    def genotype(self, sample_index: int = 0) -> Optional[str]:
        if sample_index >= len(self.samples):
            return None
        return self.samples[sample_index].get("GT")


def genotype_class(gt: Optional[str]) -> str:
    """Classify a GT string as hom_ref / het / hom_alt / missing."""
    if not gt:
        return "missing"
    alleles = gt.replace("|", "/").split("/")
    if any(a == "." for a in alleles):
        return "missing"
    if len(set(alleles)) > 1:
        return "het"
    return "hom_ref" if alleles[0] == "0" else "hom_alt"


class VCFReader:
    """Iterate over the variants of a VCF; header lines are available as ``.header``."""

    def __init__(self, path: str) -> None:
        self.path = path
        self._fh = open_text(path)
        self.meta: List[str] = []
        self.samples: List[str] = []
        self.header_line = ""
        self._first: Optional[str] = None
        for line in self._fh:
            if line.startswith("##"):
                self.meta.append(line.rstrip("\r\n"))
            elif line.startswith("#CHROM"):
                self.header_line = line.rstrip("\r\n")
                self.samples = self.header_line.split("\t")[9:]
            else:
                self._first = line
                break
        if not self.header_line:
            self._fh.close()
            raise VCFFormatError(f"{path}: missing #CHROM header line")

    @property
    def header(self) -> List[str]:
        return self.meta + [self.header_line]

    def close(self) -> None:
        self._fh.close()

    def __enter__(self) -> "VCFReader":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _parse(self, line: str) -> Variant:
        line = line.rstrip("\r\n")
        cols = line.split("\t")
        if len(cols) < 8:
            raise VCFFormatError(f"expected >= 8 columns, got {len(cols)}: {line[:60]!r}")
        fmt = cols[8].split(":") if len(cols) > 8 else []
        samples = [dict(zip(fmt, s.split(":"))) for s in cols[9:]]
        try:
            pos = int(cols[1])
        except ValueError:
            raise VCFFormatError(f"invalid POS {cols[1]!r}") from None
        return Variant(
            chrom=cols[0],
            pos=pos,
            id=cols[2],
            ref=cols[3],
            alt=cols[4].split(","),
            qual=None if cols[5] == "." else float(cols[5]),
            filter=cols[6].split(";"),
            info=parse_info(cols[7]),
            format=fmt,
            samples=samples,
            raw=line,
        )

    def __iter__(self) -> Iterator[Variant]:
        if self._first is not None and self._first.strip():
            yield self._parse(self._first)
        self._first = None
        for line in self._fh:
            if line.strip():
                yield self._parse(line)
