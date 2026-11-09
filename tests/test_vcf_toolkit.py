import gzip
import json
import shutil
from collections import Counter
from pathlib import Path

import pytest

from vcf_toolkit.analysis import VariantFilter, VariantStats, parse_region
from vcf_toolkit.cli import main
from vcf_toolkit.vcf import VCFFormatError, VCFReader, genotype_class, parse_info

EXAMPLE = str(Path(__file__).resolve().parent.parent / "examples" / "sample.vcf")


def load():
    with VCFReader(EXAMPLE) as r:
        return r.samples, list(r)


def test_header_and_samples():
    with VCFReader(EXAMPLE) as r:
        assert r.samples == ["TUMOR", "NORMAL"]
        assert r.meta[0] == "##fileformat=VCFv4.2"
        assert r.header[-1].startswith("#CHROM")


def test_parse_info_flags_and_values():
    assert parse_info("DP=10;DB;AF=0.5") == {"DP": "10", "DB": True, "AF": "0.5"}
    assert parse_info(".") == {}


def test_variant_classification():
    _, variants = load()
    assert len(variants) == 20
    by_pos = {v.pos: v for v in variants}
    assert by_pos[14653].variant_type == "SNP" and by_pos[14653].is_transition()
    assert by_pos[871215].is_transition() is False
    assert by_pos[10177].variant_type == "INDEL" and by_pos[10177].indel_lengths() == [1]
    assert by_pos[25234373].indel_lengths() == [-3]
    assert by_pos[41612].variant_type == "MNP"
    assert by_pos[73820651].variant_type == "OTHER"
    assert by_pos[15274].is_transition() is None  # multi-allelic
    assert by_pos[69511].depth() == 140
    assert by_pos[865628].allele_frequency() == 0.25


@pytest.mark.parametrize("gt,cls", [("0/1", "het"), ("1|0", "het"), ("1/1", "hom_alt"),
                                    ("0/0", "hom_ref"), ("./.", "missing"), (None, "missing"),
                                    ("1/2", "het")])
def test_genotype_class(gt, cls):
    assert genotype_class(gt) == cls


def test_stats_summary():
    samples, variants = load()
    st = VariantStats()
    for v in variants:
        st.add(v, samples)
    s = st.summary()
    assert s["total_variants"] == 20
    assert s["types"] == {"SNP": 14, "INDEL": 4, "MNP": 1, "OTHER": 1}
    assert (s["transitions"], s["transversions"]) == (10, 3)
    assert s["ti_tv_ratio"] == round(10 / 3, 3)
    assert s["multiallelic"] == 1
    assert s["filters"] == {"PASS": 17, "LowQual": 3, "LowDP": 1}
    assert s["genotypes"]["TUMOR"] == {"het": 16, "hom_alt": 4}
    assert s["genotypes"]["NORMAL"]["missing"] == 1
    assert s["indel_length_distribution"] == {"-3": 1, "-1": 1, "1": 1, "3": 1}


def test_filter_pass_snps():
    _, variants = load()
    rejected = Counter()
    kept = list(VariantFilter(pass_only=True, min_qual=30, types={"SNP"}).apply(variants, rejected))
    assert len(kept) == 12
    assert sum(rejected.values()) == 8
    assert rejected["not_pass"] == 3


def test_filter_region_and_genotype():
    _, variants = load()
    kept = list(VariantFilter(regions=[parse_region("chr2:1-1,000,000")]).apply(variants))
    assert [v.pos for v in kept] == [41612, 45895, 216283]
    hom = list(VariantFilter(genotype="hom_alt").apply(variants))
    assert len(hom) == 4


def test_parse_region():
    assert parse_region("chr1") == ("chr1", 1, float("inf"))
    assert parse_region("chr1:100-200") == ("chr1", 100, 200)
    with pytest.raises(ValueError):
        parse_region("chr1:200-100")


def test_missing_header(tmp_path):
    p = tmp_path / "bad.vcf"
    p.write_text("chr1\t1\t.\tA\tG\t.\t.\t.\n")
    with pytest.raises(VCFFormatError):
        VCFReader(str(p))


def test_cli_filter_gz_roundtrip(tmp_path, capsys):
    gz = tmp_path / "in.vcf.gz"
    with open(EXAMPLE, "rb") as src, gzip.open(gz, "wb") as dst:
        shutil.copyfileobj(src, dst)
    out = tmp_path / "out.vcf"
    assert main(["filter", str(gz), "-o", str(out), "--pass-only", "--type", "snp", "--min-depth", "50"]) == 0
    lines = out.read_text().splitlines()
    body = [l for l in lines if not l.startswith("#")]
    assert lines[0] == "##fileformat=VCFv4.2"
    assert any(l.startswith("#CHROM") for l in lines)
    assert all(int(l.split("\t")[7].split("DP=")[1].split(";")[0]) >= 50 for l in body)
    assert "kept" in capsys.readouterr().err

    js = tmp_path / "s.json"
    assert main(["stats", str(out), "--json", str(js)]) == 0
    assert json.loads(js.read_text())["total_variants"] == len(body)
