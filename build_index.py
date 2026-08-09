#!/usr/bin/env python3
"""Build a NovoMCP similarity index from the open corpus.

Reads the published NovoMCP Open Corpus (lite) parquet — a local path/glob or an
``s3://`` glob (anonymous public read) — and writes an index directory containing:

  <out>/index.h5       FPSim2 fingerprint DB (Morgan r2/2048); mol-id = PubChem CID
  <out>/meta.parquet   default-deny allowlist of display/filter fields, keyed by CID
  <out>/manifest.json  {n_molecules, fp, source, built_at, ...}

Build ONLY from the *published* open corpus (already stripped of compliance /
controlled-substance columns). The allowlist is a hard backstop: only
cid/smiles/molecular_weight/xlogp/tpsa/qed/has_pains pass through — every other
column is dropped (fail closed), so even a richer parquet can never smuggle
compliance / controlled-substance columns into a public index.
"""
import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone

# --- default-deny: the ONLY columns that may leave the corpus into the index.
#     Everything else in the source parquet is dropped (fail closed) so no
#     compliance / controlled-substance column can reach a public index. --------------
ALLOWLIST = ["cid", "smiles", "molecular_weight", "xlogp", "tpsa", "qed", "has_pains"]
META_FIELDS = [c for c in ALLOWLIST if c != "cid"]  # meta table is keyed by cid

FP_TYPE = "Morgan"
# NOTE: FPSim2 >=0.7 / RDKit >=2023 build via rdFingerprintGenerator.GetMorganGenerator,
# whose bit-count kwarg is `fpSize` (the legacy `nBits` keyword is rejected).
FP_PARAMS = {"radius": 2, "fpSize": 2048}

PRESETS = {
    "druglike": "molecular_weight < 600 AND qed > 0.3",
}


def _connect(input_glob):
    import duckdb

    con = duckdb.connect()
    con.execute("INSTALL httpfs; LOAD httpfs;")
    if input_glob.startswith("s3://"):
        # Anonymous public-bucket read — the open corpus lives in a public AWS Open Data
        # bucket, so an OSS user needs no AWS credentials. If you DO have credentials and
        # need them (e.g. a private mirror), export them and this still works via the
        # default provider chain; drop the `anon` secret to force signed requests.
        con.execute("SET s3_region='us-east-2';")
        try:
            con.execute("CREATE SECRET anon (TYPE s3, PROVIDER config, REGION 'us-east-2');")
        except Exception:
            pass
    return con


def _n_indexed(h5_path):
    """Count fingerprints actually stored in the FPSim2 db (parse successes)."""
    import tables

    with tables.open_file(h5_path, "r") as f:
        for node in f.walk_nodes("/", "Leaf"):
            if node.name == "fps":
                return int(node.nrows)
    # Fallback: load on-disk engine and measure.
    from FPSim2 import FPSim2Engine

    return int(FPSim2Engine(h5_path, in_memory_fps=False).fps.shape[0])


def build(input_glob, out_dir, limit=None, preset=None):
    import pyarrow.parquet as pq
    from FPSim2.io import create_db_file

    os.makedirs(out_dir, exist_ok=True)
    h5 = os.path.join(out_dir, "index.h5")
    meta_path = os.path.join(out_dir, "meta.parquet")
    manifest_path = os.path.join(out_dir, "manifest.json")

    con = _connect(input_glob)
    where = "smiles IS NOT NULL"
    if preset:
        if preset not in PRESETS:
            sys.exit(f"[build] unknown preset {preset!r}; choices: {list(PRESETS)}")
        where += f" AND ({PRESETS[preset]})"
    # NOTE: --limit is NOT a random sample. The corpus is sharded by molecular-weight
    # band, and this reads the first shards in order, so a --limit slice is MW-skewed
    # (fine for CI; use --preset for a usable non-CI subset).
    lim = f" LIMIT {int(limit)}" if limit else ""
    cols = ", ".join(ALLOWLIST)  # default-deny projection at the source
    query = f"SELECT {cols} FROM read_parquet('{input_glob}') WHERE {where}{lim}"

    print(f"[build] reading corpus: {input_glob}", flush=True)
    t0 = time.time()
    tbl = con.execute(query).fetch_arrow_table()
    n_rows = tbl.num_rows
    print(f"[build] pulled {n_rows:,} rows in {time.time() - t0:.0f}s", flush=True)
    if n_rows == 0:
        sys.exit("[build] ABORT: 0 rows matched — check the input path/glob and filters.")

    # Hard allowlist assertion (fail closed) even if the source had extra columns.
    extra = [c for c in tbl.column_names if c not in ALLOWLIST]
    if extra:
        sys.exit(f"[build] ABORT: non-allowlisted columns present: {extra}")

    smis = tbl.column("smiles").to_pylist()
    cids = tbl.column("cid").to_pylist()
    mols = [[s, int(c)] for s, c in zip(smis, cids) if s]

    print(f"[build] building FPSim2 db (Morgan r2/2048) for {len(mols):,} mols ...", flush=True)
    t0 = time.time()
    # FPSim2 skips unparseable SMILES internally; the db may hold fewer than len(mols).
    create_db_file(mols, h5, "smiles", FP_TYPE, FP_PARAMS)
    build_s = time.time() - t0
    print(f"[build] FP db done in {build_s:.0f}s -> {h5} ({os.path.getsize(h5) / 1e6:.0f} MB)", flush=True)

    # meta sidecar: allowlist fields, sorted by cid so the server can binary-search.
    meta = tbl.select(["cid"] + META_FIELDS).sort_by("cid")
    pq.write_table(meta, meta_path, compression="zstd")
    print(f"[build] meta sidecar -> {meta_path}", flush=True)

    n_indexed = _n_indexed(h5)
    manifest = {
        "n_molecules": n_indexed,
        "n_input_rows": n_rows,
        "n_skipped_unparseable": n_rows - n_indexed,
        "fp": "morgan-r2-2048",
        "fp_params": FP_PARAMS,
        "source": input_glob,
        "preset": preset,
        "allowlist": ALLOWLIST,
        "built_at": datetime.now(timezone.utc).isoformat(),
        "build_seconds": round(build_s, 1),
    }
    json.dump(manifest, open(manifest_path, "w"), indent=2)
    print(f"[build] manifest -> {manifest_path}\n{json.dumps(manifest, indent=2)}", flush=True)
    return manifest


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True,
                    help="Corpus parquet: local path/glob or s3:// glob "
                         "(e.g. s3://novomcp-open-corpus/novomcp-open-corpus-lite/*.parquet)")
    ap.add_argument("--out", required=True, help="Output index directory.")
    ap.add_argument("--limit", type=int, default=None,
                    help="Quick-test slice size (dev/CI only). NOT a representative sample: it "
                         "reads the first shards in order, and the corpus is partitioned by "
                         "molecular-weight band, so a --limit slice skews toward one MW range. "
                         "For a usable non-CI subset prefer --preset (e.g. druglike).")
    ap.add_argument("--preset", choices=list(PRESETS), default=None,
                    help="Subset filter, e.g. 'druglike' (MW<600 & QED>0.3).")
    args = ap.parse_args()
    build(args.input, args.out, limit=args.limit, preset=args.preset)


if __name__ == "__main__":
    main()
