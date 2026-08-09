#!/usr/bin/env python3
"""Serve a NovoMCP similarity index.

  POST /search   exact-Tanimoto nearest neighbours over the open corpus
  GET  /health   {status, n_molecules, mode, fp}

The /search response is the exact shape the NovoMCP engine forwards to the Chrome /
Word sideloads' "Similar" tab — top-level `results`, each item carrying
`cid, smiles, similarity, molecular_weight, xlogp, qed` (+ tpsa, has_pains). Do not
rename these keys; a shape drift is the "works but shows nothing" bug.

Honesty rules (this is the whole point of the OSS pass):
  * No compliance verdict. There is no compliance provider here, so no `compliance_status`
    is emitted — the sideload's compliance column stays honestly blank, never faked.
  * `exclude_controlled` is accepted but is a NO-OP (the open corpus has no
    controlled-substance flags) and we say so in `notes`.
  * A SMILES that does not parse returns 400, never an empty 200.
"""
import argparse
import json
import os
from typing import Optional

import numpy as np
import pyarrow.parquet as pq
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

app = FastAPI(title="NovoMCP Similarity Index", version="0.1.0")

# Populated at startup by load_index().
STATE = {
    "engine": None,      # FPSim2Engine
    "search": None,      # bound method: similarity | on_disk_similarity
    "mode": None,        # "in-memory" | "on-disk"
    "manifest": {},
    "cids": None,        # sorted int64 np.ndarray for binary-search lookup
    "meta": None,        # pyarrow.Table (sorted by cid), display/filter fields
}

META_FIELDS = ["smiles", "molecular_weight", "xlogp", "tpsa", "qed", "has_pains"]


class PropertyFilters(BaseModel):
    mw_min: Optional[float] = None
    mw_max: Optional[float] = None
    qed_min: Optional[float] = None
    qed_max: Optional[float] = None
    exclude_controlled: Optional[bool] = None  # accepted; no-op (see honesty rules)


class SearchRequest(BaseModel):
    smiles: str
    top_k: int = Field(10, ge=1, le=1000)
    min_similarity: float = Field(0.7, ge=0.0, le=1.0)
    property_filters: Optional[PropertyFilters] = None


def load_index(index_dir: str, mode: str):
    from FPSim2 import FPSim2Engine

    h5 = os.path.join(index_dir, "index.h5")
    meta_path = os.path.join(index_dir, "meta.parquet")
    manifest_path = os.path.join(index_dir, "manifest.json")
    for p in (h5, meta_path):
        if not os.path.exists(p):
            raise SystemExit(f"[serve] missing index file: {p} (did you run build_index.py?)")

    in_mem = mode == "in-memory"
    engine = FPSim2Engine(h5, in_memory_fps=in_mem)
    STATE["engine"] = engine
    # In-memory uses .similarity(); on-disk MUST use .on_disk_similarity()
    # (.similarity() on an on-disk engine raises "FPs not loaded into memory").
    STATE["search"] = engine.similarity if in_mem else engine.on_disk_similarity
    STATE["mode"] = mode

    meta = pq.read_table(meta_path)  # sorted by cid at build time
    STATE["meta"] = meta
    STATE["cids"] = meta.column("cid").to_numpy()
    STATE["manifest"] = json.load(open(manifest_path)) if os.path.exists(manifest_path) else {}
    n = STATE["manifest"].get("n_molecules", len(STATE["cids"]))
    print(f"[serve] loaded index: {n:,} molecules, mode={mode}", flush=True)


def _meta_for(cid: int):
    """Binary-search the sorted cid array; return the allowlist fields or None."""
    cids = STATE["cids"]
    i = int(np.searchsorted(cids, cid))
    if i >= len(cids) or int(cids[i]) != cid:
        return None
    meta = STATE["meta"]
    return {f: meta.column(f)[i].as_py() for f in META_FIELDS}


@app.get("/health")
def health():
    m = STATE["manifest"]
    return {
        "status": "ok" if STATE["engine"] is not None else "loading",
        "n_molecules": m.get("n_molecules", len(STATE["cids"]) if STATE["cids"] is not None else 0),
        "mode": STATE["mode"],
        "fp": m.get("fp", "morgan-r2-2048"),
    }


@app.post("/search")
def search(req: SearchRequest):
    if STATE["engine"] is None:
        raise HTTPException(status_code=503, detail="index not loaded")

    from rdkit import Chem

    if Chem.MolFromSmiles(req.smiles) is None:
        raise HTTPException(status_code=400, detail=f"could not parse SMILES: {req.smiles!r}")

    notes = []
    pf = req.property_filters or PropertyFilters()
    if pf.exclude_controlled:
        notes.append("exclude_controlled ignored: no compliance provider in this index")

    # FPSim2 returns (mol_id, tanimoto) sorted desc, already above threshold.
    hits = STATE["search"](req.smiles, req.min_similarity, n_workers=1)

    results = []
    for row in hits:
        cid = int(row["mol_id"])
        sim = float(row["coeff"])
        meta = _meta_for(cid)
        if meta is None:
            continue
        if pf.mw_min is not None and (meta["molecular_weight"] or 0) < pf.mw_min:
            continue
        if pf.mw_max is not None and (meta["molecular_weight"] or 1e9) > pf.mw_max:
            continue
        if pf.qed_min is not None and (meta["qed"] or 0) < pf.qed_min:
            continue
        if pf.qed_max is not None and (meta["qed"] or 1e9) > pf.qed_max:
            continue
        results.append({
            "cid": cid,
            "smiles": meta["smiles"],
            "similarity": round(sim, 4),
            "molecular_weight": meta["molecular_weight"],
            "xlogp": meta["xlogp"],
            "tpsa": meta["tpsa"],
            "qed": meta["qed"],
            "has_pains": meta["has_pains"],
        })
        if len(results) >= req.top_k:
            break

    return {
        "query_smiles": req.smiles,
        "count": len(results),
        "results": results,
        "notes": notes,
    }


def main():
    ap = argparse.ArgumentParser(description="Serve a NovoMCP similarity index.")
    ap.add_argument("--index", required=True, help="Index directory from build_index.py.")
    ap.add_argument("--mode", choices=["in-memory", "on-disk"], default="in-memory",
                    help="in-memory: fast, needs RAM for the whole FP db (~31 GB at 122M). "
                         "on-disk: low RAM, slower (~30 s/query at 122M — slice/batch only).")
    ap.add_argument("--host", default="0.0.0.0")
    ap.add_argument("--port", type=int, default=8080)
    args = ap.parse_args()

    import uvicorn

    load_index(args.index, args.mode)
    uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
