"""Hermetic tests: build a tiny local corpus, index it, exercise the server.

No network. Covers the build allowlist, the /search contract, min_similarity /
top_k / property filters, the exclude_controlled no-op note, and 400 on bad SMILES.
"""
import os
import sys

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import build_index
import server

# A small, structurally-varied set with deliberate near-neighbours (benzene family,
# aspirin/salicylic acid) so neighbour ordering is meaningful. Props are illustrative.
CORPUS = [
    # cid,    smiles,                         mw,     xlogp, tpsa,  qed,  has_pains
    (241,     "c1ccccc1",                     78.11,  1.9,   0.0,   0.40, False),
    (1140,    "Cc1ccccc1",                    92.14,  2.5,   0.0,   0.45, False),
    (996,     "Oc1ccccc1",                    94.11,  1.5,   20.2,  0.50, False),
    (6115,    "Nc1ccccc1",                    93.13,  1.1,   26.0,  0.48, False),
    (7809,    "Clc1ccccc1",                   112.6,  2.8,   0.0,   0.42, False),
    (2244,    "CC(=O)Oc1ccccc1C(=O)O",        180.16, 1.2,   63.6,  0.55, False),
    (338,     "O=C(O)c1ccccc1O",              138.12, 1.9,   57.5,  0.56, False),
    (2519,    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",   194.19, -1.0,  61.8,  0.54, False),
    (3672,    "CC(C)Cc1ccc(C(C)C(=O)O)cc1",   206.28, 3.5,   37.3,  0.60, False),
    (1983,    "CC(=O)Nc1ccc(O)cc1",           151.16, 0.5,   49.3,  0.58, False),
]
COLS = ["cid", "smiles", "molecular_weight", "xlogp", "tpsa", "qed", "has_pains"]


@pytest.fixture(scope="module")
def index_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("corpus")
    corpus_pq = str(d / "corpus.parquet")
    arrays = [pa.array([row[i] for row in CORPUS]) for i in range(len(COLS))]
    pq.write_table(pa.table(arrays, names=COLS), corpus_pq)
    out = str(d / "index")
    build_index.build(corpus_pq, out)
    return out


@pytest.fixture(scope="module")
def client(index_dir):
    server.load_index(index_dir, "in-memory")
    return TestClient(server.app)


def test_manifest_allowlist(index_dir):
    import json
    m = json.load(open(os.path.join(index_dir, "manifest.json")))
    assert m["n_molecules"] == len(CORPUS)
    assert m["fp"] == "morgan-r2-2048"
    # meta.parquet holds ONLY allowlist fields — nothing else leaked through
    meta = pq.read_table(os.path.join(index_dir, "meta.parquet"))
    assert set(meta.column_names) == set(COLS)


def test_health(client):
    h = client.get("/health").json()
    assert h["status"] == "ok"
    assert h["n_molecules"] == len(CORPUS)
    assert h["mode"] == "in-memory"
    assert h["fp"] == "morgan-r2-2048"


def test_self_similarity_and_contract(client):
    r = client.post("/search", json={"smiles": "CC(=O)Oc1ccccc1C(=O)O", "min_similarity": 0.99}).json()
    assert r["count"] >= 1
    top = r["results"][0]
    assert top["cid"] == 2244 and top["similarity"] == 1.0
    # exact contract keys the sideloads read
    for k in ("cid", "smiles", "similarity", "molecular_weight", "xlogp", "qed"):
        assert k in top
    # honesty: no compliance field is emitted
    assert "compliance_status" not in top and "status" not in top


def test_neighbours_sorted_desc(client):
    r = client.post("/search", json={"smiles": "c1ccccc1", "min_similarity": 0.1, "top_k": 5}).json()
    sims = [x["similarity"] for x in r["results"]]
    assert sims == sorted(sims, reverse=True)
    assert r["results"][0]["cid"] == 241 and r["results"][0]["similarity"] == 1.0


def test_min_similarity_and_top_k(client):
    hi = client.post("/search", json={"smiles": "c1ccccc1", "min_similarity": 0.95}).json()
    assert all(x["similarity"] >= 0.95 for x in hi["results"])
    capped = client.post("/search", json={"smiles": "c1ccccc1", "min_similarity": 0.0, "top_k": 3}).json()
    assert capped["count"] <= 3


def test_property_filter(client):
    unfiltered = client.post("/search", json={
        "smiles": "c1ccccc1", "min_similarity": 0.0, "top_k": 100}).json()
    filtered = client.post("/search", json={
        "smiles": "c1ccccc1", "min_similarity": 0.0, "top_k": 100,
        "property_filters": {"mw_max": 100}}).json()
    assert all(x["molecular_weight"] <= 100 for x in filtered["results"])
    # the filter must actually drop something (the corpus has mols > 100 Da)
    assert filtered["count"] < unfiltered["count"]


def test_exclude_controlled_is_noop_with_note(client):
    r = client.post("/search", json={
        "smiles": "c1ccccc1", "min_similarity": 0.0, "top_k": 100,
        "property_filters": {"exclude_controlled": True},
    }).json()
    assert any("exclude_controlled" in n for n in r["notes"])
    # no-op: same count as without the flag
    base = client.post("/search", json={"smiles": "c1ccccc1", "min_similarity": 0.0, "top_k": 100}).json()
    assert r["count"] == base["count"]


def test_unparseable_smiles_400(client):
    resp = client.post("/search", json={"smiles": "this-is-not-a-smiles!!!"})
    assert resp.status_code == 400
