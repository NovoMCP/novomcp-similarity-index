# NovoMCP Similarity Index

Self-hostable **exact-Tanimoto similarity search over the [NovoMCP Open Corpus](https://registry.opendata.aws/novomcp-open-corpus-lite/)** (122M PubChem compounds). It's the open backend for the "Similar" tab in the NovoMCP sideloads (Chrome, Word) and the `vector_search` tool — point the engine at it and similarity search returns real results against a local index. **No hosted service, no compliance verdicts.**

Powered by [FPSim2](https://github.com/chembl/FPSim2) (MIT, EMBL-EBI): Morgan fingerprints (radius 2, 2048 bits) + exact Tanimoto. Deterministic — same SMILES in, same neighbours out.

## Quick start (Docker)

```bash
# 1. Build an index from the open corpus (anonymous public read; no AWS creds needed)
docker run -v $PWD/index:/index novomcp/similarity-index build \
  --input 's3://novomcp-open-corpus/novomcp-open-corpus-lite/*.parquet' --out /index

# 2. Serve it
docker run -p 8080:8080 -v $PWD/index:/index novomcp/similarity-index serve --index /index

# 3. Point the NovoMCP engine at it
export NOVOMCP_MOLECULE_INDEX_URL=http://localhost:8080
```

Or run a laptop-sized slice first (dev/CI):

```bash
python build_index.py --input 's3://novomcp-open-corpus/novomcp-open-corpus-lite/*.parquet' \
  --out ./index --limit 100000
python server.py --index ./index --mode in-memory
```

> **`--limit` is not a representative sample.** The corpus is partitioned by molecular-weight
> band and `--limit` reads the first shards in order, so the slice skews toward one MW range.
> It's meant for CI / smoke tests. For a usable non-CI subset, prefer `--preset druglike`
> (`MW < 600 & QED > 0.3`), which samples across the corpus by property rather than by shard order.

## API

**`POST /search`**
```json
{ "smiles": "CC(=O)Oc1ccccc1C(=O)O", "top_k": 10, "min_similarity": 0.7,
  "property_filters": { "mw_max": 500, "qed_min": 0.3 } }
```
returns `{ query_smiles, count, results: [{ cid, smiles, similarity, molecular_weight, xlogp, tpsa, qed, has_pains }], notes }`.

**`GET /health`** → `{ status, n_molecules, mode, fp }`.

## Serving modes & resource requirements (measured)

Benchmarked on a 5M slice (Morgan r2/2048), extrapolated to the full 122M corpus:

| mode | 122M query latency | RAM | notes |
|---|---|---|---|
| **`in-memory`** (default) | **~0.25–0.65 s** | **~31 GB** | flagship path: full corpus, fast, exact. Wants a workstation/server. |
| **`on-disk`** | **~30 s/query** | modest | **slice/batch only** — too slow for interactive use at full scale. |
| **slice** (`--limit`) | ms | small | runs on a laptop; dev/CI, not the full deliverable. |

One-time full build: **~1.7 hr** single-core (parallelizable), **~11 GB** HDF5 on disk. Building the whole corpus is a workstation/server job — plan for it up front.

> **Low-RAM full-corpus users:** a `usearch`-ANN mode (approximate, fast, small RAM) is the planned fast-follow. Until then, use `in-memory` on a big-RAM box, or a `--limit`/`--preset` slice. FPSim2 `on-disk` is exact but not interactive at 122M.

## Honesty rules

- **No compliance verdict.** There's no compliance provider here, so no `compliance_status` is emitted — the sideload's compliance column stays honestly blank, never faked. `has_pains` is shipped as a raw structural-alert field only.
- **`exclude_controlled` is a no-op.** The open corpus has no controlled-substance flags; the param is accepted, not applied, and `notes` says so.
- **Unparseable SMILES → HTTP 400**, never an empty 200.
- **Column containment.** `build_index.py` uses a hard default-deny allowlist (`cid, smiles, molecular_weight, xlogp, tpsa, qed, has_pains`); every other column is dropped. Build only from the *published* open corpus.

## Hosting

We host nothing for OSS — you self-host. A managed hosted index is the paid-enterprise offering.

## License

Apache-2.0.
