# Extremal integer lifts of MDS matrices: reproducibility artifact

Code, exact finite certificates, numerical tables and plotting scripts accompanying **“Extremal integer lifts of MDS matrices and exact reduction placement.”**

**Author and contact:** Phuc-Phan Duong — duongphucphan@actvn.edu.vn.

This repository is independent of the manuscript build. **No LaTeX installation is required.** The mathematical checks use only Python's standard library. Matplotlib is needed only to regenerate figures.

The supplied data include the classification of four-dimensional integer MDS lifts, fixed-circuit normalization schedules, the restricted cancellation-free A49 catalog, cross-class witnesses, observability checks and the BN254 scalar-field calculations. Historical numbers in names such as `revision6_checks.py` identify the corresponding verification additions; they are not separate dependencies to install.

## Choose a workflow

| Goal | What to run | What changes |
|---|---|---|
| Inspect the supplied results | Open the JSON files in `data/` and figures in `figures/` | Nothing |
| Regenerate the figures from supplied verified data | Install `requirements.txt`, then run the two plotting scripts | Figure files only; no mathematical experiment is repeated |
| Check the latest BN254 and A49 cone claims | Run `code/revision6_checks.py` | Rewrites its own verification report and checks that the other data files remain unchanged |
| Recompute the mathematical experiments | Follow the complete nine-command sequence below | Recreates the associated data files, including the finite circuit catalog |

Python 3.10 or later is recommended. The packaged smoke test used Python 3.12 and Matplotlib 3.10.8. `package_manifest.json` records the actual interpreter and plotting-package versions used for that test. The supplied figures can be read without installing either package.

## Quick start in PyCharm

1. Extract the archive into a folder. Open that folder as a PyCharm project; `README.md`, `code/` and `data/` should be directly inside the project root.
2. Select a Python interpreter in the project settings. An existing interpreter or a new virtual environment is sufficient.
3. Open `code/revision6_checks.py` and choose **Run**. This check needs no third-party Python packages.
4. To draw figures, open the PyCharm terminal and install the pinned plotting dependency:

   ```bash
   python -m pip install -r requirements.txt
   ```

5. Run `code/make_figures.py`, then `code/make_supplementary_figures.py`. Both use a noninteractive plotting backend and save their output to disk; no plot window needs to appear.

Scripts locate the project root from their own file paths, so they also work when the IDE chooses a different working directory. Keep the `code/`, `data/` and `figures/` directories together. If a second interpreter is selected later, install the plotting dependency in that interpreter as well.

The comments at the top of the original plotting script mention `requirements-figures.txt`, the filename used in the manuscript project. **In this standalone repository the same pinned requirements are supplied as `requirements.txt`.**

## Regenerate plots from the supplied data

Run from the repository root:

```bash
python -m pip install -r requirements.txt
python code/make_figures.py
python code/make_supplementary_figures.py
```

The first script reads `data/extended_results.json` and `data/class_witnesses.json`. It writes:

- `figures/vertex_cover_gadget.pdf`, `.png`, `.svg` and `.eps`;
- `figures/reduction_counts.pdf`, `.png`, `.svg` and `.eps`;
- the submission-compatible copies `Fig1.pdf`, `Fig1.eps`, `Fig2.pdf` and `Fig2.eps` in the repository root.

The PNG exports use 600 dpi. PDF, SVG and EPS are vector exports. The second script reads `data/extended_results.json` and writes `figures/state_counts.pdf` and `figures/order_states.pdf`. These two supplementary plots describe optimizer state counts and traversal-order effects.

Neither plotting script writes data files. Generated PDFs and SVGs may differ byte-for-byte across runs because of metadata or internal identifiers, even when their plotted values are unchanged.

## Focused mathematical checks

```bash
python code/revision5_checks.py
python code/revision6_checks.py
```

The first command independently checks the observability determinants, all 24 row permutations of each of A7 and A17, the selected output permutation and its inherited schedules. It includes 57,122 centered-field input/schedule evaluations over F13 and writes `data/revision5_checks.json`.

The second command checks the exact BN254 **scalar** modulus

```text
21888242871839275222246405745257275088548364400416034343698204186575808495617
```

It verifies centered signed 256-bit capacity H=5; one-shot widths 256, 257 and 258 for row norms 5, 7 and 16; 414 invertible minors across six candidate matrices; four eleven-operation schedules; and nine word-capacity cases for BabyBear, Mersenne-31 and KoalaBear. It also recomputes the A49 finite certificate: 156 first-normalization forms, four successful forms, sixteen unique output representations, short-program exclusions and no-sharing predicates. The two-gate example `2*(x+y)` explicitly demonstrates why load equality alone does not imply a three-gate lower bound. The command writes only `data/revision6_checks.json`.

These focused checks use the included source code and data. They do not require re-enumerating the complete cancellation-free circuit catalog.

## Recompute the mathematical results

To rebuild the experimental data, use this dependency order:

```bash
python code/verify.py
python code/experiments.py
python code/catalog.py
python code/check_five_certificate.py
python code/revision3_checks.py
python code/class_extension.py
python code/revision4_checks.py
python code/revision5_checks.py
python code/revision6_checks.py
```

Run the plotting commands afterward if updated figures are wanted. The complete finite catalog and all reports are already supplied, so this sequence is optional when only reproducing plots. Catalog enumeration and exhaustive range checks perform substantially more work than the focused checks. No fixed execution time or hardware-independent performance claim is made.

| Script | Main purpose | Files written under `data/` |
|---|---|---|
| `verify.py` | Signed and positive classifications, minor cross-checks, DP versus subset enumeration, circuit identities and finite-field semantics | `results.json` |
| `experiments.py` | Baseline capacity sweeps, traversal-order experiments and integer widths | `extended_results.json`, `circuit_summary.json` |
| `catalog.py` | Full catalog within the documented cancellation-free grammar and the separate five-normalization certificate | `five_reduction_certificate.json`, `cancellation_free_catalog.jsonl`, `catalog_results.json` |
| `check_five_certificate.py` | Independent signed-sumset certificate check and fifteen-operation witness semantics | `independent_synthesis_check.json` |
| `revision3_checks.py` | Shared normalization lower bounds and the separate unit-cost tripling model | `revision3_checks.json` |
| `class_extension.py` | Reachability of first-normalization forms for all three signed classes | `all_class_reachability.json` |
| `revision4_checks.py` | A7/A17 witnesses, exact ranges, named-field certificates and block composition | `class_witnesses.json`, `exact_ranges.json`, `revision4_checks.json` |
| `revision5_checks.py` | Observability census, the permuted A7 witness and selected-schedule range checks | `revision5_checks.json` |
| `revision6_checks.py` | BN254 word widths and schedules; focused A49 formal-cone checks | `revision6_checks.json` |

`code/reduction_model.py` supplies the circuit model and exact placement optimizer. `code/classify.py` supplies the matrix classification and can also print the signed classification directly when run. These modules do not need a separate installation.

`data/a17_search_witness.json` records the original heuristic witness search. To rerun that optional search, use `python code/class_extension.py --search`; it rewrites both that witness file and `all_class_reachability.json`. The default command above performs the exhaustive reachability checks without rerunning the heuristic. The search is used to obtain an explicit circuit, not to establish global arithmetic optimality.

Runtime and environment fields in recomputed reports can change. The regeneration commands intentionally overwrite the corresponding supplied reports. A copy of the downloaded archive preserves the original release bytes.

## Interpreting the results

The normalization optimizer is exact for the stated compositional absolute-load certificate model and the supplied circuit DAG. It is not a technology-mapped area, delay or energy model. Raw operations must fit the capacity before a normalization can occur. Canonical outputs are part of the stated cost model.

The A49 fifteen-operation lower bound at H=4 applies to the specified binary add/subtract/double gate model and the formal-cone argument. The finite catalog is exhaustive only within its documented cancellation-free grammar; it is not an enumeration of all integer straight-line programs. A7, A17 and the selected row permutation of A7 have explicit eleven-operation witnesses, but eleven is not asserted to be their globally minimum arithmetic count.

The observability checks concern nonzero invariant subspaces contained in a coordinate hyperplane, as used for a partial SPN with one fixed nonlinear coordinate. Passing this filter does not establish complete cipher security. An output permutation preserves the supplied arithmetic DAG and its modeled costs, but changes the repeated-round primitive when a nonlinear coordinate is fixed.

The BN254 signed-word calculations use the scalar-field modulus, not the elliptic-curve base-field modulus. Exact bit widths and normalization counts do not by themselves measure a speedup on a processor or a multi-limb implementation. The F13 range-attainment results apply to the selected circuits and schedules that were exhaustively checked.

## Source provenance and licenses

`provenance/` contains the pinned Plonky3 v0.4.0 Rust baseline snapshots, their original MIT and Apache-2.0 license files, and source-location/hash records. It also contains the Poseidon2 composition source record. The complete HorizenLabs source is not redistributed here.

The baseline arithmetic counts refer to isolated four-word kernels. The source record distinguishes the reference Poseidon2 width-four convention from the pinned Plonky3 wrapper's width-four convention. This distinction matters when using a kernel inside a complete external layer.

No license grant for the original project code or data is added by this package. Third-party material retains the rights and notices stated in its supplied licenses. No publication DOI, artifact DOI or repository URL is asserted before one is assigned. `CITATION.cff` identifies the artifact and its supplied author/contact metadata.

## Integrity and package validation

`SHA256SUMS` records SHA-256 digests for the packaged files other than the checksum file itself. `package_manifest.json` records payload sizes/hashes and the packaging smoke test. That test extracts a fresh archive, runs `revision6_checks.py` and both plotting scripts, verifies the expected figure files and confirms that all fifteen supplied data files retain their original hashes. It does not claim to rerun the full nine-command mathematical reproduction during packaging.

No GitHub workflow or publishing automation is included. Upload the extracted repository files after choosing the actual repository location; the package does not create or publish a repository automatically.
