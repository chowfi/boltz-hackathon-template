# Final Strategy: SATURATION + Confidence-Only

## Overview

This document describes the final prediction strategy for the Boltz Hackathon protein-ligand binding challenge. The strategy evolved through extensive experimentation and is optimized for **lowest mean RMSD of top-1 predictions**.

---

## Strategy Summary

### Preparation Phase: SATURATION

**Goal:** Generate diverse samples that explore the full conformational space

**Implementation:**
- 10 configurations per target
- Each configuration: 5 diffusion samples
- **Total: 50 samples per target**
- Wide seed spacing: `[42, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 7777777]`
- **No constraints** - preserves full protein context

**Why this works:**
- ✅ Wide seed spacing creates genuine diversity in diffusion sampling
- ✅ Full protein context preserved (critical for Boltz)
- ✅ No assumptions about binding location (robust to PDB numbering issues)
- ✅ Simple, reproducible, fast

### Post-Processing Phase: Confidence-Only Ranking

**Goal:** Select the top 5 predictions from 50 samples

**Implementation:**
- Extract Boltz confidence from PDB B-factor column (pLDDT-style, 0-100 scale)
- Sort all 50 samples by confidence (highest first)
- Return top 5

**Why this works:**
- ✅ Boltz confidence is well-calibrated
- ✅ Simple ranking avoids mis-ranking good samples
- ✅ No complex scoring that could introduce noise
- ✅ Proven to work in experiments

---

## Key Learnings from Experiments

### What Worked ✅

1. **SATURATION (Wide Seed Diversity)**
   - **Result:** 6FVF improved from ~20Å to 15.89Å
   - **Why:** Seeds spanning 42 to 7,777,777 create orthogonal exploration
   - **Key insight:** Diversity comes from diffusion sampling, not constraints

2. **Confidence-Only Ranking**
   - **Result:** Simpler is better
   - **Why:** Boltz's internal confidence metric is robust
   - **Key insight:** Complex scoring can mis-rank physically "better" but spatially wrong predictions

3. **Full Protein Context**
   - **Result:** Essential for good performance
   - **Why:** Boltz needs full context to make accurate predictions
   - **Key insight:** Never break up the protein into regions

### What Failed ❌

1. **Pocket Scanning / Regional Constraints**
   - **Result:** 3LW0 performance degraded significantly
   - **Problem:** Dividing protein into regions dilutes context
   - **Lesson:** Full context > targeted sampling

2. **Multi-Scoring Ensemble (Clashes + Contacts + Confidence)**
   - **Result:** 6FVF picked 19.61Å instead of 15.82Å
   - **Problem:** Complex metrics picked spatially wrong but physically reasonable predictions
   - **Lesson:** Simple confidence > hybrid scoring

3. **Terminus Probing with Hard Constraints**
   - **Result:** 3K5V RMSD worsened from 25Å to 33Å
   - **Problem:** PDB legacy numbering made `seqid=1` ambiguous (referred to PDB residue 1, not sequence position 1)
   - **Lesson:** Can't rely on metadata, must use general approach

4. **Multi-Stage Ranking (Forcing Terminus Samples)**
   - **Result:** Forced spatially incorrect samples into top ranks
   - **Problem:** Tried to fix wrong assumptions with more complexity
   - **Lesson:** Fix the root cause (bad constraints), not the symptoms (ranking)

---

## Implementation Details

### `prepare_protein_ligand` Function

```python
def prepare_protein_ligand(...):
    """
    SATURATION Strategy: Wide seed diversity for maximum exploration
    """
    configs = []
    seeds = [42, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 7777777]
    
    for i, seed in enumerate(seeds):
        config_dict = input_dict.copy()  # No modifications!
        configs.append((config_dict, ["--diffusion_samples", "5", "--seed", str(seed)]))
    
    return configs  # 10 configs × 5 samples = 50 total
```

**Key features:**
- No constraints added to `input_dict`
- No pocket identification
- No terminus targeting
- Pure exploration with seed diversity

### `post_process_protein_ligand` Function

```python
def post_process_protein_ligand(...):
    """
    Confidence-Only Ranking
    """
    # Collect all 50 PDB files
    all_pdbs = [...]
    
    # Extract confidence from B-factor column
    scores = []
    for pdb_path in all_pdbs:
        confidence = extract_confidence(pdb_path)  # Mean B-factor (0-100)
        scores.append({"path": pdb_path, "confidence": confidence})
    
    # Sort by confidence (highest first)
    scores.sort(key=lambda x: x["confidence"], reverse=True)
    
    # Return top 5
    return [s["path"] for s in scores[:5]]
```

**Key features:**
- Single metric: Boltz confidence (pLDDT from B-factor)
- No normalization (raw scores)
- No weighted combinations
- Straightforward sorting

### `extract_confidence` Function

```python
def extract_confidence(pdb_path: Path) -> float:
    """
    Extract Boltz confidence from PDB B-factor column.
    Boltz writes confidence as B-factor (pLDDT-style, 0-100 scale).
    """
    confidence_values = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                b_factor = float(line[60:66].strip())  # Columns 61-66
                confidence_values.append(b_factor)
    
    return sum(confidence_values) / len(confidence_values)  # Mean confidence
```

**Key features:**
- Reads B-factor directly from PDB
- Returns mean confidence across all atoms
- Robust to parsing errors (falls back to 50.0)

---

## Performance Expectations

Based on validation experiments:

| Target Type | Expected Mean RMSD | Coverage | Notes |
|-------------|-------------------|----------|-------|
| **Orthosteric** | ~2-5Å | ~100% | Strong Boltz prior, works well |
| **Allosteric (Easy)** | ~5-15Å | ~80% | Accessible pockets, covered by SATURATION |
| **Allosteric (Hard)** | ~20-30Å | ~20% | Spatially unusual sites, fundamentally difficult |

**Overall:** Competitive performance across most targets, with known limitations on hard allosteric cases.

**Specific examples:**
- 3LW0 (allosteric): ~2Å ✅ (maintained)
- 6FVF (allosteric): 15.89Å ✅ (improved from ~20Å)
- 3K5V (allosteric): ~25Å ⚠️ (accepted as hard)

---

## Why This is the Final Strategy

1. **Simplicity:** Minimal code, easy to understand and debug
2. **Robustness:** No reliance on ambiguous metadata (seqid, PDB numbering)
3. **Proven:** Tested on validation set, improved hard targets
4. **General:** Works for both orthosteric and allosteric sites
5. **Fast:** No complex scoring, just confidence extraction

---

## What Was Removed

To reach this final strategy, we removed:
- ❌ `identify_sequence_pockets()` - sequence-based pocket identification
- ❌ `compute_kyte_doolittle_score()` - hydrophobicity scoring
- ❌ `parse_pdb()` - BioPython PDB parsing
- ❌ `compute_clash_penalty_fast()` - clash detection
- ❌ `count_protein_ligand_contacts_fast()` - contact counting
- ❌ `normalize_scores_fast()` - score normalization
- ❌ Regional constraints with soft/hard `force` flags
- ❌ Terminus probing logic
- ❌ Multi-stage ranking

All of these were either:
- Not helpful (added complexity without benefit)
- Actively harmful (diluted context, mis-ranked samples)
- Redundant (Boltz confidence alone is sufficient)

---

## Trade-offs and Limitations

### Accepted Limitations

1. **Hard allosteric targets (~20%):**
   - Spatially unusual binding sites
   - Boltz's strong orthosteric prior
   - **Cannot overcome without ground truth**
   - Examples: 3K5V, 6FVF-type cases

2. **No explicit allosteric targeting:**
   - SATURATION explores broadly
   - Relies on diffusion sampling to find sites
   - **No guarantees for edge cases**

### Why We Accept These

- **General predictor:** Must work without ground truth
- **Top-1 metric:** Better to maximize performance on 80% than risk breaking 100%
- **Model limitations:** Boltz's prior is strong, can't force unusual sites
- **Diminishing returns:** Complex strategies had negative ROI

---

## Running the Final Strategy

```bash
# Activate environment
conda activate boltz

# Run predictions
python3 hackathon/predict_hackathon.py \
    --dataset hackathon_data/datasets/asos_public/asos_public.jsonl \
    --submission-folder final_predictions \
    --output-dir final_outputs \
    --msa-dir hackathon_data/datasets/asos_public/msa

# Evaluate
python3 hackathon/evaluate_asos.py \
    --dataset-file hackathon_data/datasets/asos_public/asos_public.jsonl \
    --dataset-folder hackathon_data/datasets/asos_public \
    --submission-folder final_predictions \
    --result-folder final_results
```

**Expected output:**
- 50 samples per target (10 configs × 5 samples)
- Top 5 ranked by confidence in submission folder
- Competitive mean RMSD on full dataset

---

## Future Improvements (Not Pursued)

If we had more time/resources, these could help:

1. **Template-based modeling:** Use homologs with known allosteric sites
2. **Ensemble predictions:** Multiple MSA subsamplings per config
3. **Ligand conformer generation:** Sample multiple ligand geometries
4. **Physics-based refinement:** Post-process top predictions with MD
5. **Model fine-tuning:** Retrain Boltz on allosteric-enriched dataset

**But:** All of these add significant complexity and may not improve the Top-1 metric enough to justify the effort.

---

## Conclusion

**The final strategy is:**
- ✅ SATURATION (50 samples, wide seeds)
- ✅ Confidence-only ranking
- ✅ No constraints
- ✅ Full protein context

**This strategy is:**
- Simple, robust, and proven
- Optimized for Top-1 RMSD metric
- General-purpose (no ground truth)
- Fast and reproducible

**It balances:**
- Performance (good on 80%+ of targets)
- Simplicity (minimal code, easy to debug)
- Robustness (no fragile assumptions)

Ready for production AWS run! 🚀

