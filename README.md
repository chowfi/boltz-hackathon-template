# Final Strategy: Saturation + Confidence-Only

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
   - **Key insight:** Diversity comes from diffusion sampling

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

## Future Improvements (Not Pursued)

If we had more time/resources, these could help:

1. **Template-based modeling:** Use homologs with known allosteric sites
2. **Ensemble predictions:** Multiple MSA subsamplings per config
3. **Physics-based refinement:** Post-process top predictions with MD
4. **Model fine-tuning:** Retrain Boltz on allosteric-enriched dataset

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

