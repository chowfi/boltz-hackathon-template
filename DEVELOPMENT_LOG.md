# Development Log: Boltz Hackathon Allosteric Predictor

## Overview

This document summarizes the iterative development process for improving protein-ligand binding prediction, with a focus on allosteric sites. The winning metric is **lowest mean RMSD of top-1 predictions** on the internal test set.

---

## Experiments & Key Learnings

### 1. **Data Analysis** 📊
- **Goal:** Understand what differentiates orthosteric vs allosteric binding sites
- **Findings:** 
  - Allosteric sites are more challenging (higher RMSD)
  - ~80% of allosteric targets perform reasonably well
  - ~20% are extremely difficult (3K5V, 6FVF, etc.)
- **Takeaway:** Need specialized strategies for allosteric sites

---

### 2. **Sequence-Based Pocket Identification** 🧬
- **Approach:** Identify likely allosteric pockets using sequence heuristics
  - N/C termini (flexible regions)
  - Charged residue clusters (regulatory sites)
  - High flexibility scores (mobile regions)
- **Result:** ✅ **Worked reasonably well**
- **Performance:** Better than baseline for some allosteric targets
- **Kept as foundation** for subsequent approaches

---

### 3. **Geometric Pocket Identification** 📐
- **Approach:** Use spatial features to identify pockets
  - Surface accessibility
  - Binding depth
  - Pocket volume
- **Result:** ⚠️ **Similar to sequence-based**
- **Finding:** Sequence + geometric combined had minimal improvement
- **Takeaway:** Sequence-based was sufficient

---

### 4. **Ensemble Scoring (Multi-Scoring)** 🎯
- **Approach:** Combine multiple metrics for ranking predictions
  - Boltz confidence
  - Clash score
  - Contact score
  - Surface contacts
  - Flexibility score
  - Binding depth
- **Result:** ❌ **Failed - made things worse**
- **Problem:** Complex scoring picked physically "better" but spatially wrong predictions
- **Example:** 6FVF picked 19.61Å instead of 15.82Å (both had similar physical scores)
- **Conclusion:** **Confidence-only ranking is more robust**

---

### 5. **Systematic Pocket Scanning** 🔍
- **Approach:** Divide protein into regions, sample each systematically
  - Config 0: Baseline
  - Configs 1-N: Target specific regions with constraints
- **Result:** ❌ **Failed catastrophically**
- **Problem:** Context dilution - breaking up the protein hurt performance
- **Example:** **3LW0 got significantly worse** when protein was divided into regions
- **Key Learning:** **Full protein context is critical! Don't break it up!**
- **Reverted** to simpler sequence-based approach

---

### 6. **SATURATION Strategy (Wide Seed Diversity)** 🌊
- **Approach:** Generate many samples with widely spaced random seeds
  - 10 configs × 5 samples = 50 total samples
  - Seeds: `[42, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 7777777]`
  - No constraints, pure exploration
- **Result:** ✅ **Major success for diversity!**
- **Finding:** Wide seed spacing creates genuine diversity in diffusion sampling
- **Example:** 6FVF generated samples ranging from 15.82Å to 27Å (good diversity!)
- **Problem:** Good samples were being mis-ranked by complex scoring
- **Conclusion:** SATURATION generates diversity, but needs better ranking

---

### 7. **Simplify Scoring (Confidence-Only)** 🎖️
- **Approach:** After SATURATION failures, simplified ranking
  - Sort purely by Boltz confidence (highest first)
  - Remove all complex hybrid scoring
- **Result:** ✅ **Worked well!**
- **Example:** 6FVF improved from 19.61Å → 15.89Å
- **Why it works:** Boltz confidence is well-calibrated, complex scoring introduced noise
- **Kept as final ranking strategy**

---

### 8. **Terminus Probing (N/C-terminus Constraints)** 🎯
- **Motivation:** Hard allosteric targets like 3K5V bind at extreme N-terminus
- **Approach:** 
  - 8 SATURATION configs (general exploration)
  - 1 N-terminus hard constraint (`force: true`, residues 1-20)
  - 1 C-terminus hard constraint (`force: true`, last 20 residues)
- **Result:** ⚠️ **Mixed results**
  - 6FVF: Maintained 15.89Å ✅
  - 3K5V: Got **worse** (25Å → 33Å) ❌
- **Problem:** PDB numbering issue discovered
  - `seqid=1` in ground truth refers to **PDB residue 1** (not sequence position 1)
  - For 3K5V: PDB starts at residue 243, so `seqid=1` is actually sequence position ~0 (doesn't exist!)
  - Constraints targeted wrong spatial location
- **Conclusion:** Can't rely on `seqid` metadata due to legacy PDB numbering

---

### 9. **Multi-Stage Ranking** 🏆
- **Approach:** Force terminus-constrained samples into Top-5
  - Guarantee at least 2 terminus samples
  - Fill rest with highest confidence
- **Result:** ❌ **Made 3K5V worse**
- **Problem:** Forced spatially incorrect samples (due to PDB numbering issue) into top ranks
- **Conclusion:** Don't force samples that might be wrong

---

## Final Strategy: SATURATION + Confidence-Only

**What worked best:**

```python
# 10 configs with widely spaced seeds
seeds = [42, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 7777777]

# Each config: 5 diffusion samples, no constraints
for seed in seeds:
    cli_args = ["--diffusion_samples", "5", "--seed", str(seed)]

# Ranking: Sort by Boltz confidence (highest first), return top 5
```

**Why this works:**
- ✅ **Maximum diversity** through wide seed spacing
- ✅ **Full protein context** (no region constraints)
- ✅ **Simple, robust ranking** (confidence-only)
- ✅ **No assumptions** about binding location
- ✅ **No reliance on ambiguous metadata** (seqid)

---

## Key Insights

### What Helps:
1. ✅ **Wide seed diversity** (SATURATION) - generates genuinely different samples
2. ✅ **Full protein context** - never break up the protein into regions
3. ✅ **Simple confidence ranking** - more robust than complex hybrid scoring
4. ✅ **Many samples** (50 total) - coverage of conformational space

### What Hurts:
1. ❌ **Pocket scanning / region constraints** - dilutes context (3LW0 failure)
2. ❌ **Complex hybrid scoring** - introduces noise, mis-ranks good samples (6FVF failure)
3. ❌ **Terminus hard constraints** - PDB numbering issues make them target wrong locations (3K5V failure)
4. ❌ **Multi-stage ranking** - forces potentially wrong samples into top ranks

### Fundamental Limitations:
- **~20% of allosteric targets are extremely hard** (3K5V, 6FVF type)
- These have spatially unusual binding sites that Boltz's prior doesn't favor
- Without ground truth, can't overcome strong model bias
- **Accept this limitation** - focus on maximizing performance on the 80%

---

## Performance Summary

| Target | Type | Baseline | Final (SATURATION + Confidence) | Notes |
|--------|------|----------|--------------------------------|-------|
| 3LW0 | Allosteric | ~2Å | ~2Å | Maintained, pocket scanning hurt this |
| 6FVF | Allosteric | ~20Å | **15.89Å** | ✅ Improved via SATURATION + simple scoring |
| 3K5V | Allosteric | ~25Å | ~25Å | Hard case, PDB numbering issue |
| Orthosteric | Orthosteric | Good | Good | Maintained performance |

---

## Git Hash for Final Run

**Commit:** `6090a45efa1110838c4d72113b20b0c948f03479`

This commit represents the **SATURATION + Confidence-Only** strategy before terminus probing experiments.

---

## Lessons Learned

1. **Simple is better** - Complex scoring and multi-stage ranking introduced more problems than they solved
2. **Context is king** - Never dilute context by dividing the protein (pocket scanning failure)
3. **Diversity matters** - Wide seed spacing (10^6 - 10^7 range) creates genuine diversity
4. **Trust the model** - Boltz confidence is well-calibrated, better than hand-crafted metrics
5. **Metadata can be misleading** - Legacy PDB numbering issues make `seqid` unreliable
6. **Accept fundamental limits** - Some targets are too hard without ground truth

---

## Future Directions (Not Pursued)

- **Template-based modeling:** Use homologs with known allosteric sites
- **MSA subsampling:** Vary MSA to explore conformational diversity
- **Step scale tuning:** Lower values increase diversity (but can't overcome strong priors)
- **Potentials/steering:** Physical guidance during diffusion (unclear benefit)
- **Two-stage prediction:** Structure first, then docking (more complex pipeline)

These were deemed too complex or unlikely to significantly improve results given the fundamental limitation that Boltz's prior strongly favors orthosteric sites.

