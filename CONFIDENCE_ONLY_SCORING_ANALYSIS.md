# Confidence-Only Scoring Analysis

## 🎯 Goal
Fix the re-ranking problem where model_0 (19.61Å) is selected over better predictions (15.68Å) for 6FVF.

---

## 🔍 Root Cause Analysis

### Problem: Normalization Destroyed Useful Signal

**For 6FVF_ALLOSTERIC_503:**

```
                 Confidence  Clashes  Contacts   RMSD
model_0 (Top-1):    97.48     2866      198    19.61Å  ← WRONG!
model_1:            97.55     2866      183    15.68Å  ← Should be Top-1!
model_2:            97.66     2866      179    15.92Å
model_3:            97.55     2866      177    15.82Å
model_4:            97.60     2866      176    15.82Å
```

**Key Observations:**
1. **Confidence:** All ~97.5, but small differences exist (97.48 → 97.66)
2. **Clashes:** ALL IDENTICAL (2866) → normalization makes all = 0.5
3. **Contacts:** 176-198 → becomes the ONLY differentiator after normalization

**What Happened with Hybrid Scoring:**
```python
hybrid = 0.65 * confidence_norm - 0.25 * clash_norm + 0.10 * contacts_norm

# After normalization:
# - confidence_norm ≈ all similar (tiny range)
# - clash_norm = 0.5 for all (identical clashes)
# - contacts_norm = varies (0.0 → 1.0)

# Formula becomes:
hybrid ≈ -0.125 + 0.10 * contacts_norm  # Ranks purely by contacts!
```

**model_0 won because it had MOST contacts (198), but more contacts ≠ better binding!**
- model_0 (wrong pocket): 198 contacts (surface exposed → more protein nearby)
- model_1 (right pocket): 183 contacts (correct binding site)

---

## 🔧 Solution: Confidence-Only Scoring

### Implementation

```python
def post_process_protein_ligand(...):
    """
    Sort by raw Boltz confidence (no normalization)
    Trust the model's subtle preferences
    """
    scores = []
    for pdb_path in all_pdbs:
        boltz_conf = extract_confidence(pdb_path)  # Raw value!
        scores.append({"path": pdb_path, "confidence": boltz_conf})
    
    # Sort by raw confidence (highest first)
    scores.sort(key=lambda x: x["confidence"], reverse=True)
    
    return [s["path"] for s in scores[:5]]
```

### Why This Works (For Some Cases)

**Theory:** Boltz's confidence captures subtle structural quality differences
- Small differences (97.48 vs 97.66) are meaningful
- Don't need physics metrics if model already learned them
- Simpler = fewer failure modes

---

## 📊 Results

### Test Case: 6FVF (test_saturation predictions)

**Simulated Confidence-Only Ranking:**
```
OLD (Hybrid):
Top-1: model_0 → 19.61Å

NEW (Confidence-Only):
Rank 1: model_2 (conf=97.66) → 15.92Å  ← NEW TOP-1
Rank 2: model_4 (conf=97.60) → 15.82Å
Rank 3: model_3 (conf=97.55) → 15.82Å
Rank 4: model_1 (conf=97.55) → 15.68Å
Rank 5: model_0 (conf=97.48) → 19.61Å  ← Now ranks LAST!
```

**Result:** 19.61Å → 15.92Å (**3.7Å improvement!** ✅)

---

### Validation Set: validation_10 (re-ranked actual predictions)

**Results:**
```
Datapoint                 Type         OLD Top-1    NEW Top-1    Change    
--------------------------------------------------------------------------------
1JQH_ORTHOSTERIC_ANP      orthosteric        1.65Å        1.56Å  →  -0.09Å
2E9N_ORTHOSTERIC_76A      orthosteric        0.99Å        1.19Å  ⚠️  +0.20Å
3K5V_ALLOSTERIC_STJ       allosteric        25.52Å       25.60Å  →  +0.08Å
3LW0_ALLOSTERIC_CCX       allosteric         2.34Å        2.46Å  ⚠️  +0.12Å
5MO4_ALLOSTERIC_AY7       allosteric         1.24Å        0.94Å  ✅  -0.31Å
```

**Overall:**
- OLD Mean Top-1: 6.35Å
- NEW Mean Top-1: 6.35Å
- Change: **+0.00Å** (neutral)

**By Type:**
- Orthosteric: 1.32Å → 1.37Å (+0.05Å, slightly worse)
- Allosteric: 9.70Å → 9.67Å (-0.03Å, slightly better)

---

## 🤔 Why Mixed Results?

### When Confidence-Only Works:
✅ **6FVF:** All models have similar confidence, but ranking correlates with RMSD
✅ **5MO4_ALLOSTERIC:** 1.24Å → 0.94Å (improved)

### When Confidence-Only Fails:
❌ **2E9N_ORTHOSTERIC:** 0.99Å → 1.19Å (worse)
❌ **3LW0_ALLOSTERIC:** 2.34Å → 2.46Å (worse)

### Hypothesis:

**The problem:** Confidence doesn't always correlate with RMSD!

**Possible reasons:**
1. **Model was trained on orthosteric data** → confidence biased toward orthosteric
2. **Wrong pocket can have high confidence** → model thinks it's right even when it's not
3. **Confidence measures structural quality**, not spatial correctness
   - A ligand in the wrong pocket can be "well-formed" (high conf)
   - A ligand in the right pocket can have "steric stress" (lower conf)

**For 6FVF specifically:** The confidence happened to correlate with RMSD, but this was somewhat lucky!

---

## 🎯 Key Insight: The Fundamental Problem

**The re-ranking problem is hard because:**

1. **We can't tell which pocket is correct** without ground truth
2. **Boltz confidence doesn't distinguish pockets** (measures quality, not location)
3. **Physics metrics don't help** when all predictions are reasonable structures
4. **The best sample exists in Top-5** but picking it as Top-1 is the challenge

**For 6FVF:**
- Top-5 min RMSD: 15.68Å (good!)
- Top-1 RMSD: 19.61Å (bad)
- **The ranking problem cost us 3.9Å**

---

## 💡 Next Steps

### Option 1: Keep Confidence-Only (Current Choice)
- **Pro:** Simple, works for some cases (6FVF)
- **Con:** Neutral overall, doesn't systematically help

### Option 2: Ensemble Approach
- Generate many diverse samples (SATURATION ✅)
- Return ALL top samples to evaluation
- Let Top-5 min RMSD metric do the ranking
- **Problem:** We're judged on Top-1, not Top-5!

### Option 3: Better Diversity Strategy
- The real issue might be **lack of diversity** (3K5V: all 25Å!)
- If all 50 samples are in the wrong pocket, no ranking helps
- Need to make Boltz explore multiple pockets

### Option 4: Accept Limitations
- **For orthosteric:** Boltz works well (1.37Å average)
- **For allosteric:** Fundamentally hard (9.67Å average)
- Focus on improving diversity rather than re-ranking

---

## 📝 Current Status

**Code:** Updated to confidence-only scoring ✅

**Performance:**
- **Best case:** 6FVF improved by 3.7Å
- **Average case:** Neutral (6.35Å → 6.35Å)
- **Allosteric:** Still struggling (9.67Å)

**Main Blocker:** 3K5V has ZERO diversity (all predictions ~25Å)
- No amount of re-ranking will fix this
- Need to solve the diversity problem first

**Recommendation:** Keep confidence-only scoring, but focus effort on:
1. **Hard constraints** (`force: true`) for terminus probing
2. **Pocket ensemble** strategies (multiple explicit pocket targets)
3. **Understanding why 3K5V is stuck** (is N-terminus binding too rare in training data?)

---

## 🔬 Conclusion

**Confidence-only scoring is a reasonable default**, but it won't solve the fundamental allosteric problem. The main challenge is **generating diverse enough predictions** that cover both orthosteric and allosteric sites. Once we have that diversity, simple confidence-based ranking should work reasonably well (as it did for 6FVF).

**Current bottleneck:** Generation diversity, not re-ranking strategy.

