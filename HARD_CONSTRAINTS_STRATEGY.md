# Hard Terminus Constraints Strategy

## 🎯 Objective
Force Boltz to explore N-terminus and C-terminus regions that it would naturally avoid, specifically targeting hard cases like 3K5V where the ligand binds at seqid=1 (extreme N-terminus).

---

## 📊 Strategy Overview

### **SATURATION + HARD Terminus Probing (50 total samples)**

```
Config 0-7 (40 samples): SATURATION
├─ Seeds: 42, 1000, 5000, 10000, 50000, 100000, 500000, 7777777
├─ Unconstrained (Boltz explores naturally)
└─ 80% of total samples

Config 8 (5 samples): N-terminus HARD probe
├─ Residues: 1-20
├─ force: True (HARD constraint)
├─ Seed: 111111
└─ 10% of total samples

Config 9 (5 samples): C-terminus HARD probe
├─ Residues: last 20 residues
├─ force: True (HARD constraint)
├─ Seed: 222222
└─ 10% of total samples
```

---

## 🔄 What Changed from Previous Approach

### **Before (Soft Constraints):**
```yaml
constraints:
  - pocket:
      binder: B
      contacts: [[A, 1], [A, 2], ..., [A, 20]]
      max_distance: 6.0
      force: false  # ← "Please try N-terminus, but you can ignore me"
```

**Result:** Boltz completely ignored the constraint for 3K5V (all predictions still ~25Å)

---

### **After (Hard Constraints):**
```yaml
constraints:
  - pocket:
      binder: B
      contacts: [[A, 1], [A, 2], ..., [A, 20]]
      max_distance: 6.0
      force: true  # ← "You MUST place ligand here, no exceptions!"
```

**Expected Result:** 5/50 samples GUARANTEED to be at N-terminus, regardless of Boltz's prior

---

## 🎯 Target Cases

### **Primary Target: 3K5V_ALLOSTERIC_STJ**
- **Problem:** Ligand binds at seqid=1 (extreme N-terminus)
- **Current Performance:** All 50 samples ~25Å (zero diversity, stuck in wrong pocket)
- **Expected with Hard Constraints:** 
  - 5 samples FORCED to N-terminus → likely 1-5Å if correct pocket!
  - 40 samples still in natural pocket (~25Å)
  - 5 samples at C-terminus (likely irrelevant)
  - **If N-terminus is correct, Top-1 RMSD: 25Å → ~2-5Å** 🎯

### **Secondary Targets: C-terminus binders**
- **2YHD, 2YLO, 5A1R:** High seqid values (possibly C-terminus due to legacy PDB numbering)
- **Expected:** If any bind at C-terminus, the 5 forced samples should capture them

---

## ⚖️ Risk Assessment

### **Risks:**

1. **Forced placements might be unnatural**
   - Ligand jammed into N/C-terminus → high clashes, low confidence
   - **Mitigation:** Confidence-only scoring will naturally deprioritize these if they're bad

2. **Wasted samples if pocket is wrong**
   - If 3K5V is NOT N-terminus, we wasted 5 samples
   - **Mitigation:** 40/50 samples (80%) still explore naturally, so minimal loss

3. **Might hurt orthosteric targets**
   - If a target is clearly orthosteric, the 10 forced terminus samples are wasted
   - **Mitigation:** Orthosteric ligands rarely bind at extreme termini, so forced samples will score low and be filtered out by confidence ranking

### **Expected Impact:**

| Scenario | Probability | Impact |
|----------|-------------|--------|
| 3K5V is N-terminus (seqid=1 suggests yes) | High | **Huge win** (25Å → ~2-5Å) |
| Some targets bind C-terminus | Medium | Moderate win (coverage) |
| Most targets are internal pockets | High | Neutral (bad samples filtered by confidence) |
| Hard constraints cause crashes | Low | Need fallback strategy |

**Overall:** Low risk, high potential reward for stuck cases.

---

## 📈 Success Metrics

### **Must Have (Minimum Success):**
- 3K5V: Top-1 RMSD < 10Å (currently 25Å)
- No regression on orthosteric targets (currently ~1.37Å average)

### **Nice to Have (Stretch Goals):**
- 3K5V: Top-1 RMSD < 5Å (near-native binding)
- Other terminus-binding allosteric targets improve
- Overall allosteric average: 9.67Å → < 8Å

### **Red Flags (Abort if):**
- Orthosteric targets degrade (> 2Å average)
- Boltz crashes due to impossible constraints
- All forced samples have extremely low confidence (< 50)

---

## 🧪 Testing Plan

### **Phase 1: Quick Test on Hard Cases**
```bash
python hackathon/predict_hackathon.py \
  --input-jsonl test_hard_allosteric.jsonl \
  --msa-dir hackathon_data/datasets/asos_public/msa \
  --submission-dir hard_constraints_predictions \
  --intermediate-dir hard_constraints_intermediate \
  --result-folder hard_constraints_results
```

**Check:**
- Did 5 samples actually go to N-terminus for 3K5V?
- What's their RMSD?
- What's their confidence?

### **Phase 2: Full Validation (if Phase 1 succeeds)**
```bash
python hackathon/predict_hackathon.py \
  --input-jsonl hackathon_data/datasets/asos_public/asos_public.jsonl \
  --msa-dir hackathon_data/datasets/asos_public/msa \
  --submission-dir hackathon_data/submission/asos_public \
  --intermediate-dir hackathon_data/intermediate_files/asos_public
```

**Evaluate full dataset, compare with baseline.**

---

## 🔧 Implementation Details

### **Code Changes (predict_hackathon.py):**

**Line 100:** 
```python
"force": True  # Changed from False
```

**Line 121:**
```python
"force": True  # Changed from False
```

**Line 60-69:** Updated docstring to reflect hard constraints

**That's it!** Just 2 lines changed from `False` to `True`.

---

## 📊 Expected Console Output

```
🔥 SATURATION + Terminus Probing for 3K5V_ALLOSTERIC_STJ
  Config 0: Saturation seed=     42 (5 samples)
  Config 1: Saturation seed=   1000 (5 samples)
  ...
  Config 7: Saturation seed=7777777 (5 samples)
  Config 8: N-terminus HARD probe (residues 1-20, force=True, 5 samples)  ← NEW!
  Config 9: C-terminus HARD probe (residues 274-293, force=True, 5 samples)  ← NEW!
✅ Generated 10 configs, 50 TOTAL SAMPLES
   - 8 SATURATION configs (seeds: 42 → 7,777,777, unconstrained)
   - 1 N-terminus HARD probe (force=True, 3K5V-targeting)  ← NEW!
   - 1 C-terminus HARD probe (force=True, coverage)  ← NEW!
```

---

## 🎓 Key Learnings from Previous Attempts

### **Why Soft Constraints Failed:**
1. Boltz's prior for N-terminus binding is VERY weak (rare in training data)
2. Soft constraints (`force: false`) can be completely overridden by the model
3. For 3K5V, the prior was so strong that 0/50 samples explored N-terminus

### **Why Hard Constraints Should Work:**
1. **Guaranteed exploration:** Can't be overridden
2. **Low risk:** Only 10% of samples affected
3. **Natural filtering:** Bad forced placements will have low confidence and be filtered out
4. **Addresses root cause:** 3K5V's problem is lack of diversity, not re-ranking

### **Why This is Better Than Other Approaches:**
- ❌ **Pocket scanning:** Divides context, hurts performance
- ❌ **Multi-scoring ensemble:** Picks physically "better" but spatially wrong predictions
- ❌ **Adaptive ranking:** Misclassifies targets
- ✅ **Hard terminus constraints:** Simple, targeted, low-risk

---

## 🚀 Next Actions

1. ✅ **Code updated** with hard constraints
2. **Run Phase 1 test** on 3K5V and 6FVF
3. **Analyze results:**
   - Check if N-terminus samples actually placed at N-terminus
   - Measure RMSD improvement
   - Verify confidence scores
4. **If successful:** Run full validation
5. **If failed:** Document failure mode and try alternative (e.g., even harder distance constraints)

---

## 📝 Conclusion

**Hard terminus constraints are a low-risk, high-reward strategy to solve the "stuck" problem for cases like 3K5V.** By forcing 10% of samples to explore extreme N/C termini, we guarantee diversity in regions that Boltz naturally avoids, while keeping 80% of samples for natural (orthosteric-biased) exploration. Combined with confidence-only scoring, this should improve performance on hard allosteric targets without hurting orthosteric performance.

**Let's test it!** 🚀

