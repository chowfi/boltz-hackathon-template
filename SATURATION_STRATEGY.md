# SATURATION + Multi-Scoring Ensemble Strategy

## 🎯 Goal
Improve performance on hard allosteric cases (3K5V, 6FVF) while maintaining excellent orthosteric performance.

**Metric:** Lowest mean RMSD of **Top-1** predictions on internal test set.

---

## 📊 Current Performance (Wide Seed Diversity - 20 samples)

### ✅ EXCELLENT (8/10 targets, 80%)
**Orthosteric (all working):**
- 1AXB: 1.6Å
- 1JQH: 2.0Å
- 2AX9: 1.8Å
- 2E9N: 1.9Å
- 2PIO: 1.6Å
- 3PE2: 2.0Å

**Allosteric (2/4 working):**
- 3LW0: 2.3Å ✅
- 4YPQ: 2.4Å ✅

### ❌ FAILING (2/10 targets, 20%)
**Hard Allosteric:**
- 3K5V: 25.3Å ❌ (N-terminus site, seqid=1)
- 6FVF: 16.1Å ❌ (Unusual spatial location)

---

## 🔥 NEW STRATEGY: SATURATION + Multi-Scoring Ensemble

Based on hackathon hints image, implementing two key strategies:

### 1. **SATURATION** (Generation Phase)
**Goal:** Generate ENOUGH samples to hit rare allosteric sites

**Implementation:**
- **10 configs** × 5 samples = **50 total samples** (was 20)
- **Extreme seed spacing:** 42 → 1000 → 5000 → 10000 → 25000 → 50000 → 100000 → 500000 → 999999 → 7777777
- **Full protein context** (no pocket scanning - learned this hurts performance)

**Rationale:**
- Hard allosteric sites are RARE in diffusion sampling space
- Need to saturate the conformational space to find them
- 2.5x more samples → 2.5x better chance of hitting the site

---

### 2. **Multi-Scoring Ensemble** (Re-ranking Phase)
**Goal:** Correctly pick the BEST prediction as Top-1

**Challenge:**
- Even if we generate a perfect allosteric prediction, our scoring might rank it 20th!
- Need intelligent re-ranking for **Top-1 RMSD** metric

**Implementation:**

#### Three Scoring Schemes:

1. **ORTHOSTERIC-optimized** (trust Boltz model)
   ```python
   score = 0.80 * confidence - 0.15 * clashes + 0.05 * contacts
   ```
   - High weight on model confidence (trained on orthosteric data)
   - Model knows orthosteric sites well

2. **ALLOSTERIC-optimized** (trust physics)
   ```python
   score = 0.30 * confidence - 0.40 * clashes + 0.30 * contacts
   ```
   - LOW weight on model confidence (allosteric is hard for model)
   - HIGH weight on physics (clashes, contacts)
   - Good physics = good binding, even if model unsure

3. **BALANCED** (hybrid)
   ```python
   score = 0.60 * confidence - 0.25 * clashes + 0.15 * contacts
   ```
   - Middle ground between orthosteric and allosteric

#### Ensemble Voting Logic:

```
Get top-1 from each scheme:
  - top1_orthosteric
  - top1_allosteric
  - top1_balanced

IF all 3 agree:
  → High confidence! Pick balanced
  
ELIF orthosteric + balanced agree:
  → Likely orthosteric site, pick orthosteric
  
ELIF allosteric + balanced agree:
  → Likely allosteric site, pick allosteric
  
ELSE (all disagree):
  → Use meta-ensemble (average of all 3 schemes)
```

#### Top-5 Diversification:
1. Top-1: Ensemble winner
2-4: Best from each scheme (for diversity)
5: Next best from balanced scheme

---

## 📂 Test Files

### `test_hard_allosteric.jsonl`
Contains ONLY the 2 hard allosteric cases:
- 3K5V_ALLOSTERIC_STJ
- 6FVF_ALLOSTERIC_503

### Run Command:
```bash
python hackathon/predict_hackathon.py \
  --input-jsonl test_hard_allosteric.jsonl \
  --msa-dir hackathon_data/datasets/asos_public/msa \
  --submission-dir test_saturation_predictions \
  --intermediate-dir test_saturation_intermediate \
  --result-folder test_saturation_results
```

---

## 🎲 Expected Outcomes

### Best Case Scenario:
- **3K5V:** 25.3Å → <5Å (if we hit the N-terminus site)
- **6FVF:** 16.1Å → <5Å (if we hit the unusual site)
- **Overall:** 80% success → 100% success!

### Realistic Scenario:
- **3K5V:** 25.3Å → 10-15Å (partial improvement)
- **6FVF:** 16.1Å → 8-12Å (partial improvement)
- **Overall:** Reduce hard failures by 30-50%

### Worst Case:
- No improvement (these sites are truly impossible)
- But: Saturation shouldn't hurt existing good predictions
- Ensemble voting is robust and shouldn't degrade performance

---

## 🧪 Why This Should Work

### SATURATION:
1. ✅ Hackathon hint explicitly mentions "Saturation"
2. ✅ Math: 50 samples vs 20 → 2.5x coverage of conformational space
3. ✅ Orthogonal seeds (42, 1000, 5k, 10k, ..., 7.7M) explore different noise trajectories
4. ✅ Even if success rate is 2% for hard sites, 50 samples gives ~1 good sample

### Multi-Scoring Ensemble:
1. ✅ Hackathon hint: "Using Confidence metrics" + "Re-ranking"
2. ✅ Robust to both orthosteric AND allosteric (no guessing needed)
3. ✅ Consensus voting → more stable Top-1 selection
4. ✅ Diversified Top-5 → safety net if Top-1 is wrong

---

## 📈 Success Metrics

**After running test:**
1. Check Top-1 RMSD for 3K5V (target: <10Å, stretch: <5Å)
2. Check Top-1 RMSD for 6FVF (target: <10Å, stretch: <5Å)
3. Verify ensemble voting output (which scheme won?)
4. Check Top-5 diversity (are all 3 schemes represented?)

**If it works:**
- Run on full validation set
- Should maintain 80% performance + improve on 20% hard cases

**If it doesn't work:**
- Analyze which seeds produced best samples
- Adjust ensemble voting weights
- Consider even more extreme saturation (100 samples?)

---

## 🚀 Next Steps After Testing

1. ✅ If successful: Apply to full dataset
2. ⚠️ If partial success: Tune ensemble weights
3. ❌ If no improvement: These cases may be inherently unsolvable
   - Investigate ligand placement constraints
   - Consider alternative sampling strategies (MD, conformer generation)

