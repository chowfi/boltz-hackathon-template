# Debug Scoring Plan

## 🎯 Goal
Understand WHY 6FVF picks 19.61Å prediction instead of 15.82Å

## 📊 What We Know

**6FVF Top-5 RMSDs (from evaluation):**
```
model_0: 19.61Å  ← Selected as Top-1 (WRONG!)
model_1: 15.92Å
model_2: 15.82Å  ← Should be selected!
model_3: 15.82Å
model_4: 16.14Å
```

**Scoring Formula:**
```python
hybrid = 0.65 * confidence_norm - 0.25 * clash_norm + 0.10 * contacts_norm
```

## 🔍 Hypothesis

**Why model_0 (19.61Å) might score higher:**

### Hypothesis A: Higher Confidence
```
model_0 (19.61Å): confidence = 0.85 (high!)
model_2 (15.82Å): confidence = 0.45 (low!)

→ 0.65 * 0.85 = 0.55 (model_0)
→ 0.65 * 0.45 = 0.29 (model_2)

Even if model_2 has better geometry, the confidence gap is huge!
```

### Hypothesis B: Lower Clashes
```
model_0 (19.61Å): clashes = 5 (clean)
model_2 (15.82Å): clashes = 15 (some steric issues)

→ Model_0 looks "cleaner" but is in wrong pocket
→ Model_2 is correct pocket but tighter binding → more clashes
```

### Hypothesis C: Similar Contacts
```
Both models might have similar contact counts
→ Contacts term doesn't differentiate them
→ Decision comes down to confidence vs clashes
```

## 📝 What the Debug Output Will Show

With the new debug code, we'll see:
```
🔍 DETAILED SCORING (Top 10):
Rank  Model                           Conf    Clash   Contact  Hybrid  
-------------------------------------------------------------------------------------
1     config_X_model_Y.pdb           0.850     5.0      25      0.8234  ← Why this wins?
2     config_A_model_B.pdb           0.450    15.0      28      0.5123  ← Why this loses?
...
```

## 🎯 Next Steps Based on Results

### If Hypothesis A is True (Confidence Gap):
**Solution:** Reduce confidence weight, increase physics weight
```python
# Current:
0.65 * confidence - 0.25 * clash + 0.10 * contacts

# Try:
0.50 * confidence - 0.35 * clash + 0.15 * contacts  # More balanced
```

### If Hypothesis B is True (Clash Penalty):
**Solution:** Reduce clash penalty (tight binding ≠ bad)
```python
# Current:
-0.25 * clash

# Try:
-0.15 * clash  # Less aggressive penalty
```

### If Hypothesis C is True (Contacts Don't Help):
**Solution:** Increase contact weight
```python
# Current:
0.10 * contacts

# Try:
0.20 * contacts  # Reward binding more
```

## 🚀 Testing Plan

### Step 1: Run with Debug Output
```bash
python hackathon/predict_hackathon.py \
  --input-jsonl test_hard_allosteric.jsonl \
  --msa-dir hackathon_data/datasets/asos_public/msa \
  --submission-dir debug_scoring_predictions \
  --intermediate-dir debug_scoring_intermediate \
  --result-folder debug_scoring_results
```

### Step 2: Analyze Debug Output
Look for patterns in the scoring table:
- What's the confidence of rank 1 vs rank 2-5?
- What's the clash count of rank 1 vs rank 2-5?
- What's the contact count of rank 1 vs rank 2-5?

### Step 3: Adjust Weights
Based on findings, modify scoring formula

### Step 4: Retest
Run again with new weights and check if Top-1 improves

## 📊 Expected Outcome

**Best case:** We find the issue and fix 6FVF
- 19.61Å → 15.82Å (3.8Å improvement!)
- Same fix might help 1PZP and 3F9N (re-ranking failures)

**Realistic case:** We understand the tradeoff
- High confidence but wrong pocket vs low confidence but right pocket
- Helps us make informed decision on weight tuning

**Worst case:** The "wrong" sample has better scores on ALL metrics
- Then the RMSD evaluation is misleading
- Or we need different scoring approach entirely

---

**This is our best shot at quick improvement - debugging is key!** 🔍

