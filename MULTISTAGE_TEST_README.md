# Multi-Stage Ranking Test for 3K5V

## 🎯 Quick Test

Just run this:
```bash
./run_3k5v_multistage.sh
```

Or manually:
```bash
python hackathon/predict_hackathon.py \
  --input-jsonl test_3k5v_only.jsonl \
  --msa-dir hackathon_data/datasets/asos_public/msa \
  --submission-dir multistage_3k5v_predictions \
  --intermediate-dir multistage_3k5v_intermediate \
  --result-folder multistage_3k5v_results
```

---

## 🔍 What to Look For

### **In the Console Output:**

**1. Check if terminus samples exist:**
```
🔍 TERMINUS SAMPLES (Config 8 & 9): 10 total     ← Should see 10 samples!
Rank  Model                              Confidence  
------------------------------------------------------------
1     config_8_model_2.pdb               94.23
2     config_8_model_0.pdb               94.10
...
```

**If you see 10 samples → They were filtered out! (Hypothesis confirmed)**  
**If you see 0 samples → Generation failed (Hypothesis wrong)**

**2. Check the final ranking:**
```
✅ FINAL TOP-5:
  1.         Confidence: 96.52 - config_0_model_3.pdb
  2. 🎯 TERMINUS Confidence: 94.23 - config_8_model_2.pdb  ← Forced in!
  3.         Confidence: 96.38 - config_1_model_1.pdb
  4. 🎯 TERMINUS Confidence: 93.85 - config_9_model_1.pdb  ← Forced in!
  5.         Confidence: 96.21 - config_2_model_0.pdb
```

**Look for the 🎯 TERMINUS markers** - these are now forced into Top-5!

---

## 📊 Check the Results

```bash
cat multistage_3k5v_results/combined_results.csv
```

**Compare:**
```
Previous (confidence-only):
  Top-1 RMSD: 25.04Å
  Top-5 min RMSD: 25.04Å

New (multi-stage):
  Top-1 RMSD: ???Å
  Top-5 min RMSD: ???Å
```

---

## ✅ Success Criteria

### **Excellent:**
- Top-1 RMSD < 10Å (terminus sample becomes best)
- Shows the right answer was there, just filtered out!

### **Good:**
- Top-5 min RMSD < 10Å (even if not Top-1)
- Proves terminus samples are correct, just need better ranking

### **Okay:**
- Terminus samples exist but RMSD still high (> 20Å)
- At least we know they were generated, just not good quality

### **Bad:**
- No terminus samples found (0 total)
- Hard constraints didn't work at all

---

## 🤔 What Each Result Means

| Result | Meaning | Next Action |
|--------|---------|-------------|
| Terminus samples exist, RMSD < 10Å | ✅ **Success!** Filtering was the problem | Use multi-stage in production |
| Terminus samples exist, RMSD 10-20Å | ⚠️ Partially correct | Try tighter constraints (max_distance: 4.0) |
| Terminus samples exist, RMSD > 20Å | ❌ Wrong pocket or bad geometry | Try different constraint strategy |
| No terminus samples (0 total) | ❌ Generation failed | Check Boltz logs, try force=True with tighter distance |

---

## 📝 Quick Summary

**Problem:** 3K5V ligand should bind at residue 1 (N-terminus), but all predictions are ~27Å away

**Hypothesis:** N-terminus samples were generated but filtered out by confidence ranking

**Test:** Force terminus samples into Top-5 and check their RMSD

**Expected:** Terminus samples exist with ~94 confidence (vs ~96 for others) and have low RMSD (< 10Å)

**If true:** We found the answer! Just needed better ranking, not better sampling.

---

**Now run `./run_3k5v_multistage.sh` and let's see what happens!** 🚀

