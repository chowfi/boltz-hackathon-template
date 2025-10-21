# Multi-Stage Ranking Strategy

## 🎯 Problem Statement

**3K5V Hard Constraints Failed:**
- Used `force: True` with `max_distance: 6.0Å` for N-terminus (residues 1-20)
- Result: All Top-5 models still ~27Å away from residue 1
- **Hypothesis:** N-terminus samples WERE generated but filtered out by confidence ranking

---

## 🔍 Root Cause

**Confidence-Only Ranking Creates a Bias:**

```
Config 0-7 (40 samples): Natural exploration → High confidence (~96-97)
Config 8 (5 samples): N-terminus forced → Low confidence (~92-94?)
Config 9 (5 samples): C-terminus forced → Low confidence (~92-94?)

Confidence-only ranking → Top-5 all from Configs 0-7 → Zero diversity!
```

**Why terminus samples have lower confidence:**
- Forced placement in unusual regions
- Boltz's prior says "N-terminus binding = unlikely"
- Unnatural geometry → lower model confidence
- But might still be CORRECT spatially!

---

## 🚀 Solution: Multi-Stage Ranking

### **Strategy:**
1. **Score all 50 samples** by Boltz confidence
2. **Separate into pools:**
   - Terminus pool: Config 8 & 9 (10 samples)
   - Other pool: Config 0-7 (40 samples)
3. **Ensure diversity:**
   - Take **best 2 from terminus pool**
   - Take **best 3 from other pool**
   - Total: 5 samples with guaranteed diversity
4. **Re-rank by confidence** for Top-1 selection

---

## 📊 Expected Behavior

### **For 3K5V (N-terminus binding):**

**Before (Confidence-Only):**
```
Top-5: All from Configs 0-7, all ~27Å from residue 1
Top-1 RMSD: 25Å
```

**After (Multi-Stage):**
```
Top-5: 2 from Configs 8/9 (terminus), 3 from Configs 0-7 (other)
  - Rank 1: Config_0 sample (confidence 96.5, RMSD ~25Å)
  - Rank 2: Config_8 sample (confidence 94.2, RMSD ~3Å?) 🎯
  - Rank 3: Config_1 sample (confidence 96.2, RMSD ~25Å)
  - Rank 4: Config_9 sample (confidence 93.8, RMSD ~20Å)
  - Rank 5: Config_2 sample (confidence 96.0, RMSD ~25Å)

Top-1 RMSD: Still 25Å (if confidence dominates)
BUT: Top-5 min RMSD: ~3Å (massive improvement!)
```

**Key Insight:** Even if Top-1 doesn't improve (due to lower confidence), the **correct sample is now in Top-5**!

---

## 🔧 Implementation Details

### **Code Changes:**

```python
def post_process_protein_ligand(...):
    # 1. Score all samples
    for pdb_path in all_pdbs:
        confidence = extract_confidence(pdb_path)
        is_terminus = "config_8_" in pdb_path.name or "config_9_" in pdb_path.name
        scores.append({"path": pdb_path, "confidence": confidence, "is_terminus": is_terminus})
    
    # 2. Separate pools
    terminus_samples = [s for s in scores if s["is_terminus"]]
    other_samples = [s for s in scores if not s["is_terminus"]]
    
    # 3. Sort each pool
    terminus_samples.sort(key=lambda x: x["confidence"], reverse=True)
    other_samples.sort(key=lambda x: x["confidence"], reverse=True)
    
    # 4. Ensure diversity: 2 terminus + 3 other
    final_ranking = terminus_samples[:2] + other_samples[:3]
    
    # 5. Re-sort by confidence for Top-1
    final_ranking.sort(key=lambda x: x["confidence"], reverse=True)
    
    return [s["path"] for s in final_ranking[:5]]
```

### **Debug Output:**

The new implementation will print:
```
🔍 TERMINUS SAMPLES (Config 8 & 9): 10 total
Rank  Model                              Confidence  
------------------------------------------------------------
1     config_8_model_2.pdb               94.23
2     config_8_model_0.pdb               94.10
3     config_9_model_1.pdb               93.85
...

🔍 TOP NON-TERMINUS SAMPLES: 40 total
Rank  Model                              Confidence  
------------------------------------------------------------
1     config_0_model_3.pdb               96.52
2     config_1_model_1.pdb               96.38
3     config_2_model_0.pdb               96.21
...

✅ Multi-Stage: 2 terminus + 3 other samples

✅ FINAL TOP-5:
  1.         Confidence: 96.52 - config_0_model_3.pdb
  2. 🎯 TERMINUS Confidence: 94.23 - config_8_model_2.pdb
  3.         Confidence: 96.38 - config_1_model_1.pdb
  4. 🎯 TERMINUS Confidence: 93.85 - config_9_model_1.pdb
  5.         Confidence: 96.21 - config_2_model_0.pdb
```

This tells us:
- **Do terminus samples exist?** (Yes, 10 of them)
- **What's their confidence?** (94.23 vs 96.52 for others)
- **Are they being filtered?** (Yes, without multi-stage they'd be ranks 41-50)

---

## ⚖️ Tradeoffs

### **Pros:**
✅ **Guarantees diversity** - Always have 2 terminus samples in Top-5
✅ **Still uses confidence** - Not ignoring model's quality assessment
✅ **Low risk** - If terminus samples are terrible, they won't be Top-1
✅ **Diagnostic** - Shows us exactly what's happening

### **Cons:**
❌ **Might include bad samples** - If all terminus samples are terrible
❌ **Top-1 might not improve** - If terminus samples have lower confidence
❌ **Arbitrary 2:3 split** - Could be 1:4 or 3:2 instead

---

## 📈 Success Metrics

### **Minimum Success:**
- **Terminus samples exist and are visible** in debug output
- **Top-5 includes at least 2 terminus samples**

### **Good Success:**
- **3K5V Top-5 min RMSD < 10Å** (currently 25Å)
- **At least one terminus sample has RMSD < 15Å**

### **Excellent Success:**
- **3K5V Top-1 RMSD < 10Å** (terminus sample becomes Top-1)
- **6FVF maintains improvement** (15.89Å or better)

---

## 🧪 Testing

```bash
python hackathon/predict_hackathon.py \
  --input-jsonl test_hard_allosteric.jsonl \
  --msa-dir hackathon_data/datasets/asos_public/msa \
  --submission-dir multistage_predictions \
  --intermediate-dir multistage_intermediate \
  --result-folder multistage_results
```

**Look for in console output:**
```
🔍 TERMINUS SAMPLES (Config 8 & 9): 10 total
...
🎯 TERMINUS markers in FINAL TOP-5
```

---

## 🎓 What We're Testing

**Hypothesis:** N-terminus samples exist but are filtered out by confidence ranking

**Test:** Force them into Top-5 and measure RMSD

**Expected outcomes:**

| Result | Interpretation | Next Action |
|--------|----------------|-------------|
| Terminus samples have low conf (< 94) AND high RMSD (> 20Å) | Boltz can't place ligand at N-terminus properly | Try tighter constraints or accept limitation |
| Terminus samples have low conf (< 94) BUT low RMSD (< 10Å) | **SUCCESS!** Filtering was the problem | Use multi-stage ranking in production |
| Terminus samples have high conf (> 95) AND low RMSD (< 10Å) | Something else is wrong with ranking | Investigate further |
| No terminus samples found | Hard constraints didn't work at all | Check Boltz implementation |

---

## 💡 Key Insight

**This strategy doesn't just solve 3K5V - it solves the fundamental diversity problem:**

- **Orthosteric targets:** Naturally high confidence → dominate Top-5
- **Allosteric targets:** Lower confidence → filtered out
- **Solution:** Ensure both types are represented, then let confidence decide Top-1

**This is better than:**
- Pure confidence (current problem)
- Ignoring confidence (promotes bad samples)
- Complex scoring (hard to tune)

---

## 🚀 Ready to Test!

The code is updated. Run the command above and check the console output for:
1. Do we see terminus samples?
2. What are their confidence values?
3. Do they make it into Top-5 now?
4. What are the RMSDs?

**Let's find out if filtering was the problem!** 🔍

