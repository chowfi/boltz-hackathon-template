# Regional Sampling Strategy - Final Approach

## 🎯 Goal
Solve the PDB numbering/seqid issue by NOT relying on ground truth metadata at all!

---

## 📊 Strategy Overview

### **General Heuristics + Regional Sampling**

```
50 total samples divided into:
- 30% Pure exploration (SATURATION)
- 60% Systematic regional coverage
- 10% High diversity (edge cases)
```

### **Key Innovation:**
**Don't use `seqid` from ground truth!** Sample the entire protein systematically.

---

## 🔬 Implementation

### **Config Breakdown:**

```python
# Configs 0-2: SATURATION (15 samples)
Seeds: 42, 1000, 5000
Purpose: Pure exploration, no bias

# Configs 3-8: REGIONAL (30 samples)
Divide protein into 6 overlapping regions
Each region: ~50 residues or 25% of protein length
Soft constraints (force=False): Guide but don't force
Purpose: Systematic coverage of entire protein

# Config 9: HIGH DIVERSITY (5 samples)
Seed: 7777777
Purpose: Catch edge cases, unusual binding modes
```

### **For 3K5V (293 residues):**
```
Config 0: SATURATION seed=42
Config 1: SATURATION seed=1000
Config 2: SATURATION seed=5000
Config 3: REGION 1 residues 1-73 (soft)
Config 4: REGION 2 residues 49-122 (soft)
Config 5: REGION 3 residues 97-170 (soft)
Config 6: REGION 4 residues 145-218 (soft)
Config 7: REGION 5 residues 193-266 (soft)
Config 8: REGION 6 residues 241-293 (soft)
Config 9: HIGH DIVERSITY seed=7777777
```

**Note:** Regions overlap by ~25 residues for smooth coverage

---

## ✅ Why This Solves the PDB Numbering Problem

### **The Problem:**
```
Ground truth says: seqid=1
PDB file has: residues starting at 243
We interpreted: sequence position 1
Reality: Who knows! (ambiguous metadata)
```

### **Our Solution:**
```
Don't care about seqid!
Just sample: residues 1-73, 49-122, 97-170, ..., 241-293
One of these WILL overlap with the true binding site
Ranking picks the best → correct answer
```

### **Why It Works:**
- ✅ **No assumptions** about metadata
- ✅ **Systematic coverage** of entire protein
- ✅ **Soft constraints** let Boltz refine (not forced to wrong pockets)
- ✅ **Robust** to PDB numbering, missing residues, different constructs

---

## 📊 Expected Results

### **For 3K5V:**
- **Old (terminus probing):** 33Å (wrong target due to PDB numbering)
- **Expected (regional):** < 20Å (one region will overlap true site)
- **Best case:** < 10Å (if we hit the exact region)

### **For 6FVF:**
- **Old (confidence-only):** 15.89Å (already good!)
- **Expected (regional):** ~15Å or better (maintain or improve)

### **For Orthosteric Targets:**
- **SATURATION configs** (30%) will find them naturally
- **Expected:** No regression, stay at ~1-2Å

---

## 🎓 Key Advantages

### **1. Metadata-Independent**
- Don't rely on `seqid`, PDB numbering, or any ground truth info
- Works for any protein regardless of construct or numbering scheme

### **2. Comprehensive Coverage**
- 6 overlapping regions cover the entire protein
- No "blind spots" where we don't sample

### **3. Soft Constraints**
- `force: False` means Boltz can ignore bad regions
- Only guides exploration, doesn't force bad placements

### **4. Balanced Approach**
- 30% pure exploration (find common sites)
- 60% guided exploration (find allosteric sites)
- 10% high diversity (catch edge cases)

---

## 🔬 Comparison with Previous Approaches

| Approach | 3K5V Result | 6FVF Result | Pros | Cons |
|----------|-------------|-------------|------|------|
| **Pure SATURATION** | 25Å | 19.6Å | Simple, no bias | Limited diversity |
| **Terminus Probing (Hard)** | 33Å | 15.9Å | Forced exploration | Targeted wrong location (PDB numbering) |
| **Multi-Stage Ranking** | 33Å | 15.9Å | Ensured diversity | Still targeted wrong location |
| **Regional Sampling** | **< 20Å?** | **~15Å** | **Metadata-independent, systematic** | More configs needed |

---

## 🧪 Testing Plan

### **Test on 3K5V:**
```bash
python hackathon/predict_hackathon.py \
  --input-jsonl test_3k5v_only.jsonl \
  --msa-dir hackathon_data/datasets/asos_public/msa \
  --submission-dir regional_3k5v_predictions \
  --intermediate-dir regional_3k5v_intermediate \
  --result-folder regional_3k5v_results
```

### **What to Look For:**

**In console output:**
```
Config 3: REGION 1 residues   1- 73 (soft, 5 samples)
Config 4: REGION 2 residues  49-122 (soft, 5 samples)
...
Config 8: REGION 6 residues 241-293 (soft, 5 samples)
```

**In results:**
```
Top-1 RMSD: < 20Å  (improvement from 25Å or 33Å)
Top-5 min RMSD: < 15Å (best case)
```

### **Then test on full validation set:**
```bash
python hackathon/predict_hackathon.py \
  --input-jsonl validation_10.jsonl \
  --msa-dir hackathon_data/datasets/asos_public/msa \
  --submission-dir regional_validation_predictions \
  --intermediate-dir regional_validation_intermediate \
  --result-folder regional_validation_results
```

---

## 💡 Success Criteria

### **Minimum Success:**
- ✅ 3K5V improves from 25Å → < 20Å
- ✅ 6FVF maintains ~15-16Å
- ✅ No regression on orthosteric targets

### **Good Success:**
- ✅ 3K5V < 15Å
- ✅ 6FVF < 15Å
- ✅ Overall allosteric mean < 15Å

### **Excellent Success:**
- ✅ 3K5V < 10Å
- ✅ 6FVF < 10Å
- ✅ Competitive with state-of-the-art

---

## 🚀 Implementation Status

✅ **Code updated** in `predict_hackathon.py`
✅ **Regional sampling** implemented
✅ **Confidence-only ranking** restored (simple, robust)
✅ **Ready to test!**

---

## 📝 Final Notes

**This approach is:**
- ✅ **General:** Works for any protein
- ✅ **Robust:** No assumptions about metadata
- ✅ **Systematic:** Covers entire protein
- ✅ **Simple:** Easy to understand and debug

**It solves the fundamental problem:**
- ❌ Old: "Where does seqid=1 mean?" → Guess → Wrong
- ✅ New: "Where could the ligand bind?" → Sample everywhere → Find it

**Let's test it!** 🎯

