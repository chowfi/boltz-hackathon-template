# Stricter Constraints Plan for 3K5V

## 🔍 Problem Diagnosis

### **What We Know:**
```
3K5V Hard Constraints Results:
- All 5 final models: Distance to residue 1 = 27-28Å
- Ground truth: Ligand should be at residue 1 (seqid=1)
- Constraint: force=True, max_distance=6.0Å, residues 1-20

CONCLUSION: Hard constraints FAILED - either:
1. Boltz didn't generate N-terminus samples
2. OR samples were generated but filtered out (low confidence)
```

---

## 🎯 Root Cause Hypothesis

**Most likely: Generated but filtered out by confidence ranking**

**Evidence:**
- Boltz WITH constraints would generate 50 samples
- Config 8 (N-terminus, `force=True`) should generate 5 samples at N-terminus
- But final Top-5 has ZERO N-terminus samples
- **Conclusion:** The 5 N-terminus samples had lower confidence than the other 45 samples

**Why would N-terminus samples have low confidence?**
- N-terminus binding is RARE in training data
- Boltz learned that N-terminus = unlikely binding site
- Forced placement there → unnatural geometry → low confidence
- Confidence-only scoring → filters them out

---

## 🔧 Solution Strategies

### **Strategy 1: ULTRA-HARD Constraints (Tighter Distance)**
```python
"max_distance": 3.0,  # Changed from 6.0 - MUCH tighter
"force": True
```

**Pros:**
- Forces ligand even closer to residue 1
- Might improve fit if 6.0Å was too loose

**Cons:**
- Might be physically impossible (too tight)
- Could cause Boltz to crash or generate invalid structures

---

### **Strategy 2: Save ALL Samples (Skip Confidence Filtering)**
```python
# In post_process_protein_ligand(), return ALL samples, not just top-5 by confidence
return all_pdbs[:50]  # Return everything for analysis
```

**Pros:**
- Can see if N-terminus samples were actually generated
- Can manually inspect their quality

**Cons:**
- Not a solution for production (need Top-5 ranking)
- Just for debugging

---

### **Strategy 3: Weighted Ranking (Boost Terminus Samples)**
```python
# Give bonus to samples from config 8 and 9 (terminus configs)
for s in scores:
    if "config_8" in str(s["path"]) or "config_9" in str(s["path"]):
        s["confidence"] += 2.0  # Boost confidence by 2 points
```

**Pros:**
- Ensures terminus samples get into Top-5
- Balanced approach (not ignoring confidence completely)

**Cons:**
- Arbitrary boost amount
- Might promote bad samples

---

### **Strategy 4: Multi-Stage Ranking**
```python
# Stage 1: Force at least 1 sample from each config type into consideration
# Stage 2: Rank by confidence within each type
# Stage 3: Select Top-5 from mixed pool

terminus_samples = [s for s in scores if "config_8" or "config_9" in path]
other_samples = [s for s in scores if not terminus]

# Ensure at least 2 terminus samples in Top-5
top_terminus = sorted(terminus_samples, key=lambda x: x["confidence"])[:2]
top_other = sorted(other_samples, key=lambda x: x["confidence"])[:3]
final_top5 = sorted(top_terminus + top_other, key=lambda x: x["confidence"])[:5]
```

**Pros:**
- Guarantees diversity in Top-5
- Still uses confidence for ranking

**Cons:**
- Complex logic
- Might include terrible terminus samples

---

### **Strategy 5: Target ONLY Residue 1 (Not 1-20)**
```python
# Super focused - just residue 1!
n_term_contacts = [[protein_id, 1]]  # Only residue 1, not 1-20
```

**Pros:**
- Ultra-specific target
- Might work if issue is "residues 1-20 is too broad"

**Cons:**
- Might be TOO restrictive
- Residue 1 alone might not be enough contacts

---

## 🚀 Recommended Approach: **Multi-Pronged Attack**

### **Step 1: DEBUG - Save All Samples**
First, let's confirm the hypothesis by seeing ALL 50 samples:

```python
def post_process_protein_ligand_DEBUG(...):
    # Score all samples
    scores = [...]
    
    # Save diagnostic info
    print(f"\n🔍 ALL {len(scores)} SAMPLES:")
    for i, s in enumerate(sorted(scores, key=lambda x: x["confidence"], reverse=True)):
        config_num = "?" # extract from path
        print(f"  Rank {i+1:>2}: Config {config_num}, Confidence {s['confidence']:.2f}")
    
    # Return all 50 for manual inspection
    return [s["path"] for s in scores]
```

**Expected outcome:** We'll see that Config 8 samples have low confidence (e.g., < 90) compared to others (> 95).

---

### **Step 2: IMPLEMENT - Multi-Stage Ranking**
If Step 1 confirms filtering is the issue, use **Strategy 4**:

```python
def post_process_protein_ligand(...):
    # Separate by config type
    terminus_samples = [s for s in scores if is_terminus_config(s["path"])]
    other_samples = [s for s in scores if not is_terminus_config(s["path"])]
    
    # Get best from each pool
    top_terminus = sorted(terminus_samples, key=lambda x: x["confidence"], reverse=True)[:2]
    top_other = sorted(other_samples, key=lambda x: x["confidence"], reverse=True)[:10]
    
    # Mix and re-rank
    mixed_pool = top_terminus + top_other
    mixed_pool.sort(key=lambda x: x["confidence"], reverse=True)
    
    return [s["path"] for s in mixed_pool[:5]]
```

**Expected outcome:** At least 1-2 N-terminus samples in Top-5.

---

### **Step 3: REFINE - Tighter Constraints (if needed)**
If Step 2 gets N-terminus samples into Top-5 but RMSD is still bad, try:

```python
"max_distance": 4.0,  # Tighter than 6.0
"force": True
```

---

## 📊 Testing Plan

### **Test 1: Full Debug Run (Save All 50)**
```bash
# Modify code to return all 50 samples
python hackathon/predict_hackathon.py \
  --input-jsonl test_hard_allosteric.jsonl \
  --msa-dir hackathon_data/datasets/asos_public/msa \
  --submission-dir debug_all_samples \
  --intermediate-dir debug_all_intermediate \
  --result-folder debug_all_results
```

**Check:** Are Config 8 samples actually at N-terminus? What's their confidence?

---

### **Test 2: Multi-Stage Ranking**
```bash
# Implement Strategy 4
python hackathon/predict_hackathon.py \
  --input-jsonl test_hard_allosteric.jsonl \
  --msa-dir hackathon_data/datasets/asos_public/msa \
  --submission-dir multistage_predictions \
  --intermediate-dir multistage_intermediate \
  --result-folder multistage_results
```

**Success criteria:** 3K5V Top-1 RMSD < 15Å (currently 25Å)

---

## 🎯 Expected Outcomes

| Scenario | Probability | Action |
|----------|-------------|--------|
| N-terminus samples exist but have low confidence (< 90) | **High** | Use Multi-Stage Ranking |
| N-terminus samples don't exist at all | Medium | Investigate Boltz constraint implementation |
| N-terminus samples have good confidence but bad RMSD | Low | Try tighter max_distance |
| 3K5V is fundamentally impossible | Medium | Accept as limitation |

---

## 🔬 What We'll Learn

**From Test 1 (Debug):**
- Do N-terminus samples exist?
- What's their confidence distribution?
- Are they being filtered out?

**From Test 2 (Multi-Stage):**
- Does forcing terminus samples into Top-5 improve RMSD?
- Is the issue sampling or ranking?

---

## 💡 Key Insight

**The problem is NOT that Boltz can't generate N-terminus placements.**  
**The problem is that those placements have LOW CONFIDENCE and get filtered out.**

**Solution:** Change the ranking to ensure diverse pockets are represented, not just highest-confidence predictions.

---

**Let's start with Test 1 to confirm the hypothesis!** 🔍

