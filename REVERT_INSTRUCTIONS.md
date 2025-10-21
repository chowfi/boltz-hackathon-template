# How to Revert to SATURATION + Confidence-Only Strategy

## Target Commit
**Hash:** `6090a45efa1110838c4d72113b20b0c948f03479`  
**Message:** "upload initial prediction pipeline"

This commit contains the **SATURATION + Confidence-Only** strategy that performed best.

---

## Option 1: Hard Reset (Discard All Changes)

⚠️ **WARNING:** This will **permanently delete** all uncommitted changes!

```bash
# Navigate to project directory
cd /Users/fionachow/Documents/Projects/boltz-hackathon-template

# Hard reset to target commit (DESTRUCTIVE!)
git reset --hard 6090a45efa1110838c4d72113b20b0c948f03479

# Clean untracked files (optional)
git clean -fd
```

**Use this if:** You want a clean slate and don't need any current changes.

---

## Option 2: Create a New Branch (Keep Current Work)

✅ **RECOMMENDED:** This preserves your current work in a separate branch.

```bash
# Navigate to project directory
cd /Users/fionachow/Documents/Projects/boltz-hackathon-template

# Create a new branch from current state (save your work)
git checkout -b experimental-work

# Add and commit all current changes to this branch
git add -A
git commit -m "Save experimental work: terminus probing, multi-stage ranking, regional sampling"

# Switch back to main branch
git checkout main

# Reset main to target commit
git reset --hard 6090a45efa1110838c4d72113b20b0c948f03479
```

**Use this if:** You want to preserve your experimental work for later reference.

---

## Option 3: Checkout Specific File (Partial Revert)

This reverts only `predict_hackathon.py` to the target commit, keeping everything else.

```bash
# Navigate to project directory
cd /Users/fionachow/Documents/Projects/boltz-hackathon-template

# Checkout just the prediction script from target commit
git checkout 6090a45efa1110838c4d72113b20b0c948f03479 -- hackathon/predict_hackathon.py

# The file is now reverted, but you need to commit the change
git add hackathon/predict_hackathon.py
git commit -m "Revert predict_hackathon.py to SATURATION + Confidence-Only strategy"
```

**Use this if:** You only want to revert the prediction script, keeping other files.

---

## Option 4: Manual Cherry-Pick (Most Control)

View the target commit's changes and manually apply them.

```bash
# View what changed in that commit
git show 6090a45efa1110838c4d72113b20b0c948f03479

# View the full file as it was at that commit
git show 6090a45efa1110838c4d72113b20b0c948f03479:hackathon/predict_hackathon.py

# Copy specific sections manually to your current file
```

**Use this if:** You want fine-grained control over what to revert.

---

## Verification After Revert

After reverting, verify the strategy is correct:

```bash
# Check the current commit
git log -1 --oneline

# Check the key functions in predict_hackathon.py
grep -A 20 "def prepare_protein_ligand" hackathon/predict_hackathon.py
grep -A 20 "def post_process_protein_ligand" hackathon/predict_hackathon.py
```

**Expected behavior:**
- `prepare_protein_ligand`: Should generate 10 configs with wide seed spacing
- `post_process_protein_ligand`: Should rank by confidence only (simple sorting)

---

## Running Full Evaluation on AWS

After reverting:

```bash
# On AWS instance, pull the reverted code
git pull origin main

# Activate conda environment
conda activate boltz

# Run on full dataset
python3 hackathon/predict_hackathon.py \
    --dataset hackathon_data/datasets/asos_public/asos_public.jsonl \
    --submission-folder final_predictions \
    --output-dir final_outputs

# Evaluate results
python3 hackathon/evaluate_asos.py \
    --dataset-file hackathon_data/datasets/asos_public/asos_public.jsonl \
    --dataset-folder hackathon_data/datasets/asos_public \
    --submission-folder final_predictions \
    --result-folder final_results
```

---

## What's in the Target Commit

The commit `6090a45` contains:

### `prepare_protein_ligand`:
- 10 configurations
- Wide seed spacing: `[42, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 7777777]`
- Each config: 5 diffusion samples
- **No constraints** (pure SATURATION)
- Total: 50 samples per target

### `post_process_protein_ligand`:
- Simple confidence-only ranking
- Sort by Boltz confidence (highest first)
- Return top 5

### Why This Works:
- ✅ Maximum diversity through seed spacing
- ✅ Full protein context preserved
- ✅ Simple, robust ranking
- ✅ No reliance on ambiguous metadata

---

## Recommended: Option 2 (New Branch)

**I recommend Option 2** because:
1. ✅ Saves your experimental work
2. ✅ Clean revert for AWS run
3. ✅ Can compare results later
4. ✅ Easy to switch back if needed

```bash
# Quick copy-paste commands:
cd /Users/fionachow/Documents/Projects/boltz-hackathon-template
git checkout -b experimental-work
git add -A
git commit -m "Save experimental work: terminus probing and regional sampling"
git checkout main
git reset --hard 6090a45efa1110838c4d72113b20b0c948f03479
```

Done! Ready for AWS run. 🚀

