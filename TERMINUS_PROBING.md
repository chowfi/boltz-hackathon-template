# N/C-Terminus Probing Strategy

## 🎯 Problem Statement

**3K5V failure:** RMSD = 25.5Å (all 50 SATURATION samples)
- **Root cause:** Ligand binds at N-terminus (seqid=1, first residue!)
- **Boltz bias:** Never explores N/C termini (trained mostly on active site binding)
- **Result:** Even with 50 diverse seeds, all predictions converge to wrong pocket

---

## 💡 Solution: Explicit Terminus Probing

Use Boltz's **pocket constraint** feature to FORCE exploration of N/C termini.

### **Boltz Pocket Constraint Syntax:**

```yaml
constraints:
  - pocket:
      binder: LIGAND_CHAIN_ID      # e.g., "B"
      contacts: [[PROTEIN_CHAIN_ID, RES_IDX], ...]  # List of residues (1-indexed)
      max_distance: 6.0             # Angstroms (default)
      force: false                  # Soft constraint (conditioning, not hard constraint)
```

**Key Parameters:**
- `binder`: Which molecule to constrain (our ligand)
- `contacts`: List of protein residues to define the pocket
- `max_distance`: How close the ligand should be (6Å is standard)
- `force: false`: Uses constraint as **conditioning** (guides sampling, doesn't force it)

---

## 🔧 Implementation

### **Strategy: SATURATION + Terminus Probing**

**10 configs total:**
1. **8 SATURATION configs** (seeds: 42, 1000, 5000, 10000, 50000, 100000, 500000, 7777777)
   - Unconstrained, wide seed diversity
   - Covers "normal" allosteric and orthosteric sites

2. **1 N-terminus probe** (seed: 111111)
   - Constraint: residues 1-20
   - Targets 3K5V-like cases (N-terminus binding)

3. **1 C-terminus probe** (seed: 222222)
   - Constraint: last 20 residues
   - Coverage for C-terminus allosteric sites

**Total: 10 configs × 5 samples = 50 total samples**

---

### **Code Implementation:**

```python
def prepare_protein_ligand(datapoint_id: str, protein: Protein, 
                          ligands: list[SmallMolecule], 
                          input_dict: dict, 
                          msa_dir: Optional[Path] = None) -> List[tuple[dict, List[str]]]:
    configs = []
    
    ligand_id = ligands[0].id if ligands else "B"
    protein_id = protein.id
    protein_length = len(protein.sequence)
    
    # ========================================================================
    # Configs 0-7: SATURATION (unconstrained)
    # ========================================================================
    seeds = [42, 1000, 5000, 10000, 50000, 100000, 500000, 7777777]
    for i, seed in enumerate(seeds):
        config_dict = input_dict.copy()
        configs.append((config_dict, ["--diffusion_samples", "5", "--seed", str(seed)]))
    
    # ========================================================================
    # Config 8: N-terminus Probing
    # ========================================================================
    n_term_dict = input_dict.copy()
    n_term_residues = list(range(1, min(21, protein_length + 1)))  # 1-indexed, residues 1-20
    n_term_contacts = [[protein_id, res] for res in n_term_residues]
    
    n_term_dict["constraints"] = [{
        "pocket": {
            "binder": ligand_id,
            "contacts": n_term_contacts,
            "max_distance": 6.0,
            "force": False  # Soft constraint
        }
    }]
    
    configs.append((n_term_dict, ["--diffusion_samples", "5", "--seed", "111111"]))
    
    # ========================================================================
    # Config 9: C-terminus Probing
    # ========================================================================
    c_term_dict = input_dict.copy()
    c_term_start = max(1, protein_length - 19)  # Last 20 residues
    c_term_residues = list(range(c_term_start, protein_length + 1))
    c_term_contacts = [[protein_id, res] for res in c_term_residues]
    
    c_term_dict["constraints"] = [{
        "pocket": {
            "binder": ligand_id,
            "contacts": c_term_contacts,
            "max_distance": 6.0,
            "force": False
        }
    }]
    
    configs.append((c_term_dict, ["--diffusion_samples", "5", "--seed", "222222"]))
    
    return configs
```

---

## 📊 Expected Results

### **For 3K5V (N-terminus site at seqid=1):**

**Before (SATURATION only):**
- All 50 predictions: ~25.5Å RMSD
- No predictions near N-terminus
- Boltz never explored that region

**After (SATURATION + N-terminus probe):**
- 40 unconstrained predictions: ~25.5Å (same as before)
- **5 N-terminus constrained predictions: <10Å (hopefully!)**
- **1 C-terminus prediction: ~25Å** (wrong site, but good to check)

**Best case:** N-terminus constraint forces Boltz to try that region → finds correct site → RMSD < 5Å!

**Realistic case:** N-terminus constraint biases sampling → finds nearby pocket → RMSD ~8-12Å (improvement!)

**Worst case:** Even with constraint, can't fit ligand well → RMSD ~20Å (still better than 25Å)

---

### **For 6FVF (allosteric site at seqid=401, mid-protein):**

**Before (SATURATION + Ensemble):**
- Top-1: 19.61Å (ensemble picked wrong sample)
- Top-5 min: 15.68Å (good sample existed)

**After (SATURATION + Simple Scoring):**
- **Expected Top-1: ~15.7Å** (simple scoring should pick the good sample)
- Terminus probes won't help (site is mid-protein)
- But they won't hurt either (unconstrained samples still exist)

---

## 🧪 Testing Strategy

### **Step 1: Test on Hard Cases (3K5V, 6FVF)**

```bash
python hackathon/predict_hackathon.py \
  --input-jsonl test_hard_allosteric.jsonl \
  --msa-dir hackathon_data/datasets/asos_public/msa \
  --submission-dir test_terminus_predictions \
  --intermediate-dir test_terminus_intermediate \
  --result-folder test_terminus_results
```

**Success Metrics:**
- **3K5V:** Improve from 25.5Å → <10Å (target: <5Å)
- **6FVF:** Improve from 19.6Å → ~15.7Å
- **Mean:** Improve from 22.6Å → <13Å

---

### **Step 2: Check YAML Files**

After running, inspect the generated YAML files:

```bash
# Check N-terminus constraint was applied correctly
cat test_terminus_intermediate/input/3K5V_ALLOSTERIC_STJ_config_8.yaml
```

**Expected output:**
```yaml
sequences:
  - protein:
      id: A
      sequence: GAMDPSSPNYDKWE...
      msa: ...
  - ligand:
      id: B
      smiles: c1cc(cc(c1)C(=O)N)...
constraints:
  - pocket:
      binder: B
      contacts:
        - [A, 1]
        - [A, 2]
        - [A, 3]
        ...
        - [A, 20]
      max_distance: 6.0
      force: false
```

---

### **Step 3: Analyze Results**

```bash
# Check evaluation results
cat test_terminus_results/combined_results.csv

# Look for:
# - Did 3K5V improve?
# - Did 6FVF maintain or improve?
# - Which config produced the best sample?
```

**Analysis Questions:**
1. Did any N-terminus samples get <10Å for 3K5V?
2. If yes: SUCCESS! The constraint worked!
3. If no: Check if constraint was too strict or ligand doesn't actually fit

---

## 🎯 Why This Should Work

### **1. Constraints as "Priors"**

Boltz's pocket constraint with `force: false` acts as a **soft constraint** or **conditioning**:
- Biases the diffusion sampling toward the specified region
- Doesn't **force** the ligand there (would fail if it doesn't fit)
- Allows Boltz to still optimize for good geometry/interactions

**Analogy:** Like telling a GPS "prefer highways" vs. "only use highways"

---

### **2. Terminus Regions are Actually Reasonable**

N/C termini can be valid binding sites:
- ✅ Often flexible (good for induced fit)
- ✅ Can form pockets with nearby helices/sheets
- ✅ Allosteric regulation often happens at termini
- ❌ Just rare in training data (Boltz's prior is wrong)

**3K5V is a real PDB structure** → the ligand DOES bind there!
- Boltz just needs to be "encouraged" to try that region

---

### **3. Hedge Your Bets**

With this strategy:
- **8/10 configs are unconstrained** → maintain performance on normal cases
- **2/10 configs probe termini** → target the hard 20% failures

**Risk mitigation:**
- If terminus probing fails: Still have 40 unconstrained samples
- If terminus probing works: Get 5 good samples for hard cases
- **Worst case:** Same as before (25.5Å for 3K5V)
- **Best case:** Solve the unsolvable (3K5V < 5Å)

---

## 📈 Expected Overall Performance

### **On Full Validation Set:**

| Category | Count | Before | After | Expected |
|----------|-------|--------|-------|----------|
| **Easy Orthosteric** | 6 | <2Å | <2Å | ✅ Maintain |
| **Easy Allosteric** | 2 | ~2Å | ~2Å | ✅ Maintain |
| **Medium Allosteric** | 1 | 15.9Å | **~15.7Å** | ✅ Slight improvement |
| **Hard Allosteric** | 1 | 25.5Å | **<10Å?** | 🎯 Target |

**Overall Mean RMSD:**
- **Before:** ~8-10Å (with 3K5V dragging it up)
- **After:** ~6-8Å (if 3K5V improves to ~10Å)
- **Stretch goal:** ~5-6Å (if 3K5V improves to ~5Å)

---

## 🔬 Advanced: Other Potential Probing Strategies

If terminus probing works, could also try:

### **1. Active Site Probing** (for orthosteric)
```python
# If protein has known active site motif
active_site_residues = identify_catalytic_residues(protein.sequence)
# Use those for contacts
```

### **2. Charged Cluster Probing** (for allosteric)
```python
# Find clusters of charged residues (often regulatory sites)
charged_clusters = find_charged_clusters(protein.sequence)
# Probe each cluster
```

### **3. Hinge Region Probing** (for allosteric)
```python
# Allosteric sites often at domain interfaces
# Use MSA to identify flexible regions
hinge_residues = identify_hinge_regions(msa_data)
```

But for now, **terminus probing is the simplest and most targeted** approach for 3K5V!

---

## 🏆 Summary

**What we're doing:**
- Keep 8 SATURATION configs (works for 80% of cases)
- Add 2 terminus probes (targets the hard 20%)

**Why it should work:**
- 3K5V binds at N-terminus (seqid=1)
- Boltz never tries that region naturally
- Explicit constraint forces exploration
- Soft constraint (`force: false`) allows flexibility

**Expected outcome:**
- **3K5V:** 25.5Å → <10Å (maybe <5Å!)
- **6FVF:** 19.6Å → 15.7Å (simple scoring fix)
- **Others:** Maintain excellent performance

**Risk:** Minimal! Terminus probes are only 2/10 configs.

---

**Ready to test! This is our best shot at solving 3K5V.** 🚀

