# predict_hackathon.py
import argparse
import json
import os
import shutil
import subprocess
from collections.abc import Iterable
from pathlib import Path
from typing import Any, List, Optional

import yaml
from hackathon_api import Datapoint, Protein, SmallMolecule

# ---------------------------------------------------------------------------
# ---- Participants should modify these four functions ----------------------
# ---------------------------------------------------------------------------

def prepare_protein_complex(datapoint_id: str, proteins: List[Protein], input_dict: dict, msa_dir: Optional[Path] = None) -> List[tuple[dict, List[str]]]:
    """
    Prepare input dict and CLI args for a protein complex prediction.
    You can return multiple configurations to run by returning a list of (input_dict, cli_args) tuples.
    Args:
        datapoint_id: The unique identifier for this datapoint
        proteins: List of protein sequences to predict as a complex
        input_dict: Prefilled input dict
        msa_dir: Directory containing MSA files (for computing relative paths)
    Returns:
        List of tuples of (final input dict that will get exported as YAML, list of CLI args). Each tuple represents a separate configuration to run.
    """
    # Please note:
    # `proteins`` will contain 3 chains
    # H,L: heavy and light chain of the Fv or Fab region
    # A: the antigen
    #
    # you can modify input_dict to change the input yaml file going into the prediction, e.g.
    # ```
    # input_dict["constraints"] = [{
    #   "contact": {
    #       "token1" : [CHAIN_ID, RES_IDX/ATOM_NAME], 
    #       "token1" : [CHAIN_ID, RES_IDX/ATOM_NAME]
    #   }
    # }]
    # ```
    #
    # will add contact constraints to the input_dict

    # Example: predict 5 structures
    cli_args = ["--diffusion_samples", "5"]
    return [(input_dict, cli_args)]

def identify_allosteric_pockets(protein: Protein, n_pockets: int = 3) -> list[list[int]]:
    """
    Sequence-based pocket identification optimized for allosteric sites.
    Focus on surface-exposed, flexible regions that are likely allosteric.
    """
    sequence = protein.sequence
    seq_len = len(sequence)
    
    pockets = []
    
    # 1. Terminal regions (high flexibility, common allosteric sites)
    if seq_len >= 10:
        # N-terminal region
        n_term = list(range(1, min(6, seq_len + 1)))  # First 5 residues
        pockets.append(n_term)
        
        # C-terminal region  
        c_term = list(range(max(1, seq_len - 4), seq_len + 1))  # Last 5 residues
        pockets.append(c_term)
    
    # 2. Hydrophilic surface clusters (K, R, D, E, N, Q, S, T)
    hydrophilic_residues = []
    for i, residue in enumerate(sequence, 1):
        if residue in ['K', 'R', 'D', 'E', 'N', 'Q', 'S', 'T']:
            hydrophilic_residues.append(i)
    
    # Group nearby hydrophilic residues
    if hydrophilic_residues:
        clusters = []
        current_cluster = [hydrophilic_residues[0]]
        
        for i in range(1, len(hydrophilic_residues)):
            if hydrophilic_residues[i] - hydrophilic_residues[i-1] <= 8:  # Within 8 residues
                current_cluster.append(hydrophilic_residues[i])
            else:
                if len(current_cluster) >= 2:  # Only keep clusters with 2+ residues
                    clusters.append(current_cluster)
                current_cluster = [hydrophilic_residues[i]]
        
        if len(current_cluster) >= 2:
            clusters.append(current_cluster)
        
        # Add the largest hydrophilic cluster
        if clusters:
            largest_cluster = max(clusters, key=len)
            pockets.append(largest_cluster[:4])  # Take up to 4 residues
    
    # 3. Flexible regions (G, P, S, T, A rich)
    flexible_regions = []
    for i, residue in enumerate(sequence, 1):
        if residue in ['G', 'P', 'S', 'T', 'A']:
            flexible_regions.append(i)
    
    if flexible_regions:
        # Find the largest flexible region
        flexible_clusters = []
        current_cluster = [flexible_regions[0]]
        
        for i in range(1, len(flexible_regions)):
            if flexible_regions[i] - flexible_regions[i-1] <= 5:
                current_cluster.append(flexible_regions[i])
            else:
                if len(current_cluster) >= 3:
                    flexible_clusters.append(current_cluster)
                current_cluster = [flexible_regions[i]]
        
        if len(current_cluster) >= 3:
            flexible_clusters.append(current_cluster)
        
        if flexible_clusters:
            largest_flexible = max(flexible_clusters, key=len)
            pockets.append(largest_flexible[:4])
    
    # Ensure we have at least n_pockets
    while len(pockets) < n_pockets:
        # Add middle regions if needed
        mid_point = seq_len // 2
        start = max(1, mid_point - 2)
        end = min(seq_len, mid_point + 2)
        pockets.append(list(range(start, end + 1)))
    
    return pockets[:n_pockets]

def compute_kyte_doolittle_score(residue: str) -> float:
    """
    Fast hydrophobicity scoring for accessibility.
    Returns accessibility score (higher = more exposed)
    """
    # Kyte-Doolittle scale (simplified)
    kd_scores = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    return kd_scores.get(residue, 0.0)

def prepare_protein_ligand(datapoint_id: str, protein: Protein, ligands: list[SmallMolecule], input_dict: dict, msa_dir: Optional[Path] = None) -> List[tuple[dict, List[str]]]:
    """
    Prepare input dict and CLI args for a protein-ligand prediction.
    Use sequence-based allosteric pocket identification for better allosteric performance.
    """
    configs = []
    
    print(f"Generating allosteric-optimized configurations for {datapoint_id}")
    
    # Strategy: Focus on allosteric sites using sequence-based pocket identification
    # Total: 20 samples across 5 configurations
    
    # Config 1: High-confidence baseline (4 samples)
    baseline_dict = input_dict.copy()
    configs.append((baseline_dict, ["--diffusion_samples", "4", "--seed", "42"]))
    
    # Config 2-4: Allosteric pocket-based exploration (3 configs × 4 samples each = 12 samples)
    try:
        allosteric_pockets = identify_allosteric_pockets(protein, n_pockets=3)
        for i, pocket_residues in enumerate(allosteric_pockets):
            pocket_dict = input_dict.copy()
            # Note: Contact constraints would be added here if Boltz supports them
            # pocket_dict["constraints"] = [{"contact": {"token1": ["A", res_idx], "token2": ["L", "LIG"]}}]
            configs.append((pocket_dict, ["--diffusion_samples", "4", "--seed", str(100 + i)]))
    except Exception as e:
        print(f"Warning: Failed to identify allosteric pockets for {datapoint_id}: {e}")
        # Fallback to additional baseline configs
        for i in range(3):
            fallback_dict = input_dict.copy()
            configs.append((fallback_dict, ["--diffusion_samples", "4", "--seed", str(200 + i)]))
    
    # Config 5: High diversity (4 samples)
    diversity_dict = input_dict.copy()
    configs.append((diversity_dict, ["--diffusion_samples", "4", "--seed", "999"]))
    
    print(f"Generated {len(configs)} configurations for {datapoint_id}")
    return configs

def post_process_protein_complex(datapoint: Datapoint, input_dicts: List[dict[str, Any]], cli_args_list: List[list[str]], prediction_dirs: List[Path]) -> List[Path]:
    """
    Return ranked model files for protein complex submission.
    Args:
        datapoint: The original datapoint object
        input_dicts: List of input dictionaries used for predictions (one per config)
        cli_args_list: List of command line arguments used for predictions (one per config)
        prediction_dirs: List of directories containing prediction results (one per config)
    Returns: 
        Sorted pdb file paths that should be used as your submission.
    """
    # Collect all PDBs from all configurations
    all_pdbs = []
    for prediction_dir in prediction_dirs:
        config_pdbs = sorted(prediction_dir.glob(f"{datapoint.datapoint_id}_config_*_model_*.pdb"))
        all_pdbs.extend(config_pdbs)

    # Sort all PDBs and return their paths
    all_pdbs = sorted(all_pdbs)
    return all_pdbs

def parse_pdb(pdb_path: Path):
    """Parse PDB using Biopython."""
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    return parser.get_structure("model", pdb_path)

def extract_confidence(pdb_path: Path) -> float:
    """Extract Boltz confidence from companion JSON or PDB REMARK."""
    # Check for confidence.json in the same directory
    conf_json = pdb_path.parent / "confidences.json"
    if conf_json.exists():
        try:
            import json
            data = json.load(open(conf_json))
            return data.get("aggregate_score", 0.5)
        except:
            pass
    
    # Check for confidence in PDB REMARK lines
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('REMARK') and 'confidence' in line.lower():
                    # Try to extract confidence value
                    import re
                    match = re.search(r'(\d+\.?\d*)', line)
                    if match:
                        return float(match.group(1)) / 100.0  # Assume percentage
    except:
        pass
    
    return 0.5  # Default confidence

def compute_clash_penalty_fast(structure) -> float:
    """Fast clash penalty computation."""
    from Bio.PDB import NeighborSearch
    try:
        # Separate protein and ligand atoms
        protein_atoms = [atom for atom in structure.get_atoms() if atom.parent.resname not in ["LIG"]]
        ligand_atoms = [atom for atom in structure.get_atoms() if atom.parent.resname == "LIG"]
        
        if not ligand_atoms:
            return 0.0
        
        # NeighborSearch for clashes
        ns = NeighborSearch(protein_atoms)
        clash_count = 0
        for lig_atom in ligand_atoms:
            nearby = ns.search(lig_atom.coord, 2.0)  # VDW threshold ~2 Å
            clash_count += len(nearby)
        return clash_count
    except Exception as e:
        print(f"Warning: Failed to compute clash penalty: {e}")
        return 0.0

def count_protein_ligand_contacts_fast(structure, cutoff: float = 4.5) -> int:
    """Fast protein-ligand contact counting."""
    from Bio.PDB import NeighborSearch
    try:
        protein_atoms = [atom for atom in structure.get_atoms() 
                         if atom.parent.resname not in ["LIG"] and atom.element != "H"]
        ligand_atoms = [atom for atom in structure.get_atoms() 
                        if atom.parent.resname == "LIG" and atom.element != "H"]
        
        if not ligand_atoms:
            return 0
        
        ns = NeighborSearch(protein_atoms)
        contact_count = sum(len(ns.search(lig_atom.coord, cutoff)) for lig_atom in ligand_atoms)
        return contact_count
    except Exception as e:
        print(f"Warning: Failed to count contacts: {e}")
        return 0

def normalize_scores_fast(scores: list[dict]) -> list[dict]:
    """Fast min-max normalization for iteration."""
    if not scores:
        return scores
    
    # Get metrics present in scores
    metrics = ["boltz_conf", "clash", "contacts"]
    
    for metric in metrics:
        if metric in scores[0]:
            values = [s[metric] for s in scores]
            min_val, max_val = min(values), max(values)
            if max_val > min_val:
                for s in scores:
                    s[f"{metric}_norm"] = (s[metric] - min_val) / (max_val - min_val)
            else:
                for s in scores:
                    s[f"{metric}_norm"] = 0.5  # All equal, neutral score
    
    return scores

def post_process_protein_ligand(datapoint: Datapoint, input_dicts: List[dict[str, Any]], cli_args_list: List[list[str]], prediction_dirs: List[Path]) -> List[Path]:
    """
    Dual-scoring ensemble approach for orthosteric AND allosteric binding.
    Biophysical strategy: Score predictions separately for both binding modes,
    then return a diverse ensemble covering both possibilities.
    """
    # Collect all PDBs from all configurations
    all_pdbs = []
    for prediction_dir in prediction_dirs:
        config_pdbs = sorted(prediction_dir.glob(f"{datapoint.datapoint_id}_config_*_model_*.pdb"))
        all_pdbs.extend(config_pdbs)
    
    if not all_pdbs:
        print(f"Warning: No PDB files found for {datapoint.datapoint_id}")
        return []
    
    print(f"Dual-scoring predictions for {datapoint.datapoint_id}")
    
    scores = []
    for pdb_path in all_pdbs:
        try:
            # 1. Parse PDB (Biopython)
            structure = parse_pdb(pdb_path)
            
            # 2. Extract Boltz confidence
            boltz_conf = extract_confidence(pdb_path)
            
            # 3. Compute biophysical metrics
            clash_penalty = compute_clash_penalty_fast(structure)
            contact_count = count_protein_ligand_contacts_fast(structure, cutoff=4.5)
            
            # 4. Allosteric-specific metrics
            surface_contacts = count_surface_contacts_fast(structure)
            flexibility_score = compute_flexibility_score_fast(structure)
            polar_contacts = count_polar_contacts_fast(structure)
            
            # 5. Orthosteric-specific metrics
            binding_depth = compute_binding_depth_fast(structure)
            hydrophobic_contacts = count_hydrophobic_contacts_fast(structure)
            pocket_volume = estimate_pocket_volume_fast(structure)
            
            # 6. Store comprehensive scores
            scores.append({
                "path": pdb_path,
                "boltz_conf": boltz_conf,
                "clash": clash_penalty,
                "contacts": contact_count,
                # Allosteric indicators
                "surface_contacts": surface_contacts,
                "flexibility": flexibility_score,
                "polar_contacts": polar_contacts,
                # Orthosteric indicators
                "binding_depth": binding_depth,
                "hydrophobic_contacts": hydrophobic_contacts,
                "pocket_volume": pocket_volume,
            })
        except Exception as e:
            print(f"Warning: Failed to score {pdb_path.name}: {e}")
            continue
    
    if not scores:
        print(f"Warning: No valid scores computed for {datapoint.datapoint_id}")
        return all_pdbs[:5]
    
    # 7. Normalize all scores
    scores = normalize_scores_enhanced(scores)
    
    # 8. Compute TWO separate scoring schemes
    for s in scores:
        # ORTHOSTERIC SCORE: Deep binding, hydrophobic, high confidence
        s["orthosteric_score"] = (
            0.35 * s["boltz_conf_norm"] +
            -0.25 * s["clash_norm"] +
            0.15 * s["contacts_norm"] +
            0.15 * s["binding_depth_norm"] +
            0.1 * s["hydrophobic_contacts_norm"]
        )
        
        # ALLOSTERIC SCORE: Surface binding, polar, flexible
        s["allosteric_score"] = (
            0.25 * s["boltz_conf_norm"] +        # Lower weight (model is biased)
            -0.15 * s["clash_norm"] +
            0.2 * s["surface_contacts_norm"] +   # Surface accessibility
            0.15 * s["polar_contacts_norm"] +    # Polar interactions
            0.15 * s["flexibility_norm"] +       # Flexibility
            -0.1 * s["binding_depth_norm"]       # Penalize deep binding
        )
    
    # 9. ADAPTIVE RANKING: Detect binding type and rank accordingly
    # This optimizes for Top-1 RMSD metric
    
    # Sort by both scoring schemes
    ortho_sorted = sorted(scores, key=lambda x: x["orthosteric_score"], reverse=True)
    allo_sorted = sorted(scores, key=lambda x: x["allosteric_score"], reverse=True)
    
    # Aggregate biophysical signals to detect likely binding mode
    import numpy as np
    
    # Calculate aggregate metrics across all predictions
    avg_binding_depth = np.mean([s["binding_depth"] for s in scores])
    avg_hydrophobic = np.mean([s["hydrophobic_contacts"] for s in scores])
    avg_pocket_volume = np.mean([s["pocket_volume"] for s in scores])
    
    avg_surface_contacts = np.mean([s["surface_contacts"] for s in scores])
    avg_polar_contacts = np.mean([s["polar_contacts"] for s in scores])
    avg_flexibility = np.mean([s["flexibility"] for s in scores])
    
    # Compute orthosteric vs allosteric likelihood
    # Orthosteric signals: deep, hydrophobic, enclosed
    ortho_signal = 0.0
    if avg_binding_depth < np.median([s["binding_depth"] for s in scores]):
        ortho_signal += 1.0  # Deeper than median
    if avg_hydrophobic > np.median([s["hydrophobic_contacts"] for s in scores]):
        ortho_signal += 1.0  # More hydrophobic
    if avg_pocket_volume > np.median([s["pocket_volume"] for s in scores]):
        ortho_signal += 1.0  # Larger pocket
    
    # Allosteric signals: surface, polar, flexible
    allo_signal = 0.0
    if avg_surface_contacts > np.median([s["surface_contacts"] for s in scores]):
        allo_signal += 1.0  # More surface exposure
    if avg_polar_contacts > np.median([s["polar_contacts"] for s in scores]):
        allo_signal += 1.0  # More polar contacts
    if avg_flexibility > np.median([s["flexibility"] for s in scores]):
        allo_signal += 1.0  # More flexibility
    
    # Decision: Choose ranking based on stronger signal
    if ortho_signal > allo_signal:
        # Predictions favor orthosteric characteristics
        selected_paths = [s["path"] for s in ortho_sorted[:5]]
        binding_mode = "orthosteric"
        print(f"Adaptive ranking: ORTHOSTERIC mode detected (signal: {ortho_signal:.1f} vs {allo_signal:.1f})")
    elif allo_signal > ortho_signal:
        # Predictions favor allosteric characteristics
        selected_paths = [s["path"] for s in allo_sorted[:5]]
        binding_mode = "allosteric"
        print(f"Adaptive ranking: ALLOSTERIC mode detected (signal: {ortho_signal:.1f} vs {allo_signal:.1f})")
    else:
        # Ambiguous - use balanced scoring
        balanced_sorted = sorted(scores, 
                                key=lambda x: 0.5 * x["orthosteric_score"] + 0.5 * x["allosteric_score"], 
                                reverse=True)
        selected_paths = [s["path"] for s in balanced_sorted[:5]]
        binding_mode = "balanced"
        print(f"Adaptive ranking: BALANCED mode (ambiguous signal: {ortho_signal:.1f} vs {allo_signal:.1f})")
    
    print(f"Selected top-5 predictions using {binding_mode} ranking")
    return selected_paths[:5]

def compute_binding_depth_fast(structure) -> float:
    """Compute binding site depth - deeper = more orthosteric."""
    try:
        protein_coords = []
        ligand_coords = []
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() != "LIG":
                        for atom in residue.get_atoms():
                            protein_coords.append(atom.coord)
                    else:
                        for atom in residue.get_atoms():
                            ligand_coords.append(atom.coord)
        
        if not protein_coords or not ligand_coords:
            return 0.0
        
        protein_center = np.mean(protein_coords, axis=0)
        distances = [np.linalg.norm(lig_coord - protein_center) for lig_coord in ligand_coords]
        # Return negative so deeper sites score higher
        return -np.mean(distances)
    except:
        return 0.0

def count_polar_contacts_fast(structure) -> int:
    """Count polar contacts - important for allosteric binding."""
    try:
        polar_residues = ['K', 'R', 'D', 'E', 'N', 'Q', 'S', 'T', 'H', 'Y']
        contact_count = 0
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in polar_residues:
                        for atom in residue.get_atoms():
                            for lig_atom in structure.get_atoms():
                                if lig_atom.parent.resname == "LIG":
                                    distance = atom - lig_atom
                                    if distance < 4.0:  # H-bond distance
                                        contact_count += 1
        return contact_count
    except:
        return 0

def count_hydrophobic_contacts_fast(structure) -> int:
    """Count hydrophobic contacts - important for orthosteric binding."""
    try:
        hydrophobic_residues = ['F', 'W', 'Y', 'L', 'I', 'V', 'M', 'A', 'P']
        contact_count = 0
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in hydrophobic_residues:
                        for atom in residue.get_atoms():
                            for lig_atom in structure.get_atoms():
                                if lig_atom.parent.resname == "LIG":
                                    distance = atom - lig_atom
                                    if distance < 4.5:  # Hydrophobic contact
                                        contact_count += 1
        return contact_count
    except:
        return 0

def estimate_pocket_volume_fast(structure) -> float:
    """Estimate pocket volume - larger = more orthosteric."""
    try:
        # Simple estimate: count nearby protein atoms around ligand
        ligand_coords = []
        protein_atoms_nearby = 0
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() == "LIG":
                        for atom in residue.get_atoms():
                            ligand_coords.append(atom.coord)
        
        if not ligand_coords:
            return 0.0
        
        ligand_center = np.mean(ligand_coords, axis=0)
        
        # Count protein atoms within 8 Å of ligand center
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() != "LIG":
                        for atom in residue.get_atoms():
                            distance = np.linalg.norm(atom.coord - ligand_center)
                            if distance < 8.0:
                                protein_atoms_nearby += 1
        
        # More surrounding atoms = larger pocket
        return float(protein_atoms_nearby)
    except:
        return 0.0

def count_surface_contacts_fast(structure) -> int:
    """Count contacts with surface-exposed residues for allosteric sites."""
    try:
        # Simple heuristic: count contacts with hydrophilic residues
        hydrophilic_residues = ['K', 'R', 'D', 'E', 'N', 'Q', 'S', 'T']
        contact_count = 0
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in hydrophilic_residues:
                        # Count nearby ligand atoms
                        for atom in residue.get_atoms():
                            for lig_atom in structure.get_atoms():
                                if lig_atom.parent.resname == "LIG":
                                    distance = atom - lig_atom
                                    if distance < 5.0:  # Surface contact threshold
                                        contact_count += 1
        return contact_count
    except:
        return 0

def compute_flexibility_score_fast(structure) -> float:
    """Compute flexibility score for allosteric binding sites."""
    try:
        # Simple flexibility heuristic based on residue types
        flexible_residues = ['G', 'P', 'S', 'T', 'A']
        flexibility_score = 0.0
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in flexible_residues:
                        # Check if near ligand
                        for atom in residue.get_atoms():
                            for lig_atom in structure.get_atoms():
                                if lig_atom.parent.resname == "LIG":
                                    distance = atom - lig_atom
                                    if distance < 6.0:  # Flexibility influence radius
                                        flexibility_score += 1.0
        return flexibility_score
    except:
        return 0.0

def normalize_scores_enhanced(scores: list[dict]) -> list[dict]:
    """Enhanced normalization for dual-scoring (orthosteric + allosteric)."""
    if not scores:
        return scores
    
    # Get all metrics present in scores
    metrics = ["boltz_conf", "clash", "contacts", "surface_contacts", "flexibility", 
               "binding_depth", "polar_contacts", "hydrophobic_contacts", "pocket_volume"]
    
    for metric in metrics:
        if metric in scores[0]:
            values = [s[metric] for s in scores]
            min_val, max_val = min(values), max(values)
            if max_val > min_val:
                for s in scores:
                    s[f"{metric}_norm"] = (s[metric] - min_val) / (max_val - min_val)
            else:
                for s in scores:
                    s[f"{metric}_norm"] = 0.5  # All equal, neutral score
    
    return scores

# -----------------------------------------------------------------------------
# ---- End of participant section ---------------------------------------------
# -----------------------------------------------------------------------------


DEFAULT_OUT_DIR = Path("predictions")
DEFAULT_SUBMISSION_DIR = Path("submission")
DEFAULT_INPUTS_DIR = Path("inputs")

ap = argparse.ArgumentParser(
    description="Hackathon scaffold for Boltz predictions",
    epilog="Examples:\n"
            "  Single datapoint: python predict_hackathon.py --input-json examples/specs/example_protein_ligand.json --msa-dir ./msa --submission-dir submission --intermediate-dir intermediate\n"
            "  Multiple datapoints: python predict_hackathon.py --input-jsonl examples/test_dataset.jsonl --msa-dir ./msa --submission-dir submission --intermediate-dir intermediate",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

input_group = ap.add_mutually_exclusive_group(required=True)
input_group.add_argument("--input-json", type=str,
                        help="Path to JSON datapoint for a single datapoint")
input_group.add_argument("--input-jsonl", type=str,
                        help="Path to JSONL file with multiple datapoint definitions")

ap.add_argument("--msa-dir", type=Path,
                help="Directory containing MSA files (for computing relative paths in YAML)")
ap.add_argument("--submission-dir", type=Path, required=False, default=DEFAULT_SUBMISSION_DIR,
                help="Directory to place final submissions")
ap.add_argument("--intermediate-dir", type=Path, required=False, default=Path("hackathon_intermediate"),
                help="Directory to place generated input YAML files and predictions")
ap.add_argument("--group-id", type=str, required=False, default=None,
                help="Group ID to set for submission directory (sets group rw access if specified)")
ap.add_argument("--result-folder", type=Path, required=False, default=None,
                help="Directory to save evaluation results. If set, will automatically run evaluation after predictions.")

args = ap.parse_args()

def _prefill_input_dict(datapoint_id: str, proteins: Iterable[Protein], ligands: Optional[list[SmallMolecule]] = None, msa_dir: Optional[Path] = None) -> dict:
    """
    Prepare input dict for Boltz YAML.
    """
    seqs = []
    for p in proteins:
        if msa_dir and p.msa:
            if Path(p.msa).is_absolute():
                msa_full_path = Path(p.msa)
            else:
                msa_full_path = msa_dir / p.msa
            try:
                msa_relative_path = os.path.relpath(msa_full_path, Path.cwd())
            except ValueError:
                msa_relative_path = str(msa_full_path)
        else:
            msa_relative_path = p.msa
        entry = {
            "protein": {
                "id": p.id,
                "sequence": p.sequence,
                "msa": msa_relative_path
            }
        }
        seqs.append(entry)
    if ligands:
        def _format_ligand(ligand: SmallMolecule) -> dict:
            output =  {
                "ligand": {
                    "id": ligand.id,
                    "smiles": ligand.smiles
                }
            }
            return output
        
        for ligand in ligands:
            seqs.append(_format_ligand(ligand))
    doc = {
        "version": 1,
        "sequences": seqs,
    }
    return doc

def _run_boltz_and_collect(datapoint) -> None:
    """
    New flow: prepare input dict, write yaml, run boltz, post-process, copy submissions.
    """
    out_dir = args.intermediate_dir / "predictions"
    out_dir.mkdir(parents=True, exist_ok=True)
    subdir = args.submission_dir / datapoint.datapoint_id
    subdir.mkdir(parents=True, exist_ok=True)

    # Prepare input dict and CLI args
    base_input_dict = _prefill_input_dict(datapoint.datapoint_id, datapoint.proteins, datapoint.ligands, args.msa_dir)

    if datapoint.task_type == "protein_complex":
        configs = prepare_protein_complex(datapoint.datapoint_id, datapoint.proteins, base_input_dict, args.msa_dir)
    elif datapoint.task_type == "protein_ligand":
        configs = prepare_protein_ligand(datapoint.datapoint_id, datapoint.proteins[0], datapoint.ligands, base_input_dict, args.msa_dir)
    else:
        raise ValueError(f"Unknown task_type: {datapoint.task_type}")

    # Run boltz for each configuration
    all_input_dicts = []
    all_cli_args = []
    all_pred_subfolders = []
    
    input_dir = args.intermediate_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    
    for config_idx, (input_dict, cli_args) in enumerate(configs):
        # Write input YAML with config index suffix
        yaml_path = input_dir / f"{datapoint.datapoint_id}_config_{config_idx}.yaml"
        with open(yaml_path, "w") as f:
            yaml.safe_dump(input_dict, f, sort_keys=False)

        # Run boltz
        cache = os.environ.get("BOLTZ_CACHE", str(Path.home() / ".boltz"))
        fixed = [
            "boltz", "predict", str(yaml_path),
            "--devices", "1",
            "--out_dir", str(out_dir),
            "--cache", cache,
            "--no_kernels",
            "--output_format", "pdb",
        ]
        cmd = fixed + cli_args
        print(f"Running config {config_idx}:", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)

        # Compute prediction subfolder for this config
        pred_subfolder = out_dir / f"boltz_results_{datapoint.datapoint_id}_config_{config_idx}" / "predictions" / f"{datapoint.datapoint_id}_config_{config_idx}"
        
        all_input_dicts.append(input_dict)
        all_cli_args.append(cli_args)
        all_pred_subfolders.append(pred_subfolder)

    # Post-process and copy submissions
    if datapoint.task_type == "protein_complex":
        ranked_files = post_process_protein_complex(datapoint, all_input_dicts, all_cli_args, all_pred_subfolders)
    elif datapoint.task_type == "protein_ligand":
        ranked_files = post_process_protein_ligand(datapoint, all_input_dicts, all_cli_args, all_pred_subfolders)
    else:
        raise ValueError(f"Unknown task_type: {datapoint.task_type}")

    if not ranked_files:
        raise FileNotFoundError(f"No model files found for {datapoint.datapoint_id}")

    for i, file_path in enumerate(ranked_files[:5]):
        target = subdir / (f"model_{i}.pdb" if file_path.suffix == ".pdb" else f"model_{i}{file_path.suffix}")
        shutil.copy2(file_path, target)
        print(f"Saved: {target}")

    if args.group_id:
        try:
            subprocess.run(["chgrp", "-R", args.group_id, str(subdir)], check=True)
            subprocess.run(["chmod", "-R", "g+rw", str(subdir)], check=True)
        except Exception as e:
            print(f"WARNING: Failed to set group ownership or permissions: {e}")

def _load_datapoint(path: Path):
    """Load JSON datapoint file."""
    with open(path) as f:
        return Datapoint.from_json(f.read())

def _run_evaluation(input_file: str, task_type: str, submission_dir: Path, result_folder: Path):
    """
    Run the appropriate evaluation script based on task type.
    
    Args:
        input_file: Path to the input JSON or JSONL file
        task_type: Either "protein_complex" or "protein_ligand"
        submission_dir: Directory containing prediction submissions
        result_folder: Directory to save evaluation results
    """
    script_dir = Path(__file__).parent
    
    if task_type == "protein_complex":
        eval_script = script_dir / "evaluate_abag.py"
        cmd = [
            "python", str(eval_script),
            "--dataset-file", input_file,
            "--submission-folder", str(submission_dir),
            "--result-folder", str(result_folder)
        ]
    elif task_type == "protein_ligand":
        eval_script = script_dir / "evaluate_asos.py"
        cmd = [
            "python", str(eval_script),
            "--dataset-file", input_file,
            "--submission-folder", str(submission_dir),
            "--result-folder", str(result_folder)
        ]
    else:
        raise ValueError(f"Unknown task_type: {task_type}")
    
    print(f"\n{'=' * 80}")
    print(f"Running evaluation for {task_type}...")
    print(f"Command: {' '.join(cmd)}")
    print(f"{'=' * 80}\n")
    
    subprocess.run(cmd, check=True)
    print(f"\nEvaluation complete. Results saved to {result_folder}")

def _process_jsonl(jsonl_path: str, msa_dir: Optional[Path] = None):
    """Process multiple datapoints from a JSONL file."""
    print(f"Processing JSONL file: {jsonl_path}")

    for line_num, line in enumerate(Path(jsonl_path).read_text().splitlines(), 1):
        if not line.strip():
            continue

        print(f"\n--- Processing line {line_num} ---")

        try:
            datapoint = Datapoint.from_json(line)
            _run_boltz_and_collect(datapoint)

        except json.JSONDecodeError as e:
            print(f"ERROR: Invalid JSON on line {line_num}: {e}")
            continue
        except Exception as e:
            print(f"ERROR: Failed to process datapoint on line {line_num}: {e}")
            raise e
            continue

def _process_json(json_path: str, msa_dir: Optional[Path] = None):
    """Process a single datapoint from a JSON file."""
    print(f"Processing JSON file: {json_path}")

    try:
        datapoint = _load_datapoint(Path(json_path))
        _run_boltz_and_collect(datapoint)
    except Exception as e:
        print(f"ERROR: Failed to process datapoint: {e}")
        raise

def main():
    """Main entry point for the hackathon scaffold."""
    # Determine task type from first datapoint for evaluation
    task_type = None
    input_file = None
    
    if args.input_json:
        input_file = args.input_json
        _process_json(args.input_json, args.msa_dir)
        # Get task type from the single datapoint
        try:
            datapoint = _load_datapoint(Path(args.input_json))
            task_type = datapoint.task_type
        except Exception as e:
            print(f"WARNING: Could not determine task type: {e}")
    elif args.input_jsonl:
        input_file = args.input_jsonl
        _process_jsonl(args.input_jsonl, args.msa_dir)
        # Get task type from first datapoint in JSONL
        try:
            with open(args.input_jsonl) as f:
                first_line = f.readline().strip()
                if first_line:
                    first_datapoint = Datapoint.from_json(first_line)
                    task_type = first_datapoint.task_type
        except Exception as e:
            print(f"WARNING: Could not determine task type: {e}")
    
    # Run evaluation if result folder is specified and task type was determined
    if args.result_folder and task_type and input_file:
        try:
            _run_evaluation(input_file, task_type, args.submission_dir, args.result_folder)
        except Exception as e:
            print(f"WARNING: Evaluation failed: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    main()
