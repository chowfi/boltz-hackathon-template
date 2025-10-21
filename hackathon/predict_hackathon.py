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


def prepare_protein_ligand(datapoint_id: str, protein: Protein, ligands: list[SmallMolecule], input_dict: dict, msa_dir: Optional[Path] = None) -> List[tuple[dict, List[str]]]:
    """
    General Heuristics + Regional Sampling Strategy
    
    Key insight: Don't rely on ground truth seqid (PDB numbering issues!)
    Instead: Sample likely allosteric regions based on general biology
    
    Strategy:
    - 3 SATURATION configs (pure exploration, 30% of samples)
    - 6 REGIONAL configs (systematic coverage, 60% of samples)
    - 1 HIGH DIVERSITY config (edge cases, 10% of samples)
    
    Total: 10 configs × 5 samples = 50 total samples
    
    Why this works:
    - No assumptions about seqid or PDB numbering
    - Covers entire protein systematically
    - Soft constraints (force=False) let Boltz refine
    - Robust to metadata issues
    """
    configs = []
    
    print(f"🔥 General Heuristics + Regional Sampling for {datapoint_id}")
    
    # Get ligand ID for constraints
    ligand_id = ligands[0].id if ligands else "B"
    protein_id = protein.id
    protein_length = len(protein.sequence)
    
    print(f"   Protein length: {protein_length} residues")
    
    # ========================================================================
    # Config 0-2: SATURATION (pure exploration, 30% of samples)
    # ========================================================================
    saturation_seeds = [42, 1000, 5000]
    
    for i, seed in enumerate(saturation_seeds):
        config_dict = input_dict.copy()
        configs.append((config_dict, ["--diffusion_samples", "5", "--seed", str(seed)]))
        print(f"  Config {i}: SATURATION seed={seed:>5} (5 samples)")
    
    # ========================================================================
    # Config 3-8: REGIONAL sampling (systematic coverage, 60% of samples)
    # Divide protein into ~6 overlapping regions
    # ========================================================================
    num_regions = 6
    region_size = max(50, protein_length // 4)  # At least 50 residues, or 25% of protein
    
    # Calculate region starts to cover the protein evenly
    if protein_length <= region_size:
        # Small protein - just sample the whole thing in one region
        region_starts = [1]
        num_regions = 1
    else:
        step = (protein_length - region_size) // (num_regions - 1) if num_regions > 1 else protein_length
        step = max(step, 1)  # At least 1 residue step
        region_starts = [1 + i * step for i in range(num_regions)]
    
    for i, start in enumerate(region_starts[:6]):  # Max 6 regions
        end = min(start + region_size - 1, protein_length)
        region_residues = list(range(start, end + 1))
        
        region_dict = input_dict.copy()
        region_dict["constraints"] = [{
            "pocket": {
                "binder": ligand_id,
                "contacts": [[protein_id, res] for res in region_residues],
                "max_distance": 6.0,
                "force": False  # SOFT constraint - let Boltz decide
            }
        }]
        
        config_idx = len(configs)
        seed = 100000 + i
        configs.append((region_dict, ["--diffusion_samples", "5", "--seed", str(seed)]))
        print(f"  Config {config_idx}: REGION {i+1} residues {start:>3}-{end:>3} (soft, 5 samples)")
    
    # ========================================================================
    # Config 9: HIGH DIVERSITY (edge cases, 10% of samples)
    # ========================================================================
    high_div_dict = input_dict.copy()
    configs.append((high_div_dict, ["--diffusion_samples", "5", "--seed", "7777777"]))
    print(f"  Config {len(configs)-1}: HIGH DIVERSITY seed=7777777 (5 samples)")
    
    # Pad with SATURATION if we have < 10 configs (for small proteins)
    while len(configs) < 10:
        extra_seed = 10000 * (len(configs) + 1)
        configs.append((input_dict.copy(), ["--diffusion_samples", "5", "--seed", str(extra_seed)]))
        print(f"  Config {len(configs)-1}: Extra SATURATION seed={extra_seed} (5 samples)")
    
    # ========================================================================
    # Total: 10 configs × 5 samples = 50 samples
    # ========================================================================
    
    print(f"✅ Generated {len(configs)} configs, {len(configs)*5} TOTAL SAMPLES")
    print(f"   - {len(saturation_seeds)} SATURATION configs (pure exploration)")
    print(f"   - {min(num_regions, 6)} REGIONAL configs (systematic coverage)")
    print(f"   - 1 HIGH DIVERSITY config (edge cases)")
    print(f"   Strategy: Cover entire protein without assuming binding location")
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
    Simple Confidence-Only Ranking
    
    After testing:
    - Multi-stage ranking didn't help (targeted wrong locations due to PDB numbering)
    - Complex scoring often picks physically "better" but spatially wrong predictions
    - Simple confidence ranking is robust and works well
    
    Strategy: Sort all samples by Boltz confidence, return top 5
    """
    # Collect all PDBs from all configurations
    all_pdbs = []
    for prediction_dir in prediction_dirs:
        config_pdbs = sorted(prediction_dir.glob(f"{datapoint.datapoint_id}_config_*_model_*.pdb"))
        all_pdbs.extend(config_pdbs)
    
    if not all_pdbs:
        print(f"Warning: No PDB files found for {datapoint.datapoint_id}")
        return []
    
    print(f"🔥 Confidence-Only Ranking: {len(all_pdbs)} predictions for {datapoint.datapoint_id}")
    
    scores = []
    for pdb_path in all_pdbs:
        try:
            # Extract Boltz confidence (raw, no normalization)
            boltz_conf = extract_confidence(pdb_path)
            
            scores.append({
                "path": pdb_path,
                "confidence": boltz_conf,
            })
        except Exception as e:
            print(f"Warning: Failed to score {pdb_path.name}: {e}")
            continue
    
    if not scores:
        print(f"Warning: No valid scores computed for {datapoint.datapoint_id}")
        return all_pdbs[:5]
    
    # Sort by raw confidence (highest first)
    scores.sort(key=lambda x: x["confidence"], reverse=True)
    
    # Print top 10 for debugging
    print(f"\n🔍 TOP 10 BY CONFIDENCE:")
    print(f"{'Rank':<5} {'Model':<40} {'Confidence':<12}")
    print("-" * 60)
    for i, s in enumerate(scores[:10], 1):
        model_name = s["path"].name[-35:]
        print(f"{i:<5} {model_name:<40} {s['confidence']:.4f}")
    
    print(f"\n✅ Selected Top-5 (confidence range: {scores[0]['confidence']:.2f} - {scores[4]['confidence']:.2f})")
    
    return [s["path"] for s in scores[:5]]


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
