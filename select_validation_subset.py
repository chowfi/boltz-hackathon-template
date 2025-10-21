#!/usr/bin/env python3
"""
Select a balanced validation subset from the ASOS dataset.
This script creates a subset with 50% orthosteric and 50% allosteric datapoints.
"""

import json
import argparse
import random
from pathlib import Path

def select_validation_subset(input_file: str, output_file: str, size: int = 10):
    """
    Select a balanced subset of datapoints for validation.
    
    Args:
        input_file: Path to the input JSONL file
        output_file: Path to the output JSONL file
        size: Number of datapoints to select (default: 10)
    """
    # Read all datapoints
    datapoints = []
    with open(input_file, 'r') as f:
        for line in f:
            if line.strip():
                datapoints.append(json.loads(line))
    
    # Separate by binding type
    orthosteric = [dp for dp in datapoints if dp['ground_truth']['ligand_types'][0]['type'] == 'orthosteric']
    allosteric = [dp for dp in datapoints if dp['ground_truth']['ligand_types'][0]['type'] == 'allosteric']
    
    print(f"Found {len(orthosteric)} orthosteric and {len(allosteric)} allosteric datapoints")
    
    # Select balanced subset
    half_size = size // 2
    selected_orthosteric = random.sample(orthosteric, min(half_size, len(orthosteric)))
    selected_allosteric = random.sample(allosteric, min(half_size, len(allosteric)))
    
    # If we need more datapoints, fill from the larger group
    remaining = size - len(selected_orthosteric) - len(selected_allosteric)
    if remaining > 0:
        if len(orthosteric) > len(allosteric):
            additional = random.sample([dp for dp in orthosteric if dp not in selected_orthosteric], 
                                     min(remaining, len(orthosteric) - len(selected_orthosteric)))
            selected_orthosteric.extend(additional)
        else:
            additional = random.sample([dp for dp in allosteric if dp not in selected_allosteric], 
                                     min(remaining, len(allosteric) - len(selected_allosteric)))
            selected_allosteric.extend(additional)
    
    # Combine and shuffle
    selected = selected_orthosteric + selected_allosteric
    random.shuffle(selected)
    
    # Write to output file
    with open(output_file, 'w') as f:
        for datapoint in selected:
            f.write(json.dumps(datapoint) + '\n')
    
    print(f"Selected {len(selected)} datapoints:")
    print(f"  Orthosteric: {len(selected_orthosteric)}")
    print(f"  Allosteric: {len(selected_allosteric)}")
    print(f"Saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Select a balanced validation subset from ASOS dataset')
    parser.add_argument('--input', required=True, help='Input JSONL file')
    parser.add_argument('--output', required=True, help='Output JSONL file')
    parser.add_argument('--size', type=int, default=10, help='Number of datapoints to select (default: 10)')
    
    args = parser.parse_args()
    
    # Set random seed for reproducibility
    random.seed(42)
    
    select_validation_subset(args.input, args.output, args.size)

if __name__ == "__main__":
    main()
