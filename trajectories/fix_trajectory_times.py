#!/usr/bin/env python3
"""
Script to fix trajectory time values in concatenated trajectory file.
Each trajectory should have 29 lines with time values: 0.00, -6.00, -12.00, ..., -168.00
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "config"))
import config


def fix_trajectory_times(input_file, output_file):
    """
    Fix trajectory time values to be consistent across all trajectories.
    
    Args:
        input_file (str): Path to input file
        output_file (str): Path to output file
    """
    
    # Expected time values for each trajectory (29 time steps)
    expected_times = [0.00 + i * (-6.00) for i in range(29)]
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    fixed_lines = []
    trajectory_position = 0  # Position within current trajectory (0-28)
    
    for line in lines:
        stripped = line.strip()
        
        # Check if this is a data line (has numerical data in expected format)
        if stripped and not stripped.startswith('Reference date') and not stripped.startswith('---') and 'time' not in stripped:
            parts = stripped.split()
            
            # Check if this looks like a trajectory data line
            if len(parts) >= 11:
                try:
                    # Try to parse the first value as a float (time)
                    float(parts[0])
                    
                    # This is a data line - fix the time value
                    parts[0] = f"{expected_times[trajectory_position]:7.2f}"
                    
                    # Reconstruct the line with proper spacing
                    fixed_line = f"{parts[0]:>8s} {parts[1]:>9s} {parts[2]:>9s} {parts[3]:>6s} {parts[4]:>9s} {parts[5]:>9s} {parts[6]:>9s} {parts[7]:>9s} {parts[8]:>9s} {parts[9]:>9s} {parts[10]:>9s} {parts[11]:>9s}\n"
                    fixed_lines.append(fixed_line)
                    
                    # Move to next position in trajectory
                    trajectory_position = (trajectory_position + 1) % 29
                    
                except (ValueError, IndexError):
                    # Not a valid data line, keep as is
                    fixed_lines.append(line)
            else:
                # Not enough columns, keep as is
                fixed_lines.append(line)
        else:
            # This is a header line or empty line, keep as is
            fixed_lines.append(line)
    
    # Write the fixed data to output file
    with open(output_file, 'w') as f:
        f.writelines(fixed_lines)
    
    print(f"Fixed trajectory time values written to {output_file}")

def process_all_files():
    """
    Process all combinations of regions and timeframes.
    """
    regions = ['JS', 'NA']
    timeframes = ['present', 'future']
    airstreams = ['ALL', 'UL', 'LL']

    base_input_path = config.LAGRANTO_DIR
    
    for region in regions:
        for timeframe in timeframes:
            for airstream in airstreams:
                input_file = f"{base_input_path}/CAO_{airstream}_trajs_{region}_{timeframe}_ALL.txt"
                output_file = f"{base_input_path}/CAO_{airstream}_trajs_{region}_{timeframe}_ALL.txt.fixed"
                
                print(f"\nProcessing {region} {timeframe}...")
                print(f"Input file: {input_file}")
                print(f"Output file: {output_file}")
                
                try:
                    fix_trajectory_times(input_file, output_file)
                    
                    # Count trajectories for verification
                    with open(output_file, 'r') as f:
                        data_lines = 0
                        for line in f:
                            stripped = line.strip()
                            if stripped and not stripped.startswith('Reference date') and not stripped.startswith('---') and 'time' not in stripped:
                                parts = stripped.split()
                                if len(parts) >= 11:
                                    try:
                                        float(parts[0])
                                        data_lines += 1
                                    except ValueError:
                                        pass
                        
                        print(f"Total data lines: {data_lines}")
                        print(f"Number of trajectories: {data_lines // 29}")
                        
                except FileNotFoundError:
                    print(f"ERROR: Input file not found: {input_file}")
                except Exception as e:
                    print(f"ERROR processing {region} {timeframe}: {e}")

# Example usage
if __name__ == "__main__":
    # Process all files
    process_all_files()
    
    # Optional: Display first few lines of one fixed file to verify
    print("\n" + "="*50)
    print("Sample output from JS present:")
    print("="*50)
    
    sample_file = f"{config.LAGRANTO_DIR}/CAO_BCAR_trajs_JS_present_ALL.txt.fixed"
    try:
        with open(sample_file, 'r') as f:
            for i, line in enumerate(f):
                if i < 40:  # Show first 40 lines
                    print(line.rstrip())
                else:
                    break
    except FileNotFoundError:
        print(f"Sample file not found: {sample_file}")