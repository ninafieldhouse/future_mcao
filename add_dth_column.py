#!/usr/bin/env python3
"""
Script to add dTH column to trajectory file.
dTH = TH(x) - TH(x+1) where x+1 is the next time step in the trajectory.
At time -168.00 (first timestep of backward trajectory), dTH is always 0.00.
"""

def add_dth_column(input_file, output_file):
    """
    Add dTH and dTH_24 columns to trajectory data.
    
    Args:
        input_file (str): Path to input file
        output_file (str): Path to output file
    """
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    fixed_lines = []
    trajectory_data = []  # Store current trajectory data
    trajectory_position = 0  # Position within current trajectory (0-28)
    
    for line in lines:
        stripped = line.strip()
        
        # Check if this is a data line (has numerical data in expected format)
        if stripped and not stripped.startswith('Reference date') and not stripped.startswith('---') and 'time' not in stripped:
            parts = stripped.split()
            
            # Check if this looks like a trajectory data line
            if len(parts) >= 12:  # Now expecting 12 columns (including TH)
                try:
                    # Try to parse the first value as a float (time)
                    time_val = float(parts[0])
                    th_val = float(parts[11])  # TH is the 12th column (index 11)
                    
                    # Store this line's data
                    trajectory_data.append({
                        'line': line,
                        'parts': parts,
                        'time': time_val,
                        'th': th_val,
                        'position': trajectory_position
                    })
                    
                    trajectory_position = (trajectory_position + 1) % 29
                    
                    # If we've completed a trajectory (29 lines), process it
                    if trajectory_position == 0:
                        # Calculate dTH and dTH_24 for this trajectory and write it out
                        for i, data in enumerate(trajectory_data):
                            if i == 28:  # Last line (time = -168.00)
                                dth = 0.00
                            else:
                                # dTH = TH(current) - TH(next)
                                dth = data['th'] - trajectory_data[i + 1]['th']
                            
                            # Calculate dTH_24 = TH(x-4) - TH(x)
                            # Only calculate for positions 4, 8, 12, 16, 20, 24, 28 (every 4th line starting from 5th)
                            if i >= 4 and (i - 4) % 4 == 0:
                                # dTH_24 = TH(x-4) - TH(x)
                                dth_24 = trajectory_data[i - 4]['th'] - data['th']
                            else:
                                dth_24 = 0.00
                            
                            # Reconstruct the line with dTH and dTH_24 columns added
                            parts = data['parts']
                            new_line = f"{parts[0]:>8s} {parts[1]:>9s} {parts[2]:>9s} {parts[3]:>6s} {parts[4]:>9s} {parts[5]:>9s} {parts[6]:>9s} {parts[7]:>9s} {parts[8]:>9s} {parts[9]:>9s} {parts[10]:>9s} {parts[11]:>9s} {dth:>9.3f} {dth_24:>9.3f}\n"
                            fixed_lines.append(new_line)
                        
                        # Clear trajectory data for next trajectory
                        trajectory_data = []
                    
                except (ValueError, IndexError):
                    # Not a valid data line, keep as is
                    fixed_lines.append(line)
            else:
                # Not enough columns, keep as is
                fixed_lines.append(line)
        else:
            # This is a header line or empty line
            if 'time' in stripped and 'TH' in stripped:
                # This is the column header line - add dTH and dTH_24 to it
                if stripped.endswith('TH'):
                    new_header = stripped + '      dTH    dTH_24\n'
                else:
                    new_header = line.rstrip() + '      dTH    dTH_24\n'
                fixed_lines.append(new_header)
            elif stripped.startswith('---'):
                # This is the separator line - extend it
                new_separator = stripped + '--------------------\n'
                fixed_lines.append(new_separator)
            else:
                # Keep other lines as is
                fixed_lines.append(line)
    
    # Handle any remaining trajectory data (if file doesn't end with complete trajectory)
    if trajectory_data:
        for i, data in enumerate(trajectory_data):
            if i == len(trajectory_data) - 1:  # Last available line
                dth = 0.00
            else:
                # dTH = TH(current) - TH(next)
                dth = data['th'] - trajectory_data[i + 1]['th']
            
            # Calculate dTH_24 = TH(x-4) - TH(x)
            if i >= 4 and (i - 4) % 4 == 0:
                dth_24 = trajectory_data[i - 4]['th'] - data['th']
            else:
                dth_24 = 0.00
            
            # Reconstruct the line with dTH and dTH_24 columns added
            parts = data['parts']
            new_line = f"{parts[0]:>8s} {parts[1]:>9s} {parts[2]:>9s} {parts[3]:>6s} {parts[4]:>9s} {parts[5]:>9s} {parts[6]:>9s} {parts[7]:>9s} {parts[8]:>9s} {parts[9]:>9s} {parts[10]:>9s} {parts[11]:>9s} {dth:>9.3f} {dth_24:>9.3f}\n"
            fixed_lines.append(new_line)
    
    # Write the fixed data to output file
    with open(output_file, 'w') as f:
        f.writelines(fixed_lines)
    
    print(f"Added dTH and dTH_24 columns to {output_file}")

def process_all_files():
    """
    Process all combinations of regions and timeframes.
    """
    regions = ['JS', 'NA']
    timeframes = ['present', 'future']
    airstreams = ['ALL', 'BCAR', 'LL', 'DI', 'nonBCAR']
    airstreams = ['ALL', 'UL', 'LL']
    
    base_path = "/net/litho/atmosdyn2/nfieldhouse/data/lagranto"
    
    for region in regions:
        for timeframe in timeframes:
            for airstream in airstreams:
                input_file = f"{base_path}/CAO_{airstream}_trajs_{region}_{timeframe}_ALL.txt.fixed"
                output_file = f"{base_path}/CAO_{airstream}_trajs_{region}_{timeframe}_ALL_with_dTH.txt"
                
                print(f"\nProcessing {region} {timeframe}...")
                print(f"Input file: {input_file}")
                print(f"Output file: {output_file}")
                
                try:
                    add_dth_column(input_file, output_file)
                    
                    # Count trajectories for verification
                    with open(output_file, 'r') as f:
                        data_lines = 0
                        for line in f:
                            stripped = line.strip()
                            if stripped and not stripped.startswith('Reference date') and not stripped.startswith('---') and 'time' not in stripped:
                                parts = stripped.split()
                                if len(parts) >= 14:  # Now should have 13 columns including dTH
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
    # # Test on one file first
    # test_input = "/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_trajs_JS_present_ALL.txt.fixed"
    # test_output = "/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_trajs_JS_present_ALL_with_dTH_test.txt"
    
    # print("Testing on one file...")
    # print(f"Input: {test_input}")
    # print(f"Output: {test_output}")
    
    # try:
    #     add_dth_column(test_input, test_output)
        
    #     # Show sample output
    #     print("\nFirst 20 lines of output:")
    #     with open(test_output, 'r') as f:
    #         for i, line in enumerate(f):
    #             if i < 20:
    #                 print(line.rstrip())
    #             else:
    #                 break
                    
    # except FileNotFoundError:
    #     print(f"Test file not found: {test_input}")
    # except Exception as e:
    #     print(f"Error: {e}")


    # Process all files
    process_all_files()
    
    # Optional: Display first few lines of one output file to verify
    print("\n" + "="*50)
    print("Sample output from JS present:")
    print("="*50)
    
    sample_file = "/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_trajs_JS_present_ALL_with_dTH.txt"
    try:
        with open(sample_file, 'r') as f:
            for i, line in enumerate(f):
                if i < 40:  # Show first 40 lines
                    print(line.rstrip())
                else:
                    break
    except FileNotFoundError:
        print(f"Sample file not found: {sample_file}")