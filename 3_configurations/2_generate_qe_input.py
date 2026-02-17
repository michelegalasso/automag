#!/usr/bin/env python3
"""
Quantum ESPRESSO Magnetic Configuration Generator (Fully Generalized)

This script generates pw.x input files for different magnetic configurations
by reading a template file and a configuration list, then automatically
assigning atom labels based on distinct magnetic moments AND Hubbard U values.

Usage:
    python generate_magnetic_configs.py <template_file> <configurations_file>

Example:
    python generate_magnetic_configs.py pw_fim8.in configurations.txt
"""

import re
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Optional


def parse_configurations(config_file: str) -> Dict[str, List[float]]:
    """
    Parse the configurations file.
    
    Args:
        config_file: Path to configurations.txt file
        
    Returns:
        Dictionary mapping configuration names to lists of magnetic moments
    """
    configurations = {}
    
    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split()
            config_name = parts[0]
            magnetic_moments = [float(x) for x in parts[1:]]
            configurations[config_name] = magnetic_moments
    
    return configurations


def parse_template(template_file: str) -> Dict:
    """
    Parse the template input file preserving exact formatting.
    Automatically detects magnetic atoms from starting_magnetization.
    
    Args:
        template_file: Path to template pw.x input file
        
    Returns:
        Dictionary containing parsed template sections
    """
    with open(template_file, 'r') as f:
        lines = f.readlines()
    
    template = {'lines': lines}
    
    # Find section boundaries
    section_starts = {}
    section_ends = {}
    
    current_namelist = None
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        
        # Check for namelist starts
        if stripped.startswith('&CONTROL'):
            section_starts['CONTROL'] = i
            current_namelist = 'CONTROL'
        elif stripped.startswith('&SYSTEM'):
            if current_namelist:
                section_ends[current_namelist] = i - 1
            section_starts['SYSTEM'] = i
            current_namelist = 'SYSTEM'
        elif stripped.startswith('&ELECTRONS'):
            if current_namelist:
                section_ends[current_namelist] = i - 1
            section_starts['ELECTRONS'] = i
            current_namelist = 'ELECTRONS'
        elif stripped.startswith('ATOMIC_SPECIES'):
            if current_namelist:
                section_ends[current_namelist] = i - 1
            section_starts['ATOMIC_SPECIES'] = i
            current_namelist = None
        elif stripped.startswith('HUBBARD'):
            if 'ATOMIC_SPECIES' in section_starts and 'ATOMIC_SPECIES' not in section_ends:
                section_ends['ATOMIC_SPECIES'] = i - 1
            section_starts['HUBBARD'] = i
        elif stripped.startswith('ATOMIC_POSITIONS'):
            for sec in ['ATOMIC_SPECIES', 'HUBBARD']:
                if sec in section_starts and sec not in section_ends:
                    section_ends[sec] = i - 1
            section_starts['ATOMIC_POSITIONS'] = i
        elif stripped.startswith('K_POINTS'):
            if 'ATOMIC_POSITIONS' in section_starts:
                section_ends['ATOMIC_POSITIONS'] = i - 1
            section_starts['K_POINTS'] = i
    
    # Set end for last section
    if 'K_POINTS' in section_starts:
        section_ends['K_POINTS'] = len(lines) - 1
    
    template['section_starts'] = section_starts
    template['section_ends'] = section_ends
    
    # Parse ATOMIC_SPECIES to get element-to-pseudopotential mapping
    species_info = {}  # {label: (element, mass, pseudopotential)}
    if 'ATOMIC_SPECIES' in section_starts:
        start = section_starts['ATOMIC_SPECIES']
        end = section_ends.get('ATOMIC_SPECIES', section_starts.get('HUBBARD', len(lines)) - 1)
        
        for i in range(start + 1, end + 1):
            line = lines[i].strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) >= 3:
                label = parts[0]
                mass = parts[1]
                pseudo = parts[2]
                # Extract base element (e.g., "Tb1" -> "Tb", "Fe2" -> "Fe")
                element = re.match(r'([A-Z][a-z]?)', label).group(1)
                species_info[label] = (element, mass, pseudo)
    
    template['species_info'] = species_info
    
    # Parse HUBBARD section to get U parameters by position
    hubbard_info = {}  # {label: [(orbital, U_value), ...]}
    if 'HUBBARD' in section_starts:
        start = section_starts['HUBBARD']
        end = section_ends.get('HUBBARD', section_starts.get('ATOMIC_POSITIONS', len(lines)) - 1)
        
        for i in range(start + 1, end + 1):
            line = lines[i].strip()
            if not line or line.startswith('#') or '{' in line:
                continue
            
            # Parse lines like: "U   Tb1-4f   4.4785"
            match = re.match(r'U\s+(\w+)-(\w+)\s+([\d.]+)', line)
            if match:
                label = match.group(1)
                orbital = match.group(2)
                u_value = float(match.group(3))
                
                if label not in hubbard_info:
                    hubbard_info[label] = []
                hubbard_info[label].append((orbital, u_value))
    
    template['hubbard_info'] = hubbard_info
    
    # Determine which atom types are magnetic from starting_magnetization
    magnetic_type_indices = set()
    system_start = section_starts.get('SYSTEM', 0)
    system_end = section_ends.get('SYSTEM', 0)
    
    for i in range(system_start, system_end + 1):
        line = lines[i]
        mag_match = re.search(r'starting_magnetization\((\d+)\)\s*=\s*([-\d.]+)', line)
        if mag_match:
            type_idx = int(mag_match.group(1))
            mag_value = float(mag_match.group(2))
            if mag_value != 0.0:  # Only non-zero magnetizations
                magnetic_type_indices.add(type_idx)
    
    template['magnetic_type_indices'] = magnetic_type_indices
    
    # Parse atomic positions to identify magnetic vs non-magnetic atoms
    magnetic_atoms = []  # [(element, label, coords, orig_line, position_index)]
    non_magnetic_atoms = []  # [(label, coords, orig_line)]
    
    # First, build a mapping from type index to labels in ATOMIC_SPECIES
    type_idx_to_labels = {}  # {1: ['Tb1'], 2: ['Tb2'], ...}
    species_labels = list(species_info.keys())
    for idx, label in enumerate(species_labels, 1):
        if idx not in type_idx_to_labels:
            type_idx_to_labels[idx] = []
        type_idx_to_labels[idx].append(label)
    
    # Determine which labels are magnetic
    magnetic_labels = set()
    for type_idx in magnetic_type_indices:
        if type_idx in type_idx_to_labels:
            magnetic_labels.update(type_idx_to_labels[type_idx])
    
    if 'ATOMIC_POSITIONS' in section_starts:
        start = section_starts['ATOMIC_POSITIONS']
        end = section_ends.get('ATOMIC_POSITIONS', len(lines) - 1)
        
        position_index = 0
        for i in range(start + 1, end + 1):
            line = lines[i].strip()
            if not line:
                continue
            
            parts = line.split()
            if len(parts) >= 4:
                atom_label = parts[0]
                coords = parts[1:4]
                
                # Determine base element
                element_match = re.match(r'([A-Z][a-z]?)', atom_label)
                if element_match:
                    element = element_match.group(1)
                else:
                    element = atom_label
                
                # Check if this label is magnetic
                is_magnetic = atom_label in magnetic_labels
                
                if is_magnetic:
                    magnetic_atoms.append((element, atom_label, coords, lines[i], position_index))
                    position_index += 1
                else:
                    non_magnetic_atoms.append((atom_label, coords, lines[i]))
    
    template['magnetic_atoms'] = magnetic_atoms
    template['non_magnetic_atoms'] = non_magnetic_atoms
    
    # Extract prefix
    control_start = section_starts.get('CONTROL', 0)
    control_end = section_ends.get('CONTROL', 0)
    for i in range(control_start, control_end + 1):
        match = re.search(r"prefix\s*=\s*['\"]([^'\"]+)['\"]", lines[i])
        if match:
            template['prefix'] = match.group(1)
            break
    
    return template


def get_hubbard_params_for_position(template: Dict, element: str, 
                                     position_index: int) -> List[Tuple[str, float]]:
    """
    Get Hubbard U parameters for a specific element at a specific position.
    
    Args:
        template: Parsed template dictionary
        element: Element symbol (e.g., 'Tb', 'Fe')
        position_index: Position index of the atom (0-based)
        
    Returns:
        List of (orbital, U_value) tuples
    """
    # Find which original label was used at this position
    hubbard_info = template['hubbard_info']
    
    # Get the original label used at this position
    for elem, orig_label, coords, orig_line, pos_idx in template['magnetic_atoms']:
        if elem == element and pos_idx == position_index:
            if orig_label in hubbard_info:
                return hubbard_info[orig_label]
            break
    
    return []


def generate_input_file(config_name: str, magnetic_moments: List[float], 
                        template: Dict, output_dir: str = '.'):
    """
    Generate a pw.x input file for a specific magnetic configuration.
    
    Args:
        config_name: Name of the configuration
        magnetic_moments: List of magnetic moments for each magnetic atom
        template: Parsed template dictionary
        output_dir: Directory to save the output file
    """
    lines = template['lines']
    section_starts = template['section_starts']
    section_ends = template['section_ends']
    species_info = template['species_info']
    
    # Build atom type assignments considering BOTH magnetic moment AND Hubbard U
    # The key insight: Group atoms by (element, Hubbard U), then by magnetic moment
    
    # First, determine the Hubbard U for each position based on original label
    position_to_hubbard = {}  # {position_index: hubbard_params_tuple}
    for i, (element, orig_label, coords, orig_line, pos_idx) in enumerate(template['magnetic_atoms']):
        if orig_label in template['hubbard_info']:
            position_to_hubbard[i] = tuple(template['hubbard_info'][orig_label])
        else:
            position_to_hubbard[i] = ()
    
    # Now create type mapping based on (element, hubbard_params, moment)
    type_mapping = {}  # {(element, hubbard_params, moment): (element, new_type_num, hubbard_params)}
    atom_type_assignments = []  # [(element, type_num, hubbard_params) for each atom]
    next_type_num = {}  # {element: next_available_number}
    
    # Process each atom
    for i, (element, orig_label, coords, orig_line, pos_idx) in enumerate(template['magnetic_atoms']):
        moment = magnetic_moments[i]
        hubbard_params = position_to_hubbard.get(i, ())
        
        # Create a unique key combining element, hubbard params, and moment
        type_key = (element, hubbard_params, moment)
        
        if type_key not in type_mapping:
            # Assign a new type number for this element
            if element not in next_type_num:
                next_type_num[element] = 1
            type_mapping[type_key] = (element, next_type_num[element], hubbard_params)
            next_type_num[element] += 1
        
        atom_type_assignments.append(type_mapping[type_key])
    
    # Count total types
    num_magnetic_types = len(type_mapping)
    num_non_magnetic_types = len(set(label for label, _, _ in template['non_magnetic_atoms']))
    total_ntyp = num_magnetic_types + num_non_magnetic_types
    
    # Build reverse mapping for magnetization values
    type_to_moment = {}  # {(element, type_num): moment}
    type_to_hubbard = {}  # {(element, type_num): hubbard_params}
    for (element, hubbard_params, moment), (elem, type_num, hub_params) in type_mapping.items():
        type_to_moment[(elem, type_num)] = moment
        type_to_hubbard[(elem, type_num)] = hub_params
    
    # Sort types by element and type number for consistent ordering
    sorted_magnetic_types = sorted(set((elem, type_num) for _, (elem, type_num, _) in type_mapping.items()))
    
    output_lines = []
    
    # CONTROL section - just update prefix
    control_start = section_starts['CONTROL']
    control_end = section_ends['CONTROL']
    for i in range(control_start, control_end + 1):
        line = lines[i]
        if 'prefix' in line:
            line = re.sub(r"prefix\s*=\s*['\"][^'\"]+['\"]", f"prefix='{config_name}'", line)
        output_lines.append(line)
    
    output_lines.append('\n')
    
    # SYSTEM section - update ntyp and magnetization
    system_start = section_starts['SYSTEM']
    system_end = section_ends['SYSTEM']
    
    system_lines = []
    for i in range(system_start, system_end + 1):
        line = lines[i]
        
        # Update ntyp
        if 'ntyp' in line and not line.strip().startswith('!'):
            line = re.sub(r'ntyp\s*=\s*\d+', f'ntyp = {total_ntyp}', line)
        
        # Skip existing starting_magnetization lines
        if 'starting_magnetization' not in line:
            system_lines.append(line)
    
    # Find the closing / and insert magnetization lines before it
    for i in range(len(system_lines) - 1, -1, -1):
        if '/' in system_lines[i] and not system_lines[i].strip().startswith('!'):
            # Build magnetization lines in order of type index
            mag_lines = []
            
            # Create ordered list by type index
            type_index_map = {}  # {type_idx: moment}
            current_idx = 1
            for elem, type_num in sorted_magnetic_types:
                moment = type_to_moment[(elem, type_num)]
                type_index_map[current_idx] = moment
                current_idx += 1
            
            # Generate magnetization lines
            for type_idx in sorted(type_index_map.keys()):
                moment = type_index_map[type_idx]
                mag_lines.append(f'    starting_magnetization({type_idx}) = {moment:5.1f}')
            
            # Insert before closing
            system_lines.insert(i, ',\n'.join(mag_lines) + '\n')
            break
    
    output_lines.extend(system_lines)
    output_lines.append('\n')
    
    # ELECTRONS section
    electrons_start = section_starts['ELECTRONS']
    electrons_end = section_ends['ELECTRONS']
    for i in range(electrons_start, electrons_end + 1):
        output_lines.append(lines[i])
    
    output_lines.append('\n')
    
    # ATOMIC_SPECIES section
    output_lines.append('ATOMIC_SPECIES\n')
    output_lines.append('# the second field, atomic mass, is not actually used\n')
    output_lines.append('# except for MD calculations\n')
    
    # Add magnetic species in order
    for elem, type_num in sorted_magnetic_types:
        # Get mass and pseudopotential from original species
        # Find any label with this element in species_info
        mass = None
        pseudo = None
        for orig_label, (orig_elem, orig_mass, orig_pseudo) in species_info.items():
            if orig_elem == elem:
                mass = orig_mass
                pseudo = orig_pseudo
                break
        
        new_label = f'{elem}{type_num}'
        output_lines.append(f'    {new_label:<4s} {mass:>8s}   {pseudo}\n')
    
    # Add non-magnetic species
    non_magnetic_labels = sorted(set(label for label, _, _ in template['non_magnetic_atoms']))
    for label in non_magnetic_labels:
        # Get info from original species
        for orig_label, (elem, mass, pseudo) in species_info.items():
            if orig_label == label:
                output_lines.append(f'    {label:<4s} {mass:>8s}   {pseudo}\n')
                break
    
    output_lines.append('\n')
    
    # HUBBARD section if present
    if 'HUBBARD' in section_starts:
        hubbard_start = section_starts['HUBBARD']
        # Get the header line
        output_lines.append(lines[hubbard_start])
        
        # Generate U parameters for all magnetic types based on their Hubbard parameters
        for elem, type_num in sorted_magnetic_types:
            new_label = f'{elem}{type_num}'
            hubbard_params = type_to_hubbard.get((elem, type_num), ())
            for orbital, u_value in hubbard_params:
                output_lines.append(f'    U   {new_label}-{orbital}   {u_value}\n')
        
        output_lines.append('\n')
    
    # ATOMIC_POSITIONS section
    output_lines.append(lines[section_starts['ATOMIC_POSITIONS']])
    
    # Add magnetic atoms with appropriate labels
    for i, (element, orig_label, coords, orig_line, pos_idx) in enumerate(template['magnetic_atoms']):
        elem, type_num, hubbard_params = atom_type_assignments[i]
        new_label = f'{elem}{type_num}'
        coord_str = '    '.join(coords)
        output_lines.append(f'{new_label:<20s} {coord_str}\n')
    
    # Add non-magnetic atoms
    for label, coords, orig_line in template['non_magnetic_atoms']:
        output_lines.append(orig_line)
    
    output_lines.append('\n')
    
    # K_POINTS section
    kpoints_start = section_starts['K_POINTS']
    kpoints_end = section_ends['K_POINTS']
    for i in range(kpoints_start, kpoints_end + 1):
        output_lines.append(lines[i])
    
    # Write to file
    output_path = Path(output_dir) / f'pw_{config_name}.in'
    with open(output_path, 'w') as f:
        f.writelines(output_lines)
    
    print(f"Generated: {output_path}")


def main():
    """Main function to orchestrate the generation process."""
    if len(sys.argv) != 3:
        print("Usage: python generate_magnetic_configs.py <template_file> <configurations_file>")
        print("\nExample:")
        print("  python generate_magnetic_configs.py pw_fim8.in configurations.txt")
        sys.exit(1)
    
    template_file = sys.argv[1]
    config_file = sys.argv[2]
    
    # Validate input files
    if not Path(template_file).exists():
        print(f"Error: Template file '{template_file}' not found.")
        sys.exit(1)
    
    if not Path(config_file).exists():
        print(f"Error: Configuration file '{config_file}' not found.")
        sys.exit(1)
    
    print("=" * 70)
    print("Quantum ESPRESSO Magnetic Configuration Generator")
    print("=" * 70)
    print(f"\nTemplate file: {template_file}")
    print(f"Configuration file: {config_file}\n")
    
    # Parse inputs
    print("Parsing template file...")
    template = parse_template(template_file)
    
    print("Parsing configurations...")
    configurations = parse_configurations(config_file)
    
    print(f"\nFound {len(configurations)} configurations")
    print(f"Number of magnetic atoms: {len(template['magnetic_atoms'])}")
    print(f"Number of non-magnetic atoms: {len(template['non_magnetic_atoms'])}")
    
    # Show detected magnetic elements
    magnetic_elements = set(elem for elem, _, _, _, _ in template['magnetic_atoms'])
    print(f"Magnetic elements: {', '.join(sorted(magnetic_elements))}")
    
    # Show Hubbard info if present
    if template['hubbard_info']:
        print(f"Hubbard U corrections detected: {len(template['hubbard_info'])} atom types")
    
    print()
    
    # Generate input files
    print("Generating input files...")
    print("-" * 70)
    
    for config_name, magnetic_moments in configurations.items():
        if len(magnetic_moments) != len(template['magnetic_atoms']):
            print(f"Warning: Configuration '{config_name}' has {len(magnetic_moments)} moments "
                  f"but template has {len(template['magnetic_atoms'])} magnetic atoms. Skipping.")
            continue
        
        generate_input_file(config_name, magnetic_moments, template)
    
    print("-" * 70)
    print(f"\nSuccessfully generated {len(configurations)} input files!")
    print("\nDone!")


if __name__ == '__main__':
    main()
