#!/usr/bin/env python3
"""
Quantum ESPRESSO Magnetic Configuration Generator

This script generates pw.x input files for different magnetic configurations
by reading a template file and a configuration list, then automatically
assigning atom labels based on distinct magnetic moments.

Usage:
    python generate_magnetic_configs.py <template_file> <configurations_file>

Example:
    python generate_magnetic_configs.py pw_fm1.in configurations.txt
"""

import re
import sys
from pathlib import Path
from typing import List, Dict, Tuple


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

    in_namelist = False
    current_namelist = None

    for i, line in enumerate(lines):
        stripped = line.strip()

        # Check for namelist starts
        if stripped.startswith('&CONTROL'):
            section_starts['CONTROL'] = i
            in_namelist = True
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
            # Close previous section
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

    # Parse atomic positions to identify magnetic atoms
    magnetic_atoms = []
    non_magnetic_atoms = []

    if 'ATOMIC_POSITIONS' in section_starts:
        start = section_starts['ATOMIC_POSITIONS']
        end = section_ends.get('ATOMIC_POSITIONS', len(lines) - 1)

        for i in range(start + 1, end + 1):
            line = lines[i].strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) >= 4:
                atom_label = parts[0]
                coords = parts[1:4]

                # Determine if this is a magnetic atom (Tb in this case)
                if re.match(r'Tb\d*', atom_label):
                    magnetic_atoms.append(('Tb', coords, lines[i]))
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

    # Identify unique magnetic moments and create atom type mapping
    unique_moments = []
    atom_type_assignments = []

    for moment in magnetic_moments:
        if moment not in unique_moments:
            unique_moments.append(moment)
        type_index = unique_moments.index(moment) + 1
        atom_type_assignments.append(type_index)

    num_magnetic_types = len(unique_moments)
    num_non_magnetic_types = len(set(atom[0] for atom in template['non_magnetic_atoms']))
    total_ntyp = num_magnetic_types + num_non_magnetic_types

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
    closing_found = False
    for i in range(len(system_lines) - 1, -1, -1):
        if '/' in system_lines[i] and not system_lines[i].strip().startswith('!'):
            # Build magnetization lines
            mag_lines = []
            for j, moment in enumerate(unique_moments, 1):
                mag_lines.append(f'    starting_magnetization({j}) = {moment:5.1f}')

            # Insert magnetization lines before closing /
            system_lines.insert(i, ',\n'.join(mag_lines) + '\n')
            closing_found = True
            break

    if not closing_found:
        # If no closing / found, add magnetization lines at the end
        mag_lines = []
        for j, moment in enumerate(unique_moments, 1):
            mag_lines.append(f'    starting_magnetization({j}) = {moment:5.1f}\n')
        system_lines.extend(mag_lines)

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

    for i in range(1, num_magnetic_types + 1):
        output_lines.append(f'   Tb{i}  55.845   Tb.paw.z_19.atompaw.wentzcovitch.v1.2.upf\n')

    # Add non-magnetic species
    non_magnetic_types = set(atom[0] for atom in template['non_magnetic_atoms'])
    for atom_type in sorted(non_magnetic_types):
        output_lines.append(f'   {atom_type}    15.999   H_ONCV_PBE-1.0.oncvpsp.upf\n')

    output_lines.append('\n')

    # HUBBARD section if present
    if 'HUBBARD' in section_starts:
        hubbard_start = section_starts['HUBBARD']
        # Get the header line
        output_lines.append(lines[hubbard_start])

        # Generate U parameters for all Tb types
        for i in range(1, num_magnetic_types + 1):
            output_lines.append(f'   U   Tb{i}-3d   5.0657\n')

        output_lines.append('\n')

    # ATOMIC_POSITIONS section
    output_lines.append(lines[section_starts['ATOMIC_POSITIONS']])

    # Add magnetic atoms with appropriate labels
    for i, (base_element, coords, orig_line) in enumerate(template['magnetic_atoms']):
        type_num = atom_type_assignments[i]
        atom_label = f'Tb{type_num}'
        coord_str = '        '.join(coords)
        output_lines.append(f'{atom_label:<20s} {coord_str}\n')

    # Add non-magnetic atoms
    for atom_label, coords, orig_line in template['non_magnetic_atoms']:
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
        print("  python generate_magnetic_configs.py pw_fm1.in configurations.txt")
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
    print(f"Number of non-magnetic atoms: {len(template['non_magnetic_atoms'])}\n")

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
