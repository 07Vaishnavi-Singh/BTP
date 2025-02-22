def parse_xyz_file(filepath):
    """
    Parse XYZ file and return atom names and coordinates
    """
    atoms = []
    coordinates = []
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
        # Skip first two lines (number of atoms and comment)
        for line in lines[2:]:
            if line.strip():  # Skip empty lines
                parts = line.split()
                if len(parts) >= 4:
                    atom = parts[0]
                    x, y, z = map(float, parts[1:4])
                    atoms.append(atom)
                    coordinates.append((x, y, z))
    
    return atoms, coordinates

def find_center_atom(atoms, coordinates):
    """
    Find the atom that is most central (has smallest maximum distance to any other atom)
    """
    num_atoms = len(atoms)
    if num_atoms == 0:
        return None

    # For each atom, calculate its maximum distance to any other atom
    min_max_distance = float('inf')
    center_atom_index = None

    for i in range(num_atoms):
        max_distance = 0
        x1, y1, z1 = coordinates[i]
        
        # Calculate distances to all other atoms
        for j in range(num_atoms):
            if i != j:
                x2, y2, z2 = coordinates[j]
                distance = ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)**0.5
                max_distance = max(max_distance, distance)
        
        # If this atom has a smaller maximum distance, it's more central
        if max_distance < min_max_distance:
            min_max_distance = max_distance
            center_atom_index = i

    if center_atom_index is not None:
        return {
            'center_atom': atoms[center_atom_index],
            'center_coordinates': coordinates[center_atom_index],
            'max_distance': min_max_distance,
            'num_atoms': num_atoms
        }
    return None

def find_nearest_pb_atom(center_coords, atoms, coordinates):
    """
    Find the nearest Pb atom to the center coordinates
    """
    min_distance = float('inf')
    nearest_pb_index = None
    center_x, center_y, center_z = center_coords

    for i in range(len(atoms)):
        if atoms[i] == 'Pb':
            x, y, z = coordinates[i]
            distance = ((x - center_x)**2 + 
                       (y - center_y)**2 + 
                       (z - center_z)**2)**0.5
            if distance < min_distance:
                min_distance = distance
                nearest_pb_index = i
    
    if nearest_pb_index is not None:
        return {
            'atom': atoms[nearest_pb_index],
            'coordinates': coordinates[nearest_pb_index],
            'distance': min_distance
        }
    return None

# Usage example
try:
    # Parse the file
    atoms, coordinates = parse_xyz_file('vesta_file.xyz')
    
    # Find the center atom
    result = find_center_atom(atoms, coordinates)
    
    if result:
        print("\nCenter Atom Analysis:")
        print(f"Center atom type: {result['center_atom']}")
        print(f"Center coordinates: X: {result['center_coordinates'][0]:.6f}, "
              f"Y: {result['center_coordinates'][1]:.6f}, "
              f"Z: {result['center_coordinates'][2]:.6f}")
        print(f"Maximum distance to any other atom: {result['max_distance']:.6f}")
        print(f"Total atoms: {result['num_atoms']}")
        
        if result['center_atom'] != "Pb":
            print(f"\nFinding nearest Pb atom to the center atom...")
            nearest_pb = find_nearest_pb_atom(
                result['center_coordinates'],
                atoms,
                coordinates
            )
            
            if nearest_pb:
                print("\nNearest Pb Atom Found:")
                pb_coords = nearest_pb['coordinates']
                print(f"Coordinates: X: {pb_coords[0]:.6f}, Y: {pb_coords[1]:.6f}, Z: {pb_coords[2]:.6f}")
                print(f"Distance from center atom: {nearest_pb['distance']:.6f}")
            else:
                print("No Pb atoms found in the structure")
    else:
        print("No results found from coordinate calculation")
except Exception as e:
    print(f"An error occurred: {e}")

    