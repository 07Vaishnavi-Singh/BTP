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

# @dev function to find the 6 nearest I atoms to the centre Pb atom 
def find_nearest_I_atoms(central_pb_coords, atoms, coordinates, num_neighbors=6):
    """
    Find the nearest I (iodine) atoms to the central Pb atom
    
    Parameters:
    -----------
    central_pb_coords : tuple
        (x, y, z) coordinates of the central Pb atom
    atoms : list
        List of atom symbols
    coordinates : list
        List of (x, y, z) coordinates for each atom
    num_neighbors : int, optional
        Number of nearest I atoms to find (default: 6)
        
    Returns:
    --------
    list
        List of dictionaries containing information about the nearest I atoms
    """
    # Calculate distances from central Pb to all I atoms
    distances = []
    central_x, central_y, central_z = central_pb_coords
    
    for i in range(len(atoms)):
        if atoms[i] == 'I':  # Only consider I atoms
            x, y, z = coordinates[i]
            distance = ((x - central_x)**2 + 
                       (y - central_y)**2 + 
                       (z - central_z)**2)**0.5
            
            distances.append({
                'index': i,
                'atom': atoms[i],
                'coordinates': coordinates[i],
                'distance': distance
            })
    
    # Sort by distance
    distances.sort(key=lambda x: x['distance'])
    
    # Return the num_neighbors nearest I atoms
    return distances[:num_neighbors]


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


def main():
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
                    
                    # Find the 6 nearest I atoms to this Pb atom
                    nearest_I = find_nearest_I_atoms(pb_coords, atoms, coordinates)
                    
                    print("\n6 Nearest I Atoms to Central Pb:")
                    for i, iodine in enumerate(nearest_I, 1):
                        i_coords = iodine['coordinates']
                        print(f"{i}. I atom at X: {i_coords[0]:.6f}, Y: {i_coords[1]:.6f}, Z: {i_coords[2]:.6f}")
                        print(f"   Distance from Pb: {iodine['distance']:.6f}")
                    
                else:
                    print("No Pb atoms found in the structure")
            else:
                # The center atom is already Pb, find the 6 nearest I atoms to it
                print("\nFinding 6 nearest I atoms to the central Pb atom...")
                nearest_I = find_nearest_I_atoms(result['center_coordinates'], atoms, coordinates)
                
                print("\n6 Nearest I Atoms to Central Pb:")
                for i, iodine in enumerate(nearest_I, 1):
                    i_coords = iodine['coordinates']
                    print(f"{i}. I atom at X: {i_coords[0]:.6f}, Y: {i_coords[1]:.6f}, Z: {i_coords[2]:.6f}")
                    print(f"   Distance from Pb: {iodine['distance']:.6f}")
        else:
            print("No results found from coordinate calculation")
    except Exception as e:
        print(f"An error occurred: {e}")

# Run the main function when the script is executed
if __name__ == "__main__":
    main()
