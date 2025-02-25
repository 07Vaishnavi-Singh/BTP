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

# @dev fucntion to calcualte the volume of octahedral 
def calculate_octahedral_volume(central_atom_coords, coordinating_atoms):
    """
    Calculate the volume of an octahedron by decomposing it into tetrahedra
    
    Parameters:
    -----------
    central_atom_coords : tuple
        (x, y, z) coordinates of the central atom (Pb)
    coordinating_atoms : list
        List of dictionaries containing information about the coordinating atoms (I)
        Each dictionary should have a 'coordinates' key with (x, y, z) values
        
    Returns:
    --------
    dict
        Dictionary containing volume calculations and related information
    """
    from itertools import combinations
    
    # Extract coordinates
    center = central_atom_coords
    coords = [atom['coordinates'] for atom in coordinating_atoms]
    
    # Helper functions to replace numpy operations
    def vector_subtract(v1, v2):
        return (v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2])
    
    def vector_add(v1, v2):
        return (v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2])
    
    def vector_scale(v, scalar):
        return (v[0] * scalar, v[1] * scalar, v[2] * scalar)
    
    def vector_cross(v1, v2):
        return (
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]
        )
    
    def vector_dot(v1, v2):
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
    
    def det_3x3(matrix):
        # Calculate determinant of 3x3 matrix
        a, b, c = matrix
        return (
            a[0] * (b[1] * c[2] - b[2] * c[1]) -
            a[1] * (b[0] * c[2] - b[2] * c[0]) +
            a[2] * (b[0] * c[1] - b[1] * c[0])
        )
    
    # Calculate volume of each tetrahedron and sum
    total_volume = 0
    tetrahedra_volumes = []
    face_indices = []
    
    # Find all possible faces (combinations of 3 vertices)
    for face_combo in combinations(range(6), 3):
        i, j, k = face_combo
        
        # Get the three points that might form a face
        a, b, c = coords[i], coords[j], coords[k]
        
        # Calculate normal vector of the potential face
        ab = vector_subtract(b, a)
        ac = vector_subtract(c, a)
        normal = vector_cross(ab, ac)
        
        # Calculate centroid of the face
        sum_points = vector_add(vector_add(a, b), c)
        face_centroid = vector_scale(sum_points, 1/3)
        
        # Vector from central atom to face centroid
        center_to_face = vector_subtract(face_centroid, center)
        
        # Check if this is an outward-facing face
        if vector_dot(normal, center_to_face) > 0:
            face_indices.append(face_combo)
            
            # Calculate tetrahedron volume
            # V = (1/6) * |det(a-center, b-center, c-center)|
            matrix = [
                vector_subtract(a, center),
                vector_subtract(b, center),
                vector_subtract(c, center)
            ]
            volume = abs(det_3x3(matrix)) / 6
            
            tetrahedra_volumes.append(volume)
            total_volume += volume
    
    # For data collection and validation
    face_coordinates = []
    for face in face_indices:
        face_coordinates.append([coords[i] for i in face])
    
    return {
        'total_volume': total_volume,
        'tetrahedra_volumes': tetrahedra_volumes,
        'number_of_tetrahedra': len(tetrahedra_volumes),
        'face_indices': face_indices
    }



# @dev function to find out DI from neighbouring I atom 
def calculate_distortion_index(central_pb_coords, nearest_i_atoms):
    """
    Calculate the distortion index (D) based on Baur's formula
    
    Parameters:
    -----------
    central_pb_coords : tuple
        (x, y, z) coordinates of the central Pb atom
    nearest_i_atoms : list
        List of dictionaries containing information about the nearest I atoms
        Each dictionary should have a 'coordinates' key with (x, y, z) values
        
    Returns:
    --------
    dict
        Dictionary containing distortion index and related calculations
    """
    n = len(nearest_i_atoms)  # Number of coordinating atoms (I)
    
    # Calculate individual bond lengths (li)
    bond_lengths = []
    pb_x, pb_y, pb_z = central_pb_coords
    
    for iodine in nearest_i_atoms:
        i_x, i_y, i_z = iodine['coordinates']
        
        # Distance from Pb to I (bond length)
        li = ((i_x - pb_x)**2 + (i_y - pb_y)**2 + (i_z - pb_z)**2)**0.5
        bond_lengths.append(li)
    
    # Calculate average bond length (lav)
    lav = sum(bond_lengths) / n
    
    # Calculate individual deviations and their absolute values
    deviations = []
    abs_deviations = []
    
    for li in bond_lengths:
        dev = (li - lav) / lav  # Normalized deviation
        deviations.append(dev)
        abs_deviations.append(abs(dev))  # Take absolute value
    
    # Calculate distortion index D
    distortion_index = sum(abs_deviations) / n
    
    return {
        'bond_lengths': bond_lengths,
        'average_bond_length': lav,
        'deviations': deviations,
        'absolute_deviations': abs_deviations,
        'distortion_index': distortion_index
    }


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


def main(xyz_file_path='vesta_file.xyz'):
    try:
        # Parse the file
        print(f"Parsing XYZ file: {xyz_file_path}")
        atoms, coordinates = parse_xyz_file(xyz_file_path)
        print(f"Found {len(atoms)} atoms in the file")
        
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
            
            # Determine which Pb atom to use and find the 6 nearest I atoms
            nearest_I = None
            central_pb_coords = None
            
            if result['center_atom'] == "Pb":
                # The center atom is already Pb
                central_pb_coords = result['center_coordinates']
                print("\nCenter atom is already Pb. Finding 6 nearest I atoms...")
                nearest_I = find_nearest_I_atoms(central_pb_coords, atoms, coordinates)
            else:
                # Need to find the nearest Pb atom
                print(f"\nFinding nearest Pb atom to the center atom...")
                nearest_pb = find_nearest_pb_atom(
                    result['center_coordinates'],
                    atoms,
                    coordinates
                )
                
                if nearest_pb:
                    central_pb_coords = nearest_pb['coordinates']
                    print("\nNearest Pb Atom Found:")
                    print(f"Coordinates: X: {central_pb_coords[0]:.6f}, Y: {central_pb_coords[1]:.6f}, Z: {central_pb_coords[2]:.6f}")
                    print(f"Distance from center atom: {nearest_pb['distance']:.6f}")
                    
                    # Find the 6 nearest I atoms to this Pb atom
                    print("\nFinding 6 nearest I atoms to this Pb atom...")
                    nearest_I = find_nearest_I_atoms(central_pb_coords, atoms, coordinates)
                else:
                    print("No Pb atoms found in the structure")
                    return
            
            # Display the 6 nearest I atoms
            if nearest_I:
                print("\n6 Nearest I Atoms to Central Pb:")
                for i, iodine in enumerate(nearest_I, 1):
                    i_coords = iodine['coordinates']
                    print(f"{i}. I atom at X: {i_coords[0]:.6f}, Y: {i_coords[1]:.6f}, Z: {i_coords[2]:.6f}")
                    print(f"   Distance from Pb: {iodine['distance']:.6f}")
                
                # Calculate distortion index
                print("\nCalculating Distortion Index...")
                distortion_results = calculate_distortion_index(central_pb_coords, nearest_I)
                
                # Calculate polyhedral (octahedral) volume
                print("\nCalculating Polyhedral Volume...")
                volume_results = calculate_octahedral_volume(central_pb_coords, nearest_I)
                
                # Display results
                print("\nDistortion Index Results:")
                print(f"Bond lengths (Pb-I): {[f'{bl:.6f}' for bl in distortion_results['bond_lengths']]}")
                print(f"Average bond length (lav): {distortion_results['average_bond_length']:.6f}")
                print(f"Normalized deviations (|li-lav|/lav): {[f'{ad:.6f}' for ad in distortion_results['absolute_deviations']]}")
                print(f"\nDistortion Index (D): {distortion_results['distortion_index']:.6f}")
                
                print("\nPolyhedral Volume Results:")
                print(f"Total Octahedral Volume: {volume_results['total_volume']:.6f} Å³")
                print(f"Number of tetrahedra used: {volume_results['number_of_tetrahedra']}")
                print(f"Individual tetrahedra volumes: {[f'{v:.6f}' for v in volume_results['tetrahedra_volumes']]}")
            else:
                print("Could not find enough I atoms to calculate distortion index and polyhedral volume")
        else:
            print("No results found from coordinate calculation")
    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()

# Run the main function when the script is executed
if __name__ == "__main__":
    main()

