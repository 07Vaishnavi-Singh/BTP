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
    Calculate the volume of an octahedron by decomposing it into 8 tetrahedra
    Each tetrahedron is formed by a triangular face and the central Pb atom
    
    Parameters:
    -----------
    central_atom_coords : tuple
        (x, y, z) coordinates of the central atom (Pb)
    coordinating_atoms : list
        List of dictionaries containing information about the 6 I atoms
        Each dictionary should have a 'coordinates' key with (x, y, z) values
        
    Returns:
    --------
    dict
        Dictionary containing volume calculations and related information
    """
    from itertools import combinations
    
    # Verify we have exactly 6 coordinating atoms
    if len(coordinating_atoms) != 6:
        raise ValueError(f"Expected 6 coordinating atoms, got {len(coordinating_atoms)}")
    
    # Helper functions for vector operations
    def vector_subtract(v1, v2):
        return (v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2])
    
    def det_3x3(matrix):
        # Calculate determinant of 3x3 matrix
        a, b, c = matrix
        return (
            a[0] * (b[1] * c[2] - b[2] * c[1]) -
            a[1] * (b[0] * c[2] - b[2] * c[0]) +
            a[2] * (b[0] * c[1] - b[1] * c[0])
        )
    
    def calculate_tetrahedron_volume(p1, p2, p3, p4):
        """
        Calculate volume of a tetrahedron given 4 points
        Volume = |(a-d)·((b-d)×(c-d))| / 6
        where a,b,c are three points forming a face and d is the fourth point
        """
        # Create vectors from point 4 to other points
        v1 = vector_subtract(p1, p4)
        v2 = vector_subtract(p2, p4)
        v3 = vector_subtract(p3, p4)
        
        # Calculate determinant (which gives 6 times the volume)
        matrix = [v1, v2, v3]
        volume = abs(det_3x3(matrix)) / 6
        return volume
    
    # Get coordinates of all I atoms
    i_coords = [atom['coordinates'] for atom in coordinating_atoms]
    
    # Find all triangular faces of the octahedron
    # We'll use combinations of 3 vertices and validate each as a face
    tetrahedra_volumes = []
    face_indices = []
    total_volume = 0
    
    # Helper function to check if three points form a valid face
    def is_valid_face(p1_idx, p2_idx, p3_idx, all_points):
        """Check if three points form a valid face of the octahedron"""
        # In a valid face, all other points should be on the same side of the plane
        p1, p2, p3 = [all_points[i] for i in (p1_idx, p2_idx, p3_idx)]
        
        # Calculate normal vector of the face
        v1 = vector_subtract(p2, p1)
        v2 = vector_subtract(p3, p1)
        
        # Cross product gives normal vector
        nx = v1[1]*v2[2] - v1[2]*v2[1]
        ny = v1[2]*v2[0] - v1[0]*v2[2]
        nz = v1[0]*v2[1] - v1[1]*v2[0]
        
        # Check if all other points are on the same side of the plane
        sign = None
        for i, p in enumerate(all_points):
            if i not in (p1_idx, p2_idx, p3_idx):
                # Vector from p1 to current point
                v = vector_subtract(p, p1)
                # Dot product with normal vector
                dot = nx*v[0] + ny*v[1] + nz*v[2]
                
                if sign is None:
                    sign = dot > 0
                elif (dot > 0) != sign:
                    return False
        return True
    
    # Find all valid faces and calculate tetrahedra volumes
    for face in combinations(range(6), 3):
        if is_valid_face(face[0], face[1], face[2], i_coords):
            # Calculate volume of tetrahedron formed by this face and center point
            volume = calculate_tetrahedron_volume(
                i_coords[face[0]],
                i_coords[face[1]],
                i_coords[face[2]],
                central_atom_coords
            )
            
            tetrahedra_volumes.append(volume)
            face_indices.append(face)
            total_volume += volume
    
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
    Find the 6 nearest I atoms that form an octahedral arrangement around the central Pb atom
    
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
        List of dictionaries containing information about the nearest I atoms in octahedral arrangement
    """
    import math  # Add this import at the start of the function
    
    distances = []
    central_x, central_y, central_z = central_pb_coords
    MAX_PB_I_DISTANCE = 4.0  # Maximum reasonable Pb-I bond distance
    
    # First, collect all nearby I atoms
    for i in range(len(atoms)):
        if atoms[i] == 'I':
            x, y, z = coordinates[i]
            distance = ((x - central_x)**2 + 
                       (y - central_y)**2 + 
                       (z - central_z)**2)**0.5
            
            if distance <= MAX_PB_I_DISTANCE:
                # Calculate directional vector from Pb to I
                direction = (
                    (x - central_x) / distance,
                    (y - central_y) / distance,
                    (z - central_z) / distance
                )
                
                distances.append({
                    'index': i,
                    'atom': atoms[i],
                    'coordinates': coordinates[i],
                    'distance': distance,
                    'direction': direction,
                    'atom_index': i + 1  # Adding 1 because file starts from line 3
                })
    
    if len(distances) < 6:
        raise ValueError(f"Only found {len(distances)} I atoms within {MAX_PB_I_DISTANCE} Å. Need 6 for octahedral arrangement.")
    
    # Sort by distance initially
    distances.sort(key=lambda x: x['distance'])
    
    # Function to calculate angle between two vectors
    def angle_between(v1, v2):
        dot_product = sum(a * b for a, b in zip(v1, v2))
        # Clamp dot product to [-1, 1] to avoid floating point errors
        dot_product = max(-1.0, min(1.0, dot_product))
        return abs(180 * math.acos(dot_product) / math.pi)
    
    # Function to check if two atoms are approximately opposite (180 degrees)
    def are_opposite(atom1, atom2, tolerance=20):  # tolerance in degrees
        angle = angle_between(atom1['direction'], atom2['direction'])
        return abs(angle - 180) < tolerance
    
    # Function to check if two atoms are approximately perpendicular (90 degrees)
    def are_perpendicular(atom1, atom2, tolerance=20):  # tolerance in degrees
        angle = angle_between(atom1['direction'], atom2['direction'])
        return abs(angle - 90) < tolerance
    
    # Find best octahedral arrangement
    best_octahedral = []
    min_deviation = float('inf')
    
    # Try each of the closest atoms as the first axis
    for i in range(min(10, len(distances))):  # Check first 10 closest atoms
        for j in range(i + 1, min(11, len(distances))):
            if are_opposite(distances[i], distances[j]):
                # Found potential axis, now look for 4 atoms perpendicular to this axis
                axis_atoms = [distances[i], distances[j]]
                perpendicular_candidates = []
                
                for k in range(len(distances)):
                    if k != i and k != j:
                        atom = distances[k]
                        if all(are_perpendicular(atom, axis_atom) for axis_atom in axis_atoms):
                            perpendicular_candidates.append(atom)
                
                # If we found at least 4 perpendicular candidates
                if len(perpendicular_candidates) >= 4:
                    # Sort perpendicular candidates by distance
                    perpendicular_candidates.sort(key=lambda x: x['distance'])
                    current_arrangement = axis_atoms + perpendicular_candidates[:4]
                    
                    # Calculate deviation from ideal octahedral angles
                    total_deviation = 0
                    for a1 in range(6):
                        for a2 in range(a1 + 1, 6):
                            angle = angle_between(
                                current_arrangement[a1]['direction'],
                                current_arrangement[a2]['direction']
                            )
                            # In perfect octahedron, angles should be either 90 or 180 degrees
                            total_deviation += min(abs(angle - 90), abs(angle - 180))
                    
                    if total_deviation < min_deviation:
                        min_deviation = total_deviation
                        best_octahedral = current_arrangement
    
    if not best_octahedral:
        raise ValueError("Could not find suitable octahedral arrangement of I atoms")
    
    return best_octahedral


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
            'num_atoms': num_atoms,
            'atom_index': center_atom_index + 1  # Adding 1 because file starts from line 3
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
            'distance': min_distance,
            'atom_index': nearest_pb_index + 1  # Adding 1 because file starts from line 3
        }
    return None

def calculate_octahedral_edges(central_atom_coords, coordinating_atoms):
    """
    Calculate the average edge length of the octahedron formed by Pb and 6 I atoms
    
    Parameters:
    -----------
    central_atom_coords : tuple
        (x, y, z) coordinates of the central atom (Pb)
    coordinating_atoms : list
        List of dictionaries containing information about the 6 I atoms
        Each dictionary should have a 'coordinates' key with (x, y, z) values
        
    Returns:
    --------
    dict
        Dictionary containing edge length calculations and statistics
    """
    from itertools import combinations
    
    # Verify we have exactly 6 coordinating atoms
    if len(coordinating_atoms) != 6:
        raise ValueError(f"Expected 6 coordinating atoms, got {len(coordinating_atoms)}")
    
    def calculate_distance(p1, p2):
        """Calculate Euclidean distance between two points"""
        return ((p2[0] - p1[0])**2 + 
                (p2[1] - p1[1])**2 + 
                (p2[2] - p1[2])**2)**0.5
    
    # Get coordinates of all I atoms
    i_coords = [atom['coordinates'] for atom in coordinating_atoms]
    
    # Calculate all edge lengths
    edge_lengths = []
    edge_details = []  # Store which atoms form each edge
    
    # Calculate edges between I atoms (8 edges)
    for i, j in combinations(range(6), 2):
        length = calculate_distance(i_coords[i], i_coords[j])
        # Only consider edges that could reasonably form the octahedron
        # (exclude edges that would cross through the center)
        edge_details.append({
            'atom1_index': coordinating_atoms[i]['atom_index'],
            'atom2_index': coordinating_atoms[j]['atom_index'],
            'length': length,
            'type': 'I-I'
        })
        edge_lengths.append(length)
    
    # Calculate edges from central Pb to each I (6 edges)
    for i, coord in enumerate(i_coords):
        length = calculate_distance(central_atom_coords, coord)
        edge_details.append({
            'atom1_index': 'Pb',
            'atom2_index': coordinating_atoms[i]['atom_index'],
            'length': length,
            'type': 'Pb-I'
        })
        edge_lengths.append(length)
    
    # Calculate statistics
    avg_edge_length = sum(edge_lengths) / len(edge_lengths)
    min_edge_length = min(edge_lengths)
    max_edge_length = max(edge_lengths)
    
    return {
        'average_edge_length': avg_edge_length,
        'min_edge_length': min_edge_length,
        'max_edge_length': max_edge_length,
        'edge_details': edge_details,
        'number_of_edges': len(edge_lengths)
    }

def calculate_effective_coordination_number(central_atom_coords, coordinating_atoms):
    """
    Calculate the Effective Coordination Number (ECoN) as defined in VESTA
    
    Parameters:
    -----------
    central_atom_coords : tuple
        (x, y, z) coordinates of the central atom
    coordinating_atoms : list
        List of dictionaries containing information about the coordinating atoms
        Each dictionary should have a 'coordinates' key with (x, y, z) values
        
    Returns:
    --------
    dict
        Dictionary containing ECoN and related calculations
    """
    import math
    
    # Helper function to calculate distance between two points
    def calculate_distance(p1, p2):
        return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))
    
    # Calculate bond lengths to all coordinating atoms
    bond_lengths = []
    for atom in coordinating_atoms:
        distance = calculate_distance(central_atom_coords, atom['coordinates'])
        bond_lengths.append(distance)
    
    # Find the minimum bond length
    l_min = min(bond_lengths)
    
    # Calculate weighted average bond length (lav) using equation 11.6
    numerator = 0
    denominator = 0
    for li in bond_lengths:
        exp_factor = math.exp(1 - (li/l_min)**6)
        numerator += li * exp_factor
        denominator += exp_factor
    
    l_av = numerator / denominator if denominator != 0 else 0
    
    # Calculate bond weights and ECoN using equations 11.4 and 11.5
    bond_weights = []
    econ = 0
    
    for li in bond_lengths:
        wi = math.exp(1 - (li/l_av)**6)
        bond_weights.append(wi)
        econ += wi
    
    return {
        'bond_lengths': bond_lengths,
        'minimum_bond_length': l_min,
        'weighted_average_bond_length': l_av,
        'bond_weights': bond_weights,
        'effective_coordination_number': econ
    }

def calculate_quadratic_elongation(central_atom_coords, coordinating_atoms, polyhedron_type='octahedron'):
    """
    Calculate the quadratic elongation <λ> as defined in VESTA
    
    Parameters:
    -----------
    central_atom_coords : tuple
        (x, y, z) coordinates of the central atom
    coordinating_atoms : list
        List of dictionaries containing information about the coordinating atoms
        Each dictionary should have a 'coordinates' key with (x, y, z) values
    polyhedron_type : str
        Type of polyhedron ('tetrahedron', 'octahedron', 'cube', 'dodecahedron', 'icosahedron')
        
    Returns:
    --------
    dict
        Dictionary containing quadratic elongation and related calculations
    """
    import math
    from itertools import combinations
    
    # Helper function to calculate distance between two points
    def calculate_distance(p1, p2):
        return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))
    
    # Calculate center-to-vertex distances for actual polyhedron
    center_to_vertex_distances = []
    for atom in coordinating_atoms:
        distance = calculate_distance(central_atom_coords, atom['coordinates'])
        center_to_vertex_distances.append(distance)
    
    # Calculate mean bond length
    mean_bond_length = sum(center_to_vertex_distances) / len(center_to_vertex_distances)
    
    # For octahedron
    if polyhedron_type == 'octahedron':
        # Use direct formula for octahedron
        # In a regular octahedron, all bond lengths are equal to l0
        # So we calculate quadratic elongation using average bond length
        quad_elongation = sum((li / mean_bond_length)**2 for li in center_to_vertex_distances) / len(center_to_vertex_distances)
        
        return {
            'center_to_vertex_distances': center_to_vertex_distances,
            'ideal_center_to_vertex_distance': mean_bond_length,
            'mean_bond_length': mean_bond_length,
            'quadratic_elongation': quad_elongation
        }
    else:
        # For other polyhedra, a more complex calculation would be needed
        # This is a simplified implementation
        raise NotImplementedError(f"Quadratic elongation calculation for {polyhedron_type} not implemented")

def calculate_bond_angle_variance(central_atom_coords, coordinating_atoms, polyhedron_type='octahedron'):
    """
    Calculate the bond angle variance (σ²) as defined in VESTA
    
    Parameters:
    -----------
    central_atom_coords : tuple
        (x, y, z) coordinates of the central atom
    coordinating_atoms : list
        List of dictionaries containing information about the coordinating atoms
        Each dictionary should have a 'coordinates' key with (x, y, z) values
    polyhedron_type : str
        Type of polyhedron ('tetrahedron', 'octahedron', 'cube', 'dodecahedron', 'icosahedron')
        
    Returns:
    --------
    dict
        Dictionary containing bond angle variance and related calculations
    """
    import math
    from itertools import combinations
    
    # Helper function to calculate the angle between three points (p1-p2-p3)
    def calculate_angle(p1, p2, p3):
        # Create vectors
        v1 = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]]
        v2 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]]
        
        # Calculate dot product
        dot_product = sum(a * b for a, b in zip(v1, v2))
        
        # Calculate magnitudes
        mag1 = math.sqrt(sum(a * a for a in v1))
        mag2 = math.sqrt(sum(a * a for a in v2))
        
        # Calculate cosine of angle
        cos_angle = dot_product / (mag1 * mag2)
        cos_angle = max(-1.0, min(1.0, cos_angle))  # Clamp to avoid floating point errors
        
        # Convert to degrees
        angle_deg = math.degrees(math.acos(cos_angle))
        return angle_deg
    
    # Ideal bond angles for different polyhedra
    ideal_bond_angles = {
        'tetrahedron': 109.47,  # 109°28'
        'octahedron': 90.0,
        'cube': 90.0,
        'dodecahedron': [108.0, 116.57, 120.0],  # Depending on which face
        'icosahedron': 108.0
    }
    
    if polyhedron_type not in ideal_bond_angles:
        raise ValueError(f"Unsupported polyhedron type: {polyhedron_type}")
    
    # For simplicity, use the first ideal angle for polyhedra with multiple angles
    phi0 = ideal_bond_angles[polyhedron_type]
    if isinstance(phi0, list):
        phi0 = phi0[0]  # Just use the first angle for simplicity
    
    # Get coordinates of central atom and all vertices
    central = central_atom_coords
    vertices = [atom['coordinates'] for atom in coordinating_atoms]
    
    # Calculate all bond angles (vertex-central-vertex)
    bond_angles = []
    for i, j in combinations(range(len(vertices)), 2):
        angle = calculate_angle(vertices[i], central, vertices[j])
        bond_angles.append(angle)
    
    # For octahedron, we expect 12 bond angles (each of the 6 vertices connects to the central atom and 
    # forms angles with the other 5 vertices, giving 6*5/2 = 15 angles)
    # But in a perfect octahedron, some angles should be exactly 90° and others exactly 180°
    
    # Group angles into those that should be 90° and those that should be 180°
    angles_near_90 = []
    angles_near_180 = []
    
    for angle in bond_angles:
        if abs(angle - 90) < abs(angle - 180):
            angles_near_90.append(angle)
        else:
            angles_near_180.append(angle)
    
    # Calculate variance for the angles that should be 90°
    variance_90 = sum((angle - 90.0)**2 for angle in angles_near_90) / len(angles_near_90) if angles_near_90 else 0
    
    # For octahedron, we only consider the 90-degree angles for the variance calculation
    if polyhedron_type == 'octahedron':
        variance = variance_90
    else:
        # For other polyhedra, use all angles
        variance = sum((phi - phi0)**2 for phi in bond_angles) / len(bond_angles) if bond_angles else 0
    
    return {
        'bond_angles': bond_angles,
        'angles_near_90': angles_near_90,
        'angles_near_180': angles_near_180,
        'ideal_bond_angle': phi0,
        'number_of_angles': len(bond_angles),
        'bond_angle_variance': variance
    }

def main(xyz_file_path='vesta_file.xyz', specific_atom_index=None):
    try:
        # Parse the file
        print(f"Parsing XYZ file: {xyz_file_path}")
        atoms, coordinates = parse_xyz_file(xyz_file_path)
        print(f"Found {len(atoms)} atoms in the file")
        
        # Check if we're analyzing a specific atom
        if specific_atom_index is not None:
            # Adjust from line number to array index if needed
            # Convert from 1-based index (atom number) to 0-based array index
            atom_index = specific_atom_index - 1
            
            if atom_index < 0 or atom_index >= len(atoms):
                print(f"Error: Atom index {specific_atom_index} is out of range. The file has {len(atoms)} atoms.")
                return
            
            if atoms[atom_index] != 'Pb':
                print(f"Warning: Atom at index {specific_atom_index} is {atoms[atom_index]}, not Pb.")
                return
            
            # Use this specific Pb atom
            central_pb_coords = coordinates[atom_index]
            print(f"\nAnalyzing specific Pb atom at line number {atom_index + 3} (atom #{atom_index + 1}):")
            print(f"Coordinates: X: {central_pb_coords[0]:.6f}, Y: {central_pb_coords[1]:.6f}, Z: {central_pb_coords[2]:.6f}")
            
            # Find the 6 nearest I atoms to this specific Pb atom
            print("\nFinding 6 nearest I atoms to this Pb atom...")
            nearest_I = find_nearest_I_atoms(central_pb_coords, atoms, coordinates)
            
            # Proceed with analysis as normal
            if nearest_I:
                analyze_and_display_results(central_pb_coords, nearest_I)
            else:
                print("Could not find enough I atoms to analyze this Pb atom.")
            
            return
            
        # If no specific atom was requested, proceed with the center atom analysis as before
        result = find_center_atom(atoms, coordinates)
        
        if result:
            print("\nCenter Atom Analysis:")
            print(f"Center atom type: {result['center_atom']}")
            print(f"Center atom line number in file: {result['atom_index'] + 2}")  # +2 because atoms start from line 3
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
                    print(f"Pb atom line number in file: {nearest_pb['atom_index'] + 2}")  # +2 because atoms start from line 3
                    print(f"Coordinates: X: {central_pb_coords[0]:.6f}, Y: {central_pb_coords[1]:.6f}, Z: {central_pb_coords[2]:.6f}")
                    print(f"Distance from center atom: {nearest_pb['distance']:.6f}")
                    
                    # Find the 6 nearest I atoms to this Pb atom
                    print("\nFinding 6 nearest I atoms to this Pb atom...")
                    nearest_I = find_nearest_I_atoms(central_pb_coords, atoms, coordinates)
                else:
                    print("No Pb atoms found in the structure")
                    return
            
            # Display the 6 nearest I atoms and analyze results
            if nearest_I:
                analyze_and_display_results(central_pb_coords, nearest_I)
            else:
                print("Could not find enough I atoms to calculate distortion index and polyhedral volume")
        else:
            print("No results found from coordinate calculation")
    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()

def analyze_and_display_results(central_pb_coords, nearest_I):
    """
    Perform analysis on a Pb atom and its nearest I atoms and display the results
    
    Parameters:
    -----------
    central_pb_coords : tuple
        (x, y, z) coordinates of the central Pb atom
    nearest_I : list
        List of dictionaries containing information about the nearest I atoms
    """
    # Display the 6 nearest I atoms
    print("\n6 Nearest I Atoms to Central Pb:")
    for i, iodine in enumerate(nearest_I, 1):
        i_coords = iodine['coordinates']
        print(f"{i}. I atom at line number {iodine['atom_index'] + 2} in file")  # +2 because atoms start from line 3
        print(f"   Coordinates: X: {i_coords[0]:.6f}, Y: {i_coords[1]:.6f}, Z: {i_coords[2]:.6f}")
        print(f"   Distance from Pb: {iodine['distance']:.6f}")
    
    # Calculate distortion index
    print("\nCalculating Distortion Index...")
    distortion_results = calculate_distortion_index(central_pb_coords, nearest_I)
    
    # Calculate polyhedral (octahedral) volume
    print("\nCalculating Polyhedral Volume...")
    volume_results = calculate_octahedral_volume(central_pb_coords, nearest_I)
    
    # Calculate edge lengths
    print("\nCalculating Octahedral Edge Lengths...")
    edge_results = calculate_octahedral_edges(central_pb_coords, nearest_I)
    
    # Calculate ECoN
    print("\nCalculating Effective Coordination Number (ECoN)...")
    econ_results = calculate_effective_coordination_number(central_pb_coords, nearest_I)
    
    # Calculate quadratic elongation
    print("\nCalculating Quadratic Elongation...")
    quad_elongation_results = calculate_quadratic_elongation(central_pb_coords, nearest_I)
    
    # Calculate bond angle variance
    print("\nCalculating Bond Angle Variance...")
    bond_angle_variance_results = calculate_bond_angle_variance(central_pb_coords, nearest_I)
    
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
    
    print("\nEdge Length Results:")
    print(f"Average edge length: {edge_results['average_edge_length']:.6f} Å")
    print(f"Minimum edge length: {edge_results['min_edge_length']:.6f} Å")
    print(f"Maximum edge length: {edge_results['max_edge_length']:.6f} Å")
    print(f"Total number of edges: {edge_results['number_of_edges']}")
    
    # Display edge length results in a more readable format
    print("\nDetailed Edge Information:")
    print("\nPb-I bonds:")
    for edge in edge_results['edge_details']:
        if edge['type'] == 'Pb-I':
            print(f"Pb to I(line {edge['atom2_index'] + 2}): {edge['length']:.6f} Å")
    
    print("\nI-I edges:")
    for edge in edge_results['edge_details']:
        if edge['type'] == 'I-I':
            print(f"I(line {edge['atom1_index'] + 2}) to I(line {edge['atom2_index'] + 2}): {edge['length']:.6f} Å")
    
    print("\nECoN Results:")
    print(f"Bond lengths: {[f'{bl:.6f}' for bl in econ_results['bond_lengths']]}")
    print(f"Minimum bond length: {econ_results['minimum_bond_length']:.6f} Å")
    print(f"Weighted average bond length: {econ_results['weighted_average_bond_length']:.6f} Å")
    print(f"Bond weights: {[f'{bw:.6f}' for bw in econ_results['bond_weights']]}")
    print(f"Effective Coordination Number (ECoN): {econ_results['effective_coordination_number']:.6f}")
    
    print("\nQuadratic Elongation Results:")
    print(f"Center-to-vertex distances: {[f'{d:.6f}' for d in quad_elongation_results['center_to_vertex_distances']]}")
    print(f"Mean bond length: {quad_elongation_results['mean_bond_length']:.6f} Å")
    print(f"Quadratic elongation: {quad_elongation_results['quadratic_elongation']:.4f}")
    
    print("\nBond Angle Variance Results:")
    print(f"Number of angles: {bond_angle_variance_results['number_of_angles']}")
    print(f"Number of ~90° angles: {len(bond_angle_variance_results['angles_near_90'])}")
    print(f"Number of ~180° angles: {len(bond_angle_variance_results['angles_near_180'])}")
    print(f"Bond angle variance: {bond_angle_variance_results['bond_angle_variance']:.4f} deg.²")
    
    # Print a summary of key metrics
    print("\nSummary of Key Metrics:")
    print(f"Average bond length = {distortion_results['average_bond_length']:.4f} Å")
    print(f"Polyhedral volume = {volume_results['total_volume']:.4f} Å^3")
    print(f"Distortion index (bond length) = {distortion_results['distortion_index']:.5f}")
    print(f"Quadratic elongation = {quad_elongation_results['quadratic_elongation']:.4f}")
    print(f"Bond angle variance = {bond_angle_variance_results['bond_angle_variance']:7.4f} deg.^2")
    print(f"Effective coordination number = {econ_results['effective_coordination_number']:.4f}")
        
# Run the main function when the script is executed
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        # If atom index is provided
        try:
            atom_index = int(sys.argv[2])
            main(sys.argv[1], atom_index)
        except ValueError:
            print(f"Error: Invalid atom index: {sys.argv[2]}. Please provide a valid integer.")
    elif len(sys.argv) > 1:
        # Only file path is provided
        main(sys.argv[1])
    else:
        # No arguments, use default
        main()

