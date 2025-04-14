def parse_xyz_file(filepath):
    """
    Parse XYZ file and return atom names and coordinates
    
    Parameters:
    -----------
    filepath : str
        Path to the XYZ file to be parsed
        
    Returns:
    --------
    tuple
        A tuple containing two lists:
        1. atoms: List of atom symbols (e.g., 'C', 'N', 'Pb', 'I')
        2. coordinates: List of (x, y, z) coordinate tuples for each atom
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
        """Calculate vector difference v1 - v2"""
        return (v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2])
    
    def det_3x3(matrix):
        """Calculate determinant of 3x3 matrix"""
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
        """
        Check if three points form a valid face of the octahedron
        A valid face has all other points on the same side of the plane
        """
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
    
    # Helper function to calculate volume of the polyhedron
    def calculate_polyhedron_volume(central_coords, vertex_coords, poly_type):
        if poly_type == 'octahedron':
            # For octahedron, use the existing function from calculate_octahedral_volume
            from itertools import combinations
            
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
            tetrahedra_volumes = []
            total_volume = 0
            
            for face in combinations(range(len(vertex_coords)), 3):
                if is_valid_face(face[0], face[1], face[2], vertex_coords):
                    # Calculate volume of tetrahedron formed by this face and center point
                    volume = calculate_tetrahedron_volume(
                        vertex_coords[face[0]],
                        vertex_coords[face[1]],
                        vertex_coords[face[2]],
                        central_coords
                    )
                    
                    tetrahedra_volumes.append(volume)
                    total_volume += volume
            
            return total_volume
        
        else:
            raise NotImplementedError(f"Volume calculation for {poly_type} not implemented")
    
    # Get coordinates of all vertices
    vertices = [atom['coordinates'] for atom in coordinating_atoms]
    
    # Calculate actual distances from center to each vertex
    distances = []
    for vertex in vertices:
        dist = calculate_distance(central_atom_coords, vertex)
        distances.append(dist)
    
    # Get the number of vertices
    n = len(distances)
    
    # Check if number of vertices matches the polyhedron type
    expected_vertices = {
        'tetrahedron': 4,
        'octahedron': 6,
        'cube': 8,
        'dodecahedron': 12,
        'icosahedron': 20
    }
    
    if n != expected_vertices.get(polyhedron_type, 0):
        raise ValueError(f"Expected {expected_vertices.get(polyhedron_type)} vertices for {polyhedron_type}, got {n}")
    
    # Calculate actual volume of the polyhedron
    actual_volume = calculate_polyhedron_volume(central_atom_coords, vertices, polyhedron_type)
    
    # For octahedron: V = (√2/3) * a³
    # Therefore: a = (3V/√2)^(1/3)
    # And l₀ = a/√2
    if polyhedron_type == 'octahedron':
        # First calculate 'a' from the volume
        a = (3 * actual_volume / math.sqrt(2)) ** (1/3)
        # Then calculate l₀ = a/√2
        l0 = a / math.sqrt(2)
    else:
        raise NotImplementedError(f"l0 calculation for {polyhedron_type} not implemented")
    
    # Calculate quadratic elongation using equation 11.2
    quad_elongation = sum((li / l0)**2 for li in distances) / n
    
    return {
        'center_to_vertex_distances': distances,
        'ideal_center_to_vertex_distance': l0,
        'edge_length_a': a,  # Added to return the calculated edge length
        'actual_volume': actual_volume,
        'quadratic_elongation': quad_elongation
    }

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
        'dodecahedron': 108.0,  # Simplified
        'icosahedron': 108.0    # Simplified
    }
    
    if polyhedron_type not in ideal_bond_angles:
        raise ValueError(f"Unsupported polyhedron type: {polyhedron_type}")
    
    # Number of faces for different polyhedra
    face_counts = {
        'tetrahedron': 4,
        'octahedron': 8,
        'cube': 6,
        'dodecahedron': 12,
        'icosahedron': 20
    }
    
    # Calculate number of bond angles based on the formula: m = (number of faces) * 3 / 2
    # This is the theoretical number of angles in a regular polyhedron
    theoretical_m = int(face_counts.get(polyhedron_type, 0) * 3 / 2)
    
    # Get coordinates of central atom and all vertices
    central = central_atom_coords
    vertices = [atom['coordinates'] for atom in coordinating_atoms]
    
    # Calculate all bond angles (vertex-central-vertex)
    bond_angles = []
    for i, j in combinations(range(len(vertices)), 2):
        angle = calculate_angle(vertices[i], central, vertices[j])
        bond_angles.append(angle)
    
    # The actual number of calculated angles
    m = len(bond_angles)
    
    if polyhedron_type == 'octahedron':
        # For octahedron, we have 6 vertices and expect 12 bond angles
        # (6 choose 2 = 15 potential angles, but 3 are approximately 180°)
        
        # In a perfect octahedron, 12 angles are 90° and 3 angles are 180°
        # VESTA reference says we should only consider the 90° angles
        # Group angles into those that should be 90° and those that should be 180°
        angles_near_90 = []
        angles_near_180 = []
        
        for angle in bond_angles:
            if abs(angle - 90) < abs(angle - 180):
                angles_near_90.append(angle)
            else:
                angles_near_180.append(angle)
        
        # Calculate bond angle variance according to equation 11.3
        # Use m-1 as the divisor where m is the number of 90° angles
        if len(angles_near_90) > 1:
            variance = sum((angle - 90.0)**2 for angle in angles_near_90) / (len(angles_near_90) - 1)
        else:
            variance = 0
    else:
        # For other polyhedra, use all angles and the corresponding ideal angle
        phi0 = ideal_bond_angles[polyhedron_type]
        
        # Use m-1 as the divisor where m is the number of angles
        if m > 1:
            variance = sum((phi - phi0)**2 for phi in bond_angles) / (m - 1)
        else:
            variance = 0
    
    return {
        'bond_angles': bond_angles,
        'angles_near_90': angles_near_90 if polyhedron_type == 'octahedron' else [],
        'angles_near_180': angles_near_180 if polyhedron_type == 'octahedron' else [],
        'ideal_bond_angle': ideal_bond_angles[polyhedron_type],
        'theoretical_bond_angle_count': theoretical_m,
        'actual_bond_angle_count': m,
        'bond_angle_variance': variance
    }

def calculate_axial_and_equatorial_angles(central_pb_coords, coordinating_atoms):
    """
    Identify axial and equatorial I atoms and calculate the relevant I-Pb-I angles
    
    Parameters:
    -----------
    central_pb_coords : tuple
        (x, y, z) coordinates of the central Pb atom
    coordinating_atoms : list
        List of dictionaries containing information about the 6 I atoms
        Each dictionary should have a 'coordinates' key with (x, y, z) values
        
    Returns:
    --------
    dict
        Dictionary containing axial and equatorial angles and related information
    """
    import math
    
    # Verify we have exactly 6 coordinating atoms
    if len(coordinating_atoms) != 6:
        raise ValueError(f"Expected 6 coordinating atoms, got {len(coordinating_atoms)}")
    
    # Helper function to calculate vector between two points
    def calculate_vector(p1, p2):
        return (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    
    # Helper function to normalize a vector
    def normalize_vector(v):
        magnitude = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
        return (v[0]/magnitude, v[1]/magnitude, v[2]/magnitude)
    
    # Helper function to calculate dot product
    def dot_product(v1, v2):
        return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    
    # Helper function to calculate angle between two vectors
    def angle_between_vectors(v1, v2):
        # Normalize vectors
        v1_norm = normalize_vector(v1)
        v2_norm = normalize_vector(v2)
        
        # Calculate dot product
        dot = dot_product(v1_norm, v2_norm)
        
        # Clamp to avoid floating point errors
        dot = max(-1.0, min(1.0, dot))
        
        # Calculate angle in degrees
        angle_rad = math.acos(dot)
        angle_deg = angle_rad * 180.0 / math.pi
        
        return angle_deg
    
    # Calculate vectors from central Pb to each I atom
    vectors = []
    for i, atom in enumerate(coordinating_atoms):
        vector = calculate_vector(central_pb_coords, atom['coordinates'])
        norm_vector = normalize_vector(vector)
        vectors.append({
            'atom_index': atom['atom_index'],
            'vector': vector,
            'normalized': norm_vector,
            'position_index': i
        })
    
    # Find pairs of atoms that are approximately opposite to each other
    opposite_pairs = []
    for i in range(len(vectors)):
        for j in range(i + 1, len(vectors)):
            v1 = vectors[i]['normalized']
            v2 = vectors[j]['normalized']
            
            # Calculate dot product to find opposite vectors
            dot = dot_product(v1, v2)
            
            # If dot product is close to -1, vectors are opposite
            if dot < -0.8:  # Allow some tolerance for non-perfect octahedra
                opposite_pairs.append((i, j))
    
    # There should be 3 pairs of opposite atoms
    if len(opposite_pairs) != 3:
        print(f"Warning: Found {len(opposite_pairs)} pairs of opposite atoms, expected 3.")
    
    # Identify the axial pair: typically the axial pair is the one with the largest separation
    # We'll calculate this as the pair with the angle closest to 180 degrees
    axial_pair_index = 0
    max_angle = 0
    
    for i in range(len(opposite_pairs)):
        v1 = vectors[opposite_pairs[i][0]]['vector']
        v2 = vectors[opposite_pairs[i][1]]['vector']
        angle = angle_between_vectors(v1, v2)
        
        if angle > max_angle:
            max_angle = angle
            axial_pair_index = i
    
    # Get the axial pair
    axial_pair = opposite_pairs[axial_pair_index]
    axial_indices = [vectors[axial_pair[0]]['position_index'], vectors[axial_pair[1]]['position_index']]
    
    # The remaining atoms are equatorial
    all_indices = set(range(6))
    equatorial_indices = list(all_indices - set(axial_indices))
    
    # Calculate axial angle
    axial_angle = angle_between_vectors(
        vectors[axial_pair[0]]['vector'],
        vectors[axial_pair[1]]['vector']
    )
    
    # Calculate all possible angles between equatorial atoms
    equatorial_angles = []
    equatorial_pairs = []
    
    # Calculate angles between each pair of equatorial atoms
    for i in range(len(equatorial_indices)):
        for j in range(i + 1, len(equatorial_indices)):
            eq_idx1 = equatorial_indices[i]
            eq_idx2 = equatorial_indices[j]
            
            angle = angle_between_vectors(
                vectors[eq_idx1]['vector'],
                vectors[eq_idx2]['vector']
            )
            
            equatorial_angles.append({
                'angle': angle,
                'pair': (eq_idx1, eq_idx2),
                'atom_indices': (vectors[eq_idx1]['atom_index'], vectors[eq_idx2]['atom_index'])
            })
    
    # Sort the equatorial angles to better identify adjacent vs opposite atoms
    equatorial_angles.sort(key=lambda x: x['angle'])
    
    # In a perfect octahedron, we expect:
    # - 2 pairs of adjacent equatorial atoms with angles around 90°
    # - 1 pair of opposite equatorial atoms with angle around 180°
    
    # Find the two equatorial angles closest to 90 degrees (adjacent atoms)
    adjacent_equatorial_data = []
    opposite_equatorial_data = []
    
    for angle_data in equatorial_angles:
        if abs(angle_data['angle'] - 90) < abs(angle_data['angle'] - 180):
            # This is closer to 90°, so likely adjacent equatorial atoms
            adjacent_equatorial_data.append(angle_data)
        else:
            # This is closer to 180°, so likely opposite equatorial atoms
            opposite_equatorial_data.append(angle_data)
    
    # Take the two angles closest to 90° (there should be 4 adjacent pairs in a perfect octahedron)
    adjacent_equatorial_data = sorted(adjacent_equatorial_data, key=lambda x: abs(x['angle'] - 90))[:2]
    
    # Also identify the angle between opposite equatorial atoms (should be close to 180°)
    opposite_equatorial_angle = None
    if opposite_equatorial_data:
        opposite_equatorial_angle = opposite_equatorial_data[0]
    
    return {
        'axial_indices': [vectors[idx]['atom_index'] for idx in axial_indices],
        'axial_angle': axial_angle,
        'equatorial_indices': [vectors[idx]['atom_index'] for idx in equatorial_indices],
        'adjacent_equatorial_data': adjacent_equatorial_data,
        'opposite_equatorial_data': opposite_equatorial_angle
    }

def calculate_axial_equatorial_bond_lengths(central_pb_coords, coordinating_atoms, axial_indices, equatorial_indices):
    """
    Calculate the bond lengths for axial and equatorial I-Pb bonds separately
    
    Parameters:
    -----------
    central_pb_coords : tuple
        (x, y, z) coordinates of the central Pb atom
    coordinating_atoms : list
        List of dictionaries containing information about the 6 I atoms
    axial_indices : list
        List of indices identifying the axial I atoms
    equatorial_indices : list
        List of indices identifying the equatorial I atoms
        
    Returns:
    --------
    dict
        Dictionary containing axial and equatorial bond lengths and statistics
    """
    # Helper function to calculate distance between two points
    def calculate_distance(p1, p2):
        return ((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2)**0.5
    
    # Calculate axial bond lengths
    axial_bond_lengths = []
    for idx in axial_indices:
        atom = coordinating_atoms[idx]
        distance = calculate_distance(central_pb_coords, atom['coordinates'])
        axial_bond_lengths.append({
            'atom_index': atom['atom_index'],
            'length': distance
        })
    
    # Calculate equatorial bond lengths
    equatorial_bond_lengths = []
    for idx in equatorial_indices:
        atom = coordinating_atoms[idx]
        distance = calculate_distance(central_pb_coords, atom['coordinates'])
        equatorial_bond_lengths.append({
            'atom_index': atom['atom_index'],
            'length': distance
        })
    
    # Calculate statistics
    avg_axial_length = sum(item['length'] for item in axial_bond_lengths) / len(axial_bond_lengths) if axial_bond_lengths else 0
    avg_equatorial_length = sum(item['length'] for item in equatorial_bond_lengths) / len(equatorial_bond_lengths) if equatorial_bond_lengths else 0
    
    # Calculate the difference between axial and equatorial lengths
    axial_equatorial_difference = avg_axial_length - avg_equatorial_length
    
    return {
        'axial_bond_lengths': axial_bond_lengths,
        'equatorial_bond_lengths': equatorial_bond_lengths,
        'avg_axial_length': avg_axial_length,
        'avg_equatorial_length': avg_equatorial_length,
        'axial_equatorial_difference': axial_equatorial_difference
    }

def calculate_cn_bonds(atoms, coordinates):
    """
    Identify C-N bonds in the structure and calculate their position vectors
    
    Parameters:
    -----------
    atoms : list
        List of atom symbols
    coordinates : list
        List of (x, y, z) coordinates for each atom
        
    Returns:
    --------
    dict
        Dictionary containing C-N bond information
    """
    # Find all C atoms
    carbon_indices = [i for i, atom in enumerate(atoms) if atom == 'C']
    # Find all N atoms
    nitrogen_indices = [i for i, atom in enumerate(atoms) if atom == 'N']
    
    print(f"Found {len(carbon_indices)} C atoms and {len(nitrogen_indices)} N atoms")
    
    if len(carbon_indices) != len(nitrogen_indices):
        print("Warning: Number of C atoms doesn't match number of N atoms!")
    
    # Helper function to calculate distance between two points
    def calculate_distance(p1, p2):
        return ((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2)**0.5
    
    # Helper function to calculate vector from point 1 to point 2
    def calculate_vector(p1, p2):
        return (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    
    # First, match in order (1st C with 1st N, 2nd C with 2nd N, etc.)
    ordered_bonds = []
    
    for c_idx, n_idx in zip(carbon_indices, nitrogen_indices):
        c_coords = coordinates[c_idx]
        n_coords = coordinates[n_idx]
        distance = calculate_distance(c_coords, n_coords)
        bond_vector = calculate_vector(c_coords, n_coords)
        
        ordered_bonds.append({
            'c_atom_index': c_idx,
            'n_atom_index': n_idx,
            'c_coords': c_coords,
            'n_coords': n_coords,
            'distance': distance,
            'vector': bond_vector
        })
    
    # Now find the closest N atom for each C atom
    closest_bonds = []
    used_n_indices = set()
    
    for c_idx in carbon_indices:
        c_coords = coordinates[c_idx]
        min_distance = float('inf')
        closest_n_idx = None
        closest_vector = None
        
        for n_idx in nitrogen_indices:
            if n_idx in used_n_indices:
                continue
                
            n_coords = coordinates[n_idx]
            distance = calculate_distance(c_coords, n_coords)
            
            if distance < min_distance:
                min_distance = distance
                closest_n_idx = n_idx
                closest_vector = calculate_vector(c_coords, n_coords)
        
        if closest_n_idx is not None:
            used_n_indices.add(closest_n_idx)
            closest_bonds.append({
                'c_atom_index': c_idx,
                'n_atom_index': closest_n_idx,
                'c_coords': c_coords,
                'n_coords': coordinates[closest_n_idx],
                'distance': min_distance,
                'vector': closest_vector
            })
    
    return {
        'ordered_bonds': ordered_bonds,
        'closest_bonds': closest_bonds,
        'num_carbon': len(carbon_indices),
        'num_nitrogen': len(nitrogen_indices)
    }

def display_cn_bonds(cn_bond_results):
    """
    Display the results of C-N bond calculations
    
    Parameters:
    -----------
    cn_bond_results : dict
        Dictionary containing C-N bond information from calculate_cn_bonds
    """
    print("\n--------------------------------------------")
    print("C-N BOND ANALYSIS RESULTS:")
    print("--------------------------------------------")
    print(f"Found {cn_bond_results['num_carbon']} C atoms and {cn_bond_results['num_nitrogen']} N atoms")
    
    # Display bonds matched in order (1st C with 1st N, etc.)
    print("\nC-N Bonds (matched in order):")
    for i, bond in enumerate(cn_bond_results['ordered_bonds'], 1):
        print(f"Bond #{i}: C(line {bond['c_atom_index'] + 3}) - N(line {bond['n_atom_index'] + 3})")
        print(f"  Distance: {bond['distance']:.6f} Å")
        print(f"  Bond vector: ({bond['vector'][0]:.6f}, {bond['vector'][1]:.6f}, {bond['vector'][2]:.6f})")
    
    # Display bonds matched by closest distance
    print("\nC-N Bonds (matched by closest distance):")
    for i, bond in enumerate(cn_bond_results['closest_bonds'], 1):
        print(f"Bond #{i}: C(line {bond['c_atom_index'] + 3}) - N(line {bond['n_atom_index'] + 3})")
        print(f"  Distance: {bond['distance']:.6f} Å")
        print(f"  Bond vector: ({bond['vector'][0]:.6f}, {bond['vector'][1]:.6f}, {bond['vector'][2]:.6f})")
    
    print("--------------------------------------------")

def main(xyz_file_path='vesta_file.xyz', specific_atom_index=None):
    try:
        # Parse the file
        print(f"Parsing XYZ file: {xyz_file_path}")
        atoms, coordinates = parse_xyz_file(xyz_file_path)
        print(f"Found {len(atoms)} atoms in the file")
        
        # Analyze C-N bonds if present
        if 'C' in atoms and 'N' in atoms:
            print("\nAnalyzing C-N bonds...")
            cn_bond_results = calculate_cn_bonds(atoms, coordinates)
            display_cn_bonds(cn_bond_results)
        
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
    
    # Calculate axial and equatorial angles
    print("\nCalculating Axial and Equatorial Angles...")
    angle_results = calculate_axial_and_equatorial_angles(central_pb_coords, nearest_I)
    
    # Calculate axial and equatorial bond lengths
    print("\nCalculating Axial and Equatorial Bond Lengths...")
    bond_length_results = calculate_axial_equatorial_bond_lengths(
        central_pb_coords, 
        nearest_I, 
        [nearest_I.index(next(i for i in nearest_I if i['atom_index'] == idx)) for idx in angle_results['axial_indices']],
        [nearest_I.index(next(i for i in nearest_I if i['atom_index'] == idx)) for idx in angle_results['equatorial_indices']]
    )
    
    # Display axial and equatorial bond lengths immediately and prominently
    print("\n--------------------------------------------")
    print("AXIAL AND EQUATORIAL BOND LENGTHS RESULTS:")
    print("--------------------------------------------")
    print("2 Axial Pb-I bonds:")
    for i, bond in enumerate(bond_length_results['axial_bond_lengths'], 1):
        print(f"  Axial bond #{i}: Pb to I(line {bond['atom_index'] + 2}): {bond['length']:.6f} Å")
    
    print("\n4 Equatorial Pb-I bonds:")
    for i, bond in enumerate(bond_length_results['equatorial_bond_lengths'], 1):
        print(f"  Equatorial bond #{i}: Pb to I(line {bond['atom_index'] + 2}): {bond['length']:.6f} Å")
    print("--------------------------------------------")
    
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
    
    # Display axial and equatorial angles results 
    print("\nAxial and Equatorial Angles Results:")
    print(f"Axial I atoms (line numbers in file): {[idx + 2 for idx in angle_results['axial_indices']]}")
    print(f"Axial I-Pb-I angle: {angle_results['axial_angle']:.6f}°")
    
    print(f"\nEquatorial I atoms (line numbers in file): {[idx + 2 for idx in angle_results['equatorial_indices']]}")
    
    # Display the 2 equatorial angles
    if len(angle_results['adjacent_equatorial_data']) >= 2:
        print("Equatorial I-Pb-I angles (for adjacent I atoms):")
        for i in range(2):  # Display both angles
            angle_data = angle_results['adjacent_equatorial_data'][i]
            atom_idx1, atom_idx2 = angle_data['atom_indices']
            print(f"  I(line {atom_idx1 + 2})-Pb-I(line {atom_idx2 + 2}): {angle_data['angle']:.6f}°")
            
        # Also display the opposite equatorial angle if available
        if angle_results['opposite_equatorial_data']:
            opp_data = angle_results['opposite_equatorial_data']
            atom_idx1, atom_idx2 = opp_data['atom_indices']
            print(f"\nOpposite equatorial I-Pb-I angle:")
            print(f"  I(line {atom_idx1 + 2})-Pb-I(line {atom_idx2 + 2}): {opp_data['angle']:.6f}°")
    elif angle_results['adjacent_equatorial_data']:
        print("Found only one equatorial I-Pb-I angle (for adjacent I atoms):")
        angle_data = angle_results['adjacent_equatorial_data'][0]
        atom_idx1, atom_idx2 = angle_data['atom_indices']
        print(f"  I(line {atom_idx1 + 2})-Pb-I(line {atom_idx2 + 2}): {angle_data['angle']:.6f}°")
    else:
        print("Could not identify adjacent equatorial I atoms.")
    
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
    print(f"Ideal center-to-vertex distance (l0): {quad_elongation_results['ideal_center_to_vertex_distance']:.6f} Å")
    print(f"Actual polyhedron volume: {quad_elongation_results['actual_volume']:.6f} Å³")
    print(f"Quadratic elongation <λ>: {quad_elongation_results['quadratic_elongation']:.6f}")
    
    print("\nBond Angle Variance Results:")
    if bond_angle_variance_results['angles_near_90']:
        print(f"Number of ~90° angles: {len(bond_angle_variance_results['angles_near_90'])}")
        print(f"90° angles: {[f'{a:.6f}' for a in bond_angle_variance_results['angles_near_90']]}")
    if bond_angle_variance_results['angles_near_180']:
        print(f"Number of ~180° angles: {len(bond_angle_variance_results['angles_near_180'])}")
        print(f"180° angles: {[f'{a:.6f}' for a in bond_angle_variance_results['angles_near_180']]}")
    print(f"Theoretical number of bond angles (m): {bond_angle_variance_results['theoretical_bond_angle_count']}")
    print(f"Actual number of bond angles: {bond_angle_variance_results['actual_bond_angle_count']}")
    print(f"Ideal bond angle (φ0): {bond_angle_variance_results['ideal_bond_angle']:.6f}°")
    print(f"Bond angle variance (σ²): {bond_angle_variance_results['bond_angle_variance']:.6f} deg.²")
    
    # Print a summary of key metrics with focus on individual bond lengths
    print("\nSummary of Key Metrics:")
    print("Axial Pb-I bonds:")
    for i, bond in enumerate(bond_length_results['axial_bond_lengths'], 1):
        print(f"  Axial bond #{i}: {bond['length']:.4f} Å")
    
    print("Equatorial Pb-I bonds:")
    for i, bond in enumerate(bond_length_results['equatorial_bond_lengths'], 1):
        print(f"  Equatorial bond #{i}: {bond['length']:.4f} Å")
        
    print(f"Polyhedral volume = {volume_results['total_volume']:.4f} Å^3")
    print(f"Distortion index (bond length) = {distortion_results['distortion_index']:.5f}")
    print(f"Quadratic elongation = {quad_elongation_results['quadratic_elongation']:.4f}")
    print(f"Bond angle variance = {bond_angle_variance_results['bond_angle_variance']:7.4f} deg.^2")
    print(f"Effective coordination number = {econ_results['effective_coordination_number']:.4f}")
    print(f"Axial I-Pb-I angle = {angle_results['axial_angle']:.4f}°")
    # Display the two equatorial angles in the summary
    if len(angle_results['adjacent_equatorial_data']) >= 2:
        print(f"Equatorial I-Pb-I angles = {angle_results['adjacent_equatorial_data'][0]['angle']:.4f}°, {angle_results['adjacent_equatorial_data'][1]['angle']:.4f}°")

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

