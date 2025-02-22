# print("hello world")

def find_nearest_pb_atom(coordinates, atoms, center_coords):
    """
    Find the nearest Pb atom to the center coordinates if the centre coordinate is not Pb atom 
    
    Args:
        coordinates: List of tuples containing (x, y, z) coordinates for all atoms
        atoms: List of atom names
        center_coords: Tuple of (x, y, z) for the center point
    """
    min_distance = float('inf')
    nearest_pb_index = None
    center_x, center_y, center_z = center_coords
    
    # Find nearest Pb atom
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



def calculate_center_coordinates(filepath):
    try:
        with open(filepath, 'r') as file:
            # Skip first two lines
            next(file)  # Skip number of atoms
            next(file)  # Skip SCF Done line
            
            # Initialize lists for coordinates and atoms
            x_coords = []
            y_coords = []
            z_coords = []
            atoms = []  # List to store atom names
            
            # Process each line
            for line in file:
                if line.strip():  # Skip empty lines
                    parts = line.split()
                    if len(parts) == 4:  # Check if line has coordinates 
                        try:
                            atoms.append(parts[0])  # Store atom name
                            x_coords.append(float(parts[1]))
                            y_coords.append(float(parts[2]))
                            z_coords.append(float(parts[3]))
                        except ValueError:
                            continue  # Skip if conversion to float fails
            
            # Calculate center coordinates
            if x_coords:
                center_x = sum(x_coords) / len(x_coords)
                center_y = sum(y_coords) / len(y_coords)
                center_z = sum(z_coords) / len(z_coords)
                
                # Find the atom closest to the center
                min_distance = float('inf')
                center_atom_index = 0
                
                for i in range(len(x_coords)):
                    distance = ((x_coords[i] - center_x)**2 + 
                              (y_coords[i] - center_y)**2 + 
                              (z_coords[i] - center_z)**2)**0.5
                    if distance < min_distance:
                        min_distance = distance
                        center_atom_index = i
                
                # print("\nCenter coordinates:")
                # print(f"X: {center_x:.6f}")
                # print(f"Y: {center_y:.6f}")
                # print(f"Z: {center_z:.6f}")
                # print(f"Number of atoms processed: {len(x_coords)}")
                # print(f"Atom closest to center: {atoms[center_atom_index]}")
                # print(f"Distance from exact center: {min_distance:.6f}")
                
                return {
                    'center_coordinates': (center_x, center_y, center_z),
                    'center_atom': atoms[center_atom_index],
                    'distance_from_center': min_distance,
                    'num_atoms': len(x_coords)
                }
            
            return None
            
    except FileNotFoundError:
        print("Error: File not found")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None



# function to read the txt file 
def read_xyz_file(filepath):
    with open('vesta_file.xyz', 'r') as file:
        # Read all lines
        lines = file.readlines()
        return lines

# Usage example
try:
    result = calculate_center_coordinates('vesta_file.xyz')
    if result:
        if result['center_atom'] == "Pb":
            print("\nCenter Atom is Pb:")
            print(f"Center coordinates: {result['center_coordinates']}")
            print(f"Center atom: {result['center_atom']}")
            print(f"Distance from center: {result['distance_from_center']}")
            print(f"Total atoms: {result['num_atoms']}")
        else:
            print(f"\nCenter atom is {result['center_atom']}, finding nearest Pb atom...")
            nearest_pb = find_nearest_pb_atom(
                coordinates=result['coordinates'],
                atoms=result['atoms'],
                center_coords=result['center_coordinates']
            )
            
            if nearest_pb:
                print("\nNearest Pb Atom Found:")
                pb_coords = nearest_pb['coordinates']
                print(f"Coordinates: X: {pb_coords[0]:.6f}, Y: {pb_coords[1]:.6f}, Z: {pb_coords[2]:.6f}")
                print(f"Distance from center: {nearest_pb['distance']:.6f}")
            else:
                print("No Pb atoms found in the structure")
    else:
        print("No results found from coordinate calculation")
except Exception as e:
    print(f"An error occurred: {e}")


