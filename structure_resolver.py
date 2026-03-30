import numpy as np
import os
import sys


def read_structure(path_structure, delta_x=0, delta_y=0, delta_z=0):
    with open(path_structure, 'r') as file:
        lines = file.readlines()

        scaling_factor = float(lines[1].strip())

        lattice_vectors = [lines[i].strip().split() for i in range(2, 5)]
        lattice_vectors = scaling_factor * (np.array(lattice_vectors).astype(float))

        atomic_species = lines[5].strip().split()
        number_of_atoms = np.array(lines[6].strip().split()).astype(int)
        tot = sum(number_of_atoms)

        lattice_type = lines[7].strip()

        atomic_positions = [lines[i].strip().split() for i in range(8, 8 + tot)]
        atomic_positions = np.array(atomic_positions).astype(float)

        atoms = dict(zip(atomic_species, number_of_atoms))

    if lattice_type != "Direct":
        latticeInv = np.linalg.inv(lattice_vectors.T)
        Rd = np.array([np.dot(latticeInv, vec) for vec in atomic_positions])
        atomic_positions = Rd

    atomic_positions[:, 0] += delta_x
    atomic_positions[:, 1] += delta_y
    atomic_positions[:, 2] += delta_z
    atomic_positions = np.dot(atomic_positions, lattice_vectors)

    return atomic_positions


def write_structure(R, path_structure, path_save="contcar.vasp"):
    with open(os.path.expanduser(path_structure), "r") as f:
        header = f.readlines()[:7]

    with open(os.path.expanduser(path_save), "w") as f:
        f.writelines(header)
        f.write("Cartesian\n")
        np.savetxt(f, R)

    return None


if len(sys.argv) < 2:
    print("Usage: python script.py CONTCAR [dx dy dz]")
    sys.exit(1)

path_structure = sys.argv[1]

# Default shifts
delta_x, delta_y, delta_z = 0.0, 0.0, 0.0

# If user provides shifts
if len(sys.argv) == 5:
    delta_x = float(sys.argv[2])
    delta_y = float(sys.argv[3])
    delta_z = float(sys.argv[4])

elif len(sys.argv) != 2:
    print("Error: Provide either 0 or 3 shift values (dx dy dz)")
    sys.exit(1)

R = read_structure(path_structure, delta_x, delta_y, delta_z)
write_structure(R, path_structure)