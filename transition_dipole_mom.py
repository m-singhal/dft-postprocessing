import numpy as np
import subprocess
import os

# Define matrix size
bi = 140
bf = 146
nbands = bf - bi + 1  # Total number of bands
tdm_matrix_up = np.zeros((nbands, nbands))
tdm_matrix_down = np.zeros((nbands, nbands))

# Loop through i < j
for idxi, i in enumerate(range(bi, bf)):       # 1-based indexing (1 to 191)
    idxj = idxi  # Start from the same index for j
    for j in range(i + 1, bf + 1):
        idxj += 1  # (i+1) to 192
        

        # Create VASPKIT input
        vaspkit_input = f"""
71
713
0
{i} {j}
1
"""

        with open("vaspkit.in", "w") as f:
            f.write(vaspkit_input.strip())

        # Run VASPKIT
        subprocess.run("vaspkit < vaspkit.in", shell=True)

        # Check for the output file
        if not os.path.isfile("TDM_COMPONENTS_UP.dat"):
            print(f"Warning: Missing TDM_COMPONENTS_UP.dat for bands {i}-{j}")
            continue

        # Parse TDM total (4th column)
        with open("TDM_COMPONENTS_UP.dat", "r") as f:
            lines = f.readlines()

        # The value is in the 2nd data line (4th row, 4th column)
        try:
            total = float(lines[3].split()[3])
            tdm_matrix_up[idxi, idxj] = total
            tdm_matrix_up[idxj, idxi] = total  # Symmetric
        except Exception as e:
            print(f"Failed parsing for {i}-{j}: {e}")
            continue

        if not os.path.isfile("TDM_COMPONENTS_DW.dat"):
            print(f"Warning: Missing TDM_COMPONENTS_DW.dat for bands {i}-{j}")
            continue

        # Parse TDM total (4th column)
        with open("TDM_COMPONENTS_DW.dat", "r") as f:
            lines = f.readlines()

        # The value is in the 2nd data line (4th row, 4th column)
        try:
            total = float(lines[3].split()[3])
            tdm_matrix_down[idxi, idxj] = total
            tdm_matrix_down[idxj, idxi] = total  # Symmetric
        except Exception as e:
            print(f"Failed parsing for {i}-{j}: {e}")
            continue

# Save the matrix
np.savetxt("tdm_matrix_up.txt", tdm_matrix_up, fmt="%.6e")
np.savetxt("tdm_matrix_dw.txt", tdm_matrix_down, fmt="%.6e")
