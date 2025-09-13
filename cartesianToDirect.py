import numpy as np
def cartToDir(POSCAR):
    with open(POSCAR, "r") as f:
        lines = f.readlines()
        title = lines[0]
        scale = lines[1].strip()
        lattice = np.loadtxt(lines[2:5])
        atoms = lines[5:7]
        coordinateSystem = lines[7]
        Rc = np.loadtxt(lines[8:])
    
    newLattice = float(scale)*lattice.T 
    latticeInv = np.linalg.inv(newLattice)
    Rd = np.array([np.dot(latticeInv,vec) for vec in Rc])
    
    with open(POSCAR, "w") as fout:
        fout.write(title)
        fout.write(scale)
        fout.write("\n")
        np.savetxt(fout, lattice)
        fout.writelines(atoms)
        fout.write("Direct\n")
        np.savetxt(fout, Rd)




poscars = ["POSCAR"+str(i) for i in range(1,63)]
for i in poscars:
    cartToDir(i)