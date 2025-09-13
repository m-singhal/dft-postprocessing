import numpy as np
def generateKpointsHSE(KPATH = "KPATH.in", IBZKPT = "IBZKPT"):
    sym = []
    with open(KPATH, "r") as f:
        lines = f.readlines()
        nod = int(lines[1].strip())
        for line in lines[4:]:
            if len(line.strip().split()) > 0:
                sym.append(line.strip().split())
    high_sym = [sym[0][:3]]
    label = [sym[0][-1]]
    for i in range(2,len(sym),2):
        high_sym.append(sym[i][:3])
        label.append(sym[i][-1])
    high_sym.append(sym[-1][:3])
    label.append(sym[-1][-1])

    high_sym = np.array(high_sym, dtype=np.float32)
    lam = 1/nod
    kpoints = []
    for i in range(high_sym.shape[0]-1):
        for j in range(nod):
            kpoints.append(high_sym[i] + j*lam*(high_sym[i+1] - high_sym[i]))
    kpoints.append(high_sym[-1])
    for i in range(len(kpoints)):
        kpoints[i] = list(kpoints[i])
        kpoints[i].append(0.0)
    kpoints = np.array(kpoints, dtype=np.float32)
    

    with open(IBZKPT, "r") as f:
        lines = f.readlines()
        ibzkpt = [line.strip().split() for line in lines[3:]]
    ibzkpt = np.array(ibzkpt, dtype=np.float32)

    with open("KPOINTS", "w") as f:
        f.write("K-point mesh for HSE calculation\n")
        f.write(str(int(ibzkpt.shape[0] + kpoints.shape[0])) + "\n")
        f.write("Reciprocal lattice\n")
        np.savetxt(f, ibzkpt, fmt='%10.8f')
        np.savetxt(f, kpoints, fmt='%10.8f')
    print("KPOINTS file generated.")
generateKpointsHSE()