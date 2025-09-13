import numpy as np 
import matplotlib.pyplot as plt

class PyBands():

    def __init__(self):
        
        self.nspin = 1
        with open("INCAR", "r") as f:
            lines = f.readlines()
            for line in lines:
                if "ISPIN" in line and "2" in line:
                    self.nspin = 2
                    break

        with open("EIGENVAL", "r") as f:
            lines = f.readlines()
            self.title = lines[4].strip()
            params = lines[5].strip().split()
            self.nelec, self.nk, self.nbands = int(params[0]), int(params[1]), int(params[2])
            
            
            if self.nspin == 2:
                self.ncol = 5
            else:
                self.ncol = 3
           
            self.kpoints = np.zeros((self.nk, 4))
            self.bands = np.zeros((self.nk,self.nbands, self.ncol))

            a = 7
            self.stride = self.nbands+2
            b = (self.nk)*(self.stride)
            for index, i in enumerate(range(a,b, self.stride)):
                self.kpoints[index] = np.array(lines[i].strip().split(), dtype = np.float32)
                self.bands[index] = np.array([line.strip().split() for line in lines[i+1:i+self.stride-1]], dtype = np.float32)

        
        
        with open("OUTCAR", "r") as f:
            for line in f:
                if "E-fermi" in line:
                    self.fermi = float(line.strip().split()[2])
                    break


        if self.nspin == 2:
            
            homo_up = self.bands[:,:,1].min()
            lumo_up = self.bands[:,:,1].max()
            homo_up_k = 0
            lumo_up_k = 0
            homo_up_index = 0
            lumo_up_index = 0
            for i in range(self.nk):
                up = np.where(self.bands[i,:,3] == 0.0)[0][0]
                if self.bands[i,up-1,1] > homo_up:
                    homo_up = self.bands[i,up-1,1]
                    homo_up_k = i
                    homo_up_index = up-1
                if self.bands[i,up,1] < lumo_up:
                    lumo_up = self.bands[i,up,1]
                    lumo_up_k = i
                    lumo_up_index = up
        
            homo_down = self.bands[:,:,2].min()
            lumo_down = self.bands[:,:,2].max()
            homo_down_k = 0
            lumo_down_k = 0
            homo_down_index = 0
            lumo_down_index = 0
            for i in range(self.nk):
                down = np.where(self.bands[i,:,4] == 0.0)[0][0]
                if self.bands[i,down-1,2] > homo_down:
                    homo_down = self.bands[i,down-1,2]
                    homo_down_k = i
                    homo_down_index = down-1
                if self.bands[i,down,2] < lumo_down:
                    lumo_down = self.bands[i,down,2]
                    lumo_down_k = i
                    lumo_down_index = down

            min_fermi = max(homo_up, homo_down, self.fermi)
            max_fermi = min(lumo_up, lumo_down)

            if min_fermi != self.fermi:
                self.fermi = min_fermi + ((max_fermi - min_fermi)/10)
            self.homo_up = homo_up
            self.homo_up_index = homo_up_index
            self.homo_up_k = homo_up_k
            self.lumo_up = lumo_up
            self.lumo_up_index = lumo_up_index
            self.lumo_up_k = lumo_up_k
            self.homo_down = homo_down
            self.homo_down_index = homo_down_index
            self.homo_down_k = homo_down_k
            self.lumo_down = lumo_down
            self.lumo_down_index = lumo_down_index
            self.lumo_down_k = lumo_down_k
        
        else:
            homo = self.bands[:,:,1].min()
            lumo = self.bands[:,:,1].max()
            homo_k = 0
            lumo_k = 0
            homo_index = 0
            lumo_index = 0
            for i in range(self.nk):
                up = np.where(self.bands[i,:,2] == 0)[0][0]
                if self.bands[i,up-1,1] > homo:
                    homo = self.bands[i,up-1,1]
                    homo_k = i
                    homo_index = up-1
                if self.bands[i,up,1] < lumo:
                    lumo = self.bands[i,up,1]
                    lumo_k = i
                    lumo_index = up
            
            min_fermi = max(homo, self.fermi)
            max_fermi = lumo

            if min_fermi != self.fermi:
                self.fermi = min_fermi + ((max_fermi - min_fermi)/10)
            self.homo = homo
            self.homo_index = homo_index
            self.homo_k = homo_k
            self.lumo = lumo
            self.lumo_index = lumo_index
            self.lumo_k = lumo_k
        
        
                
    def bandStructure(self, Emin=None, Emax=None, KPOINTS="KPOINTS", markerSize=1, mode = "line", HSE=False):
        
        if HSE:
            weight0 = np.where(self.kpoints[:,-1]==0)[0][0]
            self.kpoints = self.kpoints[weight0:]
            self.bands = self.bands[weight0:]
            self.nk = self.kpoints.shape[0]
        
        

        labels = []
        high_sym_dist = []
        with open(KPOINTS, "r") as f:
            lines = f.readlines()
            self.nk_bet = int(lines[1].strip().split()[0])
            for index,line in enumerate(lines[1:]):
                line = line.strip().split()
                if len(line) > 2:
                    if len(labels) == 0:
                        labels.append(line[-1])
                    elif line[-1] != labels[-1]:
                        labels.append(line[-1])
                    
                    if len(line[-1]) > 1:
                        high_sym_dist.append(line[:3])
            high_sym_dist = np.array(high_sym_dist, dtype=np.float32)
            high_sym_dist = np.reshape(high_sym_dist, (int(high_sym_dist.shape[0]/2), 2, 3))
                    
        self.kpath = [0.0]
        for i in range(1,self.nk):
            if any(
                np.allclose(self.kpoints[i-1][:3], seg[0], atol=1e-08) and
                np.allclose(self.kpoints[i][:3], seg[1], atol=1e-08)
                for seg in high_sym_dist
            ):                
                dist = 0
            else:
                dist = np.linalg.norm(self.kpoints[i][:3] - self.kpoints[i-1][:3])
            self.kpath.append(self.kpath[-1] + dist)
        self.kpath = np.array(self.kpath)

        self.klabels = []
        for index, i in enumerate(range(0, len(self.kpath), self.nk_bet)):
            self.klabels.append((float(self.kpath[i]), labels[index]))
        self.klabels.append((float(self.kpath[-1]), labels[-1]))

        xticks, xlabels = zip(*self.klabels)

       
       

        if self.nspin == 2:   
            print(f"""SPIN UP: 
            Band Gap (eV) = {self.lumo_up - self.homo_up}
            Homo band index = {self.homo_up_index+1}, k-points = {self.kpoints[self.homo_up_k,:3]}, Energy (eV) = {self.homo_up}
            Lumo band index = {self.lumo_up_index+1}, k-points = {self.kpoints[self.lumo_up_k,:3]}, , Energy (eV) = {self.lumo_up}
            
            Fermi Energy (ev) = {self.fermi}
            """)
            
            print(f"""SPIN DOWN: 
            Band Gap (eV) = {self.lumo_down - self.homo_down}
            Homo band index = {self.homo_down_index+1}, k-points = {self.kpoints[self.homo_down_k,:3]}, , Energy (eV) = {self.homo_down}
            Lumo band index = {self.lumo_down_index+1}, k-points = {self.kpoints[self.lumo_down_k,:3]}, , Energy (eV) = {self.lumo_down}
            """)
            
            
                        
            fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
            for i in range(self.nbands):    
                if mode=="line":
                    axs[0].plot(self.kpath, self.bands[:,i,1], color='blue', label="Spin Up")
                else:
                    axs[0].scatter(self.kpath, self.bands[:,i,1], s=markerSize, color='blue', label="Spin Up")
            axs[0].axhline(y=self.fermi, color='black', linestyle='--')
            axs[0].set_title(self.title + " (SPIN UP)")
            axs[0].set_xticks(xticks)
            axs[0].text(x=plt.xlim()[-1]+0.05, y=self.fermi, s=r"$E_f$", va='bottom', ha='right')
            axs[0].set_xticklabels(xlabels)
            axs[0].set_ylabel("Energy (eV)")
            axs[0].set_ylim(Emin,Emax)
            for i in xticks:
                axs[0].axvline(x=i, color='grey', linestyle='--', alpha=0.2)
            
            
            
            for i in range(self.nbands):    
                if mode=="line":
                    axs[1].plot(self.kpath, self.bands[:,i,2], color='red', label="Spin Down")
                else:
                    axs[1].scatter(self.kpath, self.bands[:,i,2], s=markerSize, color='red', label="Spin Down")
            axs[1].axhline(y=self.fermi, color='black', linestyle='--')
            axs[1].set_title(self.title + " (SPIN DOWN)")
            axs[1].set_xticks(xticks)
            axs[1].text(x=plt.xlim()[-1]+0.05, y=self.fermi, s=r"$E_f$", va='bottom', ha='right')
            axs[1].set_xticklabels(xlabels)
            axs[1].set_ylabel("Energy (eV)")
            axs[1].set_ylim(Emin,Emax)
            for i in xticks:
                axs[1].axvline(x=i, color='grey', linestyle='--', alpha=0.2)
            
            plt.tight_layout()
            plt.show()

            
        
        else:
            print(f"""
            Band Gap (eV) = {self.lumo - self.homo}
            Homo band index = {self.homo_index+1}, k-points = {self.kpoints[self.homo_k,:3]}, Energy (eV) = {self.homo}
            Lumo band index = {self.lumo_index+1}, k-points = {self.kpoints[self.lumo_k,:3]}, Energy (eV) = {self.lumo}
            Fermi Energy (ev) = {self.fermi}
            """)
            for i in range(self.nbands):    
                if mode=="line":
                    plt.plot(self.kpath, self.bands[:,i,1], color='red')
                else:
                    plt.scatter(self.kpath, self.bands[:,i,1], s=markerSize, color='red')
            plt.axhline(y=self.fermi, color='black', linestyle='--')
            plt.title(self.title)
            plt.xticks(xticks, xlabels)
            plt.text(x=plt.xlim()[-1]+0.05, y=self.fermi, s=r"$E_f$", va='bottom', ha='right')
            plt.ylabel("Energy (eV)")
            plt.ylim(Emin,Emax)
            for i in xticks:
                plt.axvline(x=i, color='grey', linestyle='--', alpha=0.2)
            plt.show()
            
    

    def DOS(self, Emin=None, Emax= None,
            polarization = 4 # 0 - spin up, 1 - spin down, 2 - both, 3 - total dos, 4 - all three
            ):
        
        with open("DOSCAR", "r") as f:
            lines = f.readlines()
            title = lines[4].strip()
            nedos = int(lines[5].strip().split()[2])
            data = np.loadtxt(lines[6:6+nedos])
        

        E = data[:,0]
        print(f"Fermi energy (eV) = {self.fermi}")
        if self.nspin == 2:
            print(f"""
                SPIN UP:   
                Band Gap (eV) = {self.lumo_up - self.homo_up}
                SPIN DOWN:   
                Band Gap (eV) = {self.lumo_down - self.homo_down}
""")
            self.dos_up = data[:,1]
            self.dos_down = data[:,2]
            self.dos = self.dos_up + self.dos_down
            if polarization == 0:
                plt.plot(self.dos_up,E, label = "SPIN UP")
            elif polarization == 1:
                plt.plot(self.dos_down,E, label = "SPIN DOWN")
            elif polarization == 2:
                plt.plot(self.dos_up,E, label = "SPIN UP")
                plt.plot(self.dos_down,E, label = "SPIN DOWN", alpha=0.5)
            elif polarization == 3:
                plt.plot(self.dos,E, label = "TOTAL")
            else:
                plt.plot(self.dos_up,E, label = "SPIN UP")
                plt.plot(self.dos_down,E, label = "SPIN DOWN", alpha=0.5)
                plt.plot(self.dos,E, label = "TOTAL", alpha=0.2)
            plt.xlabel("DOS")
            plt.axhline(y=self.fermi, color='black', linestyle='--')
            plt.ylabel("Energy (eV)")
            plt.text(x=plt.xlim()[-1]+0.5, y=self.fermi, s=r"$E_f$", va='bottom', ha='right')
            plt.ylim(Emin, Emax)
            plt.title(title)
            plt.legend()
            plt.show()
        
        else:
            print(f"Band Gap (eV) = {self.lumo - self.homo}")
            self.dos = data[:,1]
            plt.plot(self.dos,E)
            plt.xlabel("DOS")
            plt.axhline(y=self.fermi, color='black', linestyle='--')
            plt.ylabel("Energy (eV)")
            plt.text(x=plt.xlim()[-1]+0.5, y=self.fermi, s=r"$E_f$", va='bottom', ha='right')
            plt.ylim(Emin, Emax)
            plt.title(title)
            plt.show()


pb = PyBands()

pb.DOS(2,8,polarization=4)
pb.bandStructure(2,8,mode="line")
