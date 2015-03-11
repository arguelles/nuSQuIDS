import numpy as np
import matplotlib.pyplot as plt
import nuSQUIDSpy as nsq

# we will use the boost injector metaclass

BoostPythonMetaclass = nsq.nuSQUIDS.__class__

class injector(object):
    class __metaclass__(BoostPythonMetaclass):
        def __init__(self, name, bases, dict):
            for b in bases:
                if type(b) not in (self, type):
                    for k,v in dict.items():
                        setattr(b,k,v)
            return type.__init__(self, name, bases, dict)

# functions to construct SU(N) algebra representation matrix

def __ExtShape(A):
    Ar = A.shape[0];
    if len(A.shape) == 1:
        Ac = 1;
    else:
        Ac = A.shape[1]
    return (Ar,Ac);

def __DirectSum(A,B):
    SA = __ExtShape(A);
    SB = __ExtShape(B);
    # direct sum of numpy arrays
    Z = np.zeros(shape = (SA[0]+SB[0],
                          SA[1]+SB[1]));

    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if i < SA[0] and j < SA[1]:
                if SA[1] == 1 and j == 0:
                    Z[i,j] = A[i];
                elif SA[1] == 1 and j!= 0 :
                    Z[i,j] = 0;
                else :
                    Z[i,j] = A[i,j];
            elif i >= SA[0] and j >= SA[0]:
                if SB[1] == 1 and j == SA[1]:
                    Z[i,j] = B[i - SA[0]];
                elif SB[1] == 1 and j!= SA[1] :
                    Z[i,j] = 0;
                else :
                    Z[i,j] = B[i - SA[0] ,j - SA[1] ];
    return Z;

def __EArray(j,k,d):
    # auxiliary E arrays
    E = np.zeros(shape = (d,d));
    E[j,k] = 1.0;
    return E;

def _FArray(j,k,d):
    # rules to generate SU(N) matrixes
    if d == 1 :
        return np.array([1.0]);
    elif j == k and j != (d-1) and j != 0 :
        return __DirectSum(_FArray(j,j,d-1), np.array([0]));
    elif j == k == (d-1) :
        return np.sqrt(2.0/(d*(d-1))) * __DirectSum( _FArray(0,0,d-1), np.array([1-d]) );
    elif j == k == 0 :
        return np.identity(d);
    elif k < j :
        return __EArray(j,k,d) + __EArray(k,j,d);
    elif k > j :
        return np.complex(0,-1) * ( __EArray(j,k,d) - __EArray(k,j,d) );

# just a simple matrix rotation thingy
    def R(self,i,j,cp,param):
        """ Rotation Matrix
        Calculates the R_ij rotations. Also incorporates CP-phases when necesary.
        @type   i       :   int
        @param  i       :   i-column.
        @type   j       :   int
        @param  j       :   j-row.
        @type   cp      :   int
        @param  cp      :   if cp = 0 : no CP-phase. else CP-phase = CP_array[cp]

        @rtype          :   numpy.array
        @return         :   returns the R_ij rotation matrix.
        """
        # if cp = 0 -> no complex phase
        # R_ij, i<j
        if(j<i):
            # no funny business
            l = i
            i = j
            j = l

        # rotation matrix generator
        R = np.zeros([param.numneu,param.numneu],complex)
        # diagonal terms
        for k in np.arange(0,param.numneu,1):
            if(k != i-1 and k != j-1):
                R[k,k] = 1.0
            else :
                R[k,k] = param.c[i,j]
        # non-diagonal terms
        if(cp != 0):
            sd = np.sin(param.dcp[cp])
            cd = np.cos(param.dcp[cp])
            faseCP = complex(cd,sd)
        else:
            faseCP = complex(1.0,0.0)

        R[i-1,j-1] = param.s[i,j]*faseCP.conjugate()
        R[j-1,i-1] = -param.s[i,j]*faseCP
        return R

# nuSQUIDS python extensions

class ExtNuSQUIDS(injector,nsq.nuSQUIDS):
    """ nuSQUIDS class extension """

    def SUBasis(self):
        nsun = self.GetNumNeu();
        basis = [];
        for i in range(nsun):
            for j in range(nsun):
                basis.append(_FArray(j,i,nsun));
        return np.array(basis);

    def PrintStatus(self,filepath):
        """ Print the nuSQUIDS status to a file"""
        file = open(filepath,"w")
        file.write('% -- nuSQuIDS Automatically Generated Report -- \n')
        file.write('\\documentclass{article}\n')
        file.write('\\usepackage{fullpage}\n')

        file.write('\\begin{document}\n')

        file.write('\\section{Basic}\n')
        file.write('Number of neutrinos:\t{0}'.format(self.GetNumNeu()))
        file.write('\\section{Physics Parameters}\n')
        file.write('\\subsection{Mixing Angles}\n')
        file.write('\\subsection{CP Phases}\n')
        file.write('\\subsection{Square Mass Differences}\n')
        file.write('\\subsection{Mixing Matrix}\n')
        file.write('\\section{Hamiltonian}\n')

        file.write('\\end{document}\n')
        file.close()

    def PrintStatus(self):
        return 0

    def PlotFlavor(self):
        return 0

    def PlotMass(self):
        return 0

    def PlotComposition(self):
        return 0

    def _PrintSUVector(self,suv):
        basis = self.SUBasis();
        suv_array = np.zeros((self.GetNumNeu(),self.GetNumNeu()));
        for i,c in enumerate(suv.GetComponents()):
            suv_array += c*basis[i];
        return suv_array;

    def PrintH0(self,E):
        """
            Returns the H0 hamiltonian at a given energy
        """
        return self._PrintSUVector(self.H0(E));

    def PrintHI(self,ie):
        """
            Returns the HI hamiltonian at a given energy node
            in the interaction-massbasis.
        """
        if ie < 0 or ie > self.GetNumE():
            print "Error::Invalid energy node number.";
            exit();
        return self._PrintSUVector(self.HI(ie));

#    def PrintHamiltonian(self,ie):
#        """
#            Returns full Hamiltonian at a given energy node
#            in the mass basis.
#        """
#        if ie < 0 or ie > self.GetNumE():
#            print "Error::Invalid energy node number.";
#        return self._PrintSUVector(self.GetHamiltonian(ie,rho));

    def PrintMixingMatrix(self):
        return 0;

    def PrintState(self,ie,rho = 0):
        """
            Returns full Hamiltonian at a given energy node
            in the mass basis.
        """
        if ie < 0 or ie > self.GetNumE():
            print "Error::Invalid energy node number.";

        return self._PrintSUVector(self.GetState(ie,rho));

class ExtNuSQUIDSAtm(injector,nsq.nuSQUIDSAtm):
    """ nuSQUIDSAtm class extension """
    def PlotFlavor(self,flavor,neutype,
            cosmin = None, cosmax = None, enumin = None, enumax = None,
            enu_step = 400.0, costh_step = 100.0, contour_divisions = 100, colorscale = "log",
            fontsize = 16, fontname = "Arial", colormap = "jet" ):

        if colorscale == "log":
            get_flux = lambda costh,log_enu : np.log10(self.EvalFlavor(flavor,costh,10**log_enu,neutype))
        else:
            get_flux = lambda costh,log_enu : self.EvalFlavor(flavor,costh,10**log_enu,neutype)
        vgfnumu = np.vectorize(get_flux)

        erange = self.GetERange()
        enu_min = np.log10(erange[0]+0.1)
        enu_max = np.log10(erange[-1]-0.1)

        cosrange = self.GetCosthRange()
        costh_min = cosrange[0]
        costh_max = cosrange[-1]

        if cosmin != None:
            costh_min = cosmin
        if cosmax != None:
            costh_max = cosmax
        if enumin != None:
            enu_min = enumin
        if enumax != None:
            enu_max = enumax

        step_enu = (enu_max - enu_min)/enu_step
        step_costh = (costh_max - costh_min)/costh_step

        nenu = np.arange(enu_min,enu_max,step_enu)
        ncosth = np.arange(costh_min,costh_max,step_costh)

        NCOSTH, NENU = np.meshgrid(ncosth, nenu)
        PHINUMU = vgfnumu(NCOSTH,NENU)

        fig = plt.figure(figsize = (10,8))

        ax = plt.subplot()


        plt.ylim(enu_min,enu_max)
        plt.xlim(costh_min,costh_max)

        color_map = plt.get_cmap(colormap)

        plt.contourf(NCOSTH,NENU,PHINUMU,contour_divisions, cmap = color_map)

        plt.xlabel(r"$\cos\ \theta_z$", fontsize = fontsize, fontname = fontname)
        plt.ylabel(r"$\mathrm{log}_{10}(E/\mathrm{GeV})$", fontsize = fontsize, fontname = fontname)

        cbar = plt.colorbar(cmap = color_map)
        if colorscale == "log":
            cbar.set_label(r"${\rm log}_{10}({\rm Weight})$",fontsize = fontsize, fontname = fontname)
        else:
            cbar.set_label(r"{\rm Weight}", fontsize = fontsize, fontname = fontname)

        for label in (ax.get_xticklabels() + ax.get_yticklabels() + \
                      cbar.ax.get_yticklabels() + cbar.ax.get_xticklabels()):
            label.set_fontname(fontname)
            label.set_fontsize(fontsize)
        #return fig
        #plt.close()

