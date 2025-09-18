import numpy as np
import matplotlib.pyplot as plt

# from STT_Tool import Tools
# color = Tools.ColorList()
print("=====")
print("We have the following methods:")
methodList = ("DataLoad()", "DataPlot()")
for name in methodList:
    print(f"    {name}") 
print("=====")
print("If want ama please use Plot_adv")

class DataLoad():
    def __init__(self, dataPath, Border=2, 
                 BasisParameter=[('t', 1), ('m', 3), ('Bext', 3)],
                 FurthParameter=[('Region-m1', 3), ('Region-m2', 3)]
                 ):
        
        data = np.loadtxt(dataPath)
        print(np.shape(data))

        dataDict = {}

        StartNum = 0
        for name, Dnum in BasisParameter:
            EndNum = StartNum+Dnum
            dataDict[name] = data[:,(StartNum):(EndNum)]
            StartNum = EndNum

        for name, Dnum in FurthParameter:
            EndNum = StartNum+Dnum
            dataDict[name] = data[:,(StartNum):(EndNum)]
            StartNum = EndNum
        
        self.BasisParameter = BasisParameter
        self.FurthParameter = FurthParameter
        self.dataDict = dataDict
        # ==========
        print("We have the following methods:")
        methodList = ("ALLDataLoad()", "AngleLoad()")
        for name in methodList:
            print(f"    {name}") 

    def ALLDataLoad(self):
        dataDict = self.dataDict
        return dataDict
    
    def _AngleCal(self, mx, my, mz, zeroAxis):
        if zeroAxis == 'z':
            r = np.sqrt(mx**2+my**2+mz**2)
            # -----
            theta  = np.arctan((np.sqrt(mx**2+my**2))/abs(mz))*(180/np.pi)
            if mx>0 and mz>0:
                theta = theta
            elif mx>0 and mz<0:
                theta = 180 - theta
            elif mx<0 and mz>0:
                theta = -theta
            elif mx<0 and mz<0:
                theta = -180 + theta
            elif mx==0 and mz<0:
                theta = 180
            elif mx==0 and mz>0:
                theta = 0
            # -----    
            if mx == 0:
                phi =0
            else:
                phi = np.arctan(my/mx)*180/np.pi
        else:
            pass
        
        return r, theta, phi
    
    def AngleLoad(self, BextAxis='z', LoadDir='F', zeroAxis='z',
                        CalParameter=['m', 'Region-m1', 'Region-m2']):
        
        dataDict = self.dataDict
        BextDict = {}
        BextArray = dataDict['Bext']
        match BextAxis:
            case 'x':
                BextInd = 0
            case 'y':
                BextInd = 1
            case 'z':
                BextInd = 2
            case _:
                print("Wrong Input, plz check the BextAxis")
        
        if LoadDir.upper()=='F':
            start = 0
            end   = (len(BextArray)//2)+1
        elif LoadDir.upper()=='B':
            start = (len(BextArray)//2)
            end   = (len(BextArray))+1
        else:
            print("!! No This direction !!")
        
        for i, Barray in enumerate(BextArray[start:end]):
            Car_r_array  = np.array([])
            theta_array  = np.array([])
            phi_array    = np.array([])
            for Cal in CalParameter:
                data = dataDict[Cal][start+i]
                mx, my, mz = data

                r, theta, phi = self._AngleCal(mx, my, mz, zeroAxis)
                
                Car_r_array = np.append(Car_r_array, r)
                theta_array = np.append(theta_array, theta)
                phi_array   = np.append(phi_array, phi)

            Car_r_array = Car_r_array.reshape(len(Car_r_array)//3, 3)
            BextDict[str(Barray[BextInd])] = {'Car_r': Car_r_array,
                                              'theta': theta_array, 
                                              'phi': phi_array}
            self.BextDict = BextDict
            
        return BextDict
            

class DataPlot():        
    def __init__(self, dataPath, Border=2, 
                BasisParameter=[('t', 1), ('m', 3), ('Bext', 3)],
                FurthParameter=[('Region-m1', 3), ('Region-m2', 3)]
                ):
        
        data = np.loadtxt(dataPath)
        print(np.shape(data))

        dataDict = {}

        StartNum = 0
        for name, Dnum in BasisParameter:
            EndNum = StartNum+Dnum
            dataDict[name] = data[:,(StartNum):(EndNum)]
            StartNum = EndNum

        for name, Dnum in FurthParameter:
            EndNum = StartNum+Dnum
            dataDict[name] = data[:,(StartNum):(EndNum)]
            StartNum = EndNum
        
        self.BasisParameter = BasisParameter
        self.FurthParameter = FurthParameter
        self.dataDict = dataDict
        # ==========
        print("We have the following methods:")
        methodList = ("MHLoop()", "Trajectory()", "Trajectory_Plot3D()")
        for name in methodList:
            print(f"    {name}") 

    def MHLoop(self, Bext="Bext", Bext_dir=2,
                     minput="m", m_dir=2, rescale=1,
                     plotpos=111, 
                     line=("black", 1, 3, "mz"), line2=("-", "o"),
                     xylabel=(False, False, 14), 
                     xytick=((False, False), (False, False), 12)):
        
        # ===== Data Load =====
        dataDict = self.dataDict
        Bext = dataDict[Bext][:,Bext_dir]
        m = dataDict[minput][:,m_dir]*rescale
        # ===== Basis Plot =====
        ax = plt.subplot(plotpos)
        ax.plot(Bext, m, 
                color=line[0], alpha=line[1], linewidth=line[2], label=line[3],
                linestyle=line2[0], marker=line2[1]
                )
        ax.grid("--")
        # ===== Label =====
        if xylabel[0]:
            ax.set_xlabel(xylabel[0], fontsize=xylabel[2])
        else:
            pass
    
        if xylabel[1]:
            ax.set_ylabel(xylabel[1], fontsize=xylabel[2])
        else:
            pass
        # ===== Tict =====
        if xytick[0][0] and xytick[0][1]:
            ax.set_xticks(ticks=xytick[0][0], labels=xytick[0][1], fontsize=xytick[2])
        elif xytick[0][1]: # Only ticks no labels
            ax.set_xticklabels("")
        else:
            pass

        if xytick[1][0] and xytick[1][1]:
            ax.set_yticks(ticks=xytick[1][0], labels=xytick[1][1], fontsize=xytick[2])
        elif xytick[1][1]: # Only ticks no labels
            ax.set_yticklabels("")
        else:
            pass

        return ax
    
    def Trajectory(self, time="t",
                        minput="m", m_dir=2, rescale=1,
                        plotpos=111, 
                        line=("black", 1, 3, "mz"), line2=("-", "o"),
                        xylabel=(False, False, 14), 
                        xytick=((False, False), (False, False), 12)):
        # ===== Data Load =====
        dataDict = self.dataDict
        time = dataDict[time]
        m = dataDict[minput][:,m_dir]*rescale
        # ===== Basis Plot =====
        ax = plt.subplot(plotpos)
        ax.plot(time, m, 
                color=line[0], alpha=line[1], linewidth=line[2], label=line[3],
                linestyle=line2[0], marker=line2[1]
                )
        ax.grid("--")
        # ===== Label =====
        if xylabel[0]:
            ax.set_xlabel(xylabel[0], fontsize=xylabel[2])
        else:
            pass
    
        if xylabel[1]:
            ax.set_ylabel(xylabel[1], fontsize=xylabel[2])
        else:
            pass
        # ===== Tict =====
        if xytick[0][0] and xytick[0][1]:
            ax.set_xticks(ticks=xytick[0][0], labels=xytick[0][1], fontsize=xytick[2])
        elif xytick[0][1]: # Only ticks no labels
            ax.set_xticklabels("")
        else:
            pass

        if xytick[1][0] and xytick[1][1]:
            ax.set_yticks(ticks=xytick[1][0], labels=xytick[1][1], fontsize=xytick[2])
        elif xytick[1][1]: # Only ticks no labels
            ax.set_yticklabels("")
        else:
            pass

        return ax

    def Trajectory_3D(self, input_ax=False, 
                            minput="m", mrange=False, 
                            color="red"):
        '''
        These codes were created by Lise
        '''
        #===========================
        def set_axes_equal(ax: plt.Axes):
            """Set 3D plot axes to equal scale.

            Make axes of 3D plot have equal scale so that spheres appear as
            spheres and cubes as cubes.  Required since `ax.axis('equal')`
            and `ax.set_aspect('equal')` don't work on 3D.
            """
            limits = np.array([
                ax.get_xlim3d(),
                ax.get_ylim3d(),
                ax.get_zlim3d(),
            ])
            origin = np.mean(limits, axis=1)
            radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
            _set_axes_radius(ax, origin, radius)
            ax.set_box_aspect([1,1,1])
        def _set_axes_radius(ax, origin, radius):
            x, y, z = origin
            ax.set_xlim3d([x - radius, x + radius])
            ax.set_ylim3d([y - radius, y + radius])
            ax.set_zlim3d([z - radius, z + radius])
        def create_sphere_mesh(ntheta, nphi):
            u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:60j]
            x = np.cos(u)*np.sin(v)
            y = np.sin(u)*np.sin(v)
            z = np.cos(v)
            return x, y, z
        #===========================
        def plot_sphere_frame(ax, num=10, **kwargs):
            theta = np.linspace(0, 2*np.pi, 100)
            for angle in np.linspace(-np.pi, np.pi, num):
                xlat = np.cos(theta)*np.cos(angle)
                ylat = np.sin(theta)*np.cos(angle)
                zlat = np.ones_like(theta)*np.sin(angle)
                xlon = np.sin(angle)*np.sin(theta)
                ylon = np.cos(angle)*np.sin(theta)
                zlon = np.cos(theta)
                ax.plot(xlat, ylat, zlat, **kwargs)
                ax.plot(xlon, ylon, zlon, **kwargs)
        def plot_sphere_surface(ax, ntheta=10, nphi=20, **kwargs):
            sphere_mesh = create_sphere_mesh(ntheta, nphi)
            ax.plot_surface(*sphere_mesh, **kwargs)
        def plot_sphere_axis(ax):
            ax.plot([-1,1], [ 0,0], [ 0,0], color='gray')
            ax.plot([ 0,0], [-1,1], [ 0,0], color='gray')
            ax.plot([ 0,0], [ 0,0], [-1,1], color='gray')
            ax.text(1.1, 0, 0, 'x')
            ax.text(0, 1.1, 0, 'y')
            ax.text(0, 0, 1.1, 'z')
        #===========================
        dataDict = self.dataDict

        if input_ax:
            ax = input_ax
        else:
            ax = plt.subplot2grid((1,1), (0,0), projection='3d')

            plot_sphere_frame(ax, num=10, color='black', alpha=0.05)

            plot_sphere_surface(ax, ntheta=30, nphi=60, color='#808080', alpha=0.05)
            
            plot_sphere_axis(ax)

            ax.view_init(elev=30, azim=30)

            ax.set_xbound(-0.6, 0.6)
            ax.set_ybound(-0.6, 0.6)
            ax.set_zbound(-0.6, 0.6)
            ax.set_axis_off()
            set_axes_equal(ax)
    #-----------------------------------------------
        if mrange:
            mx = dataDict[minput][:,0][mrange[0]:mrange[1]]
            my = dataDict[minput][:,1][mrange[0]:mrange[1]]
            mz = dataDict[minput][:,2][mrange[0]:mrange[1]]
        else:
            mx = dataDict[minput][:,0]
            my = dataDict[minput][:,1]
            mz = dataDict[minput][:,2]
        
        ax.plot(mx, my, mz, color=color)
    #-----------------------------------------------
        
        return ax
