import numpy as np
import matplotlib.pyplot as plt

from STT_Tool import Tools
color = Tools.ColorList()

class TwoDomain():
    def __init__(self, dataPath, dir, FirstNum, SecondNum, mType='m'):
        
        data = np.loadtxt(dataPath)
        if mType=="m":
            data = data
        elif mType=="M":
            data = np.delete(data, [1, 2, 3], axis=1)
        
        if dir == "x":
            dirInd = 0
        elif dir == "y":
            dirInd = 1
        elif dir == "z":
            dirInd = 2

        self.dir, self.dirInd = dir, dirInd
        self.time  = data[:,0]
        self.m_tot = data[:,1:4]
        self.B_ext = data[:,4:7]

        self.data = data

        self.FirstNum = FirstNum
        self.SecondNum  = SecondNum
        self.N_num = FirstNum + SecondNum

    def Plot_totMH(self, xlim, ylim, fig=True,
                   posGrid=((1,1), (0, 0), 1, 1), 
                   width=3, 
                   wordsize=12,
                   save=[False, ]):
        m_tot = self.m_tot
        B_ext = self.B_ext
        # =====
        dir    = self.dir
        dirInd = self.dirInd
        # =====
        if fig:
            plt.figure(figsize=(4, 3))
        else:
            pass
        ax = plt.subplot2grid(posGrid[0], posGrid[1], 
                              rowspan=posGrid[2], colspan=posGrid[3]
                              )
        ax.plot(B_ext[:,dirInd], m_tot[:,dirInd], 
                label="Mtot", 
                color='k',
                linewidth=width
                )
        ax.set_xlabel(r"$B_{ext}$ (T)", fontsize=wordsize)
        ax.set_ylabel("Magnetization (A/m)", fontsize=wordsize)
        ax.set_xlim(xlim[0], xlim[1])
        if ylim==False:
            ylimArray = ax.set_ylim()
        else:
            ylimArray = ax.set_ylim(ylim[0], ylim[1])
        ax.legend(ncol=2, 
                bbox_to_anchor=(0.985, 0.05), loc=4, borderaxespad=0,
                fontsize=wordsize-4)
        ax.tick_params(axis='both', labelsize=(wordsize-2))
        ax.grid("--")
        plt.tight_layout()
        
        if save[0]:
            plt.savefig(save[1], 
                        transparent=True
                        )
        else:
            pass
        
        if fig:
            plt.show()
        else:
            pass
        
        return ax, ylimArray[0], ylimArray[1]
    
    def Plot_Arrow(self, datax, datay, dataz, FirstIgnore,
                   fig=True, 
                   posGrid=((1,1), (0, 0), 1, 1), 
                   arrowCol=("red", "green", "blue", "gray"),
                   ):

        FirstNum  = self.FirstNum
        SecondNum = self.SecondNum

        datax = np.array(datax)
        datay = np.array(datay)
        dataz = np.array(dataz)

        data_len = np.sqrt(datax**2+datay**2+dataz**2)

        dx = (datax/data_len)
        dy = (dataz/data_len)

        xpoint = dx/2*(-1)
        ypoint = np.linspace(1, len(xpoint)*2, len(xpoint))

        ypoint = ypoint-dy[::-1]/2
        ypoint = ypoint[::-1]

        if fig:
            plt.figure(figsize=(1,12))
        else:
            pass
        ax = plt.subplot2grid(posGrid[0], posGrid[1], 
                              rowspan=posGrid[2], colspan=posGrid[3]
                              )

        for i in range(len(xpoint)):
            
            if i < FirstNum:
                c = arrowCol[0]
            elif i < FirstNum+SecondNum:
                c = arrowCol[1]

            if i < FirstIgnore:
                pass
            else:
                ax.arrow(xpoint[i], ypoint[i], dx[i], dy[i], 
                        width=0.2, 
                        head_width=0.5, head_length=0.5,
                        length_includes_head=True,
                        color=c)
            
        ax.set_yticks(ticks=[], labels=[])
        ax.set_xticks(ticks=[], labels=[])
        ax.set_xlim(-1.5, 1.5)

        plt.tight_layout()
        
        if fig:
            plt.show()
        else:
            pass

        return ax

    def Plot_subMH(self, cutB, xlim, label, FirstIgnore,
                subPLot=[False,], forward=True, 
                arrowCol=("red", "green", "blue", "gray"),
                save=[False, ]
                ):

        data = self.data
        
        B_ext = self.B_ext
        N_num = self.N_num

        FirstNum  = self.FirstNum
        SecondNum = self.SecondNum
        
        # =====
        plt.figure(figsize=(16, 12))
        ax = self.Plot_totMH(xlim=xlim,
                             posGrid=((2,8), (0, 0), 1, 7), fig=False, 
                             width=5, 
                             wordsize=18, save=[False, ]
                             )
        # =====
        cutPointx  = []
        cutPointy  = []
        cutPointz  = []
        label_list = []
        for i in range(0, N_num):
            m_sub = data[:, (7+3*i):(7+3*(i+1))]

            if "all" in subPLot:
                ax.plot(B_ext[:,0], m_sub[:,0], label=f"SRO{i}", linewidth=3)
            elif i in subPLot:
                if len(subPLot)==1:
                    pass
                else:
                    ax.plot(B_ext[:,0], m_sub[:,0], label=f"SRO{i}", linewidth=3)

            # -----
            B_ext = np.round(B_ext, 3)
            cutB  = np.round(cutB, 3)
            if forward:
                specBPoint = np.min(np.where(B_ext[:,0]<=cutB)[0]) 
                # min => The first time reach cut_B
            else:
                specBPoint = np.max(np.where(B_ext[:,0]<=cutB)[0]) 
                # min => The last time reach cut_B
            
            cutPointx.append(m_sub[:,0][specBPoint])
            cutPointy.append(m_sub[:,1][specBPoint])
            cutPointz.append(m_sub[:,2][specBPoint])
            # -----
            label_list.append(f"{i}")
        
        # =========================
        cutPoint = cutPointx
        ax2 = plt.subplot2grid((2,8), (1, 0), rowspan=1, colspan=7)
        xtick = np.linspace(1, len(cutPoint), len(cutPoint))
        #  -----
        start = FirstIgnore
        end   = FirstNum
        ax2.scatter(xtick[start:end:], cutPoint[start:end:], 
                    label=label[0],color=arrowCol[0],
                    s=500
                    )
        ax2.vlines(x=end+0.5, ymin=-1.5, ymax=1.5, 
                    color='k',
                    linestyles='--')
        #  -----
        start = FirstNum
        end   = FirstNum + SecondNum
        ax2.scatter(xtick[start:end:], cutPoint[start:end:], 
                    label=label[1],color=arrowCol[1],
                    s=500
                    )
        ax2.vlines(x=end+0.5, ymin=-1.5, ymax=1.5, 
                    color='k',
                    linestyles='--')
        
        ax2.set_ylabel("Magnetization (A/m)", fontsize=18)
        
        ax2.set_xticks(ticks  = xtick[FirstIgnore:(FirstNum+SecondNum)], 
                       labels = label_list[FirstIgnore:(FirstNum+SecondNum)]
                    )
        ax2.set_xlabel("Layer Number", fontsize=18)
        ax2.tick_params(axis='both', labelsize=16)
        ax2.legend(loc="best", fontsize=14)
        ax2.grid("--")
        
        # =========================

        ax3 = self.Plot_Arrow(datax=cutPointx, datay=cutPointy, dataz=cutPointz,
                               fig=False, FirstIgnore=FirstIgnore,
                               posGrid=((2,8), (0, 7), 2, 1), arrowCol=arrowCol
                   )

        plt.tight_layout()
        if save[0]:
            plt.savefig(save[1], 
                        transparent=True
                        )
        else:
            plt.show()

class ThreeDomain():
    def __init__(self, dataPath, dir, FirstNum, SecondNum, ThirdNum, mType='m'):
        
        data = np.loadtxt(dataPath)
        if mType=="m":
            data = data
        elif mType=="M":
            data = np.delete(data, [1, 2, 3], axis=1)
        
        if dir == "x":
            dirInd = 0
        elif dir == "y":
            dirInd = 1
        elif dir == "z":
            dirInd = 2

        self.dir, self.dirInd = dir, dirInd
        self.time  = data[:,0]
        self.m_tot = data[:,1:4]
        self.B_ext = data[:,4:7]

        self.data = data

        self.FirstNum = FirstNum
        self.SecondNum  = SecondNum
        self.ThirdNum = ThirdNum
        self.N_num = FirstNum + SecondNum + ThirdNum

    def Plot_totMH(self, xlim, ylim, fig=True,
                   posGrid=((1,1), (0, 0), 1, 1), 
                   width=3, 
                   wordsize=12,
                   save=[False, ]):
        
        m_tot = self.m_tot
        B_ext = self.B_ext
        # =====
        dir    = self.dir
        dirInd = self.dirInd
        # =====
        if fig:
            plt.figure(figsize=(4, 3))
        else:
            pass
        ax = plt.subplot2grid(posGrid[0], posGrid[1], 
                              rowspan=posGrid[2], colspan=posGrid[3]
                              )
        ax.plot(B_ext[:,dirInd], m_tot[:,dirInd], 
                label="Mtot", 
                color='k',
                linewidth=width
                )
        ax.set_xlabel(r"$B_{ext}$ (T)", fontsize=wordsize)
        ax.set_ylabel("Magnetization (A/m)", fontsize=wordsize)
        ax.set_xlim(xlim[0], xlim[1])
        if ylim==False:
            ylimArray = ax.set_ylim()
        else:
            ylimArray = ax.set_ylim(ylim[0], ylim[1])

        ax.legend(ncol=2, 
                bbox_to_anchor=(0.985, 0.05), loc=4, borderaxespad=0,
                fontsize=wordsize-4)
        ax.tick_params(axis='both', labelsize=(wordsize-2))
        ax.grid("--")
        plt.tight_layout()
        
        if save[0]:
            plt.savefig(save[1], 
                        transparent=True
                        )
        else:
            pass
        
        if fig:
            plt.show()
        else:
            pass
        
        return ax, ylimArray[0], ylimArray[1]
    
    def Plot_Arrow(self, datax, datay, dataz, FirstIgnore,
                   fig=True, 
                   posGrid=((1,1), (0, 0), 1, 1), 
                   arrowCol=("red", "green", "blue", "gray"),
                   ):

        FirstNum   = self.FirstNum
        SecondNum  = self.SecondNum
        ThirdNum   = self.ThirdNum

        datax = np.array(datax)
        datay = np.array(datay)
        dataz = np.array(dataz)

        data_len = np.sqrt(datax**2+datay**2+dataz**2)

        dx = (datax/data_len)
        dy = (dataz/data_len)

        xpoint = dx/2*(-1)
        ypoint = np.linspace(1, len(xpoint)*2, len(xpoint))

        ypoint = ypoint-dy[::-1]/2
        ypoint = ypoint[::-1]

        if fig:
            plt.figure(figsize=(1,12))
        else:
            pass
        ax = plt.subplot2grid(posGrid[0], posGrid[1], 
                              rowspan=posGrid[2], colspan=posGrid[3]
                              )

        for i in range(len(xpoint)):
            
            if i < FirstNum:
                c = arrowCol[0]
            elif i < FirstNum+SecondNum:
                c = arrowCol[1]
            elif i < FirstNum+SecondNum+ThirdNum:
                c = arrowCol[2]

            if i < FirstIgnore:
                pass
            else:
                ax.arrow(xpoint[i], ypoint[i], dx[i], dy[i], 
                        width=0.2, 
                        head_width=0.5, head_length=0.5,
                        length_includes_head=True,
                        color=c)
            
        ax.set_yticks(ticks=[], labels=[])
        ax.set_xticks(ticks=[], labels=[])
        ax.set_xlim(-1.5, 1.5)

        plt.tight_layout()
        
        if fig:
            plt.show()
        else:
            pass

        return ax
    
    def Plot_subMH(self, cutB, xlim, label, FirstIgnore,
                subPLot=[False,], forward=True, 
                arrowCol=("red", "green", "blue", "gray"),
                save=[False, ]
                ):

        data = self.data
        
        B_ext = self.B_ext
        N_num = self.N_num

        FirstNum  = self.FirstNum
        SecondNum = self.SecondNum
        ThirdNum  = self.ThirdNum
        
        # =====
        plt.figure(figsize=(16, 12))
        ax = self.Plot_totMH(xlim=xlim,
                             posGrid=((2,8), (0, 0), 1, 7), fig=False, 
                             width=5, 
                             wordsize=18, save=[False,]
                             )
        # =====
        cutPointx  = []
        cutPointy  = []
        cutPointz  = []
        label_list = []
        for i in range(0, N_num):
            m_sub = data[:, (7+3*i):(7+3*(i+1))]

            if "all" in subPLot:
                ax.plot(B_ext[:,0], m_sub[:,0], label=f"SRO{i}", linewidth=3)
            elif i in subPLot:
                if len(subPLot)==1:
                    pass
                else:
                    ax.plot(B_ext[:,0], m_sub[:,0], label=f"SRO{i}", linewidth=3)

            # -----
            B_ext = np.round(B_ext, 3)
            cutB  = np.round(cutB, 3)
            if forward:
                specBPoint = np.min(np.where(B_ext[:,0]<=cutB)[0]) 
                # min => The first time reach cut_B
            else:
                specBPoint = np.max(np.where(B_ext[:,0]<=cutB)[0]) 
                # min => The last time reach cut_B
            
            cutPointx.append(m_sub[:,0][specBPoint])
            cutPointy.append(m_sub[:,1][specBPoint])
            cutPointz.append(m_sub[:,2][specBPoint])
            # -----
            label_list.append(f"{i}")
        
        # =========================
        cutPoint = cutPointx
        ax2 = plt.subplot2grid((2,8), (1, 0), rowspan=1, colspan=7)
        xtick = np.linspace(1, len(cutPoint), len(cutPoint))
        #  -----
        start = FirstIgnore
        end   = FirstNum
        ax2.scatter(xtick[start:end:], cutPoint[start:end:], 
                    label=label[0],color=arrowCol[0],
                    s=500
                    )
        ax2.vlines(x=end+0.5, ymin=-1.5, ymax=1.5, 
                    color='k',
                    linestyles='--')
        #  -----
        start = FirstNum
        end   = FirstNum + SecondNum
        ax2.scatter(xtick[start:end:], cutPoint[start:end:], 
                    label=label[1],color=arrowCol[1],
                    s=500
                    )
        ax2.vlines(x=end+0.5, ymin=-1.5, ymax=1.5, 
                    color='k',
                    linestyles='--')
        #  -----
        start = FirstNum + SecondNum
        end   = FirstNum + SecondNum + ThirdNum
        ax2.scatter(xtick[start:end:], cutPoint[start:end:], 
                    label=label[2], color=arrowCol[2],
                    s=500
                    )
        ax2.vlines(x=end+0.5, ymin=-1.5, ymax=1.5, 
                    color='k',
                    linestyles='--')
        
        ax2.set_ylabel("Magnetization (MA/m)", fontsize=22)
        
        ax2.set_xticks(ticks  = xtick[FirstIgnore:(FirstNum+SecondNum)], 
                       labels = label_list[FirstIgnore:(FirstNum+SecondNum)]
                    )
        # ax2.set_xlabel("Layer Number", fontsize=18)
        ax2.tick_params(axis='both', labelsize=20)
        ax2.legend(loc="best", fontsize=18)
        ax2.grid("--")
        
        # =========================

        ax3 = self.Plot_Arrow(datax=cutPointx, datay=cutPointy, dataz=cutPointz,
                               fig=False, FirstIgnore=FirstIgnore,
                               posGrid=((2,8), (0, 7), 2, 1), arrowCol=arrowCol
                   )

        plt.tight_layout()
        if save[0]:
            plt.savefig(save[1], 
                        transparent=True
                        )
        else:
            plt.show()

class FourDomain():
    def __init__(self, dataPath, dir, FirstNum, SecondNum, ThirdNum, FourNum, mType='m'):
        
        data = np.loadtxt(dataPath)
        if mType=="m":
            data = data
        elif mType=="M":
            data = np.delete(data, [1, 2, 3], axis=1)
        
        if dir == "x":
            dirInd = 0
        elif dir == "y":
            dirInd = 1
        elif dir == "z":
            dirInd = 2

        self.dir, self.dirInd = dir, dirInd
        self.time  = data[:,0]
        self.m_tot = data[:,1:4]
        self.B_ext = data[:,4:7]

        self.data = data

        self.FirstNum = FirstNum
        self.SecondNum  = SecondNum
        self.ThirdNum = ThirdNum
        self.FourNum = FourNum
        self.N_num = FirstNum + SecondNum + ThirdNum + FourNum

    def Plot_totMH(self, xlim, ylim, fig=True,
                   posGrid=((1,1), (0, 0), 1, 1), 
                   width=3, 
                   wordsize=12,
                   save=[False, ]):
        
        m_tot = self.m_tot
        B_ext = self.B_ext
        # =====
        dir    = self.dir
        dirInd = self.dirInd
        # =====
        if fig:
            plt.figure(figsize=(4, 3))
        else:
            pass
        ax = plt.subplot2grid(posGrid[0], posGrid[1], 
                            rowspan=posGrid[2], colspan=posGrid[3]
                            )
        ax.plot(B_ext[:,dirInd], m_tot[:,dirInd], 
                label="Mtot", 
                color='k',
                linewidth=width
                )
        ax.set_xlabel(r"$B_{ext}$ (T)", fontsize=wordsize)
        ax.set_ylabel("Magnetization (A/m)", fontsize=wordsize)
        ax.set_xlim(xlim[0], xlim[1])
        if ylim==False:
            ylimArray = ax.set_ylim()
        else:
            ylimArray = ax.set_ylim(ylim[0], ylim[1])

        ax.legend(ncol=2, 
                bbox_to_anchor=(0.985, 0.05), loc=4, borderaxespad=0,
                fontsize=wordsize-4)
        ax.tick_params(axis='both', labelsize=(wordsize-2))
        ax.grid("--")
        plt.tight_layout()
        if save[0]:
            plt.savefig(save[1], 
                        transparent=True
                        )
        else:
            pass
        if fig:
            plt.show()
        else:
            pass
        
        return ax, ylimArray[0], ylimArray[1]
    
    def Plot_Arrow(self, datax, datay, dataz, FirstIgnore,
                fig=True, 
                posGrid=((1,1), (0, 0), 1, 1), 
                arrowCol=("red", "green", "blue", "gray"),
                ):

        FirstNum   = self.FirstNum
        SecondNum  = self.SecondNum
        ThirdNum   = self.ThirdNum
        FourNum    = self.FourNum

        datax = np.array(datax)
        datay = np.array(datay)
        dataz = np.array(dataz)

        data_len = np.sqrt(datax**2+datay**2+dataz**2)

        dx = (datax/data_len)
        dy = (dataz/data_len)

        xpoint = dx/2*(-1)
        ypoint = np.linspace(1, len(xpoint)*2, len(xpoint))

        ypoint = ypoint-dy[::-1]/2
        ypoint = ypoint[::-1]

        if fig:
            plt.figure(figsize=(1,12))
        else:
            pass
        ax = plt.subplot2grid(posGrid[0], posGrid[1], 
                            rowspan=posGrid[2], colspan=posGrid[3]
                            )

        for i in range(len(xpoint)):
            
            if i < FirstNum:
                c = arrowCol[0]
            elif i < FirstNum+SecondNum:
                c = arrowCol[1]
            elif i < FirstNum+SecondNum+ThirdNum:
                c = arrowCol[2]
            elif i < FirstNum+SecondNum+ThirdNum+FourNum:
                c = arrowCol[3]

            if i < FirstIgnore:
                pass
            else:
                ax.arrow(xpoint[i], ypoint[i], dx[i], dy[i], 
                        width=0.2, 
                        head_width=0.5, head_length=0.5,
                        length_includes_head=True,
                        color=c)
            
        ax.set_yticks(ticks=[], labels=[])
        ax.set_xticks(ticks=[], labels=[])
        ax.set_xlim(-1.5, 1.5)

        plt.tight_layout()
        
        if fig:
            plt.show()
        else:
            pass

        return ax
    
    def Plot_subMH(self, cutB, xlim, label, FirstIgnore,
                subPLot=[False,], forward=True, 
                arrowCol=("red", "green", "blue", "gray"),
                save=[False, ]
                ):

        data = self.data
        
        B_ext = self.B_ext
        N_num = self.N_num

        FirstNum  = self.FirstNum
        SecondNum = self.SecondNum
        ThirdNum  = self.ThirdNum
        FourNum   = self.FourNum
        
        # =====
        plt.figure(figsize=(16, 12))
        ax = self.Plot_totMH(xlim=xlim,
                            posGrid=((2,8), (0, 0), 1, 7), fig=False, 
                            width=5, 
                            wordsize=18, save=[False,]
                            )
        # =====
        cutPointx  = []
        cutPointy  = []
        cutPointz  = []
        label_list = []
        for i in range(0, N_num):
            m_sub = data[:, (7+3*i):(7+3*(i+1))]

            if "all" in subPLot:
                ax.plot(B_ext[:,0], m_sub[:,0], label=f"SRO{i}", linewidth=3)
            elif i in subPLot:
                if len(subPLot)==1:
                    pass
                else:
                    ax.plot(B_ext[:,0], m_sub[:,0], label=f"SRO{i}", linewidth=3)

            # -----
            B_ext = np.round(B_ext, 3)
            cutB  = np.round(cutB, 3)
            if forward:
                specBPoint = np.min(np.where(B_ext[:,0]<=cutB)[0]) 
                # min => The first time reach cut_B
            else:
                specBPoint = np.max(np.where(B_ext[:,0]<=cutB)[0]) 
                # min => The last time reach cut_B
            
            cutPointx.append(m_sub[:,0][specBPoint])
            cutPointy.append(m_sub[:,1][specBPoint])
            cutPointz.append(m_sub[:,2][specBPoint])
            # -----
            label_list.append(f"{i}")
        
        # =========================
        cutPoint = cutPointx
        ax2 = plt.subplot2grid((2,8), (1, 0), rowspan=1, colspan=7)
        xtick = np.linspace(1, len(cutPoint), len(cutPoint))
        #  -----
        start = FirstIgnore
        end   = FirstNum
        ax2.scatter(xtick[start:end:], cutPoint[start:end:], 
                    label=label[0],color=arrowCol[0],
                    s=500
                    )
        ax2.vlines(x=end+0.5, ymin=-1.5, ymax=1.5, 
                    color='k',
                    linestyles='--')
        #  -----
        start = FirstNum
        end   = FirstNum + SecondNum
        ax2.scatter(xtick[start:end:], cutPoint[start:end:], 
                    label=label[1],color=arrowCol[1],
                    s=500
                    )
        ax2.vlines(x=end+0.5, ymin=-1.5, ymax=1.5, 
                    color='k',
                    linestyles='--')
        #  -----
        start = FirstNum + SecondNum
        end   = FirstNum + SecondNum + ThirdNum
        ax2.scatter(xtick[start:end:], cutPoint[start:end:], 
                    label=label[2], color=arrowCol[2],
                    s=500
                    )
        ax2.vlines(x=end+0.5, ymin=-1.5, ymax=1.5, 
                    color='k',
                    linestyles='--')
        #  -----
        start = FirstNum + SecondNum + ThirdNum
        end   = FirstNum + SecondNum + ThirdNum + FourNum
        ax2.scatter(xtick[start:end:], cutPoint[start:end:], 
                    label=label[2], color=arrowCol[2],
                    s=500
                    )
        ax2.vlines(x=end+0.5, ymin=-1.5, ymax=1.5, 
                    color='k',
                    linestyles='--')
        # -----
        
        ax2.set_ylabel("Magnetization (A/m)", fontsize=18)
        
        ax2.set_xticks(ticks  = xtick[FirstIgnore:(FirstNum+SecondNum)], 
                    labels = label_list[FirstIgnore:(FirstNum+SecondNum)]
                    )
        ax2.set_xlabel("Layer Number", fontsize=18)
        ax2.tick_params(axis='both', labelsize=16)
        ax2.legend(loc="best", fontsize=14)
        ax2.grid("--")
        
        # =========================

        ax3 = self.Plot_Arrow(datax=cutPointx, datay=cutPointy, dataz=cutPointz,
                            fig=False, FirstIgnore=FirstIgnore,
                            posGrid=((2,8), (0, 7), 2, 1), arrowCol=arrowCol
                )

        plt.tight_layout()
        if save[0]:
            plt.savefig(save[1], 
                        transparent=True
                        )
        else:
            plt.show()