import numpy as np
import matplotlib.pyplot as plt

class TB_DataLoad():
    def __init__(self, dataPath):
        
        self.Rdata = np.loadtxt(dataPath)
    
    def ALLDataLoad(self):
        Rdata = self.Rdata
        return Rdata

class TB_DataPlot():
    def __init__(self, dataPath):
        
        self.Rdata = np.loadtxt(dataPath)
    
    def Plot(self, xinput,
            plotpos=111, 
            line=("black", 1, 3, "mz"), line2=("-", "o"),
            xylabel=(False, False, 14), 
            xytick=((False, False), (False, False), 12)):
        # ===== Data Load =====
        Rdata = self.Rdata
        # ===== Basis Plot =====
        ax = plt.subplot(plotpos)
        ax.plot(xinput, Rdata, 
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