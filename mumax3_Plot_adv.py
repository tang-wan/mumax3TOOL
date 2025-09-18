import numpy as np
import matplotlib.pyplot as plt

import matplotlib.animation as ama
from tqdm.autonotebook import tqdm
import matplotlib.colors as mcolors


# from STT_Tool import Tools
# color = Tools.ColorList()
print("=====")
print("We have the following methods:")
methodList = ("DataPlot_ama()",)
for name in methodList:
    print(f"    {name}") 

class DataPlot_ama():
    def __init__(self, dataPath,
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
            
    def ALLDataLoad(self):
        dataDict = self.dataDict
        return dataDict

    def MHLoop(self):
        pass

    def Trajectory(self):
        pass

    def Trajectory_3D(self, time="t", minput="m", mrange=False, 
                            sep=10, file_name=None, 
                            color="red", method="line", taillen=5):
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

        fig = plt.figure(figsize=(4, 3))
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
        dataDict = self.dataDict
        if mrange:
            datax  = dataDict[time][mrange[0]:mrange[1]][::sep]
        else:
            datax  = dataDict[time][::sep]
        n = len(datax)
        x = []
    #-----------------------------------------------
        
        if mrange:
            datay_x = dataDict[minput][:,0][mrange[0]:mrange[1]][::sep]
            datay_y = dataDict[minput][:,1][mrange[0]:mrange[1]][::sep]
            datay_z = dataDict[minput][:,2][mrange[0]:mrange[1]][::sep]
        else:
            datay_x = dataDict[minput][:,0][::sep]
            datay_y = dataDict[minput][:,1][::sep]
            datay_z = dataDict[minput][:,2][::sep]
            # print(len(datay_x))
    #----------------------------------------------
        match method.lower():
            case "line":
                yx, yy, yz = [], [], []
                line1, = ax.plot(yx, yy, yz, color=color)
                #----------------------------------------------
                def run(fram_ind):
                    
                    x.append(datax[fram_ind])

                    yx.append(datay_x[fram_ind])
                    yy.append(datay_y[fram_ind])
                    yz.append(datay_z[fram_ind])
                    line1.set_data_3d(yx, yy, yz)
            case "scatter":
                point, = ax.plot([], [], [], marker="o", color=color, markersize=6)
                #----------------------------------------------
                def run(fram_ind):
                    point.set_data_3d([datay_x[fram_ind]],
                                    [datay_y[fram_ind]],
                                    [datay_z[fram_ind]])
                    return point,
            case "scatter tail":
                base_rgb = mcolors.to_rgba(color)[:3]   # 取得 RGB，alpha 我們會在下面設定

                # main
                point, = ax.plot([], [], [], marker="o", color=color, markersize=6, zorder=3)

                # tail
                trail_length = taillen   # 你可以改這個數字來調整尾巴長短
                trail_lines = [ax.plot([], [], [], linewidth=2, solid_capstyle='round')[0]
                            for _ in range(trail_length)]

                def run(fram_ind):
                    # 更新主角點位置
                    px = datay_x[fram_ind]
                    py = datay_y[fram_ind]
                    pz = datay_z[fram_ind]
                    point.set_data_3d([px], [py], [pz])

                    # 取出要顯示在尾巴上的點（最多 trail_length+1 個點）
                    start = max(0, fram_ind - trail_length)
                    xs = datay_x[start:fram_ind+1]
                    ys = datay_y[start:fram_ind+1]
                    zs = datay_z[start:fram_ind+1]

                    seg_count = max(0, len(xs) - 1)  # 有幾條線段要畫
                    if seg_count == 0:
                        # 尚無足夠點，全部隱藏
                        for ln in trail_lines:
                            ln.set_data_3d([], [], [])
                    else:
                        # alpha 從最舊 (頭) 透明 -> 最靠近點 (尾) 實心
                        alphas = np.linspace(0.1, 1.0, seg_count)
                        for i, ln in enumerate(trail_lines):
                            if i < seg_count:
                                # 每個 line 代表一段連線 (p_i -> p_{i+1})
                                xseg = [xs[i], xs[i+1]]
                                yseg = [ys[i], ys[i+1]]
                                zseg = [zs[i], zs[i+1]]
                                ln.set_data_3d(xseg, yseg, zseg)
                                rgba = (base_rgb[0], base_rgb[1], base_rgb[2], float(alphas[i]))
                                ln.set_color(rgba)
                            else:
                                # 超出目前尾巴長度的 line 隱藏
                                ln.set_data_3d([], [], [])

                    # 回傳被更新的 artists（FuncAnimation 不設定 blit 也可）
                    return (point, *trail_lines)
            case "arrow":
                arrow = [ax.quiver(0, 0, 0, 
                        datay_x[0], datay_y[0], datay_z[0], 
                        color=color, arrow_length_ratio=0.2)]

                def run(fram_ind):
                    # 清掉舊的箭頭
                    arrow[0].remove()
                    # 畫新的箭頭：尾巴在 (0,0,0)，頭在 (x,y,z)
                    arrow[0] = ax.quiver(0, 0, 0,
                                        datay_x[fram_ind], 
                                        datay_y[fram_ind], 
                                        datay_z[fram_ind],
                                        color=color,
                                        linewidth=2, arrow_length_ratio=0.3
                                        )
                    return arrow
            case "arrow_3d":
                def create_arrow_mesh(length=1.0, radius=0.05, head_length=0.2, head_radius=0.1, n=30):
                    """ 建立一個沿 z 軸的箭頭 (圓柱+圓錐)，回傳 (X,Y,Z) 座標網格 """
                    # 圓柱部分
                    z_cyl = np.linspace(0, length-head_length, 2)
                    theta = np.linspace(0, 2*np.pi, n)
                    Zc, Tc = np.meshgrid(z_cyl, theta)
                    Xc = radius * np.cos(Tc)
                    Yc = radius * np.sin(Tc)

                    # 圓錐部分
                    z_cone = np.linspace(length-head_length, length, 2)
                    Zcone, Tcone = np.meshgrid(z_cone, theta)
                    r_cone = np.linspace(head_radius, 0, 2)
                    Rcone, _ = np.meshgrid(r_cone, theta)
                    Xcone = Rcone * np.cos(Tcone)
                    Ycone = Rcone * np.sin(Tcone)

                    return (Xc, Yc, Zc), (Xcone, Ycone, Zcone)

                def rotation_matrix_from_vectors(vec1, vec2):
                    """ 計算把 vec1 轉成 vec2 的旋轉矩陣 """
                    a, b = (vec1 / np.linalg.norm(vec1)), (vec2 / np.linalg.norm(vec2))
                    v = np.cross(a, b)
                    c = np.dot(a, b)
                    s = np.linalg.norm(v)
                    if s == 0:
                        return np.eye(3)
                    kmat = np.array([[0, -v[2], v[1]],
                                    [v[2], 0, -v[0]],
                                    [-v[1], v[0], 0]])
                    return np.eye(3) + kmat + kmat.dot(kmat) * ((1-c)/(s**2))

                def transform_mesh(X, Y, Z, R):
                    """ 將 mesh (X,Y,Z) 用旋轉矩陣 R 轉換 """
                    shape = X.shape
                    coords = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
                    coords_rot = R @ coords
                    Xr, Yr, Zr = coords_rot
                    return Xr.reshape(shape), Yr.reshape(shape), Zr.reshape(shape)

                # 建立沿 z 軸的箭頭 mesh
                cyl, cone = create_arrow_mesh(length=1.0)

                # 預設方向 = (0,0,1)
                arrow_parts = []
                Xc, Yc, Zc = cyl
                Xcone, Ycone, Zcone = cone
                arrow_parts.append(ax.plot_surface(Xc, Yc, Zc, color="r", alpha=0.8))
                arrow_parts.append(ax.plot_surface(Xcone, Ycone, Zcone, color="r", alpha=0.8))

                def run(fram_ind):
                    # 刪掉舊箭頭
                    for surf in arrow_parts:
                        surf.remove()
                    arrow_parts.clear()

                    # 目標方向
                    vec = np.array([datay_x[fram_ind], datay_y[fram_ind], datay_z[fram_ind]])

                    # 箭頭長度 = 向量大小
                    length = np.linalg.norm(vec)
                    if length < 1e-6:
                        return []

                    cyl, cone = create_arrow_mesh(length=length)
                    R = rotation_matrix_from_vectors(np.array([0,0,1]), vec)

                    # 旋轉 mesh
                    Xc, Yc, Zc = transform_mesh(*cyl, R)
                    Xcone, Ycone, Zcone = transform_mesh(*cone, R)

                    # 畫新的箭頭
                    arrow_parts.append(ax.plot_surface(Xc, Yc, Zc, color="r", alpha=0.8))
                    arrow_parts.append(ax.plot_surface(Xcone, Ycone, Zcone, color="r", alpha=0.8))

                    return arrow_parts
            case _:
                run = None   # 先預設
                print("!!!! No this method !!!!")

        if run is not None:
            fig.tight_layout()
            bar = tqdm(total=n)
            ani = ama.FuncAnimation(fig, run, frames=n, interval=1)
            ani.save(f'{file_name}', fps=n, progress_callback=lambda i, n: bar.update(1))
            bar.close()
            plt.show()

