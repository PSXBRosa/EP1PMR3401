import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import sys

# plotagem com plotly
def plot_with_plotly(x,y,z):
    fig = go.Figure(data=[go.Surface(x=x,y=y,z=z,colorscale='Cividis')])
    fig.update_traces(contours_z=dict(show=True, usecolormap=True,
                                    highlightcolor="skyblue", project_z=True))
    fig.update_layout(title='Distribuição de Tensões', autosize=False,
                    scene_camera_eye=dict(x=1.87, y=0.88, z=-0.64),
                    width=1200, height=700,
                    margin=dict(l=65, r=50, b=65, t=90)
    )

    fig.show()

# plotagem com matplotlib
def plot_with_matplotlib(x,y,z, figname):
    fig = plt.figure(figsize=(21,7))
    ax = fig.add_subplot(projection='3d')
    ax.contour3D(x,y,z,100)
    # ax.plot_surface(x,y,z)
    plt.savefig(fname=figname)
    plt.show()

def prep(arr):
    arr = np.hstack((np.flip(arr,axis=1),arr))
    m,n = arr.shape

    # manipulações com o numpy
    r,p = np.linspace(R1,R2,m), np.linspace(-T,T,n)
    R,P = np.meshgrid(r,p)
    X = R*np.cos(P)
    Y = R*np.sin(P)
    return X,Y,arr

if __name__=='__main__':
    # preparações
    r1,r2,R1,R2,t,T = 0.05, 0.08, 0.03, 0.11, np.deg2rad(18), np.deg2rad(40)

    if len(sys.argv) == 1:
        file = 'mesh.txt'
        fig_name = file[:-4] + ' - plot'
        arr = np.loadtxt(file)
        X,Y,arr = prep(arr)
        plot_with_matplotlib(X,Y,arr.T, fig_name)
    else:
        for file in sys.argv[1:]:
            fig_name = file[:-4] + ' - plot'
            arr = np.loadtxt(file)
            X,Y,arr = prep(arr)
            plot_with_matplotlib(X,Y,arr.T, fig_name) 

    

