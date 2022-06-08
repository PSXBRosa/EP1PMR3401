import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

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
def plot_with_matplotlib(x,y,z):
    fig = plt.figure(figsize=(21,7))
    ax = fig.add_subplot(projection='3d')
    ax.contour3D(x,y,z,100)
    plt.show()

if __name__=='__main__':
    # preparações
    r1,r2,R1,R2,t,T = 0.05, 0.08, 0.03, 0.11, np.deg2rad(18), np.deg2rad(40)
    arr = np.loadtxt('T.txt')
    arr = np.hstack((np.flip(arr,axis=1),arr))
    m,n = arr.shape

    # manipulações com o numpy
    r,p = np.linspace(R1,R2,m), np.linspace(-T,T,n)
    R,P = np.meshgrid(r,p)
    X = R*np.cos(P)
    Y = R*np.sin(P)

    plot_with_matplotlib(X,Y,arr.T)
