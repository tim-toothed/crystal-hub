import numpy as np
import matplotlib.pyplot as plt

from ..calculation import calculation
from ..Magnetics import Magnetics
from .tables_legends_conventions import color_atoms, COLORS_list
from .crystal_field import Freeion_charge_dist

#######FIGURES###########

def plot_energy_levels(eigenvalues, ax=None, color='b', label=None, tolerance=0.05, offset=0, delta=0, ylabel=None):
    """
    Plot energy levels (e.g., crystal field levels), optionally labeling them on the x-axis.

    Parameters:
    - eigenvalues: list or array of energy levels to plot.
    - ax: matplotlib Axes object. If None, a new figure and axis will be created.
    - color: color used to draw energy levels.
    - label: label to associate with this set of energy levels, shown on the x-axis at the given offset.
    - tolerance: max difference to consider levels as degenerate.
    - offset: x-position to plot this set of energy levels.
    - delta: extra shift applied to horizontal line positions (for fine adjustment).
    - ylabel: optional label for the y-axis.
    """

    if ax is None:
        fig, ax = plt.subplots()

    # Initialize storage for custom ticks if not present
    if not hasattr(ax, "_custom_xticks"):
        ax._custom_xticks = []
        ax._custom_xticklabels = []

    # Group nearly degenerate levels
    unique_levels = []
    grouped_levels = []

    for ev in sorted(eigenvalues):
        if not unique_levels or abs(ev - unique_levels[-1]) > tolerance:
            unique_levels.append(ev)
            grouped_levels.append([ev])
        else:
            grouped_levels[-1].append(ev)

    x_offset = 0.15

    for level_group in grouped_levels:
        energy = level_group[0]
        n_deg = len(level_group)
        x_positions = np.linspace(-x_offset * (n_deg - 1) / 2,
                                  x_offset * (n_deg - 1) / 2, n_deg) + offset
        for x in x_positions:
            ax.hlines(y=energy, xmin=x - 0.05 + delta, xmax=x + 0.05 + delta,
                      color=color, linewidth=2)

    # Add custom label
    if label and offset not in ax._custom_xticks:
        ax._custom_xticks.append(offset)
        ax._custom_xticklabels.append(label)

    # Apply the full list of custom ticks and labels
    ax.set_xticks(ax._custom_xticks)
    ax.set_xticklabels(ax._custom_xticklabels)
    ax.tick_params(axis='x', which='both', bottom=False, top=False)

    ax.get_xaxis().set_visible(True)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.grid(axis='y', linestyle='--', alpha=0.5)

    return ax

def fig_tensor_rep_1(tensor, n_points=40):

    u = np.linspace(0, 2 * np.pi, n_points)
    v = np.linspace(0, np.pi, n_points)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))

    points = np.stack((x.ravel(), y.ravel(), z.ravel()), axis=1)  # Shape: (N, 3)
    points = points / np.linalg.norm(points, axis=1)[:, np.newaxis]

    N = points.shape[0]
    magnitudes = np.zeros(N)

    # Compute M(n) for each vector n
    for idx in range(N):
        n = points[idx]  # Extract the vector
        magnitudes[idx] = n.T @ tensor @ n  # Matrix-vector operations

    x_scaled = (x.ravel() * magnitudes).reshape(x.shape)
    y_scaled = (y.ravel() * magnitudes).reshape(y.shape)
    z_scaled = (z.ravel() * magnitudes).reshape(z.shape)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    norm = plt.Normalize(magnitudes.min(), magnitudes.max())
    colors = plt.cm.turbo(norm(magnitudes).reshape(x.shape))

    ax.plot_surface(x_scaled, y_scaled, z_scaled, facecolors=colors, edgecolor='k', alpha=1, linewidth=0.5)

    # Set transparent background
    fig.patch.set_alpha(0.0)
    ax.patch.set_alpha(0.0)

    # Hide the axes and labels
    # ax.set_axis_off()

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    x_limits = [x_scaled.min(), x_scaled.max()]
    y_limits = [y_scaled.min(), y_scaled.max()]
    z_limits = [z_scaled.min(), z_scaled.max()]

    all_limits = np.array([x_limits, y_limits, z_limits])
    widest_range = all_limits.max() - all_limits.min()

    ax.set_xlim(all_limits.min(), all_limits.min() + widest_range)
    ax.set_ylim(all_limits.min(), all_limits.min() + widest_range)
    ax.set_zlim(all_limits.min(), all_limits.min() + widest_range)

    mappable = plt.cm.ScalarMappable(cmap='turbo', norm=norm)
    mappable.set_array(magnitudes)
    plt.colorbar(mappable, ax=ax, shrink=0.5, aspect=8)

    plt.show()

def fig_tensor_rep_2(tensor, n_points=100):

    from mpl_toolkits.mplot3d import Axes3D

    eigenvalues, eigenvectors = np.linalg.eigh(tensor)

    u = np.linspace(0, 2 * np.pi, n_points)
    v = np.linspace(0, np.pi, n_points)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))

    sphere = np.stack((x.ravel(), y.ravel(), z.ravel()))

    ellipsoid = eigenvectors @ np.diag(np.sqrt(eigenvalues)) @ sphere
    x_ellipsoid, y_ellipsoid, z_ellipsoid = ellipsoid.reshape(3, *x.shape)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(x_ellipsoid, y_ellipsoid, z_ellipsoid, color='b', alpha=0.6, edgecolor='k')

    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')

    x_limits = [x_ellipsoid.min(), x_ellipsoid.max()]
    y_limits = [y_ellipsoid.min(), y_ellipsoid.max()]
    z_limits = [z_ellipsoid.min(), z_ellipsoid.max()]

    all_limits = np.array([x_limits, y_limits, z_limits])
    widest_range = all_limits.max() - all_limits.min()

    ax.set_xlim(all_limits.min(), all_limits.min() + widest_range)
    ax.set_ylim(all_limits.min(), all_limits.min() + widest_range)
    ax.set_zlim(all_limits.min(), all_limits.min() + widest_range)

    plt.show()

def fig_susc_field(conf, dic_Bkq, temp=2., n_points=20, delta=0.01):

    from mpl_toolkits.mplot3d import Axes3D

    def use_nja_(conf, dic_Bkq, field_vecs, wordy=False):

        calc = calculation(conf, ground_only=True, TAB=True, wordy=wordy)
        
        dic = {}
        dic['dic_bkq'] = dic_Bkq

        Magn = Magnetics(calc, ['Hcf','Hz'], dic)
        _, susc_field = Magn.susceptibility_field(fields=field_vecs, temp=temp, delta = delta, wordy=wordy)

        return susc_field

    u = np.linspace(0, 2 * np.pi, n_points)
    v = np.linspace(0, np.pi, n_points)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))

    points = np.stack((x.ravel(), y.ravel(), z.ravel()), axis=1)  # Shape: (N, 3)

    points = points / np.linalg.norm(points, axis=1)[:, np.newaxis]

    magnitudes = use_nja_(conf, dic_Bkq, points)

    #subtract the average
    #magnitudes -= np.average(magnitudes)

    x_scaled = (x.ravel() * magnitudes).reshape(x.shape)
    y_scaled = (y.ravel() * magnitudes).reshape(y.shape)
    z_scaled = (z.ravel() * magnitudes).reshape(z.shape)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x_scaled, y_scaled, z_scaled, cmap='viridis', edgecolor='k', alpha=0.8)

    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')

    x_limits = [x_scaled.min(), x_scaled.max()]
    y_limits = [y_scaled.min(), y_scaled.max()]
    z_limits = [z_scaled.min(), z_scaled.max()]

    all_limits = np.array([x_limits, y_limits, z_limits])
    widest_range = all_limits.max() - all_limits.min()

    ax.set_xlim(all_limits.min(), all_limits.min() + widest_range)
    ax.set_ylim(all_limits.min(), all_limits.min() + widest_range)
    ax.set_zlim(all_limits.min(), all_limits.min() + widest_range)

    plt.show()

    return magnitudes, points

def fig_rep_magnfield(Mvec, xyz, data=None):

    import matplotlib
    import matplotlib.cm as cm

    def cmap2list(cmap, N=10, start=0, end=1):
        x = np.linspace(start, end, N)
        colors = cmap(x)
        return colors

    ### plot magnetization surface
    fig = plt.figure()
    ax = fig.add_subplot(121, projection='3d', facecolor='white')
    Mvec = np.reshape(Mvec,(len(Mvec),1))
    data_p = np.hstack((Mvec,xyz))
    data_p = data_p[data_p[:, 0].argsort()]
    vectors = data_p[:,1:]
    norm_or = data_p[:,0]
    
    colorlist = cmap2list(cm.coolwarm, N=vectors.shape[0])

    ax.scatter(vectors[:,0],vectors[:,1],vectors[:,2], color=colorlist)

    box = plt.axes([0.75, 0.3, 0.02, 0.45])
    norm = matplotlib.colors.Normalize(norm_or[0],norm_or[-1])
    plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.coolwarm), cax=box)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    if data is not None:
        for i in range(data.shape[0]):
            vector = data[i,1:-1]
            ax.plot([0.,vector[0]],[0.,vector[1]],[0.,vector[2]],'--',lw=0.2,c='k')
            if data[i,0] in color_atoms().keys():
                ax.scatter(vector[0],vector[1],vector[2],'o',c = color_atoms()[data[i,0]],lw=3)
            else:
                ax.scatter(vector[0],vector[1],vector[2],'o',c = color_atoms()['_'],lw=3)
            ax.text(vector[0]+0.4*np.sign(vector[0]),vector[1]+0.4*np.sign(vector[1]),vector[2]+0.4*np.sign(vector[2]),data[i,-1], size=8)

    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.set_zlim(-4, 4)

    # Remove the grid
    ax.grid(False)

    plt.show()

def calc_segm(E_val, x=1, spanx=0.5):
    """ Costruisce il segmento (x-spanx/x, E_val), (x+spanx/2, E_val) """

    x1 = x - spanx/2
    x2 = x + spanx/2

    segment = (x1, x2), (E_val, E_val) #tupla di tuple
    return segment

def plot_segm(ax, segment, lw=0.5, c='k', ls='-'):
    """ Plotta il segmentino, ritorna l'oggetto 'linea' """
    line, = ax.plot(segment[0], segment[1],
            lw=lw,      # linewidth
            c=c,        # color
            ls=ls,      # linestyle
            )
    return line

def text_on_segm(ax, segment, text, text_kwargs={}):
    """ Scrive text sul lato sinistro del segment che gli passi """
    x_t = segment[0][0]
    y_t = segment[1][0]
    text += ' '     # Per mettere text non appiccicato alla barra
    text_plot = ax.text(x_t+0.11, y_t+80, text,
            horizontalalignment='right',
            verticalalignment='center',
            **text_kwargs)
    return text_plot

def plot_lines(ax, S1, S2, e1, e2, l, c='k'):

    for i,e in enumerate(e2):
        for ii, item in enumerate(l[int(i+1)].items()):
            key, value = item
            if str(e1) == key:
                line, = ax.plot([S1[0][1], S2[ii][0][0]], [e1, e], lw=value/100, c=c)
            else:
                line=None
    return line

def level_fig_tot(E_matrix, theories, proj_LS_dict, proj_prev_dict):

    COLORS = COLORS_list()
    levels = [str(int(w)) for w in range(E_matrix.shape[1])]  #number of levels
    deg = {theory:np.ones(len(set(E_matrix[k,:])), dtype='int32') for k,theory in enumerate(theories)}
    segm = {}
    spanx = 0.5
    for k, theory in enumerate(theories):
        x = k + 1   # La scala delle x deve partire da 1 altrimenti fa schifo
        segm[theory] = [calc_segm(E, x, spanx) for E in np.sort(list(set(E_matrix[k,:])))]    # Costruisco i segmentini e li salvo nel dizionario di prima
                                                                                              # solo uno per valore di energia
        count = 0
        for i in range(len(deg[theory])):
            prec = E_matrix[k,count]
            for j in range(int(count),len(E_matrix[k,:])):
                #print(E_matrix[k,j],prec)
                if E_matrix[k,j]==prec:
                    deg[theory][i] += 1
                    #print(i,deg[theory][i])
                else:
                    break
            deg[theory][i] -= 1
            count += deg[theory][i]

    fig = plt.figure()
    fig.set_size_inches(11,8)  #larghezza altezza
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.10, top=0.95)
    ax = fig.add_subplot()

    for k, theory in enumerate(theories):
        [plot_segm(ax, S,
            lw=2.0,
            c=C,
            ls='-') for S, C in zip(segm[theory], COLORS[:len(segm[theory])])]

        if k==0:
            keys_ee = []
            for kk in range(len(deg[theory])):
                [keys_ee.append(key) for key in proj_LS_dict[theory][sum(deg[theory][:kk+1])].keys()]
            [text_on_segm(ax, S,
                keys_ee[kk]+' ({})'.format(deg[theory][kk]),
                text_kwargs={'fontsize':12, 'color':COLORS[kk]})
                for kk, S in enumerate(segm[theory])]

        else:
            [text_on_segm(ax, S,
                '({})'.format(deg[theory][kk]),
                text_kwargs={'fontsize':12, 'color':COLORS[kk]})
                for kk, S in enumerate(segm[theory]) if deg[theory][kk] > 1]

        if k>0:
            [plot_lines(ax, S1 = segm[theories[k-1]][kk], S2 = segm[theory],
                e1 = np.sort(list(set(E_matrix[k-1,:])))[kk], e2 = E_matrix[k,:],
                l = proj_prev_dict[theory], c = C)
                for kk,C in enumerate(COLORS[:len(segm[theories[k-1]])])]

    ax.tick_params(labelsize=12)
    ax.ticklabel_format(axis='y', style='scientific', useMathText=True)#, scilimits=(0,0))
    ax.set_xticks(np.arange(len(theories))+1)   # Scala finta coi numeri (che parte da 1)
    ax.set_xticklabels(theories, fontsize=12)                # Rimpiazzo i numeri con la lista delle teorie
    ax.set_ylabel('Energy (cm$^{-1})$', fontsize=12)

    plt.show()

def E_vs_field(field, Epert, Etrue, name=None):

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    fig.set_size_inches(12,8)
    plt.subplots_adjust(left=0.1, bottom=0.1, top=0.9, right=0.875)

    ax.set_xlim(1,11)
    ax.set_xticks(np.arange(int(min(field)),int(max(field))+2,2))
    ax.set_xlabel('B$_0$ (T)')
    ax.set_ylabel('Energy (cm$^{-1}$)')
    ax.plot(0,0,'^', label='2° order', c='grey', transform=fig.transFigure)
    ax.plot(0,0,'.', label='exact', c='grey', transform=fig.transFigure)

    for i in range(Etrue.shape[1]):
        dots, = ax.plot(field, Etrue[:,i], '.')
        ax.plot(field, Etrue[:,i], '-', c=dots.get_color(), label=str(i+1))
        ax.plot(field, Epert[:,i], '--', c=dots.get_color())
        ax.plot(field, Epert[:,i], '^', c=dots.get_color())

    ax.legend(loc='upper left', bbox_to_anchor=(0.89,0.9), bbox_transform=fig.transFigure)
    if name is None:
        plt.show()
    else:
        plt.savefig(name, dpi=300)

def plot_charge_density(A2, A4, A6):
    #Reproduce the plots in fig 8 of Jeffrey D. Rinehart and Jeffrey R. Long Chem.Sci., 2011, 2, 2078-2085

    theta = np.linspace(0, np.pi, 50)
    phi = np.linspace(0, 2*np.pi, 50)
    theta, phi = np.meshgrid(theta, phi)
    r = np.zeros_like(theta)
    for i in range(theta.shape[0]):
        for j in range(theta.shape[1]):
            r[i,j] = Freeion_charge_dist(theta[i,j], phi[i,j], A2, A4, A6)
    xyz = np.array([r*np.sin(theta) * np.sin(phi),
                    r*np.sin(theta) * np.cos(phi),
                    r*np.cos(theta)])
    X, Y, Z = xyz

    fig = plt.figure(figsize=plt.figaspect(1.))
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    asse = 1
    ax.set_xlim(-asse, asse)
    ax.set_ylim(-asse, asse)
    ax.set_zlim(-asse, asse)
    ax.plot_surface(X, Y, Z, edgecolor='royalblue', lw=0.5, rstride=2, cstride=2,
                alpha=0.3)
    plt.show()

def plot_charge_density_data(A2, A4, A6, data):
    #Reproduce the plots in fig 8 of Jeffrey D. Rinehart and Jeffrey R. Long Chem.Sci., 2011, 2, 2078-2085

    theta = np.linspace(0, np.pi, 50)
    phi = np.linspace(0, 2*np.pi, 50)
    theta, phi = np.meshgrid(theta, phi)
    r = np.zeros_like(theta)
    for i in range(theta.shape[0]):
        for j in range(theta.shape[1]):
            r[i,j] = Freeion_charge_dist(theta[i,j], phi[i,j], A2, A4, A6)
    xyz = np.array([r*np.sin(theta) * np.sin(phi),
                    r*np.sin(theta) * np.cos(phi),
                    r*np.cos(theta)])
    X, Y, Z = xyz

    fig = plt.figure(figsize=plt.figaspect(1.))
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    asse = 2
    ax.set_xlim(-asse, asse)
    ax.set_ylim(-asse, asse)
    ax.set_zlim(-asse, asse)
    ax.plot_surface(X, Y, Z, edgecolor='royalblue', lw=0.5, rstride=2, cstride=2,
                alpha=0.3)
    for i in range(data.shape[0]):
        vector = data[i,1:-1]
        ax.quiver(0.,0.,0.,vector[0],vector[1],vector[2],color='b')
        ax.text(vector[0],vector[1],vector[2],data[i,-1])
    plt.show()
