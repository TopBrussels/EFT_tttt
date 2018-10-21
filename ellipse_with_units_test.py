import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches


fig = plt.figure()
x = np.linspace(-1.0, 1.0, 100)
y = x

x_marg = 0.45
x_marginal = np.array([x_marg]*100)
y_marginal = np.linspace(-1.0, 1.0, 100)

x_indep = 0.29
x_indep = np.array([x_indep]*100)
y_indep = np.linspace(-1.0, 1.0, 100)

ax = fig.add_subplot(1, 1, 1)
ax.set_title('Illustration of marginal and independent limits')
ax.set_xlabel(r'$c_{R}$')
ax.set_ylabel(r'$c_{L}$')
ax.xaxis.set_label_coords(0.95, 0.55)
ax.yaxis.set_label_coords(0.55, 0.95)
ax.plot(x, y, alpha=0.0)
marg_line, = ax.plot(x_marginal,y_marginal,color='red')
ax.plot(-x_marginal,y_marginal,color='red',label=r'Marginal')
indep_line, = ax.plot(x_indep,y_indep,color='green', label='Independent\n'+r'$C_{L}=0$')
ax.plot(-x_indep,y_indep,color='green')

plt.legend()

ax.spines['left'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('center')
ax.spines['top'].set_color('none')
ax.spines['left'].set_smart_bounds(True)
ax.spines['bottom'].set_smart_bounds(True)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

xcenter, ycenter = 0., 0.
width, height = 1.0, 0.35
angle = -30
theta = np.deg2rad(np.arange(0.0, 360.0, 1.0))
x = 0.5 * width * np.cos(theta)
y = 0.5 * height * np.sin(theta)
rtheta = np.radians(angle)
R = np.array([
    [np.cos(rtheta), -np.sin(rtheta)],
    [np.sin(rtheta),  np.cos(rtheta)],
    ])
x, y = np.dot(R, np.array([x, y]))
x += xcenter
y += ycenter
e1 = patches.Ellipse((xcenter, ycenter), width, height,
                     angle=angle, linewidth=2, fill=False, zorder=2)
ax.add_patch(e1)

# ----------------------------------------------------


def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])

plt.savefig('marginal_and_independent_limits.png')
plt.show()

