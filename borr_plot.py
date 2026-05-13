from mpl_toolkits import mplot3d
import numpy as np
from matplotlib import pyplot as plt
from concurrent.futures import ProcessPoolExecutor
import time

from JONES import *

def endpoints(C):
    return C[-1], C[0]

def intermediate(C, t):
    r1 = np.array(endpoints(C)[0])
    r2 = np.array(endpoints(C)[1])
    return list(r1 + t * (r2 - r1))

def build_segments(t):
    u1 = [[0,0,0],[1,1,0],[2,2,0.5],[3,3,0.5],[4,4,0],[5,5,0],
          [6,6,0.5],[7,7,0.5],[8,7,0.5],[9,5,0.2],[9,3,0.2],
          [8,0,0.2],[8,-1,0.2],[6,-1.5,0],[4,-2,0],[2,-1.5,0]]

    u2 = [[1,0,0.5],[4,0,0],[5,1,0],[5,4,0.5],[4,5,0.5],
          [3,6,0],[2,7,0],[-1,6,0],[-1,3,0.5]]

    u3 = [[6,0,0.5],[7,6,0],[6,7,0],[3,7,0.5],
          [2,6,0.5],[2,3,0],[3,2,0],[4,1,0.5]]

    if t > 0:
        for l in (u1, u2, u3):
            l.append(intermediate(l, t))

    return [u1, u2, u3]

# -----------------------
# plotting 
# -----------------------

Col = ['r', 'b', 'k']

def plot_curve_set(BM2, t):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.set_box_aspect([1, 1, 1])
    ax.set_xlim((-4, 12))
    ax.set_ylim((-5, 15))
    ax.set_zlim((-2, 2))
    ax.view_init(elev=82, azim=-90)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    for i, C in enumerate(BM2):
        ex = np.array([p[0] for p in C])
        ey = np.array([p[1] for p in C])
        ez = np.array([p[2] for p in C])
        ax.plot(ex, ey, ez, color=Col[i])

    fig.savefig(f"borr_{round(t,2)}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

# -----------------------
# worker
# -----------------------

def compute_for_t(t):

    BM2 = build_segments(t)

    # ---- PLOT curve
    #plot_curve_set(BM2, t)

    # ---- JONES ----
    input_knot = [list(arr) for arr in BM2]

    if t == 1.0:
        return None
    # ----projections
    n = 1000
    start = time.time()
    result = expected_jones(input_knot, n)
    end = time.time()

    print(f"t={t:.2f} done in {end-start:.2f}s")

    return result

# -----------------------
# main
# -----------------------

if __name__ == "__main__":

    K = [0, 0.22, 0.44, 0.67, 0.68, 0.70, 0.89, 1]
    print("K:", K)

    with ProcessPoolExecutor() as executor:
        results = list(executor.map(compute_for_t, K))

    results = [r for r in results if r is not None]

    with open("borr_lists.txt", "w") as f:
        print(results, file=f)
        print("", file=f)
        print(len(results), file=f)

    print("DONE. Total cases:", len(results))