import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
import pandas as pd
import trimesh

def analytical_Y(ell, m, theta, phi):
    if ell == 0 and m == 0:
        return np.sqrt(1 / (4 * np.pi)) * np.ones_like(theta, dtype=complex)
    elif ell == 1:
        if m == -1:
            return np.sqrt(3 / (8 * np.pi)) * np.sin(theta) * np.exp(-1j * phi)
        elif m == 0:
            return np.sqrt(3 / (4 * np.pi)) * np.cos(theta)
        elif m == 1:
            return -np.sqrt(3 / (8 * np.pi)) * np.sin(theta) * np.exp(1j * phi)

def Psi_vec(ell, m, k, points, x0):
    diff = points - x0
    r = np.linalg.norm(diff, axis=1)
    r[r==0] = 1e-16  # avoid singularity
    theta = np.arccos(diff[:,1]/r)
    phi = np.arctan2(diff[:,2], diff[:,0])
    h2 = sp.spherical_jn(ell, k*r) - 1j*sp.spherical_yn(ell, k*r)
    # Y = sp.sph_harm_y(m, ell, phi, theta)
    Y = analytical_Y(ell, m, theta, phi)
    return h2 * Y

def pressure_vec(k, points, ell_m, sources, coeffs):
    p = np.zeros(points.shape[0], dtype=complex)
    for i, source in enumerate(sources):
        ell, m = ell_m[i]
        p += coeffs[i] * Psi_vec(ell, m, k, points, source)
    return p

import os

def plot_heatmap(filename, grid_resolution, plane="xy", z_fixed=0, x_fixed=0, y_fixed=0, outdir="plots"):
    # --- 1. Read metadata ---
    with open(filename, "r") as f:
        first_line = f.readline().strip().split(",")  # CSV
    meta = {}
    i = 0
    while i < len(first_line)-1:
        meta[first_line[i]] = first_line[i+1]
        i += 2

    geo_name = meta.get('geo_name', 'unknown')
    mode = int(meta['mode'])
    k = float(meta['k'])
    max_r = float(meta['max_r'])
    print("geo_name:", geo_name, "mode:", mode, "k:", k, "max_r:", max_r)

    # --- 2. Read table ---
    df = pd.read_csv(filename, skiprows=1)
    df["coeff"] = df["re"] + 1j * df["im"]
    sources = df[["source_x","source_y","source_z"]].to_numpy()
    ell_m = df[["n","m"]].to_numpy()
    coeffs = df["coeff"].to_numpy()

    # --- 3. Define grid depending on plane ---
    grid_min, grid_max = -5*max_r, 5*max_r
    num_points = grid_resolution
    vals = np.linspace(grid_min, grid_max, num_points)

    if plane == "xy":
        X, Y = np.meshgrid(vals, vals)
        points = np.stack([X.ravel(), Y.ravel(), np.full(X.size, z_fixed)], axis=1)
        xlabel, ylabel = "x", "y"
        fixed_info = f"z={z_fixed}"
    elif plane == "yz":
        Y, Z = np.meshgrid(vals, vals)
        points = np.stack([np.full(Y.size, x_fixed), Y.ravel(), Z.ravel()], axis=1)
        X, Y = Y, Z  # for plotting
        xlabel, ylabel = "y", "z"
        fixed_info = f"x={x_fixed}"
    elif plane == "xz":
        X, Z = np.meshgrid(vals, vals)
        points = np.stack([X.ravel(), np.full(X.size, y_fixed), Z.ravel()], axis=1)
        X, Y = X, Z  # for plotting
        xlabel, ylabel = "x", "z"
        fixed_info = f"y={y_fixed}"
    else:
        raise ValueError("plane must be 'xy', 'yz', or 'xz'")

    # --- 4. Compute pressure magnitude ---
    P_flat = pressure_vec(k, points, ell_m, sources, coeffs)
    P_flat_abs = np.abs(P_flat)
    P = P_flat_abs.reshape(X.shape)

    # --- 5. Load mesh ---
    mesh = trimesh.load(f"../../asset/{geo_name}/{geo_name}.obj")
    inside_mask = mesh.contains(points)
    outside_mask = ~inside_mask
    vmax_outside = P_flat_abs[outside_mask].max()

    # --- 6. Plot heatmap ---
    plt.figure(figsize=(6,5))
    plt.pcolormesh(X, Y, P, shading='auto', cmap='viridis', vmax=vmax_outside)
    plt.colorbar(label='|p|')

    # Overlay mesh edges
    vertices = mesh.vertices
    if plane == "xy":
        vertices_proj = vertices[:, :2]
    elif plane == "yz":
        vertices_proj = vertices[:, 1:3]
    else:  # xz
        vertices_proj = vertices[:, [0,2]]

    for face in mesh.faces:
        pts = vertices_proj[face]
        pts = np.vstack([pts, pts[0]])
        plt.plot(pts[:,0], pts[:,1], color='black', linewidth=1)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f'acoustic transfer in {plane.upper()} plane ({fixed_info}) \n {geo_name} mode {mode}')

    # --- Save instead of show ---
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{geo_name}_mode{mode}.png")
    plt.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close()  # free memory
    print(f"Saved: {outfile}")



def worker(geo, i):
    outdir = os.path.join("../../acoustic_transfer_plots/", geo)
    return plot_heatmap(f"../../data/{geo}_mode{i}.csv", 800, plane="xy", outdir=outdir)

