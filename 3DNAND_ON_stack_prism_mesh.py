# This file follows CFD-GEOM journaling style.
# 3D NAND Oxide/Nitride (ON) stack with a circular HAR trench,
# built directly as an EXTRUDED PRISM MESH (semi-structured domains).
#
# Method
#   1. Make the square footprint at z = 0, split it once with the trench
#      circle (annulus + disk), and mesh those two surfaces with triangles.
#   2. Sweep that 2-D mesh upward layer by layer with
#      GMesh.CreateExtrudedVolumeGrids. Each sweep starts from the top
#      surfaces of the previous one, so every interface is one shared,
#      perfectly flat surface with identical nodes on both sides.
#   3. Cells in z are set per layer (default: 1 cell per oxide/nitride
#      layer, which already captures a flat layer exactly).
#
# Cell count = (triangles in footprint) x (cells in z). The footprint is fine
# only along the trench circle and coarsens toward the box edges.

import GUtils
import GGeometry
import GMesh
import GBCVC
import GFileIO

# =====================================================================
# USER INPUTS (lengths in nm)
# =====================================================================
OXIDE_THICKNESS_NM = 20.0
NITRIDE_THICKNESS_NM = 30.0
NUM_ON_PAIRS = 10
TRENCH_RADIUS_NM = 60.0
PROFILE_WIDTH_NM = 400.0

SUBSTRATE_THICKNESS_NM = 100.0
TOP_REGION_THICKNESS_NM = 500.0

# ---- mesh controls ----
CIRCLE_POINTS = 64            # points on the trench circle (roundness of hole)
EDGE_POINTS = 21              # points on each outer box edge
FOOTPRINT_MAX_CELL_NM = 20.0  # largest triangle size (box edges / far field)
FOOTPRINT_MIN_CELL_NM = 5.0   # smallest triangle size (near the circle)
TRANSITION_FACTOR = 1.2       # triangle growth rate away from the circle

CELLS_PER_OXIDE = 1           # cells through each oxide layer
CELLS_PER_NITRIDE = 1         # cells through each nitride layer
CELLS_SUBSTRATE = 3
CELLS_TOP_REGION = 10

EXPORT_DTF = False
EXPORT_PATH = "C:/Tests/Initial_Geo/HAR/HAR_3DNAND_ON_prism.DTF"
# =====================================================================

NM = 1.0 / 1000.0             # model units are microns
OXIDE_THICKNESS = OXIDE_THICKNESS_NM * NM
NITRIDE_THICKNESS = NITRIDE_THICKNESS_NM * NM
TRENCH_RADIUS = TRENCH_RADIUS_NM * NM
WIDTH = PROFILE_WIDTH_NM * NM
SUBSTRATE_THICKNESS = SUBSTRATE_THICKNESS_NM * NM
TOP_REGION_THICKNESS = TOP_REGION_THICKNESS_NM * NM

X_MIN, X_MAX = 0.0, WIDTH
Y_MIN, Y_MAX = 0.0, WIDTH
TRENCH_CENTER_X = WIDTH / 2.0
TRENCH_CENTER_Y = WIDTH / 2.0

if NUM_ON_PAIRS < 1:
    raise ValueError("NUM_ON_PAIRS must be >= 1")
if not (0.0 < TRENCH_RADIUS < WIDTH / 2.0):
    raise ValueError("Trench radius must be > 0 and < width/2")

GUtils.SetUnits(5)   # microns


# ---------------------------------------------------------------------
# Call CFD-GEOM API functions by parameter NAME, not position.
# (The argument order in the 2026.0 user guide does not always match the
#  installed GMesh.py, so we read the real signature at runtime.)
# ---------------------------------------------------------------------
import inspect


def call_by_name(func, first_arg, patterns, positional_fallback):
    """patterns: ordered list of (substring, value); the first substring
    found in a parameter name supplies that parameter's value."""
    try:
        spec = inspect.getargspec(func)
        names = list(spec.args)
        n_default = len(spec.defaults or ())
    except Exception:
        names, n_default = [], 0

    if not names:   # wrapper is f(*args): use the known internal order
        print("%s: no named args, calling positionally" % func.__name__)
        return func(first_arg, *positional_fallback)

    required = names[:len(names) - n_default]
    kwargs = {names[0]: first_arg}
    for name in names[1:]:
        low = name.lower()
        for sub, val in patterns:
            if sub in low:
                kwargs[name] = val
                break
        else:
            if name in required:
                raise RuntimeError("Don't know a value for '%s' in %s(%s)"
                                   % (name, func.__name__, ", ".join(names)))
    print("%s(%s)" % (func.__name__,
                      ", ".join("%s=%r" % (n, kwargs[n]) for n in names[1:] if n in kwargs)))
    return func(**kwargs)


# ---------------------------------------------------------------------
# Footprint: square split by the trench circle, at z = 0
# ---------------------------------------------------------------------
p1 = GGeometry.CreatePoint(X_MIN, Y_MIN, 0.0)
p2 = GGeometry.CreatePoint(X_MIN, Y_MAX, 0.0)
p3 = GGeometry.CreatePoint(X_MAX, Y_MAX, 0.0)
p4 = GGeometry.CreatePoint(X_MAX, Y_MIN, 0.0)
box_lines = GGeometry.CreateLines([p1, p4, p3, p2, p1], 0)
base_surface = GGeometry.FitSurfaceThroughCurves(
    [box_lines["lines"][1], box_lines["lines"][0],
     box_lines["lines"][2], box_lines["lines"][3]]
)

c_center = GGeometry.CreatePoint(TRENCH_CENTER_X, TRENCH_CENTER_Y, 0.0)
c_radius = GGeometry.CreatePoint(TRENCH_CENTER_X - TRENCH_RADIUS, TRENCH_CENTER_Y, 0.0)
trench_circle = GGeometry.Create2PtCircle(c_center, c_radius)

split = GGeometry.Make2Surfaces(base_surface, trench_circle)
annulus_surf = base_surface            # square with hole (same handle, as in your scripts)
disk_surf = split["surfs"][0]          # trench disk (same index as your working script)

# ---------------------------------------------------------------------
# 2-D triangle mesh on the footprint
# ---------------------------------------------------------------------
# Uniform point distributions on the boundary curves (power = 1.0)
GMesh.DirectMeshLinesCurvesPowerLaw(box_lines["lines"], EDGE_POINTS, 1.0, 1)
GMesh.DirectMeshLinesCurvesPowerLaw([trench_circle], CIRCLE_POINTS, 1.0, 1)

GMesh.SetGlobalSurfaceMeshType(1)      # 1 = triangular
MAX_CELL = FOOTPRINT_MAX_CELL_NM * NM
MIN_CELL = FOOTPRINT_MIN_CELL_NM * NM
call_by_name(
    GMesh.MeshSurfaces,
    [annulus_surf, disk_surf],
    [   # most specific substrings first
        ("curve_mesh_max_ratio", 1.0),
        ("force_remesh_surf", 1),
        ("force_remesh_curv", 0),
        ("remesh_surf", 1),
        ("remesh_curv", 0),        # keep the curve meshes made above
        ("override", 1),
        ("max_cell", MAX_CELL),
        ("min_cell", MIN_CELL),
        ("background", 0),         # use_background_mesh: off
        ("octree", 0),
        ("transition", TRANSITION_FACTOR),
        ("tfactor", TRANSITION_FACTOR),
        ("expansion", 1.1),
        ("efactor", 1.1),
        ("adapt", 0),              # no curvature adaption (int)
        ("angle", 30.0),
        ("type", 1),               # 1 = triangular
    ],
    # installed 2026.0 order: mesh_type, use_background_mesh, max_cell_size,
    # min_cell_size, transition_factor, adapt_to_curvature, max_angle,
    # remesh_curves, force_remesh_curves, force_remesh_surfaces
    [1, 0, MAX_CELL, MIN_CELL, TRANSITION_FACTOR, 0, 30.0, 0, 0, 1],
)


# ---------------------------------------------------------------------
# Layer-by-layer sweep of the 2-D mesh
# ---------------------------------------------------------------------
def sweep_layer(bottom_surfs, thickness, n_cells, annulus_vc, annulus_type,
                disk_vc, disk_type):
    """Extrude [annulus, disk] meshes by `thickness` with n_cells in z.

    Returns (annulus_domain, disk_domain, [top_annulus, top_disk]).
    """
    res = call_by_name(
        GMesh.CreateExtrudedVolumeGrids,
        bottom_surfs,
        [
            ("dx", 0.0), ("dy", 0.0), ("dz", 1.0),
            ("distance", thickness),
            ("num_grid", int(n_cells) + 1),   # grid points in z
            ("numtimes", 1),
            ("multiple", 1),
            ("link_extru", 0),
            ("link_trans", 0),
        ],
        [0.0, 0.0, 1.0, thickness, int(n_cells) + 1, 1, 1, 0, 0],
    )
    doms = res["ssdomains"]
    tops = res["surfs"]
    if len(tops) != 2:
        print("Warning: sweep returned %d surfs (expected 2); "
              "using the last two as top surfaces." % len(tops))
        tops = tops[-2:]
    GBCVC.SetVC(doms[0], annulus_vc, annulus_type)
    GBCVC.SetVC(doms[1], disk_vc, disk_type)
    return doms[0], doms[1], [tops[0], tops[1]]


all_domains = []

# 1. Silicon substrate (annulus and disk are both silicon)
a, d, top = sweep_layer([annulus_surf, disk_surf], SUBSTRATE_THICKNESS,
                        CELLS_SUBSTRATE, "Silicon_Sub", "Solid",
                        "Silicon_Sub_Core", "Solid")
all_domains += [a, d]

# 2. ON pairs with the trench (fluid) in the disk
for i in range(1, NUM_ON_PAIRS + 1):
    a, d, top = sweep_layer(top, OXIDE_THICKNESS, CELLS_PER_OXIDE,
                            "Oxide_%02d" % i, "Solid",
                            "HAR_Fluid_Ox_%02d" % i, "Fluid")
    all_domains += [a, d]

    a, d, top = sweep_layer(top, NITRIDE_THICKNESS, CELLS_PER_NITRIDE,
                            "Nitride_%02d" % i, "Solid",
                            "HAR_Fluid_Nit_%02d" % i, "Fluid")
    all_domains += [a, d]

# 3. Top fluid region
a, d, top = sweep_layer(top, TOP_REGION_THICKNESS, CELLS_TOP_REGION,
                        "Top_Fluid", "Fluid", "HAR_Fluid_Top", "Fluid")
all_domains += [a, d]

GBCVC.SetBC(top[0], "NoName", "Inlet")
GBCVC.SetBC(top[1], "NoName", "Inlet")

try:
    print("Domains: %d, total cells: %s" % (len(all_domains), GMesh.GetNumCells()))
except Exception:
    print("Domains: %d" % len(all_domains))

if EXPORT_DTF:
    GFileIO.ExportDTF(EXPORT_PATH, 3, 1, 5, 0)