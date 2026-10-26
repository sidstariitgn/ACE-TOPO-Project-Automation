# This file follows CFD-GEOM journaling style.
# Single-material HAR: an Oxide block with a circular trench 500 nm deep,
# and a fluid region above it.
#
# Same conformal method as the ON-stack script:
#   1. Oxide base block (below the trench bottom).
#   2. Split its top face ONCE with the trench circle (annulus + disk).
#   3. Extrude [annulus, disk] by the trench depth:
#        annulus -> Oxide (trench wall), disk -> HAR_Fluid (trench).
#   4. Extrude the new top faces again for the top fluid region.
#   5. UniteDomains: the 2 oxide pieces -> ONE Oxide domain, the 3 fluid
#      pieces -> ONE fluid domain. Internal faces used only by the united
#      pieces are discarded; faces shared with the other material are kept,
#      so the oxide/fluid interface stays a single shared surface.
#   6. VCs/BCs are assigned AFTER the unions (Booleans can drop them).
#
# Geometry only - no meshing.

import GUtils
import GGeometry
import GBCVC

# =====================================================================
# USER INPUTS (nm)
# =====================================================================
TRENCH_DEPTH_NM = 600.0        # depth of the trench
TRENCH_RADIUS_NM = 60.0        # radius of the trench
PROFILE_WIDTH_NM = 400.0       # width of the whole (square) profile
BASE_THICKNESS_NM = 100.0      # oxide below the trench bottom
TOP_REGION_THICKNESS_NM = 500.0  # flux plane height above the oxide top;
                                 # 200 nm keeps the TOPO domain short
                                 # (finer cells at the same refinement level)

MATERIAL_VC = "Oxide"
FLUID_VC = "Fluid"
MERGE_REGIONS = True             # False -> keep the separate pieces
# =====================================================================

NM = 1.0 / 1000.0              # model units are microns
TRENCH_DEPTH = TRENCH_DEPTH_NM * NM
TRENCH_RADIUS = TRENCH_RADIUS_NM * NM
WIDTH = PROFILE_WIDTH_NM * NM
BASE_THICKNESS = BASE_THICKNESS_NM * NM
TOP_REGION_THICKNESS = TOP_REGION_THICKNESS_NM * NM

X_MIN, X_MAX = 0.0, WIDTH
Y_MIN, Y_MAX = 0.0, WIDTH
TRENCH_CENTER_X = WIDTH / 2.0
TRENCH_CENTER_Y = WIDTH / 2.0

if not (0.0 < TRENCH_RADIUS < WIDTH / 2.0):
    raise ValueError("Trench radius must be > 0 and < width/2")
if TRENCH_DEPTH <= 0.0 or BASE_THICKNESS <= 0.0:
    raise ValueError("Trench depth and base thickness must be > 0")


def set_vc(domain, name, region_type="Solid"):
    GBCVC.SetVC(domain, name, region_type)


def create_rectangle_surface():
    p1 = GGeometry.CreatePoint(X_MIN, Y_MIN, 0.0)
    p2 = GGeometry.CreatePoint(X_MIN, Y_MAX, 0.0)
    p3 = GGeometry.CreatePoint(X_MAX, Y_MAX, 0.0)
    p4 = GGeometry.CreatePoint(X_MAX, Y_MIN, 0.0)

    curves = GGeometry.CreateLines([p1, p4, p3, p2, p1], 0)
    return GGeometry.FitSurfaceThroughCurves(
        [
            curves["lines"][1],
            curves["lines"][0],
            curves["lines"][2],
            curves["lines"][3],
        ]
    )


def create_circle_at_z(z_value):
    center = GGeometry.CreatePoint(TRENCH_CENTER_X, TRENCH_CENTER_Y, z_value)
    radius_point = GGeometry.CreatePoint(
        TRENCH_CENTER_X - TRENCH_RADIUS,
        TRENCH_CENTER_Y,
        z_value,
    )
    return GGeometry.Create2PtCircle(center, radius_point)


def get_next_top_surfaces(split_solid):
    # CreateExtrudedSolids([annulus, disk], ...) -> next top is [surfs][5], [surfs][6]
    return [split_solid["surfs"][5], split_solid["surfs"][6]]


GUtils.SetUnits(5)   # microns

# -----------------------------
# 1. Oxide base block
# -----------------------------
base_surface = create_rectangle_surface()
base = GGeometry.CreateExtrudedSolids(base_surface, 0, 0, 1, BASE_THICKNESS, 1, 1)

# -----------------------------
# 2. Split base top surface (once)
# -----------------------------
circle = create_circle_at_z(BASE_THICKNESS)
split_top = GGeometry.Make2Surfaces(base["surfs"][4], circle)
current_top_surfaces = [base["surfs"][4], split_top["surfs"][0]]

# -----------------------------
# 3. Trench layer (500 nm): Oxide wall + fluid trench
# -----------------------------
trench = GGeometry.CreateExtrudedSolids(current_top_surfaces, 0, 0, 1, TRENCH_DEPTH, 1, 1)
current_top_surfaces = get_next_top_surfaces(trench)

# -----------------------------
# 4. Top fluid region
# -----------------------------
top_region = GGeometry.CreateExtrudedSolids(current_top_surfaces, 0, 0, 1, TOP_REGION_THICKNESS, 1, 1)
top_faces = get_next_top_surfaces(top_region)

oxide_pieces = [base["udomains"][0], trench["udomains"][0]]
fluid_pieces = [trench["udomains"][1], top_region["udomains"][0], top_region["udomains"][1]]

# -----------------------------
# 5. Merge into one Oxide domain and one Fluid domain
# -----------------------------
if MERGE_REGIONS:
    oxide_domains = GGeometry.UniteDomains(oxide_pieces[0], oxide_pieces[1:])["udomains"]
    fluid_domains = GGeometry.UniteDomains(fluid_pieces[0], fluid_pieces[1:])["udomains"]
    print("After merging: %d oxide domain(s), %d fluid domain(s) (expected 1 and 1)"
          % (len(oxide_domains), len(fluid_domains)))
else:
    oxide_domains = oxide_pieces
    fluid_domains = fluid_pieces

# -----------------------------
# 6. VCs and BCs (after the Booleans)
# -----------------------------
for dom in oxide_domains:
    set_vc(dom, MATERIAL_VC, "Solid")
for dom in fluid_domains:
    set_vc(dom, FLUID_VC, "Fluid")

try:
    GBCVC.SetBC(top_faces[0], "NoName", "Inlet")
    GBCVC.SetBC(top_faces[1], "NoName", "Inlet")
except Exception as exc:
    # Only a label - the Flux BC is set in TOPO anyway.
    print("Could not label the top faces as Inlet (%s); set it in the GUI/TOPO." % exc)