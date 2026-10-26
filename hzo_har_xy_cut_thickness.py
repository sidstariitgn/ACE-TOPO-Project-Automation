# hzo_har_xy_cut_thickness.py   (v2: with HZO node index)
# Run inside CFD-VIEW:  File > Run Python Script
#
# HAR (hole) version of hzo_halfcycle_thickness.py, with materials HZO on Oxide
# (Oxide replaces SUB of the planar-wafer case).
#
# For every cycle n and half cycle (Hf, Zr) the script
#   1. switches the reader to the LAST time-step DTF of that half cycle
#   2. cuts the HZO surface with an X-Slice at X = CUT_X  -> profile in the Y-Z plane
#      and with a Y-Slice at Y = CUT_Y                  -> profile in the X-Z plane
#      (both planes pass through the hole axis, i.e. along the diameter)
#   3. evaluates EXPRESSION ('"Oxide"-"HZO"') on each cut curve = HZO thickness per point
#   4. finds for every cut point the nearest node of the HZO surface (hzo_node, node_dist)
#   5. writes the results to  OUT_DIR/x-cut/  and  OUT_DIR/y-cut/
#
# HOW TO USE
#   1. File > New, File > Import START_FILE with your usual settings.
#      Check that the Oxide and HZO objects are visible.
#   2. Check the USER SETTINGS below (folders, CUT_X / CUT_Y, HOLE_RADIUS, HOLE_DEPTH).
#   3. Run the script. Afterwards run har_thickness_analyzer.py (normal Python) on OUT_DIR.
#
# OUTPUT  (same layout in x-cut and y-cut)
#   <n>-<Hf|Zr>.txt   one line per cut point:
#       X Y Z        coordinates of the point (m)
#       s            position along the cut (Y for x-cut, X for y-cut)
#       offset       s - hole axis (negative = left wall / left side)
#       depth        Z_TOP - Z (positive inside the hole)
#       side         left / right
#       region       field (top) / sidewall / bottom / corner
#       thickness    calculator value (EXPRESSION)
#       thk_geo      distance to the original hole outline (cross-check)
#       hzo_node     0-based index of the nearest node of the HZO surface object
#                    (exact node if method = fallback; nearest node if a Slice was used,
#                     because slice points lie between HZO nodes)
#       node_dist    distance from the cut point to that HZO node (m)
#   summary.csv       one line per half cycle: field, sidewall at depths (left/right), bottom

import os
import re
import glob
import math

# ----------------------------- USER SETTINGS -----------------------------
DTF_DIR     = r"J:\Saketh_TOPO\Calibration_w_o_O-recom\HAR"        # <-- HAR DTF folder
OUT_DIR     = r"J:\Saketh_TOPO\Calibration_w_o_O-recom\HAR_out"    # x-cut / y-cut go in here
N_CYCLES    = 100
HALF_CYCLES = ["Hf", "Zr"]
FILE_PATTERN = "Cycle-{n}-{hc}-cycle.*.DTF"
START_FILE  = "Cycle-1-Hf-cycle.0000175.DTF"   # file you imported by hand

FILM_OBJECT = "HZO"          # film surface object, as named in the Object Explorer
BASE_OBJECT = "Oxide"        # underlying material object (was SUB for the planar wafer)

# Calculator expression = thickness (same as typed by hand in the calculator)
EXPRESSION  = '"Oxide"-"HZO"'

# Cut positions. 2E-07 lies on the hole diameter, so both cuts pass through the axis.
CUT_X       = 2e-07          # X-Slice plane  X = CUT_X   -> saved in x-cut
CUT_Y       = 2e-07          # Y-Slice plane  Y = CUT_Y   -> saved in y-cut
CUT_TOL     = 2e-09          # fallback only (if slicing fails): HZO nodes within +/- CUT_TOL
NODE_SLAB   = 1e-08          # nearest-node search only looks at HZO nodes within this distance
                             # of the cut plane (speed); widened automatically if none found

# Hole geometry (model units, metres) - REQUIRED
HOLE_CENTER = (CUT_X, CUT_Y) # (X, Y) of the hole axis
HOLE_RADIUS = None           # e.g. 5e-08  radius of the ORIGINAL (uncoated) hole
HOLE_DEPTH  = None           # e.g. 1e-06  depth of the ORIGINAL hole
Z_TOP       = None           # Z of the original top surface; None = max Z of BASE_OBJECT

SIDEWALL_DEPTH_FRACTIONS = [0.1, 0.25, 0.5, 0.75, 0.9]   # depths as fraction of HOLE_DEPTH
DEPTH_WINDOW_FRACTION    = 0.02     # +/- window (fraction of HOLE_DEPTH) around each depth
CORNER_MARGIN_FRACTION   = 0.05     # of HOLE_RADIUS; points this close to a corner = "corner"
BOTTOM_R_FRACTION        = 0.5      # bottom average uses points with r < this * HOLE_RADIUS
# Thickness is written in metres; the analyzer converts units.
# --------------------------------------------------------------------------

if HOLE_RADIUS is None or HOLE_DEPTH is None:
    raise Exception("Set HOLE_RADIUS and HOLE_DEPTH (metres) in the USER SETTINGS first.")


# --------------------------- file handling --------------------------------
def last_timestep_file(n, hc):
    pattern = os.path.join(DTF_DIR, FILE_PATTERN.format(n=n, hc=hc))
    best, best_num = None, -1
    for p in glob.glob(pattern):
        m = re.search(r"-cycle\.(\d+)\.dtf$", os.path.basename(p), re.IGNORECASE)
        if m and int(m.group(1)) > best_num:
            best, best_num = p, int(m.group(1))
    return best


def load_file(reader, path):
    try:
        reader.switchFile(path, reader.getSwitches(), reader.getSimulation(), True)
        Reflow()
    except Exception as e:
        print("       switchFile failed (%s) -> trying SetFilenames/Reread" % str(e))
        SetFilenames(reader, [path])
        Reread()


def find_reader():
    names = [START_FILE]
    for p in sorted(glob.glob(os.path.join(DTF_DIR, "*.DTF"))):
        b = os.path.basename(p)
        if b not in names:
            names.append(b)
    for b in names:
        for full in (False, True):
            try:
                arg = os.path.join(DTF_DIR, b) if full else b
                r = GetFileReader(arg, full)
            except Exception:
                r = None
            if r is not None:
                print("       found file reader via '%s' (fullpath=%s)" % (arg, str(full)))
                return r
    return None


def get_object(name):
    objs = GetObjectsByName(name)
    return objs[0] if objs else None


# ------------------------------ cutting -----------------------------------
def make_slice(surface, axis, value):
    """X-Slice / Y-Slice of a surface object -> curve object (or None)."""
    try:
        try:
            fn = CreateXSlice if axis == "X" else CreateYSlice
        except NameError:
            import VCurve
            fn = VCurve.CreateXSlice if axis == "X" else VCurve.CreateYSlice
        sl = fn([surface])
        if isinstance(sl, (list, tuple)):
            sl = sl[0]
        SetValue(sl, value)
        try:
            ForceUpdate()
        except Exception:
            pass
        return sl
    except Exception as e:
        print("       %s-Slice failed: %s" % (axis, str(e)))
        return None


def thickness_on(obj):
    """Evaluate EXPRESSION on obj; fallback: Oxide - HZO from nodal data."""
    for attempt in (0, 1):
        try:
            if attempt == 0:
                SetExpression(obj, EXPRESSION)
            else:
                ForceUpdate()
            r = GetExpressionResult(obj)
            if hasattr(r, "__len__"):
                return list(r)
        except Exception:
            pass
    try:
        a = GetNodalData(obj, BASE_OBJECT)
        b = GetNodalData(obj, FILM_OBJECT)
        return [a[i] - b[i] for i in range(len(a))]
    except Exception as e:
        print("       could not evaluate thickness: " + str(e))
        return None


def nearest_nodes(pts, FX, FY, FZ, axis, value):
    """For every (x, y, z) in pts: (index, distance) of the nearest HZO surface node."""
    c = FX if axis == "X" else FY
    slab = NODE_SLAB
    cand = []
    while not cand and slab < 1e3:
        cand = [i for i in range(len(c)) if abs(c[i] - value) <= slab]
        slab *= 4.0
    out = []
    for (x, y, z) in pts:
        best, bd = -1, 1e300
        for i in cand:
            d = (FX[i] - x) ** 2 + (FY[i] - y) ** 2 + (FZ[i] - z) ** 2
            if d < bd:
                best, bd = i, d
        out.append((best, math.sqrt(bd)))
    return out


def cut_points(film, FX, FY, FZ, axis, value):
    """List of (x, y, z, thickness, hzo_node, node_dist) on the cut, plus method text."""
    sl = make_slice(film, axis, value)
    if sl is not None:
        pts = None
        try:
            X = GetNodalData(sl, "X"); Y = GetNodalData(sl, "Y"); Z = GetNodalData(sl, "Z")
            T = thickness_on(sl)
            if X and T is not None and len(T) == len(X):
                nn = nearest_nodes([(X[i], Y[i], Z[i]) for i in range(len(X))], FX, FY, FZ, axis, value)
                pts = [(X[i], Y[i], Z[i], T[i], nn[i][0], nn[i][1]) for i in range(len(X))]
        except Exception as e:
            print("       reading %s-Slice failed: %s" % (axis, str(e)))
        try:
            Delete([sl])
        except Exception:
            pass
        if pts:
            return pts, "%s-Slice" % axis
    # fallback: HZO surface nodes close to the plane (exact node indices)
    T = thickness_on(film)
    if T is None:
        return None, "failed"
    c = FX if axis == "X" else FY
    pts = [(FX[i], FY[i], FZ[i], T[i], i, 0.0) for i in range(len(FX)) if abs(c[i] - value) <= CUT_TOL]
    return pts, "HZO nodes within %g of plane (fallback)" % CUT_TOL


# ------------------------------ geometry ----------------------------------
def classify(x, y, z):
    """(r, depth, region, thk_geo) for one point; region based on the original hole."""
    r = math.hypot(x - HOLE_CENTER[0], y - HOLE_CENTER[1])
    d = Z_TOP - z
    R, D = HOLE_RADIUS, HOLE_DEPTH
    cm = CORNER_MARGIN_FRACTION * R
    d_top = math.hypot(max(R - r, 0.0), d)
    d_wall = math.hypot(R - r, max(-d, 0.0, d - D))
    d_bot = math.hypot(max(r - R, 0.0), d - D)
    t = min(d_top, d_wall, d_bot)
    if t == d_top and r >= R:
        region = "field"
    elif t == d_bot and r <= R:
        region = "bottom"
    elif t == d_wall and 0.0 <= d <= D:
        region = "sidewall"
    else:
        region = "corner"
    if math.hypot(R - r, d) < cm + t and d < cm:                 # top corner
        region = "corner"
    if (R - r) < cm + t and (D - d) < cm + t and r < R:          # bottom corner
        region = "corner"
    return r, d, region, t


def mean(v):
    return sum(v) / len(v) if v else float("nan")


def fmt(v):
    return "nan" if v != v else "%.8g" % v


# ------------------------------- main -------------------------------------
if GetActiveModel() is None:
    raise Exception("No active model. Import " + START_FILE + " by hand first.")
reader = find_reader()
if reader is None:
    raise Exception("Could not find the file reader for a DTF in " + DTF_DIR)

if Z_TOP is None:
    base = get_object(BASE_OBJECT)
    if base is None:
        raise Exception("Z_TOP is None and object '%s' not found - set Z_TOP by hand." % BASE_OBJECT)
    Z_TOP = max(GetNodalData(base, "Z"))
    print("       Z_TOP taken from %s: %.8g" % (BASE_OBJECT, Z_TOP))

CUTS = [("x-cut", "X", CUT_X, 1),   # X-Slice: position along the cut = Y (index 1)
        ("y-cut", "Y", CUT_Y, 0)]   # Y-Slice: position along the cut = X (index 0)

summaries = {}
for folder, axis, value, s_idx in CUTS:
    d = os.path.join(OUT_DIR, folder)
    if not os.path.isdir(d):
        os.makedirs(d)
    f = open(os.path.join(d, "summary.csv"), "w")
    cols = ["label", "source_file", "method", "n_points", "min", "max", "mean", "field", "bottom"]
    for fr in SIDEWALL_DEPTH_FRACTIONS:
        cols += ["side@%.3gD_left" % fr, "side@%.3gD_right" % fr]
    cols += ["stepcov_bottom", "max_abs_thk_minus_geo", "max_node_dist"]
    f.write(",".join(cols) + "\n")
    summaries[folder] = f

for n in range(1, N_CYCLES + 1):
    for hc in HALF_CYCLES:
        label = "%d-%s" % (n, hc)
        dtf = last_timestep_file(n, hc)
        if dtf is None:
            print("[skip] no numbered DTF file for " + label)
            continue
        load_file(reader, dtf)
        film = get_object(FILM_OBJECT)
        if film is None:
            print("[skip] object '%s' not found in %s" % (FILM_OBJECT, dtf))
            continue
        FX = GetNodalData(film, "X"); FY = GetNodalData(film, "Y"); FZ = GetNodalData(film, "Z")

        for folder, axis, value, s_idx in CUTS:
            pts, method = cut_points(film, FX, FY, FZ, axis, value)
            if not pts:
                print("[skip] %s %s: no points (%s)" % (label, folder, method))
                continue

            centre = HOLE_CENTER[s_idx]
            rows = []
            for (x, y, z, t, node, nd) in pts:
                s = (x, y)[s_idx]
                r, dep, region, tg = classify(x, y, z)
                side = "left" if s < centre else "right"
                rows.append(dict(X=x, Y=y, Z=z, s=s, offset=s - centre, depth=dep, side=side,
                                 region=region, thickness=t, thk_geo=tg, hzo_node=node, node_dist=nd))
            rows.sort(key=lambda rw: rw["offset"])

            outd = os.path.join(OUT_DIR, folder)
            fo = open(os.path.join(outd, label + ".txt"), "w")
            fo.write("# source: %s\n# cut: %s = %.10g (%s)\n# expression: %s\n"
                     % (dtf, axis, value, method, EXPRESSION))
            fo.write("# geometry: R=%.10g D=%.10g Z_TOP=%.10g XC=%.10g YC=%.10g cut_axis=%s cut_value=%.10g\n"
                     % (HOLE_RADIUS, HOLE_DEPTH, Z_TOP, HOLE_CENTER[0], HOLE_CENTER[1], axis, value))
            fo.write("X\tY\tZ\ts\toffset\tdepth\tside\tregion\tthickness\tthk_geo\thzo_node\tnode_dist\n")
            for rw in rows:
                fo.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\n"
                         % (fmt(rw["X"]), fmt(rw["Y"]), fmt(rw["Z"]), fmt(rw["s"]), fmt(rw["offset"]),
                            fmt(rw["depth"]), rw["side"], rw["region"], fmt(rw["thickness"]),
                            fmt(rw["thk_geo"]), rw["hzo_node"], fmt(rw["node_dist"])))
            fo.close()

            T = [rw["thickness"] for rw in rows]
            field = [rw["thickness"] for rw in rows if rw["region"] == "field"]
            bottom = [rw["thickness"] for rw in rows if rw["region"] == "bottom"
                      and abs(rw["offset"]) < BOTTOM_R_FRACTION * HOLE_RADIUS]
            side_rows = [rw for rw in rows if rw["region"] == "sidewall"]
            win = DEPTH_WINDOW_FRACTION * HOLE_DEPTH
            vals = [label, os.path.basename(dtf), method, str(len(T)),
                    fmt(min(T)), fmt(max(T)), fmt(mean(T)), fmt(mean(field)), fmt(mean(bottom))]
            for fr in SIDEWALL_DEPTH_FRACTIONS:
                dp = fr * HOLE_DEPTH
                for sd in ("left", "right"):
                    vals.append(fmt(mean([rw["thickness"] for rw in side_rows
                                          if rw["side"] == sd and abs(rw["depth"] - dp) <= win])))
            mf, mb = mean(field), mean(bottom)
            vals.append(fmt(mb / mf) if mf == mf and mf != 0 and mb == mb else "nan")
            diffs = [abs(rw["thickness"] - rw["thk_geo"]) for rw in rows if rw["region"] != "corner"]
            vals.append(fmt(max(diffs) if diffs else float("nan")))
            vals.append(fmt(max(rw["node_dist"] for rw in rows)))
            summaries[folder].write(",".join(vals) + "\n")
            summaries[folder].flush()
            print("[ok]   %s %s  points=%d  mean=%s  (%s)" % (label, folder, len(T), fmt(mean(T)), method))

for f in summaries.values():
    f.close()
print("Done. Now run har_thickness_analyzer.py on " + OUT_DIR)
