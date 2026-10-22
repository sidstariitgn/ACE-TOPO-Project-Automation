# hzo_halfcycle_thickness.py
# Run inside CFD-VIEW:  File > Run Python Script
#
# HOW TO USE
#   1. In CFD-VIEW: File > New, then File > Import with your usual settings
#      (DTF, Original Layout, VC and BC Groups, Virtual Data, Visualization Objects ...)
#      and import the FIRST file this script will use, i.e. the last time step of
#      cycle 1 Hf:   Cycle-1-Hf-cycle.0000182.DTF   (set START_FILE below to the same name).
#   2. Check that SUB and HZO are visible, then run this script.
#
# The script keeps that model and its import settings, and for every cycle n and half
# cycle (Hf, Zr) it
#   - switches the reader to the LAST time-step DTF of that half cycle
#   - evaluates the calculator expression on OBJECT_NAME
#   - writes  <OUT_DIR>/<n>-<Hf|Zr>.txt   (node index, X, Y, Z, thickness)
# and appends min / max / mean to <OUT_DIR>/summary.csv

import os
import re
import glob

# ----------------------------- USER SETTINGS -----------------------------
DTF_DIR     = r"J:\Saketh_TOPO\Test"
OUT_DIR     = r"J:\Saketh_TOPO\Test\thickness"
N_CYCLES    = 1
HALF_CYCLES = ["Hf", "Zr"]

# File pattern: Cycle-<n>-<Hf|Zr>-cycle.<time step>.DTF
# Only files whose suffix is a number are used; the highest number = last time step.
FILE_PATTERN = "Cycle-{n}-{hc}-cycle.*.DTF"

# File you imported by hand in step 1 (file name only, no folder)
START_FILE  = "Cycle-1-Hf-cycle.0000182.DTF"

# Object you select in the right-hand column (name as in the Object Explorer)
OBJECT_NAME = "HZO"

# Calculator expression, exactly what you type in the calculator
EXPRESSION  = '"SUB"-"HZO"'      # same as typed by hand (no 'Thk =' assignment)

UNIT_SCALE  = 1.0       # e.g. 1e9 if the data is in metres and you want nm
# --------------------------------------------------------------------------


def last_timestep_file(n, hc):
    """Among Cycle-n-hc-cycle.<number>.DTF return the file with the largest number.
    Files like '...cycle.IC.DTF' or '...cycle.DTF' (no number) are ignored."""
    pattern = os.path.join(DTF_DIR, FILE_PATTERN.format(n=n, hc=hc))
    best, best_num = None, -1
    for p in glob.glob(pattern):
        m = re.search(r"-cycle\.(\d+)\.dtf$", os.path.basename(p), re.IGNORECASE)
        if m and int(m.group(1)) > best_num:
            best, best_num = p, int(m.group(1))
    return best


def get_coords(obj):
    """Return (X, Y, Z) nodal lists, or (None, None, None) if unavailable."""
    try:
        return (GetNodalData(obj, "X"), GetNodalData(obj, "Y"), GetNodalData(obj, "Z"))
    except Exception:
        return (None, None, None)


def get_thickness(obj):
    """Evaluate EXPRESSION on obj and return a list with one value per node.
    Tries several ways, because GetExpressionResult can fail even when the GUI
    calculator shows the result."""
    errors = []

    # 1) normal route: set expression, read its result
    try:
        SetExpression(obj, EXPRESSION)
        r = GetExpressionResult(obj)
        return list(r) if hasattr(r, "__len__") else [r]
    except Exception as e:
        errors.append("GetExpressionResult: " + str(e))

    # 2) force the data-flow network to update, then try again
    try:
        ForceUpdate()
        r = GetExpressionResult(obj)
        return list(r) if hasattr(r, "__len__") else [r]
    except Exception as e:
        errors.append("after ForceUpdate: " + str(e))

    # 3) compute SUB - HZO ourselves from the nodal data of the two variables
    try:
        a = GetNodalData(obj, "SUB")
        b = GetNodalData(obj, "HZO")
        return [a[i] - b[i] for i in range(len(a))]
    except Exception as e:
        errors.append("GetNodalData SUB/HZO: " + str(e))

    for line in errors:
        print("       " + line)
    return None


def load_file(reader, path):
    """Point the existing reader at another DTF file, keeping its import settings."""
    try:
        reader.switchFile(path, reader.getSwitches(), reader.getSimulation(), True)
        Reflow()
    except Exception as e:
        print("       switchFile failed (%s) -> trying SetFilenames/Reread" % str(e))
        SetFilenames(reader, [path])
        Reread()


# ------------------------------- main -------------------------------------
if GetActiveModel() is None:
    raise Exception("No active model. Import " + START_FILE + " by hand first (see top of script).")

def find_reader():
    """Find the file reader of the model you imported by hand. Tries START_FILE first,
    then every DTF file in DTF_DIR (with and without the folder in the name)."""
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


reader = find_reader()
if reader is None:
    print("No file reader found. Objects in model: %d" % len(GetAllObjects()))
    raise Exception("Could not find the file reader. Make sure a DTF from " + DTF_DIR +
                    " is imported and selected as the active model, then run again.")

if not os.path.isdir(OUT_DIR):
    os.makedirs(OUT_DIR)

summary = open(os.path.join(OUT_DIR, "summary.csv"), "w")
summary.write("label,source_file,nodes,min,max,mean\n")

for n in range(1, N_CYCLES + 1):
    for hc in HALF_CYCLES:
        label = "%d-%s" % (n, hc)
        dtf = last_timestep_file(n, hc)
        if dtf is None:
            print("[skip] no numbered DTF file found for " + label)
            continue

        load_file(reader, dtf)

        objs = GetObjectsByName(OBJECT_NAME)
        if not objs:
            print("[skip] object '" + OBJECT_NAME + "' not found after loading " + dtf)
            print("       %d surfaces, %d volumes in the model" % (len(GetAllSurfaces()), len(GetAllVolumes())))
            continue
        obj = objs[0]

        thk = get_thickness(obj)
        if thk is None:
            print("[skip] could not get the thickness values for " + label)
            continue
        thk = [v * UNIT_SCALE for v in thk]
        X, Y, Z = get_coords(obj)

        out_path = os.path.join(OUT_DIR, label + ".txt")
        f = open(out_path, "w")
        f.write("# source: %s\n" % dtf)
        f.write("# expression: %s\n" % EXPRESSION)
        f.write("# node index is 0-based, same ordering as the object's nodes in CFD-VIEW\n")
        f.write("node\tX\tY\tZ\tthickness\n")
        for i in range(len(thk)):
            if X is None:
                f.write("%d\t\t\t\t%.8g\n" % (i, thk[i]))
            else:
                f.write("%d\t%.8g\t%.8g\t%.8g\t%.8g\n" % (i, X[i], Y[i], Z[i], thk[i]))
        f.close()

        mean = sum(thk) / len(thk)
        summary.write("%s,%s,%d,%.8g,%.8g,%.8g\n"
                      % (label, os.path.basename(dtf), len(thk), min(thk), max(thk), mean))
        summary.flush()
        print("[ok]   %s  nodes=%d  mean=%.5g  -> %s" % (label, len(thk), mean, out_path))

summary.close()
print("Done.")
