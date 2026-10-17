# -*- coding: utf-8 -*-
"""CFD-TOPO+ HZO controller: Hf-Cycle -> Zr-Cycle, repeated N times."""

import glob
import os
import re
import time

import GuiVersion
import GuiBC
import GuiFILE
import GuiIC
import GuiMO
import GuiRun


N_CYCLES = 50

MASTER_DIRECTORY = r"J:\Saketh_TOPO\Calibration_w_o_O-recom"
OUTPUT_DIRECTORY = r"J:\Saketh_TOPO\Calibration_w_o_O-recom\0.125"

HF_CYCLE_TEMPLATE_DTF = os.path.join(MASTER_DIRECTORY, "Hf-125.DTF")
ZR_CYCLE_TEMPLATE_DTF = os.path.join(MASTER_DIRECTORY, "Zr-125.DTF")

DRY_RUN = False

HALF_CYCLE_SECONDS = 5.125
O_PULSE_SECONDS = 0.125
O_PULSE_EPSILON = 0.0000001

FLUX_BOUNDARY = [67]
TIME_PROFILE_SPECIES = ("O", "O2", "TDMAHF", "TDMAZ")

MODEL_OUTPUT_FILE_PATTERN = "{stem}.out"
COMPLETION_TEXT = "Normal Termination"
POLL_SECONDS = 10
MAX_SOLVER_WAIT_SECONDS = 86400
OUTPUT_RULE = "newest_dtf"


def _require_file(path, description):
    if not os.path.isfile(path):
        raise RuntimeError("%s does not exist: %s" % (description, path))


def _ensure_directory(path):
    if os.path.isdir(path):
        return

    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise


def _working_dtf_for_template(template_dtf, cycle, pulse_name):
    return os.path.join(
        OUTPUT_DIRECTORY,
        "Cycle-%d-%s-cycle.DTF" % (cycle, pulse_name),
    )


def _model_output_file(working_dtf):
    stem = os.path.splitext(os.path.basename(working_dtf))[0]

    return os.path.join(
        os.path.dirname(working_dtf),
        MODEL_OUTPUT_FILE_PATTERN.format(stem=stem),
    )


def wait_for_solver_completion(working_dtf):
    output_file = _model_output_file(working_dtf)
    started_at = time.time()
    completion_marker = COMPLETION_TEXT.lower().encode("ascii")

    print("  waiting for completion marker in: %s" % output_file)

    while time.time() - started_at <= MAX_SOLVER_WAIT_SECONDS:
        if os.path.isfile(output_file):
            with open(output_file, "rb") as stream:
                if completion_marker in stream.read().lower():
                    print("  solver completion marker found")
                    return

        time.sleep(POLL_SECONDS)

    raise RuntimeError(
        "Timed out waiting for '%s' in %s. Do not start the next half-cycle."
        % (COMPLETION_TEXT, output_file)
    )


def find_final_output(working_dtf, submitted_at):
    if OUTPUT_RULE != "newest_dtf":
        raise RuntimeError("No verified output-selection rule is configured.")

    working_stem = os.path.splitext(os.path.basename(working_dtf))[0]

    timestep_name = re.compile(
        r"^%s\.\d+\.DTF$" % re.escape(working_stem),
        re.IGNORECASE,
    )

    candidates = [
        path
        for path in glob.glob(
            os.path.join(os.path.dirname(working_dtf), "*.DTF")
        )
        if os.path.getmtime(path) >= submitted_at
        and timestep_name.match(os.path.basename(path))
    ]

    if not candidates:
        raise RuntimeError(
            "No final timestep DTF appeared after Submit. "
            "Do not start the next half-cycle."
        )

    return max(candidates, key=os.path.getmtime)


def _set_flux_time_profile(end_time):
    start_time = end_time - HALF_CYCLE_SECONDS

    profile_times = (
        start_time,
        start_time + O_PULSE_SECONDS,
        start_time + O_PULSE_SECONDS + O_PULSE_EPSILON,
        end_time,
    )

    for species in TIME_PROFILE_SPECIES:
        for index, profile_time in enumerate(profile_times, 1):
            GuiBC.SetArray(
                FLUX_BOUNDARY,
                "Chem/Flux_t",
                [species, index],
                profile_time,
            )


def run_pulse(cycle, half_index, pulse_name, template_dtf, restart_dtf):
    half_cycle_number = (cycle - 1) * 2 + half_index
    end_time = round(half_cycle_number * HALF_CYCLE_SECONDS, 6)

    working_dtf = _working_dtf_for_template(
        template_dtf,
        cycle,
        pulse_name,
    )

    print("Cycle %d / %s" % (cycle, pulse_name))
    print("  template: %s" % template_dtf)
    print("  EndTime:  %.6f" % end_time)
    print("  restart:  %s" % (
        restart_dtf or "template initial condition"
    ))
    print("  working:  %s" % working_dtf)

    if DRY_RUN:
        return "<final timestep DTF from Cycle-%d-%s-cycle.DTF>" % (
            cycle,
            pulse_name,
        )

    _require_file(template_dtf, "%s template" % pulse_name)

    if restart_dtf is not None:
        _require_file(restart_dtf, "%s restart DTF" % pulse_name)

    if os.path.exists(working_dtf):
        raise RuntimeError(
            "Refusing to overwrite existing working DTF: %s" % working_dtf
        )

    GuiFILE.Open(template_dtf)

    GuiMO.Set("Shared/EndTime", end_time)

    _set_flux_time_profile(end_time)

    if restart_dtf is not None:
        GuiIC.Set("Shared/Restart", 1)
        GuiIC.Set("Shared/Continuation_Run", 1)
        GuiIC.Set("Shared/RestartFile_Name", restart_dtf)

    GuiFILE.SaveAs(working_dtf)

    submitted_at = time.time()

    GuiRun.Submit()

    wait_for_solver_completion(working_dtf)

    final_dtf = find_final_output(working_dtf, submitted_at)

    print("  selected final timestep: %s" % final_dtf)

    GuiFILE.CloseCurrentSim()

    return final_dtf


def run_cycles(n_cycles):
    if n_cycles < 1:
        raise ValueError("N_CYCLES must be at least 1")

    if not DRY_RUN:
        _ensure_directory(OUTPUT_DIRECTORY)

    previous_final_dtf = None

    for cycle in range(1, n_cycles + 1):
        hf_final_dtf = run_pulse(
            cycle,
            1,
            "Hf",
            HF_CYCLE_TEMPLATE_DTF,
            previous_final_dtf,
        )

        zr_final_dtf = run_pulse(
            cycle,
            2,
            "Zr",
            ZR_CYCLE_TEMPLATE_DTF,
            hf_final_dtf,
        )

        previous_final_dtf = zr_final_dtf

    print("Completed %d Hf -> Zr cycle(s)." % n_cycles)


GuiVersion.RecordVersion("2026.0.0.14417")
GuiFILE.SetMode("TOPO")

run_cycles(N_CYCLES)