"""Which screenshot of the tutorial shows what: one row per picture in docs/tutorial/img (plain data, used by update_tutorial.py).

    shot name: (keys to pass to capture_shots.py, run shown, dashboard source files that draw it, sweeps shown ("" = none))

Keep it in step with capture_shots.py: when a shot is added there, add a row here (update_tutorial.py status says if one is missing).
A shot is reported STALE when one of its source files, its run's results or its sweeps changed after it was captured.
"""
RUNS = {                                            # the same names as RUN_A ... RUN_E in capture_shots.py
    "A": "loci_demo_seed11_24bba2",                 # showcase: Vutara, bursty gene, MS2, probes, calibration
    "B": "tracking_demo_seed3_e15523",              # tracking errors
    "C": "loci_demo_seed8_f5c88d",                  # measured error effects
    "D": "loci_3d_seed23_2163f4",                   # 3D z-stacks + piezo + beads
    "E": "loci_library_seed12_87a16d",              # polymer library, two copies, stray spots
}
DASH = "app/dashboard.py"
SHOTS = {
    "home": (["home"], "A", [DASH, "app/style.py"], ""),
    "guide": (["guide"], "A", [DASH], ""),
    "glossary": (["guide"], "A", [DASH], ""),
    "sb_sim": (["sb_sim"], "A", [DASH], ""),
    "sb_loci": (["sb_loci"], "A", [DASH], ""),
    "pg_movie": (["movie"], "B", [DASH, "app/fastview.py"], ""),
    "pg_movie_player": (["movie_player"], "B", [DASH, "app/player.py", "app/fastview.py"], ""),
    "pg_metrics": (["metrics"], "B", [DASH], ""),
    "pg_traj": (["traj"], "B", [DASH], ""),
    "pg_scans": (["scans"], "B", [DASH, "app/sweep_plot.py"], "speed*"),
    "pg_cells": (["cells"], "A", ["app/cells_tab.py"], ""),
    "pg_cells_3d": (["cells_3d"], "D", ["app/cells_tab.py"], ""),
    "pg_cells_copies": (["cells_copies"], "E", ["app/cells_tab.py"], ""),
    "pg_spots": (["spots"], "A", ["app/spots_tab.py"], "ms2_brightness"),
    "pg_spots_copies": (["spots_copies"], "E", ["app/spots_tab.py"], "copies*"),
    "pg_metrics_prop": (["metrics_prop"], "E", [DASH], ""),
    "pg_images": (["images"], "A", ["app/images_tab.py"], ""),
    "pg_images_build": (["images", "images_build"], "A", ["app/images_tab.py"], ""),
    "pg_images_err": (["images_err"], "C", ["app/images_tab.py", "app/errors_view.py"], ""),
    "pg_cal": (["cal"], "A", ["app/calibration_tab.py"], ""),
    "pg_stage": (["stage"], "D", ["app/stage_tab.py"], ""),
    "pg_tx": (["tx"], "A", ["app/transcription_tab.py"], ""),
    "pg_ms2": (["ms2"], "A", ["app/ms2_tab.py"], ""),
    "pg_ep": (["ep", "ep_run"], "A", ["app/ep_tab.py", "app/sweep_plot.py"], "ep_*"),
    "pg_probes": (["probes"], "A", ["app/probe_tab.py"], ""),
    "pg_lib": (["lib"], "A", ["app/library_tab.py"], ""),
    "pg_chrom": (["chrom"], "A", ["app/chromatin_tab.py", "app/chromatin_view.py"], ""),
    "pg_cfg": (["cfg"], "A", [DASH], ""),
}
