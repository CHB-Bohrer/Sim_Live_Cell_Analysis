"""Tall page captures + element rectangles for the tutorial.  python grab.py <outdir> [shot ...]"""
import json, pathlib, sys, time
from playwright.sync_api import sync_playwright

URL = "http://localhost:8502"
W, TALL = 1600, 9000
RUN_A = "loci_demo_seed11_24bba2"          # showcase: Vutara, bursty gene, MS2, probes, calibration
RUN_B = "tracking_demo_seed3_e15523"       # tracking errors (33 identity switches)
RUN_C = "loci_demo_seed8_f5c88d"           # has measured error effects
RUN_D = "loci_3d_seed23_2163f4"            # 3D z-stacks + piezo + beads
RUN_E = "loci_library_seed12_87a16d"       # polymer library, two promoter-enhancer copies in some cells, stray spots, MS2 (2026-10-19)
OUT = pathlib.Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
ONLY = set(sys.argv[2:])

EC = "[data-testid=stElementContainer]"
HB = "[data-testid=stHorizontalBlock]"
EXP = "[data-testid=stExpander]"
ALERT = "[data-testid=stAlert]"


def idle(pg, extra=1.2, timeout=180):
    t0 = time.time(); time.sleep(0.7)
    while time.time() - t0 < timeout:
        if pg.locator("[data-testid=stStatusWidget]").count() == 0:
            break
        time.sleep(0.4)
    time.sleep(extra)


def wait_gone(pg, text, timeout=240):
    t0 = time.time()
    while time.time() - t0 < timeout:
        if pg.get_by_text(text).locator("visible=true").count() == 0:
            return
        time.sleep(0.5)


def goto_page(pg, name):
    pg.locator("[data-testid=stButtonGroup] button", has_text=name).first.click(); idle(pg, 2.5)


def select_run(pg, run):
    sel = pg.locator("[data-testid=stSidebar] [data-testid=stSelectbox]").last
    if run in sel.inner_text():
        return
    sel.click(); time.sleep(0.4); pg.keyboard.type(run); time.sleep(0.5); pg.keyboard.press("Enter"); idle(pg, 3)


def rect(pg, spec, ox, oy):
    """spec = (text, up_css) | ('css', selector); the smallest visible `up` element whose text contains `text`; [x, y, w, h] relative to the crop origin."""
    kind, up = spec
    if kind in ("css", "css_last"):
        loc = pg.locator(up).locator("visible=true")
        loc = loc.first if kind == "css" else loc.last
        if loc.count() == 0:
            return None
        r = loc.evaluate("el => { const b = el.getBoundingClientRect(); return [b.x, b.y, b.width, b.height]; }")
    else:
        r = pg.evaluate("""([text, up]) => {
            const sel = up === '[data-testid=stExpander]' ? up + ' summary' : up;
            const els = [...document.querySelectorAll(sel)].filter(e => (e.innerText || '').toLowerCase().includes(text.toLowerCase()) && e.getBoundingClientRect().height > 0);
            if (!els.length) return null;
            els.sort((a, b) => a.getBoundingClientRect().width * a.getBoundingClientRect().height - b.getBoundingClientRect().width * b.getBoundingClientRect().height);
            let el = els[0]; if (up === '[data-testid=stExpander]') el = el.closest(up) || el;
            const b = el.getBoundingClientRect(); return [b.x, b.y, b.width, b.height]; }""", [kind, up])
        if r is None:
            return None
    return [round(r[0] - ox), round(r[1] - oy), round(r[2]), round(r[3])]


def exact(pg, text, up, ox, oy, nth=0):
    loc = pg.get_by_text(text, exact=True).locator("visible=true")
    if loc.count() <= nth:
        return None
    r = loc.nth(nth).evaluate("(el, up) => { const t = el.closest(up) || el; const b = t.getBoundingClientRect(); return [b.x, b.y, b.width, b.height]; }", up)
    return [round(r[0] - ox), round(r[1] - oy), round(r[2]), round(r[3])]


def save_page(pg, name, anchors, exacts=None, top_css="[data-testid=stButtonGroup]", crop_x=300):
    """Tall capture of the main column from the navigation bar down; anchors give named rectangles."""
    nav = pg.locator(top_css).first.bounding_box()
    oy = max(nav["y"] - 10, 0)
    bottom = pg.evaluate("(() => { const b=document.querySelector('[data-testid=stMainBlockContainer]').getBoundingClientRect(); return b.bottom; })()")
    h = int(min(bottom + 30, TALL) - oy)
    pg.screenshot(path=str(OUT / f"{name}.jpg"), type="jpeg", quality=84, clip={"x": crop_x, "y": oy, "width": W - crop_x, "height": h})
    rects = {}
    for k, spec in anchors.items():
        r = rect(pg, spec, crop_x, oy)
        if r is None:
            print("   (missing)", name, k, flush=True)
        else:
            rects[k] = r
    for k, (text, up) in (exacts or {}).items():
        r = exact(pg, text, up, crop_x, oy)
        if r is not None:
            rects[k] = r
    (OUT / f"{name}.json").write_text(json.dumps({"w": W - crop_x, "h": h, "rects": rects}, ensure_ascii=False, indent=0), encoding="utf-8")
    print("saved", name, W - crop_x, h, len(rects), flush=True)


def click_text(pg, text, exact_=False, nth=0):
    """Open/close an expander by its label, or press a button by its name."""
    btn = pg.get_by_role("button", name=text).locator("visible=true")
    ex = pg.locator("[data-testid=stExpander] summary", has_text=text).locator("visible=true")
    if btn.count() and not ex.count():
        btn.first.click()
    elif ex.count():
        ex.first.click()
    else:
        btn.first.click()
    idle(pg, 1.5)


INTRO = ("WHAT THIS SHOWS", HB)


def want(n):
    return not ONLY or n in ONLY


def main():
    with sync_playwright() as p:
        b = p.chromium.launch(channel="msedge", headless=True)
        pg = b.new_page(viewport={"width": W, "height": 1000}, device_scale_factor=int(__import__("os").environ.get("DSF","1")))
        pg.goto(URL); pg.wait_for_selector("text=Sim Live Cell Analysis", timeout=90000); idle(pg, 3)
        select_run(pg, RUN_A)

        # ---- home (normal viewport, with the sidebar)
        if want("home"):
            pg.evaluate("document.querySelector('[data-testid=stMain]').scrollTop = 0"); time.sleep(0.5)
            pg.screenshot(path=str(OUT / "home.jpg"), type="jpeg", quality=86)
            A = {"sidebar": ("css", "[data-testid=stSidebar]"), "hero": ("Cell tracking: how well did we do?", "[data-testid=stMarkdownContainer]"),
                 "guide": ("New here?", EXP), "glossary": ("Glossary of terms", EXP), "summary": ("What was simulated in this run", EXP),
                 "errglance": ("Image errors at a glance", EXP), "nav": ("css", "[data-testid=stButtonGroup]"),
                 "runsel": ("css", "[data-testid=stSidebar] [data-testid=stSelectbox]"), "newsim": ("New simulation", EXP), "newloci": ("New single-cell loci run", EXP),
                 "scores_auto": ("Automatic segmentation (realistic)", EC), "scores_gt": ("Ground-truth masks (linking error only)", EC)}
            rr = {k: rect(pg, s, 0, 0) for k, s in A.items()}
            rr = {k: v for k, v in rr.items() if v}
            (OUT / "home.json").write_text(json.dumps({"w": W, "h": 1000, "rects": rr}), encoding="utf-8")
            print("saved home", flush=True)

        # ---- tall viewport for everything else
        pg.set_viewport_size({"width": W, "height": TALL}); time.sleep(1.5); idle(pg, 2)
        pg.add_style_tag(content=".stApp{background-image:none !important}")   # the fixed page background stretches over a 9000 px viewport

        if want("guide"):
            click_text(pg, "New here?")
            pg.evaluate("document.querySelector('[data-testid=stMain]').scrollTop = 0")
            nav = pg.locator("[data-testid=stButtonGroup]").first.bounding_box()
            exp = pg.locator(EXP, has_text="New here?").first.bounding_box()
            pg.screenshot(path=str(OUT / "guide.jpg"), type="jpeg", quality=86, clip={"x": 300, "y": exp["y"] - 10, "width": W - 300, "height": exp["height"] + 20})
            (OUT / "guide.json").write_text(json.dumps({"w": W - 300, "h": int(exp["height"] + 20), "rects": {}}), encoding="utf-8")
            click_text(pg, "New here?")
            click_text(pg, "Glossary of terms")
            exp = pg.locator(EXP, has_text="Glossary of terms").first.bounding_box()
            pg.screenshot(path=str(OUT / "glossary.jpg"), type="jpeg", quality=86, clip={"x": 300, "y": exp["y"] - 10, "width": W - 300, "height": exp["height"] + 20})
            (OUT / "glossary.json").write_text(json.dumps({"w": W - 300, "h": int(exp["height"] + 20), "rects": {}}), encoding="utf-8")
            click_text(pg, "Glossary of terms")
            print("saved guide, glossary", flush=True)

        # ---- sidebar forms (whole form as one tall image)
        for nm, label, anchors in [
            ("sb_sim", "New simulation", ["Random seed", "Acquisition", "Cells and motion", "Nuclear shape", "Nuclear texture", "Imaging", "Segmentation and tracker", "Image errors", "Run simulation"]),
            ("sb_loci", "New single-cell loci run", ["Random seed", "Cells in the field", "Frame interval", "Cell speed scale", "Photons per locus per second", "Microscope",
                                                      "Field of view", "Room outside the field", "Bleach each cell", "3D z-stacks", "Record 3D z-stacks", "Calibration",
                                                      "Where do the locus movements come from?", "Model probe binding", "Model distance-dependent transcription",
                                                      "Promoter-enhancer copies and stray spots", "Share of cells with two copies", "Stray spots: new ones per nucleus per minute (each colour)",
                                                      "Spot finding (stage 6)", "Copies per cell the spot finder looks for", "Pair gate",
                                                      "Image errors", "Run single-cell loci simulation"])]:
            if not want(nm):
                continue
            pg.locator("[data-testid=stSidebar]").get_by_text(label).first.click(); idle(pg, 1.5)
            sb = pg.locator("[data-testid=stSidebar]").first
            box = sb.bounding_box()
            bottom = pg.locator('[data-testid=stSidebar] [data-testid=stSelectbox]').locator('visible=true').last.bounding_box()
            bottom = bottom['y'] + bottom['height'] + 10
            h = int(min(bottom + 20, TALL))
            pg2 = pg
            sbimg = OUT / f"{nm}.jpg"
            pg.screenshot(path=str(sbimg), type="jpeg", quality=86, clip={"x": 0, "y": 0, "width": 300, "height": h})
            rr = {}
            for a in anchors:
                sbl = pg.locator("[data-testid=stSidebar]")
                r = sbl.get_by_text(a, exact=True).locator("visible=true")
                if r.count() == 0:
                    r = sbl.get_by_text(a, exact=False).locator("visible=true")
                if r.count():
                    bb = r.first.evaluate("el => { const t = el.closest('[data-testid=stElementContainer]') || el; const b = t.getBoundingClientRect(); return [b.x, b.y, b.width, b.height]; }")
                    rr[a] = [round(v) for v in bb]
            (OUT / f"{nm}.json").write_text(json.dumps({"w": 300, "h": h, "dpr": int(__import__("os").environ.get("DSF","1")), "rects": rr}, ensure_ascii=False), encoding="utf-8")
            pg.locator("[data-testid=stSidebar]").get_by_text(label).first.click(); idle(pg, 1)
            print("saved", nm, h, len(rr), flush=True)

        # ---- movie (run B, with identity switches)
        if want("movie") or want("movie_player"):
            select_run(pg, RUN_B); goto_page(pg, "Movie")
            if want("movie"):
                click_text(pg, "Next error", False); idle(pg, 2)
                save_page(pg, "pg_movie", {"intro": INTRO, "player": ("Play the movie", EXP), "tracked": ("Tracked on:", EC), "prev": ("Previous error", EC), "next": ("Next error", EC),
                                           "raw": ("Also show raw image panel", EC), "slider": ("css", "[data-testid=stSlider]"), "errcount": ("frames contain an error", EC),
                                           "images": ("css", "[data-testid=stImage]"), "errbox": ("Errors in frame", ALERT), "howto": ("How to read the pictures", EXP),
                                           "bars": ("Identity switches over time", EC)})
            if want("movie_player"):
                click_text(pg, "Play the movie")
                pg.locator("[data-testid=stToggle], [data-testid=stCheckbox]", has_text="Prepare player").first.click(); idle(pg, 6, 240); time.sleep(6)
                save_page(pg, "pg_movie_player", {"intro": INTRO, "player": ("Play the movie", EXP), "prepare": ("Prepare player", EC)})
                click_text(pg, "Play the movie")

        if want("metrics"):
            goto_page(pg, "Metrics")
            click_text(pg, "What do these numbers mean?")
            save_page(pg, "pg_metrics", {"intro": INTRO, "inwords": ("In words.", ALERT), "table": ("css", "[data-testid=stDataFrame]"), "meaning": ("What do these numbers mean?", EXP),
                                         "percell": ("Per-cell identity table for:", EC)})
        if want("traj"):
            goto_page(pg, "Trajectories")
            save_page(pg, "pg_traj", {"intro": INTRO, "tracked": ("Tracked on:", EC), "plot": ("css", "[data-testid=stImage]")})
        if want("scans"):
            goto_page(pg, "Scans")
            # choose the speed scan
            sel = pg.locator("[data-testid=stMain] [data-testid=stSelectbox]").first
            sel.click(); time.sleep(0.4); pg.keyboard.type("speed_scale"); time.sleep(0.5); pg.keyboard.press("Enter"); idle(pg, 3)
            save_page(pg, "pg_scans", {"intro": INTRO, "scan": ("css", "[data-testid=stMain] [data-testid=stSelectbox]"), "caption": ("runs (", EC), "metric": ("Metric", EC),
                                       "masks": ("Masks", EC), "roles": ("X axis", EC), "plot": ("css", "[data-testid=stImage]"), "table": ("Table (mean and standard", EXP)})

        # ---- run A pages
        if any(want(n) for n in ("cells", "spots", "images", "cal", "tx", "ms2", "ep", "probes", "lib", "chrom", "cfg")):
            select_run(pg, RUN_A)
        if want("cells"):
            goto_page(pg, "Cells")
            save_page(pg, "pg_cells", {"intro": INTRO, "ident": ("Cells identified by", EC), "table": ("Cell summary table", EXP), "playall": ("Play all cells together", EXP),
                                       "allframe": ("All cells at one frame", EC), "grid": ("css", "[data-testid=stImage]"), "onecell": ("One cell through time", EC),
                                       "cellsel": ("Clean track", ALERT), "legend": ("White + = spot found by the analysis", EC), "traces": ("Solid = copy 0", EC)})
        if want("spots"):
            goto_page(pg, "Spot finding"); wait_gone(pg, "Running the spot finder")
            idle(pg, 2)
            save_page(pg, "pg_spots", {"intro": INTRO, "method": ("1 · The method and where it comes from", EC), "table": ("dog_viterbi: the MS2 spot pipeline", EXP),
                                       "step2": ("2 · Every step on one cell", EC), "settings": ("Settings of the spot finder", EXP), "cellsel": ("Cell (true identity)", EC),
                                       "spotsel": ("locus 0 (promoter)", EC), "frames": ("Frame (this cell)", EC),
                                       "s3": ("3 · This run's result", EC), "s3b": ("3b · Photobleaching correction", EC), "s3c": ("3c · Uneven illumination", EC), "s3d": ("3d · How bright must the MS2 spot be", EC), "s4": ("4 · Compare methods on this run", EC)})
        if want("cells_copies") or want("spots_copies"):          # the two-copy showcase: copy labels, stray spots, the linker's paths
            select_run(pg, RUN_E)
            if want("cells_copies"):
                goto_page(pg, "Cells")
                save_page(pg, "pg_cells_copies", {"intro": INTRO, "onecell": ("One cell through time", EC), "legend": ("White + = spot found by the analysis", EC),
                                                  "fig": ("css_last", "[data-testid=stImage]"), "traces": ("Solid = copy 0", EC)})
            if want("spots_copies"):
                goto_page(pg, "Spot finding"); wait_gone(pg, "Running the spot finder"); wait_gone(pg, "Running the multi-copy linker"); idle(pg, 3)
                save_page(pg, "pg_spots_copies", {"intro": INTRO, "step2": ("2 · Every step on one cell", EC), "copysel": ("Copy", EC),
                                                  "s2b": ("2b · Copies and stray spots", EC), "s3": ("3 · This run's result", EC),
                                                  "s3e": ("3e · Several copies and stray spots", EC)})
            select_run(pg, RUN_A)
        if want("images"):
            goto_page(pg, "Images")
            save_page(pg, "pg_images", {"intro": INTRO, "p1": ("1. The recorded image and the truth behind it", EC), "frame": ("Frame", EC), "channels": ("Channels", EC),
                                        "view": ("Whole field", EC), "draw": ("Draw the ground truth on the images", EC), "errfree": ("Also show the error-free version", EC),
                                        "recorded": ("Recorded: nucleus", EC), "stats": ("Stage drift this frame", EC), "truth": ("Ground truth for cell 1 in this frame", EC),
                                        "p2": ("2 · One error source", EC), "p3": ("3. How this image was built", EC), "build": ("Build the image up", EC),
                                        "p4": ("4. How much each error matters on its own", EC), "p5": ("5. Every source of image error, explained", EC)})
            # build-up variant
            if want("images_build"):
                click_text(pg, "Build the image up"); idle(pg, 15, 240); time.sleep(8)
                save_page(pg, "pg_images_build", {"p3": ("3. How this image was built", EC), "build": ("Build the image up", EC)})
        if want("images_err"):
            select_run(pg, RUN_C); goto_page(pg, "Images")
            sel = pg.get_by_text("Error source", exact=True).locator("visible=true").first
            box = pg.locator("[data-testid=stMain] [data-testid=stSelectbox]", has_text="Error source").first
            box.wait_for(timeout=180000)
            box.click(); time.sleep(0.4); pg.keyboard.type("Chromatic"); time.sleep(0.5); pg.keyboard.press("Enter"); idle(pg, 6)
            save_page(pg, "pg_images_err", {"p2": ("2 · One error source", EC), "source": ("Error source", EC), "what": ("What it is", EC), "does": ("What it does to the picture", EC),
                                            "howfix": ("How it is corrected", EC), "results": ("What it does to the results", EC)})
            select_run(pg, RUN_A)
        if want("cal"):
            goto_page(pg, "Calibration")
            save_page(pg, "pg_cal", {"intro": INTRO, "steps": ("1Acquire", HB), "table": ("Estimated versus true", EC), "dark": ("Dark frames: the camera with the shutter closed", EC),
                                     "flat": ("Flat-field slide: a uniform light source", EC), "bead": ("Bead slide: tiny sources", EC), "help": ("Does it help?", EC)})
        if want("stage"):
            select_run(pg, RUN_D); goto_page(pg, "Stage & fiducials")
            save_page(pg, "pg_stage", {"intro": INTRO, "cards": ("The piezo is asked", HB), "s1": ("1. Where the planes really were", EC), "s2": ("2. The fiducial beads on the slide", EC), "s3": ("3. Does correcting help?", EC)})
            select_run(pg, RUN_A)
        if want("tx"):
            goto_page(pg, "Transcription")
            save_page(pg, "pg_tx", {"intro": INTRO, "diagram": ("From chromosome distance to the MS2 spot", EC), "s1": ("1. How the movie samples time", EC), "s2": ("2. Distance between promoter and enhancer", EC),
                                    "s3": ("3. How the distance sets a rate", EC), "s4": ("4. Promoter switching", EC), "s5": ("5. Polymerases travel along the gene", EC),
                                    "s6": ("6. The MS2 signal", EC), "s7": ("7. What the microscope sees", EC), "s8": ("8. Parameter sweep", EC), "s9": ("9.", EC)})
        if want("ms2"):
            goto_page(pg, "MS2")
            save_page(pg, "pg_ms2", {"intro": INTRO, "run": ("Run (any run with MS2 spots", EC), "s1": ("1 · What goes in", EC), "s2": ("2 · What one Pol II looks like", EC), "s3": ("3 · How noisy the trace is", EC),
                                     "s4": ("4 · One cell: the trace", EC), "burst": ("Burst window w (frames)", EC), "s5": ("5 · The whole run against the truth", EC),
                                     "s5b": ("5b · Has the sampler converged?", EC), "s6": ("6 · Choosing the burst rule", EC)})
        if want("ep"):
            goto_page(pg, "E–P curve")
            if want("ep_run"):
                pg.get_by_role("button", name="Simulate and fit").first.click(); idle(pg, 5, 300)
            save_page(pg, "pg_ep", {"intro": INTRO, "s1": ("1 · Distances from the polymer simulations", EC), "controls": ("Microscope", EC), "s1b": ("1b · The earlier version", EC),
                                    "s2": ("2 · Saved scans", EC), "scan": ("Scan", EC)})
        if want("probes"):
            goto_page(pg, "Probes")
            save_page(pg, "pg_probes", {"intro": INTRO, "cycle": ("How probes bind a locus", EC), "promoter": ("Promoter (locus 0)", EC), "enhancer": ("Enhancer (locus 1)", EC),
                                        "sweep": ("Parameter sweep", EC), "use": ("Use these settings in a movie", EXP)})
        if want("lib"):
            goto_page(pg, "Library"); time.sleep(3)
            save_page(pg, "pg_lib", {"intro": INTRO, "how": ("How simulations reach a movie", EC), "live": ("Library generation: live status", EC), "choose": ("Choose which simulations a movie uses", EC),
                                     "radio": ("How should cells get their simulation?", EC), "preview": ("Preview: what each cell will get", EC), "use": ("Use this in a movie", EXP)})
        if want("chrom"):
            goto_page(pg, "Chromatin")
            save_page(pg, "pg_chrom", {"intro": INTRO, "validation": ("Validation report", EXP), "lib": ("Library", EC), "sim": ("Simulation", EC), "player": ("Watch the chromosome move", EXP),
                                       "contact": ("Show contact map", EC), "two": ("Two loci: distance over time", EC), "locusA": ("Locus A position (kb)", EC), "sep": ("Separation to locus B (kb)", EC),
                                       "loops": ("Loop-extruding factors over time", EC), "zoom": ("Zoom: individual factors", EC)})
        if want("cfg"):
            goto_page(pg, "Config")
            click_text(pg, "loci · the two labelled loci")
            save_page(pg, "pg_cfg", {"intro": INTRO, "sections": ("Settings by section", EC), "acq": ("acquisition · frame count", EXP), "loci": ("loci · the two labelled loci", EXP),
                                     "whole": ("The whole file", EXP), "files": ("Files for this run", EC)})
        b.close()


main()
