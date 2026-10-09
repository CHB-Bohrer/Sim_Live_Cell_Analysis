"""Write docs/IMAGE_ERRORS.md from the table in stage3_microscopy/errors.py (so the doc never disagrees with the code).

    scripts\\run.cmd python scripts\\make_image_error_docs.py
"""
from pathlib import Path

from simlive.stage3_microscopy.errors import GROUPS, ORDER, SOURCES

HEAD = """# Sources of image error (stage 3)

Generated from `src/simlive/stage3_microscopy/errors.py` by `scripts/make_image_error_docs.py`; do not edit by hand.
All numbers below are **placeholders** until measured on the real microscope. Switch them in the config block `imaging_errors:`
(each source has `enabled: true|false` plus its parameters), or in the dashboard sidebar ("Image errors"). With
`imaging_errors.enabled: false` only the original model remains: blur, bleaching, constant background, quantum efficiency,
shot noise and read noise.

How a frame is made, in physical order: objects at their true positions -> stage drift and chromatic shift move them -> optical blur
(diffraction, defocus, exposure motion, out-of-focus haze) -> photobleaching -> uneven illumination and flicker -> colour cross-talk ->
background (constant, autofluorescence) -> detector (quantum efficiency, shot noise, dark current) -> camera electronics (hot pixels,
cosmic rays, pixel gain, read noise, row noise, offset pattern, digitisation and saturation).

Ground truth saved with every run (`stage3_microscopy/`): `labels.tif` (true masks, moved by the drift), `imaging_truth.csv` (per frame: drift,
focus error, blur per channel, light level per channel), `cells_image.csv` and `loci_truth_image.csv` (truth in image coordinates; the latter
also has where each colour's spot is really drawn), `imaging_errors.json` (which sources were on, with all parameters).
The dashboard's **Images** tab re-renders any frame from this truth, with any subset of sources, and checks it equals the saved image.

Reproducibility: every source draws from its own random stream seeded by (run, frame, source, channel), so switching one source never changes
another source's noise, and any frame can be rendered alone.
"""


def main() -> None:
    out = [HEAD]
    for g in GROUPS:
        out.append(f"\n## {g}\n")
        for sid in ORDER:
            s = SOURCES[sid]
            if s["group"] != g:
                continue
            out.append(f"### {s['title']}  (`{sid}`{', original model' if s['baseline'] else ''})\n")
            out.append(f"{s['what']}\n\n**Effect:** {s['effect']}\n")
            if s["params"]:
                out.append("| parameter | default | unit | meaning |\n|---|---|---|---|")
                for k, (d, u, t) in s["params"].items():
                    out.append(f"| `{k}` | `{d}` | {u} | {t} |")
                out.append("")
    p = Path(__file__).resolve().parents[1] / "docs" / "IMAGE_ERRORS.md"
    p.write_text("\n".join(out), encoding="utf-8")
    print("wrote", p)


if __name__ == "__main__":
    main()
