"""Constants of the FCC-ee LCC_107 Z LOCO example (4D, cavities off, no radiation)."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results"

# scripts/00_fetch_lattice.sh clones acc-models-fcc-ee (LCC_107) here
LATTICE_DIR = ROOT / "data" / "acc-models-fcc-ee" / "lattices" / "z"
MADX_NAME = "fccee_z.madx"
XSUITE_JSON = LATTICE_DIR / "fccee_z.json"
XSUITE_LINE = "fccee_p_ring"
# scripts/00_make_model.py: lattice with corrector kicks deferred to MAD-X variables k_<name>
SEQUENCE_FILE = ROOT / "data" / "fccee_z_loco.madx"

# Arc quadrupoles only: the final-focus and straight-section quads are too sensitive
ARC_QUAD = r"^q[df][123]a(m)?\.\d+$"
QUAD_REL_DK1 = 0.5e-2  # relative k1L error, rms
TRIM_STEP = 0.5e-6  # corrector kick step, rad; every corrector is scanned at +-step
BPM_NOISE = 0.1e-6  # m, rms per BPM reading
