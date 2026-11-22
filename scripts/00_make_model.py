# /// script
# requires-python = ">=3.11"
# dependencies = ["fcc-loco"]
# [tool.uv.sources]
# fcc-loco = { path = "..", editable = true }
# ///
"""Write data/fccee_z_loco.madx: the Z lattice with deferred corrector kicks, cavities off."""

from pymadng_utils.model_creator.create_fcc_model import create_fcc_model

from fcc_loco.config import LATTICE_DIR, MADX_NAME, SEQUENCE_FILE

create_fcc_model(LATTICE_DIR, SEQUENCE_FILE, madx_name=MADX_NAME)
text = SEQUENCE_FILE.read_text()
assert "rf_on = 1.0;" in text
SEQUENCE_FILE.write_text(text.replace("rf_on = 1.0;", "rf_on = 0.0;"))
