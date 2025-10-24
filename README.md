# SPIR: Structure Prediction Input Representation

Convert between AF3, AF3 Server, Boltz, Chai, and Protenix input formats with an intermediate Python IR that preserves glycans and covalent bonds.

Quick start:

```bash
pip install -e .
spir convert --src path/to/input.json --to af3 --out out.json
spir validate --src path/to/input.json
spir validate --src path/to/input.json --format af3 --explain
```

Python API:

```python
from spir import validate

# Boolean result
is_valid = validate("path/to/input.json")

# Force a format and get explanation string (empty string means valid)
msg = validate("path/to/input.json", format="af3", explain=True)
print("valid?", msg == "", "issues:", msg or "(none)")
```

## Glycans and bonds

- AF3‑Server inputs with `proteinChain.glycans[].residues` are parsed into a bonded glycan graph. We preserve intra‑glycan bonds (parent O[n] → child C1) and add a protein↔glycan anchor at the specified residue. Conversions to AF3, Boltz, Chai and Protenix include these bonds.
- Boltz outputs normalize multi‑CCD glycans into one ligand per sugar and express all connectivity via YAML `constraints: bond:` pairs.
