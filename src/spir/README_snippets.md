# Quick usage

High-level Python API:

```python
from spir import read, convert_af3_server_glycan

ir = read("input.yaml")         # auto-detects format
ir.write_alphafold3("out/")    # writes out/input.af3.json
ir.write_boltz("out/")         # writes out/input.boltz.yaml
ir.write_chai("out/")          # writes out/input.fasta and maybe .restraints.csv
```

Inspect an AF3 Server glycan across formats via API:

```python
# Minimal snippet for AF3 (AlphaFold3) format
print(convert_af3_server_glycan("NAG(NAG)(BMA)", to="af3"))

# Full sample Protenix JSON
print(convert_af3_server_glycan("NAG(NAG)(BMA)", to="protenix", detail="full"))
```

Convert AF3 Server -> AF3:

```bash
spir convert --src job_server.json --to af3 --out job_af3.json
```

Convert Chai FASTA+CSV -> Protenix:

```bash
spir convert --src input.fasta --src bonds.restraints.csv --to protenix --out out.protenix.json
```
