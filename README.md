# Spir

Intermediate representation and converters for protein folding model inputs.

## Quick start

```bash
pip install -e .
spir convert --from alphafoldserver input.json --to alphafold3 output
```

## Layout

- `src/spir/ir`: IR models and normalization
- `src/spir/ir/glycans`: glycan parsing/rendering/resolution
- `src/spir/dialects`: format adapters
- `src/spir/convert.py`: conversion orchestrator
