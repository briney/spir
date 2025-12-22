# FoldIR

Intermediate representation and converters for protein folding model inputs.

## Quick start

```bash
pip install -e .
foldir convert --from alphafoldserver input.json --to alphafold3 output
```

## Layout

- `src/foldir/ir`: IR models and normalization
- `src/foldir/ir/glycans`: glycan parsing/rendering/resolution
- `src/foldir/dialects`: format adapters
- `src/foldir/convert.py`: conversion orchestrator
