![](https://img.shields.io/pypi/v/spir.svg?colorB=blue)
[![tests](https://github.com/briney/spir/actions/workflows/pytest.yaml/badge.svg)](https://github.com/briney/spir/actions/workflows/pytest.yaml)
![](https://img.shields.io/badge/license-MIT-blue.svg)

# SPIR

<!-- ![SPIR Logo](./src/spir/img/logo.png) -->

SPIR (**S**tructure **P**rediction **I**ntermediate **R**epresentation) exists to make it practical to compare and iterate across multiple structure prediction models without constantly rewriting inputs by hand. Different predictors (AlphaFold3 Server/non-Server, Chai-1, Boltz-2, Protenix) can yield meaningfully different structures, confidence metrics, and binding/interface hypotheses on the same biological system; being able to run the *same* job across models is essential for validating conclusions, spotting model-specific artifacts, and choosing the best tool for a given target or constraint set.  
  
In practice, such comparisons are hindered by the format fragmentation across different models, especially for glycans, where representations range from compact tree strings with implicit chemistry (e.g., AF3 Server) to fully specified multi-component ligands with explicitly specified bonded atom pairs. Reliably converting between formats requires more than renaming fields: it necessitates an intermediate graph-like representation that preserves residue identity, connectivity, attachment sites, and (when needed) explicit linkage atoms/positions, while also handling cases where a target format omits or infers chemistry. SPIR provides that IR together with model-specific converters so scientific questions, not input wrangling, drive the workflow.

## Installation

The easiest way to install SPIR is to use `pip`:

```bash
pip install spir
```

If you want to build from source, you can clone the repository and run:

```bash
git clone https://github.com/briney/spir
cd spir
pip install -e .
```

## Usage

SPIR provides a CLI for converting between different structure prediction model inputs.

```bash
spir convert --from DIALECT INPUT_FILE --to DIALECT OUTPUT_PREFIX
```

For example, to convert an AlphaFold3 Server input to an AlphaFold3 (non-Server) output, you can run:

```bash
spir convert --from alphafold3server path/to/input.json --to alphafold3 path/to/output
```

> [!NOTE]
> The `output_prefix` should only contain the prefix for the output files (no extension). The correct extension will be added automatically.

If your input is Chai-1 formatted and includes restraints, you can specify the restraints file with the `--restraints` option:

```bash
spir convert --from chai1 input.fasta --to protenix output --restraints restraints.csv
```

## Supported Formats

SPIR supports the following formats:

- AlphaFold3 Server: `alphafold3server`
- AlphaFold3 (non-Server): `alphafold3`
- Boltz-2: `boltz2`
- Chai-1: `chai1`
- Protenix: `protenix`

## License

SPIR is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
