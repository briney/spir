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

## Supported Models

SPIR supports the following structure prediction models:

| model | dialect |
|-------|---------|
| AlphaFold3 Server | `alphafold3server` |
| AlphaFold3 (non-Server) | `alphafold3` |
| Boltz-2 | `boltz2` |
| Chai-1 | `chai1` |
| Protenix | `protenix` |


## Custom MSAs

AlphaFold3 (non-Server) and Boltz-2 support custom MSA paths as part of their respective input formats. We anticipate many users will want to convert from the AlphaFold3 Server format to one of these dialects, since the AlphaFold3 Server format is particularly user-friendly with respect to glycans. While the official [AlphaFold3 Server input format](https://github.com/google-deepmind/alphafold/blob/main/server/README.md) does not support custom MSA paths, SPIR allows users to supply custom MSA using an unofficial `msa_path` field for any `proteinChain`, `dnaSequence`, or `rnaSequence` in an AlphaFold3 Server input, like so:

```json
{
  "name": "msa_test",
  "modelSeeds": [42],
  "sequences": [
    {
      "proteinChain": {
        "id": "A",
        "sequence": "MVLSPADKTN",
        "msa_path": "/path/to/msa/protein_a.a3m"
      }
    }
  ]
}
```

SPIR will then add the custom MSA path to the appropriate format for the target output dialect. For example, if you convert to AlphaFold3 (non-Server), the supplied MSA path will be added to the `unpairedMsaPath` field. For Boltz-2, the supplied MSA path will be added to the `msa` field.

> [!NOTE]
> The unofficial `msa_path` field in AlphaFold3 Server is only supported for input files. If an AlphaFold3 (non-Server) or Boltz-2 input file containing an MSA path is converted to AlphaFold3 Server format, the MSA path will be ignored.


## License

SPIR is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
