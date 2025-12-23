# SPIR

![SPIR Logo](./src/spir/img/logo.png)

SPIR (**S**tructure **P**rediction **I**ntermediate **R**epresentation) exists to make it practical to compare and iterate across multiple structure prediction models without constantly rewriting inputs by hand. Different predictors (AlphaFold3 Server/non-Server, Chai-1, Boltz-2, Protenix) can yield meaningfully different structures, confidence metrics, and binding/interface hypotheses on the same biological system; being able to run the same job across models is essential for validating conclusions, spotting model-specific artifacts, and choosing the best tool for a given target or constraint set.  
  
In practice, such comparisons are hindered by the format fragmentation across model ecosystems, especially for glycans, where representations range from compact tree strings with implicit chemistry (e.g., AF3 Server) to fully specified multi-component ligands plus explicit bonded atom pairs. Reliably converting between formats requires more than renaming fields: it requires an intermediate graph-like representation that preserves residue identity, connectivity, attachment sites, and (when needed) explicit linkage atoms/positions, while also handling cases where a target format omits or infers chemistry. SPIR provides that IR together with model-specific converters so scientific questions, not input wrangling, drive the workflow.

## Installation

The easiest way to install SPIR is to use `pip`:

```bash
pip install spir
```

If you want to build from source, you can clone the repository and run:

```bash
git clone https://github.com/foldir/spir.git
cd spir
pip install -e .
```

## Usage

SPIR provides a CLI for converting between different structure prediction model inputs.

```bash
spir convert --from alphafold3server input.json --to alphafold3 output
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
