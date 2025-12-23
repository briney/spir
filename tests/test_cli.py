"""Tests for SPIR CLI commands."""

import json

from typer.testing import CliRunner

from spir.cli import app

runner = CliRunner()


class TestConvertCommand:
    def test_convert_help(self):
        result = runner.invoke(app, ["convert", "--help"], color=False)
        assert result.exit_code == 0
        assert "--from" in result.output
        assert "--to" in result.output

    def test_convert_basic(self, tmp_path):
        payload = [
            {
                "name": "test",
                "modelSeeds": [],
                "sequences": [{"proteinChain": {"sequence": "MLKK", "count": 1}}],
            }
        ]
        in_path = tmp_path / "input.json"
        in_path.write_text(json.dumps(payload))
        out_prefix = tmp_path / "output"

        result = runner.invoke(
            app,
            [
                "convert",
                str(in_path),
                "--from",
                "alphafoldserver",
                str(out_prefix),
                "--to",
                "alphafold3",
            ],
            color=False,
        )
        assert result.exit_code == 0
        assert (tmp_path / "output.json").exists()


class TestValidateCommand:
    def test_validate_help(self):
        result = runner.invoke(app, ["validate", "--help"], color=False)
        assert result.exit_code == 0
        assert "--dialect" in result.output

    def test_validate_valid_file(self, tmp_path):
        payload = {
            "name": "test",
            "modelSeeds": [1],
            "sequences": [{"protein": {"id": "A", "sequence": "MLKK"}}],
            "dialect": "alphafold3",
            "version": 4,
        }
        path = tmp_path / "valid.json"
        path.write_text(json.dumps(payload))

        result = runner.invoke(
            app, ["validate", str(path), "--dialect", "alphafold3"], color=False
        )
        assert result.exit_code == 0
        assert "passed" in result.output.lower()

    def test_validate_invalid_file(self, tmp_path):
        payload = {
            "name": "test",
            # Missing modelSeeds
            "sequences": [{"protein": {"id": "A", "sequence": "MLKK"}}],
            "dialect": "alphafold3",
            "version": 4,
        }
        path = tmp_path / "invalid.json"
        path.write_text(json.dumps(payload))

        result = runner.invoke(
            app, ["validate", str(path), "--dialect", "alphafold3"], color=False
        )
        assert result.exit_code == 1
        assert "failed" in result.output.lower() or "ERROR" in result.output

    def test_validate_nonexistent_file(self, tmp_path):
        path = tmp_path / "nonexistent.json"

        result = runner.invoke(
            app, ["validate", str(path), "--dialect", "alphafold3"], color=False
        )
        assert result.exit_code == 1

    def test_validate_with_short_option(self, tmp_path):
        payload = {
            "name": "test",
            "modelSeeds": [1],
            "sequences": [{"protein": {"id": "A", "sequence": "MLKK"}}],
            "dialect": "alphafold3",
            "version": 4,
        }
        path = tmp_path / "valid.json"
        path.write_text(json.dumps(payload))

        # Using -d short option
        result = runner.invoke(
            app, ["validate", str(path), "-d", "alphafold3"], color=False
        )
        assert result.exit_code == 0


class TestMainHelp:
    def test_main_help(self):
        result = runner.invoke(app, ["--help"], color=False)
        assert result.exit_code == 0
        assert "convert" in result.output.lower()
        assert "validate" in result.output.lower()
