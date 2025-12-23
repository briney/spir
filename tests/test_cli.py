"""Tests for SPIR CLI commands."""

import json
import re

from typer.testing import CliRunner

from spir.cli import app

runner = CliRunner()


def strip_ansi(text: str) -> str:
    """Remove ANSI escape codes from text."""
    ansi_pattern = re.compile(r'\x1b\[[0-9;]*m')
    return ansi_pattern.sub('', text)


class TestConvertCommand:
    def test_convert_help(self):
        result = runner.invoke(app, ["convert", "--help"])
        assert result.exit_code == 0
        output = strip_ansi(result.output)
        assert "--from" in output
        assert "--to" in output

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
        )
        assert result.exit_code == 0
        assert (tmp_path / "output.json").exists()


class TestValidateCommand:
    def test_validate_help(self):
        result = runner.invoke(app, ["validate", "--help"])
        assert result.exit_code == 0
        output = strip_ansi(result.output)
        assert "--dialect" in output

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
            app, ["validate", str(path), "--dialect", "alphafold3"]
        )
        assert result.exit_code == 0
        output = strip_ansi(result.output)
        assert "passed" in output.lower()

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
            app, ["validate", str(path), "--dialect", "alphafold3"]
        )
        assert result.exit_code == 1
        output = strip_ansi(result.output)
        assert "failed" in output.lower() or "ERROR" in output

    def test_validate_nonexistent_file(self, tmp_path):
        path = tmp_path / "nonexistent.json"

        result = runner.invoke(
            app, ["validate", str(path), "--dialect", "alphafold3"]
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
            app, ["validate", str(path), "-d", "alphafold3"]
        )
        assert result.exit_code == 0


class TestMainHelp:
    def test_main_help(self):
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        output = strip_ansi(result.output)
        assert "convert" in output.lower()
        assert "validate" in output.lower()
