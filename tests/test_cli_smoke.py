from __future__ import annotations

import json
from pathlib import Path

from validation.run_pyexocross import main


def test_end_to_end_pyexocross_smoke(tmp_path: Path) -> None:
    fixtures = Path(__file__).parent / "fixtures"
    output = tmp_path / "report"
    result = main([
        "--ion", "C-III",
        "--charge", "2",
        "--states", str(fixtures / "C_III__Kurucz.states"),
        "--trans", str(fixtures / "C_III__Kurucz.trans"),
        "--pf", str(fixtures / "C_III__Kurucz.pf"),
        "--temperatures", "100", "1000", "10000",
        "--range", "1", "2500",
        "--output", str(output),
    ])
    assert result == 0
    summary = json.loads((output / "validation_summary.json").read_text())
    assert summary["overall_status"] in {"PASS", "PASS_WITH_WARNINGS"}
    assert summary["pyexocross"]["version"] == "1.1.9"
    assert summary["spectrum"]["generated"]
    for name in [
        "pf_comparison.csv",
        "lifetime_comparison.csv",
        "strongest_lines.csv",
        "partition_function.png",
        "stick_spectrum_wavenumber.png",
        "validation_report.md",
    ]:
        assert (output / name).exists()
