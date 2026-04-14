import json
import subprocess
import sys


def test_package_import_does_not_eagerly_import_cli_modules():
    code = """
import json
import sys
import cdhit_reader

print(json.dumps({
    "has_read_cdhit": hasattr(cdhit_reader, "read_cdhit"),
    "cli_loaded": "cdhit_reader._cli" in sys.modules,
    "compare_loaded": "cdhit_reader._compare" in sys.modules,
    "test_loaded": "cdhit_reader._testit" in sys.modules,
}))
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=True,
    )

    payload = json.loads(result.stdout)
    assert payload == {
        "has_read_cdhit": True,
        "cli_loaded": False,
        "compare_loaded": False,
        "test_loaded": False,
    }


def test_lazy_cli_export_imports_on_demand():
    code = """
import json
import sys
import cdhit_reader

before = "cdhit_reader._cli" in sys.modules
_ = cdhit_reader.cli
after = "cdhit_reader._cli" in sys.modules

print(json.dumps({"before": before, "after": after}))
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=True,
    )

    payload = json.loads(result.stdout)
    assert payload == {"before": False, "after": True}
