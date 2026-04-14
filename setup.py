from setuptools import setup

if __name__ == "__main__":
    console_scripts = [
        "cdhit-parser = cdhit_reader._cli:cli",
        "cdhit-compare = cdhit_reader._compare:compare",
        "cdhit-compare-otus = cdhit_reader._compare_otus:compare_otus",
    ]
    setup(entry_points=dict(console_scripts=console_scripts))
