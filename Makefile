test: lint
	pytest

setup:
	python -m pip install --upgrade pip
	python -m pip install flake8 pytest
	if [ -f requirements.txt ]; then pip install -r requirements.txt; fi

install:
	pip install .

lint:
	flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
	flake8 . --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics

build: test
	python setup.py sdist bdist_wheel

.PHONY: install lint test build
