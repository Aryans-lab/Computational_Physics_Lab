# Computational Physics & Scientific Computing — task runner
#
# Targets:
#   make bundle  — copy the canonical mylib.py into every assignment folder
#   make run     — bundle + run every question script (regenerates outputs)
#   make test    — bundle + run the full pytest suite (unit + golden outputs)
#   make lint    — ruff over the whole repository
#   make all     — lint + test + run (what CI does)

PYTHON ?= python3

.PHONY: bundle run test lint all clean

bundle:
	$(PYTHON) scripts/bundle.py

run: bundle
	$(PYTHON) scripts/run_all.py

test: bundle
	$(PYTHON) -m pytest tests -q

lint:
	$(PYTHON) -m ruff check .

all: lint test run

clean:
	rm -f assignments/*/mylib.py
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
	rm -rf .pytest_cache .ruff_cache
