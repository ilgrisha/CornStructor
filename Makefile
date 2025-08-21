# File: Makefile
# Version: v0.1.3

.PHONY: venv install test lint fmt bench

VENV?=backend/.venv
PY?=$(VENV)/bin/python
PIP?=$(VENV)/bin/pip
PYTEST?=$(VENV)/bin/pytest
RUFF?=$(VENV)/bin/ruff
BLACK?=$(VENV)/bin/black

venv:
	python3 -m venv $(VENV)

install: venv
	$(PIP) install -U pip
	$(PIP) install -r backend/requirements.txt
	$(PIP) install black ruff mypy

test:
	PYTHONPATH=$(PWD) $(PYTEST) -q

lint:
	$(RUFF) backend
	mypy backend

fmt:
	$(BLACK) backend

bench:
	PYTHONPATH=$(PWD) $(PY) backend/app/bench/bench_fitness.py
