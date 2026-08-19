# ---------------------------------------------------------------------------
# Loblaw Bio cell-count assignment
#   make setup      create the virtualenv and install dependencies
#   make pipeline   build the database and run every analysis
#   make dashboard  launch the interactive dashboard
# ---------------------------------------------------------------------------
PYTHON ?= python3
VENV   := .venv
BIN    := $(VENV)/bin
PY     := $(BIN)/python
DB     := db/cell_counts.db
CSV    ?= data/cell-count.csv

.DEFAULT_GOAL := help
.PHONY: help setup pipeline dashboard db analysis test clean distclean

help:
	@echo "make setup      - create $(VENV) and install requirements"
	@echo "make pipeline   - load $(CSV) into $(DB) and run parts 2-4"
	@echo "make dashboard  - launch the Streamlit dashboard"
	@echo "make test       - run the test suite"
	@echo "make clean      - remove generated database and outputs"

$(BIN)/activate:
	$(PYTHON) -m venv $(VENV)
	$(BIN)/python -m pip install --upgrade pip
	$(BIN)/python -m pip install -r requirements.txt

setup: $(BIN)/activate
	@echo "Environment ready. Next: make pipeline"

db: setup
	$(PY) load_data.py --csv $(CSV) --db $(DB)

analysis:
	$(PY) -m analysis.run_all

pipeline: db
	$(PY) -m analysis.run_all

dashboard: setup
	@test -f $(DB) || $(PY) load_data.py --csv $(CSV) --db $(DB)
	$(BIN)/streamlit run dashboard/app.py

test: setup
	$(PY) -m pytest -q tests

clean:
	rm -rf $(DB) outputs/*.csv outputs/figures/*.png
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +

distclean: clean
	rm -rf $(VENV)
