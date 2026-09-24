#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
check_models.py

Simple summary script for the final Erythroxylum novogranatense
leaf and root genome-scale metabolic models.

Reports:
1) Model dimensions
2) Bopt for Leaf and Root
3) Structural Jaccard similarity
4) Functional Jaccard similarity using pFBA
5) Production potential at >=90% Bopt for the 8 metabolites used in Figure 7

Expected default files in the current directory:
  model_leaf_E_novogranatense.xml
  model_root_E_novogranatense.xml
  medium_leaf_E_novogranatense.tsv
  medium_root_E_novogranatense.tsv

The script is safe to run from:
  - a terminal
  - Jupyter
  - Google Colab
  - `%run check_models.py`

COBRApy is installed automatically if it is missing.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path


# =============================================================================
# DEPENDENCY SETUP
# =============================================================================

def ensure_package(import_name, pip_name=None):
    """Install a missing package automatically."""
    try:
        __import__(import_name)
    except ImportError:
        package = pip_name or import_name
        print(f"[SETUP] Installing missing package: {package}")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "-q", package]
        )


ensure_package("cobra", "cobra>=0.29")

from cobra import Reaction
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba


# =============================================================================
# ANALYSIS SETTINGS
# =============================================================================

OBJECTIVE = "bio2_biomass"
ACTIVE_EPS = 1e-9
BIOMASS_FRACTION = 0.90

TARGETS = [
    ("Nicotinamide", "cpd00133_c0"),
    ("Gynesine", "cpd00737_c0"),
    ("S-Adenosyl-homocysteine", "cpd00019_c0"),
    ("Sinapaldehyde", "cpd03333_c0"),
    ("Niacin", "cpd00218_c0"),
    ("Cyanidin", "cpd05905_c0"),
    ("Kaempferol", "cpd05903_c0"),
    ("Naringenin", "cpd27666_c0"),
]


# =============================================================================
# HELPERS
# =============================================================================

def section(title):
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def first_present(fieldnames, candidates):
    lower = {str(x).strip().lower(): x for x in fieldnames}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def read_medium(path):
    """Read a medium TSV with reaction ID, lower bound and upper bound."""
    reaction_cols = [
        "reaction_id", "reaction", "rxn_id", "rxn",
        "exchange", "exchange_id", "id",
    ]
    lower_cols = [
        "lower_bound", "lowerbound", "lb", "lower",
        "min_flux", "minimum",
    ]
    upper_cols = [
        "upper_bound", "upperbound", "ub", "upper",
        "max_flux", "maximum",
    ]

    path = Path(path)

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        rid_col = first_present(reader.fieldnames or [], reaction_cols)
        lb_col = first_present(reader.fieldnames or [], lower_cols)
        ub_col = first_present(reader.fieldnames or [], upper_cols)

        if not all([rid_col, lb_col, ub_col]):
            raise ValueError(
                f"{path}: could not identify reaction/lower/upper bound columns."
            )

        rows = []
        for row in reader:
            rid = str(row[rid_col]).strip()
            if not rid or rid.startswith("#"):
                continue
            rows.append(
                (rid, float(row[lb_col]), float(row[ub_col]))
            )

    return rows


def apply_medium(model, medium_rows, label):
    """
    Close only uptake/source directions of existing boundary reactions,
    preserve secretion/drain directions, then apply the tissue medium.
    """
    m = model.copy()

    for rxn in m.boundary:
        if len(rxn.metabolites) != 1:
            continue

        _, coeff = next(iter(rxn.metabolites.items()))
        lb, ub = map(float, rxn.bounds)

        if coeff < 0:
            lb = max(lb, 0.0)
        elif coeff > 0:
            ub = min(ub, 0.0)

        rxn.bounds = (lb, ub)

    missing = []

    for rid, lb, ub in medium_rows:
        if rid not in m.reactions:
            missing.append(rid)
        else:
            m.reactions.get_by_id(rid).bounds = (lb, ub)

    if missing:
        raise KeyError(
            f"{label}: medium reactions absent from model: "
            + ", ".join(missing)
        )

    return m


def model_counts(model):
    return {
        "reactions": len(model.reactions),
        "metabolites": len(model.metabolites),
        "genes": len(model.genes),
        "exchanges": len(model.exchanges),
    }


def calculate_bopt(model, label):
    m = model.copy()

    if OBJECTIVE not in m.reactions:
        raise KeyError(f"{label}: objective {OBJECTIVE!r} not found.")

    m.objective = m.reactions.get_by_id(OBJECTIVE)
    m.objective_direction = "max"

    sol = m.optimize()

    if sol.status != "optimal":
        raise RuntimeError(
            f"{label}: biomass optimization failed: {sol.status}"
        )

    return float(sol.objective_value)


def internal_ids(model):
    boundary = {rxn.id for rxn in model.boundary}
    return {
        rxn.id
        for rxn in model.reactions
        if rxn.id not in boundary
    }


def stoich_dict(rxn):
    return {
        met.id: float(coeff)
        for met, coeff in rxn.metabolites.items()
        if abs(float(coeff)) > 1e-12
    }


def same_stoich(rxn_a, rxn_b, tol=1e-12):
    a = stoich_dict(rxn_a)
    b = stoich_dict(rxn_b)

    if set(a) != set(b):
        return False

    return all(
        abs(a[mid] - b[mid]) <= tol
        for mid in a
    )


def chemistry_keys(reaction_ids, tissue, mismatch_ids):
    """
    Biomass is excluded from Jaccard sets.
    Same reaction ID with different stoichiometry is treated as tissue-specific.
    """
    keys = set()

    for rid in reaction_ids:
        if rid == OBJECTIVE:
            continue

        if rid in mismatch_ids:
            keys.add(f"{tissue}::{rid}")
        else:
            keys.add(rid)

    return keys


def jaccard(set_a, set_b):
    intersection = set_a & set_b
    union = set_a | set_b

    value = (
        len(intersection) / len(union)
        if union
        else float("nan")
    )

    return len(intersection), len(union), value


# =============================================================================
# JACCARD
# =============================================================================

def structural_jaccard(leaf, root):
    leaf_internal = internal_ids(leaf)
    root_internal = internal_ids(root)

    shared_ids = leaf_internal & root_internal

    mismatch_ids = {
        rid
        for rid in shared_ids
        if not same_stoich(
            leaf.reactions.get_by_id(rid),
            root.reactions.get_by_id(rid),
        )
    }

    leaf_keys = chemistry_keys(
        leaf_internal, "LEAF", mismatch_ids
    )
    root_keys = chemistry_keys(
        root_internal, "ROOT", mismatch_ids
    )

    intersection, union, value = jaccard(
        leaf_keys, root_keys
    )

    return {
        "leaf_internal": len(leaf_keys),
        "root_internal": len(root_keys),
        "shared": intersection,
        "union": union,
        "jaccard": value,
        "stoich_mismatch_ids": sorted(mismatch_ids),
    }


def functional_jaccard(leaf, root, mismatch_ids):
    """
    Active internal reactions from pFBA solutions at the biomass optimum.
    """
    leaf_m = leaf.copy()
    root_m = root.copy()

    leaf_m.objective = leaf_m.reactions.get_by_id(OBJECTIVE)
    root_m.objective = root_m.reactions.get_by_id(OBJECTIVE)

    leaf_m.objective_direction = "max"
    root_m.objective_direction = "max"

    leaf_sol = pfba(
        leaf_m,
        fraction_of_optimum=1.0,
    )
    root_sol = pfba(
        root_m,
        fraction_of_optimum=1.0,
    )

    leaf_internal = internal_ids(leaf_m)
    root_internal = internal_ids(root_m)

    leaf_active = {
        rid
        for rid in leaf_internal
        if rid != OBJECTIVE
        and abs(float(leaf_sol.fluxes[rid])) > ACTIVE_EPS
    }

    root_active = {
        rid
        for rid in root_internal
        if rid != OBJECTIVE
        and abs(float(root_sol.fluxes[rid])) > ACTIVE_EPS
    }

    leaf_keys = chemistry_keys(
        leaf_active, "LEAF", mismatch_ids
    )
    root_keys = chemistry_keys(
        root_active, "ROOT", mismatch_ids
    )

    intersection, union, value = jaccard(
        leaf_keys, root_keys
    )

    return {
        "leaf_active": len(leaf_keys),
        "root_active": len(root_keys),
        "shared": intersection,
        "union": union,
        "jaccard": value,
    }


# =============================================================================
# PRODUCTION POTENTIAL
# =============================================================================

def production_potential(model, metabolite_id, bopt):
    """
    Maximize an in-memory temporary demand reaction while requiring
    biomass >= 90% of the tissue-specific Bopt.

    The SBML file is never modified.
    """
    m = model.copy()

    if metabolite_id not in m.metabolites:
        return float("nan")

    biomass = m.reactions.get_by_id(OBJECTIVE)
    biomass.lower_bound = BIOMASS_FRACTION * bopt

    met = m.metabolites.get_by_id(metabolite_id)

    dm = Reaction(f"TMP_PROD_{metabolite_id}")
    dm.name = f"Temporary demand for {met.name}"
    dm.bounds = (0.0, 1000.0)
    dm.add_metabolites({met: -1.0})

    m.add_reactions([dm])

    m.objective = dm
    m.objective_direction = "max"

    sol = m.optimize()

    if sol.status != "optimal":
        return float("nan")

    return float(sol.objective_value)


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_analysis(
    leaf_model="model_leaf_E_novogranatense.xml",
    root_model="model_root_E_novogranatense.xml",
    leaf_medium="medium_leaf_E_novogranatense.tsv",
    root_medium="medium_root_E_novogranatense.tsv",
):
    paths = {
        "leaf_model": Path(leaf_model),
        "root_model": Path(root_model),
        "leaf_medium": Path(leaf_medium),
        "root_medium": Path(root_medium),
    }

    missing = [
        str(path)
        for path in paths.values()
        if not path.exists()
    ]

    if missing:
        raise FileNotFoundError(
            "Required file(s) not found:\n  - "
            + "\n  - ".join(missing)
            + "\n\nPlace the two XML models and two medium TSV files "
              "in the current working directory, or pass their paths explicitly."
        )

    section("1. MODEL DIMENSIONS")

    leaf_raw = read_sbml_model(
        str(paths["leaf_model"])
    )
    root_raw = read_sbml_model(
        str(paths["root_model"])
    )

    leaf_counts = model_counts(leaf_raw)
    root_counts = model_counts(root_raw)

    print(
        f"Leaf : reactions={leaf_counts['reactions']} | "
        f"metabolites={leaf_counts['metabolites']} | "
        f"genes={leaf_counts['genes']} | "
        f"exchanges={leaf_counts['exchanges']}"
    )

    print(
        f"Root : reactions={root_counts['reactions']} | "
        f"metabolites={root_counts['metabolites']} | "
        f"genes={root_counts['genes']} | "
        f"exchanges={root_counts['exchanges']}"
    )

    leaf = apply_medium(
        leaf_raw,
        read_medium(paths["leaf_medium"]),
        "LEAF",
    )

    root = apply_medium(
        root_raw,
        read_medium(paths["root_medium"]),
        "ROOT",
    )

    section("2. BIOMASS OPTIMUM (Bopt)")

    leaf_bopt = calculate_bopt(
        leaf, "LEAF"
    )
    root_bopt = calculate_bopt(
        root, "ROOT"
    )

    print(f"Leaf Bopt = {leaf_bopt:.12f}")
    print(f"Root Bopt = {root_bopt:.12f}")

    section("3. STRUCTURAL JACCARD")

    structural = structural_jaccard(
        leaf_raw, root_raw
    )

    print(
        f"Leaf internal reactions = "
        f"{structural['leaf_internal']}"
    )
    print(
        f"Root internal reactions = "
        f"{structural['root_internal']}"
    )
    print(
        f"Shared reactions        = "
        f"{structural['shared']}"
    )
    print(
        f"Union reactions         = "
        f"{structural['union']}"
    )
    print(
        f"Structural Jaccard      = "
        f"{structural['jaccard']:.6f}"
    )

    section("4. FUNCTIONAL JACCARD (pFBA)")

    functional = functional_jaccard(
        leaf,
        root,
        set(structural["stoich_mismatch_ids"]),
    )

    print(
        f"Leaf active reactions = "
        f"{functional['leaf_active']}"
    )
    print(
        f"Root active reactions = "
        f"{functional['root_active']}"
    )
    print(
        f"Shared active         = "
        f"{functional['shared']}"
    )
    print(
        f"Union active          = "
        f"{functional['union']}"
    )
    print(
        f"Functional Jaccard    = "
        f"{functional['jaccard']:.6f}"
    )

    section("5. PRODUCTION POTENTIAL AT >= 90% Bopt")

    print(
        f"{'Metabolite':30s} "
        f"{'Leaf':>12s} "
        f"{'Root':>12s}"
    )
    print("-" * 58)

    production = []

    for name, metabolite_id in TARGETS:
        leaf_value = production_potential(
            leaf,
            metabolite_id,
            leaf_bopt,
        )

        root_value = production_potential(
            root,
            metabolite_id,
            root_bopt,
        )

        production.append({
            "metabolite": name,
            "metabolite_id": metabolite_id,
            "leaf": leaf_value,
            "root": root_value,
        })

        leaf_text = (
            "NA"
            if math.isnan(leaf_value)
            else f"{leaf_value:.6f}"
        )

        root_text = (
            "NA"
            if math.isnan(root_value)
            else f"{root_value:.6f}"
        )

        print(
            f"{name:30s} "
            f"{leaf_text:>12s} "
            f"{root_text:>12s}"
        )

    section("6. SUMMARY")

    print(f"Leaf Bopt          = {leaf_bopt:.6f}")
    print(f"Root Bopt          = {root_bopt:.6f}")
    print(
        f"Structural Jaccard = "
        f"{structural['jaccard']:.6f}"
    )
    print(
        f"Functional Jaccard = "
        f"{functional['jaccard']:.6f}"
    )

    return {
        "leaf_counts": leaf_counts,
        "root_counts": root_counts,
        "leaf_bopt": leaf_bopt,
        "root_bopt": root_bopt,
        "structural_jaccard": structural,
        "functional_jaccard": functional,
        "production_potential": production,
    }


# =============================================================================
# COMMAND-LINE ENTRY POINT
# =============================================================================

def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Check the main results of the final "
            "E. novogranatense leaf/root GEMs."
        )
    )

    parser.add_argument(
        "--leaf-model",
        default="model_leaf_E_novogranatense.xml",
    )
    parser.add_argument(
        "--root-model",
        default="model_root_E_novogranatense.xml",
    )
    parser.add_argument(
        "--leaf-medium",
        default="medium_leaf_E_novogranatense.tsv",
    )
    parser.add_argument(
        "--root-medium",
        default="medium_root_E_novogranatense.tsv",
    )

    return parser


def main(argv=None):
    """
    Colab/Jupyter-safe entry point.

    parse_known_args() deliberately ignores notebook kernel arguments
    such as:
        -f /root/.local/share/jupyter/runtime/kernel-....json
    """
    parser = build_parser()

    args, unknown = parser.parse_known_args(argv)

    # Jupyter/Colab injects arguments unrelated to this script.
    # They are intentionally ignored.
    if unknown and not any(
        str(x).startswith("-f")
        for x in unknown
    ):
        print(
            "[INFO] Ignoring unrecognized environment arguments: "
            + " ".join(unknown)
        )

    return run_analysis(
        leaf_model=args.leaf_model,
        root_model=args.root_model,
        leaf_medium=args.leaf_medium,
        root_medium=args.root_medium,
    )


if __name__ == "__main__":
    main()
