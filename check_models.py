#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import numpy as np
import pandas as pd

from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba
from cobra.sampling import ACHRSampler


EPS = 1e-9
BIOMASS_ID = "bio2_biomass"
GROWTH_FRACTION = 0.90
N_SAMPLES = 500
THINNING = 100
SEED = 123

LEAF_MODEL_FILE = "model_leaf_E_novogranatense.xml"
ROOT_MODEL_FILE = "model_root_E_novogranatense.xml"
LEAF_MEDIUM_FILE = "medium_leaf_E_novogranatense.tsv"
ROOT_MEDIUM_FILE = "medium_root_E_novogranatense.tsv"

TARGET_REACTIONS = ["rxn02402_c0", "rxn00126"]

TARGET_METABOLITES = [
    ("Gynesine", "cpd00737_c0"),
    ("Nicotinamide", "cpd00133_c0"),
]


def first_existing(df, candidates):
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


def load_medium_table(path):
    df = pd.read_csv(path, sep="\t")

    rxn_col = first_existing(df, ["reaction_id", "reaction", "rxn_id", "rxn", "exchange", "id"])
    lb_col = first_existing(df, ["lower_bound", "lb", "min_flux", "lower", "uptake_lb"])
    ub_col = first_existing(df, ["upper_bound", "ub", "max_flux", "upper"])
    uptake_col = first_existing(df, ["uptake", "max_uptake", "import", "value", "flux"])

    if rxn_col is None:
        raise RuntimeError(f"Could not identify reaction ID column in {path}")

    parsed = []
    for _, row in df.iterrows():
        rxn_id = str(row[rxn_col]).strip()
        if not rxn_id or rxn_id.lower() == "nan":
            continue

        if lb_col is not None or ub_col is not None:
            lb = float(row[lb_col]) if lb_col is not None and pd.notna(row[lb_col]) else None
            ub = float(row[ub_col]) if ub_col is not None and pd.notna(row[ub_col]) else None
        elif uptake_col is not None and pd.notna(row[uptake_col]):
            uptake = float(row[uptake_col])
            lb = -abs(uptake)
            ub = None
        else:
            raise RuntimeError(f"Could not interpret medium bounds in {path}")

        parsed.append((rxn_id, lb, ub))
    return parsed


def close_all_imports(model):
    for rxn in model.exchanges:
        if rxn.lower_bound < 0:
            rxn.lower_bound = 0.0


def apply_medium(model, medium_path):
    parsed = load_medium_table(medium_path)
    close_all_imports(model)
    for rxn_id, lb, ub in parsed:
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            if lb is not None:
                rxn.lower_bound = lb
            if ub is not None:
                rxn.upper_bound = ub


def is_transport_like(rxn):
    if rxn.boundary:
        return False

    comps = set()
    for met in rxn.metabolites:
        comp = getattr(met, "compartment", None)
        if comp:
            comps.add(comp)
        elif "_" in met.id:
            comps.add(met.id.rsplit("_", 1)[1])

    rid = rxn.id.lower()
    return (len(comps) > 1) or rid.startswith("t_") or rid.startswith("tx_") or "_to_" in rid


def summarize_model(label, model):
    return {
        "label": label,
        "n_reactions": len(model.reactions),
        "n_metabolites": len(model.metabolites),
        "n_genes": len(model.genes),
    }


def optimize_biomass(model, biomass_id=BIOMASS_ID):
    model.objective = biomass_id
    sol = model.optimize()
    if sol.status != "optimal":
        raise RuntimeError(f"Biomass optimization failed: {sol.status}")
    return float(sol.objective_value), sol


def active_reaction_counts(model, biomass_id=BIOMASS_ID):
    model.objective = biomass_id
    psol = pfba(model)

    active_all = {rid for rid, v in psol.fluxes.items() if abs(v) > EPS}
    active_internal = {
        rxn.id for rxn in model.reactions
        if (rxn.id in active_all) and (not rxn.boundary) and (not is_transport_like(rxn))
    }

    return {
        "n_active_all_pfba": len(active_all),
        "n_active_internal_pfba": len(active_internal),
        "pfba_solution": psol,
    }


def prepare_sampling_model(model_file, medium_file, biomass_id=BIOMASS_ID, growth_fraction=GROWTH_FRACTION):
    model = read_sbml_model(model_file)
    apply_medium(model, medium_file)
    model.objective = biomass_id

    opt_sol = model.optimize()
    if opt_sol.status != "optimal":
        raise RuntimeError(f"Initial optimization failed for {model_file}: {opt_sol.status}")

    optimum = float(opt_sol.objective_value)
    biomass_rxn = model.reactions.get_by_id(biomass_id)
    biomass_rxn.lower_bound = growth_fraction * optimum

    return model, optimum


def sample_shared_fluxes(leaf_model, root_model, n_samples=N_SAMPLES, thinning=THINNING, seed=SEED):
    leaf_rxns = set(r.id for r in leaf_model.reactions if not r.boundary)
    root_rxns = set(r.id for r in root_model.reactions if not r.boundary)
    common = sorted(leaf_rxns & root_rxns)

    leaf_sampler = ACHRSampler(leaf_model, thinning=thinning, seed=seed)
    root_sampler = ACHRSampler(root_model, thinning=thinning, seed=seed + 1)

    leaf_df = leaf_sampler.sample(n_samples)
    root_df = root_sampler.sample(n_samples)

    common = [r for r in common if r in leaf_df.columns and r in root_df.columns]

    all_df = pd.concat([leaf_df[common], root_df[common]], axis=0, ignore_index=True)
    std = all_df.std(axis=0)
    keep = std[std > 1e-12].index.tolist()

    return leaf_df[keep].copy(), root_df[keep].copy()


def resolve_reaction_id(df_columns, query):
    if query in df_columns:
        return query
    matches = [c for c in df_columns if c.startswith(query)]
    if len(matches) == 1:
        return matches[0]
    return None


def compute_reaction_zscores(leaf_df, root_df, query_ids):
    rows = []
    all_df = pd.concat([leaf_df, root_df], axis=0, ignore_index=True)
    mu = all_df.mean(axis=0)
    sd = all_df.std(axis=0).replace(0, np.nan)

    for q in query_ids:
        rid = resolve_reaction_id(leaf_df.columns, q)
        if rid is None:
            rows.append({
                "query": q,
                "reaction_id_used": None,
                "leaf_mean_flux": None,
                "root_mean_flux": None,
                "leaf_mean_z": None,
                "root_mean_z": None,
                "abs_delta_z": None,
                "status": "not_found",
            })
            continue

        leaf_z = ((leaf_df[rid] - mu[rid]) / sd[rid]).mean()
        root_z = ((root_df[rid] - mu[rid]) / sd[rid]).mean()

        rows.append({
            "query": q,
            "reaction_id_used": rid,
            "leaf_mean_flux": float(leaf_df[rid].mean()),
            "root_mean_flux": float(root_df[rid].mean()),
            "leaf_mean_z": float(leaf_z),
            "root_mean_z": float(root_z),
            "abs_delta_z": float(abs(leaf_z - root_z)),
            "status": "ok",
        })

    return pd.DataFrame(rows)


def max_production_for_metabolite(model, metabolite):
    with model:
        dm_id = f"DM_TMP_{metabolite.id}"
        if dm_id in model.reactions:
            rxn = model.reactions.get_by_id(dm_id)
        else:
            rxn = model.add_boundary(metabolite, type="demand", reaction_id=dm_id)

        model.objective = rxn
        sol = model.optimize()

        value = float(sol.objective_value) if sol.status == "optimal" and sol.objective_value is not None else None
        return {
            "metabolite_id": metabolite.id,
            "metabolite_name": metabolite.name,
            "status": sol.status,
            "max_flux": value,
        }


def ensure_files_exist():
    required = [
        LEAF_MODEL_FILE,
        ROOT_MODEL_FILE,
        LEAF_MEDIUM_FILE,
        ROOT_MEDIUM_FILE,
    ]
    missing = [f for f in required if not Path(f).exists()]
    if missing:
        raise FileNotFoundError(f"Missing required files: {missing}")


def main():
    ensure_files_exist()

    leaf_model = read_sbml_model(LEAF_MODEL_FILE)
    root_model = read_sbml_model(ROOT_MODEL_FILE)

    apply_medium(leaf_model, LEAF_MEDIUM_FILE)
    apply_medium(root_model, ROOT_MEDIUM_FILE)

    summary_df = pd.DataFrame([
        summarize_model("Leaf", leaf_model),
        summarize_model("Root", root_model),
    ])

    leaf_obj, _ = optimize_biomass(leaf_model, BIOMASS_ID)
    root_obj, _ = optimize_biomass(root_model, BIOMASS_ID)

    leaf_active = active_reaction_counts(leaf_model, BIOMASS_ID)
    root_active = active_reaction_counts(root_model, BIOMASS_ID)

    biomass_df = pd.DataFrame([
        {
            "label": "Leaf",
            "objective_value": leaf_obj,
            "n_active_all_pfba": leaf_active["n_active_all_pfba"],
            "n_active_internal_pfba": leaf_active["n_active_internal_pfba"],
        },
        {
            "label": "Root",
            "objective_value": root_obj,
            "n_active_all_pfba": root_active["n_active_all_pfba"],
            "n_active_internal_pfba": root_active["n_active_internal_pfba"],
        },
    ])

    leaf_sampling_model, leaf_optimum = prepare_sampling_model(
        LEAF_MODEL_FILE, LEAF_MEDIUM_FILE, BIOMASS_ID, GROWTH_FRACTION
    )
    root_sampling_model, root_optimum = prepare_sampling_model(
        ROOT_MODEL_FILE, ROOT_MEDIUM_FILE, BIOMASS_ID, GROWTH_FRACTION
    )

    leaf_samples, root_samples = sample_shared_fluxes(
        leaf_sampling_model,
        root_sampling_model,
        n_samples=N_SAMPLES,
        thinning=THINNING,
        seed=SEED
    )

    zscore_df = compute_reaction_zscores(leaf_samples, root_samples, TARGET_REACTIONS)

    prod_rows = []
    for label, model_file, medium_file in [
        ("Leaf", LEAF_MODEL_FILE, LEAF_MEDIUM_FILE),
        ("Root", ROOT_MODEL_FILE, ROOT_MEDIUM_FILE),
    ]:
        model_gc, optimum = prepare_sampling_model(model_file, medium_file, BIOMASS_ID, GROWTH_FRACTION)

        for metabolite_name, metabolite_id in TARGET_METABOLITES:
            if metabolite_id not in model_gc.metabolites:
                prod_rows.append({
                    "label": label,
                    "query": metabolite_name,
                    "metabolite_id": metabolite_id,
                    "metabolite_name": None,
                    "status": "not_found",
                    "max_flux": None,
                    "growth_fraction": GROWTH_FRACTION,
                    "growth_lb": GROWTH_FRACTION * optimum,
                })
                continue

            metabolite = model_gc.metabolites.get_by_id(metabolite_id)
            result = max_production_for_metabolite(model_gc, metabolite)

            prod_rows.append({
                "label": label,
                "query": metabolite_name,
                "metabolite_id": result["metabolite_id"],
                "metabolite_name": result["metabolite_name"],
                "status": result["status"],
                "max_flux": result["max_flux"],
                "growth_fraction": GROWTH_FRACTION,
                "growth_lb": GROWTH_FRACTION * optimum,
            })

    production_df = pd.DataFrame(prod_rows)

    print("\n=== BASIC MODEL SUMMARY ===")
    print(summary_df.to_string(index=False))

    print("\n=== OBJECTIVE REACTION ===")
    print("Leaf objective ID:", BIOMASS_ID)
    print("Leaf objective equation:", leaf_model.reactions.get_by_id(BIOMASS_ID).reaction)
    print()
    print("Root objective ID:", BIOMASS_ID)
    print("Root objective equation:", root_model.reactions.get_by_id(BIOMASS_ID).reaction)

    print("\n=== BIOMASS OPTIMIZATION AND ACTIVE REACTIONS ===")
    print(biomass_df.to_string(index=False))

    print("\n=== REACTION COMPARISON Z-SCORES ===")
    print(zscore_df.to_string(index=False))

    print("\n=== MAXIMUM PRODUCTION UNDER 90% OPTIMAL BIOMASS ===")
    print(production_df.to_string(index=False))

    summary_df.to_csv("repo_model_summary.csv", index=False)
    biomass_df.to_csv("repo_biomass_and_active_reactions.csv", index=False)
    zscore_df.to_csv("repo_reaction_zscores.csv", index=False)
    production_df.to_csv("repo_metabolite_production.csv", index=False)

    print("\nSaved files:")
    print("- repo_model_summary.csv")
    print("- repo_biomass_and_active_reactions.csv")
    print("- repo_reaction_zscores.csv")
    print("- repo_metabolite_production.csv")


if __name__ == "__main__":
    main()