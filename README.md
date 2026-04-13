# Tissue-Specific Genome-Scale Metabolic Models of *Erythroxylum novogranatense*

This repository contains two tissue-specific genome-scale metabolic models (GEMs) of *Erythroxylum novogranatense*, one for leaf and one for root, together with the medium files used for the constraint-based simulations described in the associated manuscript.

## Repository contents

### Models
- `model leaf E novogranatense.xml`
- `model root E novogranatense.xml`

### Medium files used in simulations
- `medium leaf E novogranatense.tsv`
- `medium root E novogranatense.tsv`

### Objective function
- `bio2_biomass`

## Project overview

*Erythroxylum novogranatense* is a chemically rich non-model plant that produces specialized metabolites in a tissue-dependent manner. To provide a systems-level framework for studying this metabolic specialization, we reconstructed and curated two tissue-specific GEMs, one for leaf and one for root.

The models were built from:
- a chromosome-level genome assembly,
- functional genome annotation,
- targeted curation of central and specialized metabolism,
- and untargeted LC–QTOF–MS metabolomics used for tissue-specific contextualization.

The final models were used to study:
- structural and functional similarities between tissues,
- tissue-specific flux behavior,
- and differences in the production potential of selected metabolites under comparable growth constraints.

## Reconstruction strategy

The workflow followed the general strategy described in the manuscript:

1. Retrieval of a chromosome-level genome assembly of *E. novogranatense*.
2. Functional annotation using complementary resources, including:
   - Plant RAST / PlantSEED
   - eggNOG-mapper
   - plantiSMASH
3. Generation of a draft metabolic reconstruction.
4. Manual and semi-automatic curation of primary and specialized metabolism.
5. Integration of untargeted LC–QTOF–MS metabolomics data from leaves and roots.
6. Tissue-specific contextualization of the models.
7. Constraint-based simulations under tissue-specific media.

## Biological scope

These models were developed to provide a tissue-resolved metabolic framework for exploring how *E. novogranatense* distributes metabolic capacity between growth and the synthesis of specialized metabolites.

The manuscript focuses on:
- alkaloid, flavonoid, terpenoid, and lipid-related metabolism,
- tissue-specific metabolite incorporation from metabolomics,
- flux-space separation between leaf and root,
- and tissue-dependent production potential for selected metabolites.

## Graphical abstract

<img width="1063" height="911" alt="image" src="https://github.com/user-attachments/assets/1c49e066-e7fe-47aa-bc2e-7f50010a1e37" />


## File description

### Models
The SBML files contain the leaf and root tissue-specific GEMs.

### Media
The TSV files define the simulation media used for each tissue.

### Objective
All basic growth simulations in this repository use:
- `bio2_biomass`

## Minimal local usage with COBRApy

```python
from cobra.io import read_sbml_model

leaf = read_sbml_model("model leaf E novogranatense.xml")
root = read_sbml_model("model root E novogranatense.xml")

print("Leaf:", len(leaf.reactions), len(leaf.metabolites), len(leaf.genes))
print("Root:", len(root.reactions), len(root.metabolites), len(root.genes))
