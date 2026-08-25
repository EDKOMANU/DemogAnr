---
title: 'DemogAnr: Readable, Textbook-Style Demographic Analysis in R'
tags:
  - R
  - demography
  - population projection
  - life tables
  - official statistics
authors:
  - name: Edward Owusu Manu
    orcid: 0009-0006-8808-0018
    affiliation: 1
  - name: Emmanuel Tetteh Amponsah
    affiliation: 2
affiliations:
  - index: 1
    name: Ghana Statistical Service, Accra, Ghana
  - index: 2
    name: Educational Assessment and Research Centre, Accra, Ghana
date: 25 August 2026
bibliography: paper.bib
---

# Summary

DemogAnr is an R package that performs the standard calculations of population demography: birth and death rates, life tables, life expectancy, estimates of fertility and child mortality when registration data are incomplete, and population projections down to the level of individual districts. These calculations underpin official population statistics used for national and regional planning, resource allocation, and public health, and have traditionally been carried out largely in spreadsheets. DemogAnr brings the same calculations into R while keeping every function's inputs and outputs labelled the way they are already taught and reported: a demographer supplies a data table and the names of its columns, and gets back a result that reads like the table they would write into an official report, rather than a nested object that must be further processed before it is usable. The package spans the everyday measures (crude and age-specific rates, complete and abridged life tables), methods for incomplete or defective data (indirect estimation of fertility and child mortality, model life table fitting), and population projection, from simple trend extrapolation to simulation-based methods that carry uncertainty about future rates through to the projected population.

# Statement of need

National statistical offices and demography courses need to move a national population estimate through a specific, well-established sequence: split published age groups into the intervals a district plan requires, build life tables from partial or grouped vital registration data, estimate fertility and child mortality where registration is too incomplete for direct measurement, project each region or district forward, and reconcile the sum of the regional projections back to the agreed national total. Each step is standard, textbook demography, documented in works such as @preston2001 and the United Nations' Manual X [@un1983]. But the R packages that implement pieces of this sequence do so under naming conventions and object structures built for downstream statistical modelling rather than for the way a practising demographer reads a life table or a fertility schedule, and no single package follows the whole production workflow from a census through to a district-level report. DemogAnr was built to close that gap: a small, dependency-light package whose functions and outputs use the same notation, structure, and labelling a demographer already uses on paper, spanning the everyday rate calculations, life tables, indirect estimation, and both simple and simulation-based subnational projection.

# State of the field

R already has a mature ecosystem of demographic software, and DemogAnr is intended to complement it rather than duplicate it. The demography package [@demography-pkg] is built around forecasting: Lee-Carter and functional-data models for mortality and fertility, and stochastic population forecasting; it assumes the user is working within that modelling framework rather than seeking a general calculator for textbook measures. DemoTools [@demotools-pkg], developed for the United Nations Population Division, is the most comprehensive toolkit for evaluating, adjusting, and standardizing demographic data, offering many alternative methods for age graduation, smoothing, and life-table construction; its breadth is its strength, but for an analyst who wants one clearly documented default rather than a menu of methods, that same breadth raises the entry cost. MortalityLaws [@mortalitylaws-pkg] fits parametric mortality laws and builds life tables from fitted models, but is narrower than a general demographic workflow. demogR [@jones2007] focuses on the mathematical analysis of age-structured (Leslie-matrix) population models. bayesPop [@sevcikova2016], together with its companion packages, implements the Bayesian hierarchical probabilistic projection methodology used by the United Nations for the World Population Prospects [@raftery2012]; it is the authoritative tool for fully probabilistic national projections but is a substantial modelling apparatus, not a lightweight calculator. Across this landscape there is no small package that computes the everyday rate and life-table measures, indirect estimation, and subnational projection together, under textbook notation, returning report-ready tables rather than objects built for further modelling. DemogAnr targets exactly that gap and is deliberately complementary to, not a replacement for, the packages above.

This difference is concrete rather than a matter of taste. Given the identical published probabilities of dying from @preston2001 (Box 4.1, United States females, 1991), DemoTools's `lt_abridged()` returns an eleven-column data frame — `Age`, `AgeInt`, `nMx`, `nAx`, `nqx`, `lx`, `ndx`, `nLx`, `Sx`, `Tx`, `ex` — mixing life-table quantities with the age-interval bookkeeping and survivor-ratio (`Sx`) columns its documentation notes are "required for population projections" downstream. DemogAnr's `lifetable_nqx()`, run on the same input, returns the eight columns a demographer would write into a report, named exactly as in @preston2001 (`Age`, `n`, `nqx`, `lx`, `dx`, `Lx`, `Tx`, `ex`), with life expectancy at birth read directly off the `ex` column rather than extracted from a larger computational structure.

# Software design

Three design commitments guided the implementation, written entirely in R [@rcore]. Function and argument names follow standard demographic textbook vocabulary: rates are requested by conventional abbreviations (CBR, TFR, ASFR, CDR), and life-table columns are named exactly as in @preston2001: nMx, nax, nqx, lx, dx, Lx, Tx, ex. Every core function returns a lightweight S3-classed object with a custom print method that lays results out as a labelled table, rather than printing diagnostic output or returning an opaque nested structure; the full numeric detail remains available in the object for further computation, so presentation and computation are separated rather than traded off against each other. Every high-level function takes a data frame plus character column names and returns either a data frame or a classed list, so the same calling convention spans the whole workflow, and the dependency footprint stays small (only dplyr, tidyr, and base stats), keeping installation light and the package easy to maintain.

The corresponding trade-off is scope: DemogAnr is deliberately a calculator for established methods, not a modelling framework, and does not attempt the parametric mortality laws, functional-data forecasting, or Bayesian hierarchical projection that the packages above provide. Where more than one textbook procedure is defensible for the same problem — for example, two accepted assumptions for constructing a cause-deleted life table — the package exposes the choice as an explicit argument rather than silently picking one, so the analyst, not the software, makes the methodological decision.

# Research impact statement

The core of DemogAnr's projection and calibration logic was built for, and used in, producing Ghana's 2021 Population and Housing Census-based subnational population projections, covering all 16 regions and 261 districts, as part of the author's work at the Ghana Statistical Service; the district- and region-level age processing, age-group splitting, and subnational projection functions were developed to carry out that production task before being generalized into the present package. The decision to generalize the code into a public, versioned, documented package followed that production use rather than preceding it. The regional projection reports it supported are published by the Ghana Statistical Service as part of its 2021 Population and Housing Census thematic report series (for example, the *Population Projections 2021-2050* reports issued for each of the 16 regions) [@gss2021projections]; as an internal analytical tool at the time, those reports do not cite the package by name, but its design was shaped directly by that production requirement, and formal, citable use is expected going forward now that the package is public.

# AI usage disclosure

The original DemogAnr codebase, developed for the 2021 Population and Housing Census subnational projection work and the package's initial CRAN release, was written without generative AI assistance. From the subsequent development cycle onward — including this paper — the author used Claude (Anthropic), accessed through Cowork and Claude Code, for code generation, documentation and test authoring, and manuscript drafting. All AI-assisted code and text were reviewed, checked against published worked examples from the demographic literature [@preston2001; @un1983], and edited by the author, who made all methodological and design decisions.

# Acknowledgements

We thank the Ghana Statistical Service for the professional context in which this work originated. No specific funding supported this work.

# References
