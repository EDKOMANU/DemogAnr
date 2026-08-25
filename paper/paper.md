---
title: "DemogAnr: Readable, Textbook-Style Demographic Analysis in R"
tags:
- R
- demography
- population projection
- life tables
- official statistics
date: "25 August 2026"
output:
  word_document: default
  pdf_document: default
authors:
- name: Edward Owusu Manu
  orcid: "0009-0006-8808-0018"
  affiliation: 1
- name: Emmanuel Tetteh Amponsah
  affiliation: 2
bibliography: paper.bib
affiliations:
- index: 1
  name: Ghana Statistical Service, Accra, Ghana
- index: 2
  name: Educational Assessment and Research Centre, Accra, Ghana
---

# Summary

DemogAnr is an R package that performs the standard calculations of population demography: birth and death rates, life tables, life expectancy, estimates of fertility and child mortality when registration data are incomplete, and population projections down to the level of individual districts. These calculations underpin official population statistics used for planning, resource allocation and public health, and have traditionally been carried out largely in spreadsheets. DemogAnr brings them into R while keeping every function's inputs and outputs labelled the way they are already taught and reported: a demographer supplies a data table and the names of its columns, and gets back a result that reads like the table they would write into an official report, rather than a nested object needing further processing. The package spans the everyday measures (crude and age-specific rates, complete and abridged life tables); the classical methods for incomplete or defective data, comprising indirect estimation of fertility and of child and adult mortality, the Coale-Demeny and United Nations model life table systems [@coale1983; @un1982], the Brass relational logit model [@brass1975], and the death distribution methods by which the completeness of death registration is itself estimated [@hill1987; @preston1982]; and population projection, from trend extrapolation to simulation-based methods that carry uncertainty about future rates through to the projected population.

# Statement of need

National statistical offices and demography courses must move a national population estimate through a specific, well-established sequence: establish how completely deaths are registered before any life table is built from them, split published age groups into the intervals a district plan requires, build life tables from partial or grouped vital registration data, estimate fertility and child and adult mortality where registration is too incomplete for direct measurement, convert those fragmentary estimates into a complete age pattern of mortality, project each region or district forward, and reconcile the sum of the regional projections back to the agreed national total. Each step is standard textbook demography, documented in @preston2001 and the United Nations' Manual X [@un1983]. But the R packages that implement pieces of this sequence do so under naming conventions and object structures built for downstream statistical modelling rather than for the way a practising demographer reads a life table or a fertility schedule, and no single package follows the whole production workflow from a census through to a district-level report. DemogAnr was built to close that gap: a small, dependency-light package whose functions and outputs use the same notation, structure, and labelling a demographer already uses on paper, spanning the everyday rate calculations, life tables, the assessment of registration completeness, indirect estimation of fertility and of child and adult mortality, model life tables, and both simple and simulation-based subnational projection.

# State of the field

R has a mature ecosystem of demographic software, which DemogAnr complements rather than duplicates. The demography package [@demography-pkg] is built around forecasting — Lee-Carter and functional-data models — and assumes the user is working within that framework rather than seeking a calculator for textbook measures. DemoTools [@demotools-pkg], developed for the United Nations Population Division, is the most comprehensive toolkit for evaluating and adjusting demographic data, offering many alternative methods for graduation, smoothing and life-table construction; its breadth is its strength, but for an analyst wanting one documented default rather than a menu, that breadth raises the entry cost. MortalityLaws [@mortalitylaws-pkg] fits parametric mortality laws and demogR [@jones2007] analyses Leslie-matrix models; both are narrower than a general workflow. bayesPop [@sevcikova2016] implements the Bayesian hierarchical projection methodology the United Nations uses for the World Population Prospects [@raftery2012], and is authoritative for probabilistic national projections, but is a modelling apparatus rather than a calculator. Across this landscape no small package computes the everyday rate and life-table measures, indirect estimation and subnational projection together, under textbook notation, returning report-ready tables. No individual method DemogAnr implements is novel; all are established procedures from @preston2001, @brass1975 and @un1983. Its contribution is their assembly under a single calling convention, and the verification described below.

Why build rather than contribute? Several of these methods already exist elsewhere, and were the contribution a method, upstreaming it would have been the right course. It is not. The contribution is a convention — applied uniformly to how arguments are named, how results are structured, and how they print — and a convention cannot be added to a mature package as a patch, because it is not additive: adopting it would mean redesigning an interface that package's users depend on. Building separately, and stating that the two are complementary, imposes no such cost.

The difference is concrete rather than a matter of taste. On the identical published probabilities of dying from @preston2001 (Box 4.1, United States females, 1991), DemoTools's `lt_abridged()` returns eleven columns — `Age`, `AgeInt`, `nMx`, `nAx`, `nqx`, `lx`, `ndx`, `nLx`, `Sx`, `Tx`, `ex` — mixing life-table quantities with age-interval bookkeeping and the survivor ratios its documentation notes are required for projections downstream. DemogAnr's `lifetable_nqx()` returns the eight columns a demographer would write into a report, named as in @preston2001 (`Age`, `n`, `nqx`, `lx`, `dx`, `Lx`, `Tx`, `ex`), with life expectancy read off the `ex` column rather than extracted from a larger structure.

# Software design

Three design commitments guided the implementation, written entirely in R [@rcore]. Names follow textbook vocabulary: rates are requested by their conventional abbreviations (CBR, TFR, ASFR, CDR), and life-table columns are named as in @preston2001. Every core function returns a lightweight S3-classed object with a custom print method that lays results out as a labelled table, rather than printing diagnostic output or returning an opaque nested structure; the full numeric detail remains available in the object for further computation, so presentation and computation are separated rather than traded off against each other. Every high-level function takes a data frame plus character column names and returns either a data frame or a classed list, so the same calling convention spans the whole workflow, and the dependency footprint stays small (dplyr, tidyr, ggplot2, and base stats), keeping installation light and the package easy to maintain.

The corresponding trade-off is scope: DemogAnr is deliberately a calculator for established methods, not a modelling framework, and does not attempt the parametric mortality laws, functional-data forecasting or Bayesian hierarchical projection the packages above provide. Where more than one textbook procedure is defensible for the same problem; for example, two accepted assumptions for constructing a cause-deleted life table — the package exposes the choice as an argument, so the analyst makes the methodological decision.

# Validation

Because the package is intended for the production of official statistics rather than for exploratory modelling, each method is verified against an external standard rather than against its own output. Where a published worked example exists, the test suite reproduces it: the Brass P/F ratio, indirect child mortality and maternal orphanhood on the Bangladesh 1974, Panama 1976 and Bolivia 1975 examples of @un1983 (Chapters II-IV); the crude death rate, cohort-component projection, singulate mean age at marriage and direct age standardization on Boxes 1.2, 6.1, 4.4 and 2.1 of @preston2001; and the United Nations age-sex accuracy index on the Ghana 1960 example of @kpedekpo1982. The orphanhood implementation reproduces its example exactly: all seven weighting factors to the four decimal places published, and the seven survivorship ratios and six time-reference periods to the precision *Manual X* reports [@brasshill1973].

Methods for which the literature supplies a construction rather than numbers are verified against the identities they must satisfy: the Kitagawa and Arriaga decompositions, for instance, must sum exactly to the difference in crude rates and in life expectancy each decomposes.

Where neither a published example nor a closed identity exists, verification uses a synthetic population whose true value is known by construction. The death distribution methods are tested against a stable population obtained by integrating a Siler mortality schedule on a fine age grid, from which a known fraction of deaths is withheld to simulate incomplete registration; growth balance, extinct generations and the combined procedure each recover the withheld fraction to within one per cent for completeness between 0.3 and 1.0.

All 37 exported functions are covered by the test suite, which comprises more than 500 assertions and runs on each commit against macOS, Windows and Linux across four versions of R, including the oldest the package declares.

# Research impact statement

DemogAnr's projection and calibration logic was built for, and used in, producing Ghana's 2021 Population and Housing Census-based subnational projections covering all 16 regions and 261 districts, as part of the author's work at the Ghana Statistical Service; the age processing, age-group splitting and subnational projection functions were developed for that production task before being tested and generalized into the present package. The decision to generalize the code into a public, versioned, documented package followed that production use. The regional projection reports it supported are published by the Ghana Statistical Service as part of its 2021 Population and Housing Census thematic report series (for example, the *Population Projections 2021-2050* reports issued for each of the 16 regions) [@gss2021projections]; as an internal analytical tool at the time, those reports do not cite the package by name, but its design was shaped directly by that production requirement. Beyond that production history, the package is released on CRAN, and the published worked examples reproduced above are executable by any reader directly from its test files, so the correctness claims made here can be checked.

# AI usage disclosure

The original DemogAnr codebase, developed for the 2021 Population and Housing Census subnational projection work and the package's initial CRAN release, was written without generative AI assistance. From the subsequent development cycle onward — including this paper — the author used Claude and (Anthropic) and Gemini (Google), accessed through Cowork,Claude Code and Google Antigravity, for code generation, documentation and test authoring. All AI-assisted code and text were reviewed and edited by the author, who made all methodological and design decisions.

# Acknowledgements

We thank the Ghana Statistical Service for the professional context in which this work originated. No specific funding supported it.

# References
