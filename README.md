# HYCOM Spatial Bottom-Pressure Correlation (Trench Reference-Station Study)

This repository contains scripts and workflows to quantify **spatial coherence of bottom pressure** using **HYCOM** model outputs.  
We sample bottom-pressure time series at points **across-trench** and **along-trench**, compute pairwise correlations, and use the resulting correlation-length scales to guide:

- **How to choose a reference station**
- **How far a reference station can be from a target site** while remaining effective at suppressing common-mode oceanographic variability

The main deliverable is a set of correlation maps/curves that summarize the **most efficient reference-site range** for differenced pressure analyses.

---

## Scientific motivation

Differential bottom pressure (e.g., `ΔP = P_target − P_ref`) can suppress shared oceanographic variability, but its effectiveness depends on how **correlated** bottom-pressure signals are over distance and direction.

This project evaluates correlation structure by:
- selecting points **along-trench** (directional coherence)
- selecting points **across-trench** (cross-slope coherence)
- comparing correlation decay and identifying the distance window where a reference station performs best

---

## Data

**Model:** HYCOM (HYbrid Coordinate Ocean Model)  
**Variable(s):** bottom pressure (or sea-floor pressure / pressure anomaly depending on product)  
**Spatial sampling:** user-defined points along and across the trench  
**Time window:** configurable in the scripts

> Note: HYCOM files are not included in this repo. See `config/` and “How to run” for expected inputs.

---



