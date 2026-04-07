# VEN Gene Signature Analysis - Frontal Cortex Spaceflight

**Do spaceflight conditions dysregulate VEN-associated genes in frontal cortex?**

> Companion analysis to [The Fast Lane Hypothesis](https://github.com/esila-keskin/fast-lane-hypothesis)  
> and [Spaceflight VEN Behavioral Analysis](https://github.com/esila-keskin/spaceflight-ven-analysis)

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://python.org)
[![Data](https://img.shields.io/badge/data-NASA%20OSD--698-orange.svg)](https://osdr.nasa.gov/bio/repo/data/study/OSD-698)
[![GEO](https://img.shields.io/badge/GEO-GSE239336-green.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239336)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![AWG](https://img.shields.io/badge/NASA-Brain%20AWG-blue.svg)](https://awg.osdr.space)

---

## Overview

The **Fast Lane Hypothesis** proposes that Von Economo Neurons (VENs) implement a biological speed-accuracy tradeoff in social decision circuits via fast projection from frontal cortex. A key prediction is that spaceflight stressors should affect the molecular machinery underlying this fast pathway.

This repository tests that prediction using **GeoMx Digital Spatial Profiling** of mouse frontal cortex tissue from real ISS spaceflight (SpaceX-24 mission, 35 days).

We define a **VEN gene signature** from the literature and ask: are these genes enriched among spaceflight-dysregulated genes in frontal cortex?

---

## Key Findings

| Gene | Category | Log2FC | p-value | Direction |
|------|----------|--------|---------|-----------|
| **Cnp** | Myelination | +0.618 | 0.011* | ↑ spaceflight |
| **Sod2** | Oxidative Stress | +0.438 | 0.041* | ↑ spaceflight |
| **Snap25** | Fast Signalling | +0.422 | 0.044* | ↑ spaceflight |
| Nos1 | VEN Identity | -0.354 | 0.086 | ↓ trending |
| Mag | Myelination | +0.545 | 0.158 | ↑ trending |
| Mbp | Myelination | +0.531 | 0.174 | ↑ trending |

**VEN gene enrichment:** 11.1% of VEN genes significant vs 4.6% background  
Fisher exact OR=2.62, p=0.12 - enriched but underpowered at n=12 samples

**Interpretation:**
- Myelination genes (CNP, MAG, MBP) are upregulated - compensatory response to spaceflight-induced stress in fast-conducting axons
- SNAP25 upregulation suggests increased synaptic vesicle machinery under spaceflight stress
- NOS1 (nNOS - a direct VEN biochemical marker) trends downward - consistent with VEN pathway disruption

---

## Dataset

**NASA OSD-698 / GEO GSE239336** - Kremsky et al. (2023)  
*Spaceflight-Induced Gene Expression Profiles in the Mouse Brain Are Attenuated by Treatment with the Antioxidant BuOE*  
DOI: [10.3390/ijms241713569](https://doi.org/10.3390/ijms241713569)

- Mice flown on SpaceX-24 (ISS, 35 days, Dec 2021–Jan 2022)
- GeoMx Digital Spatial Profiling - frontal cortex tissue
- This analysis: Ground Control vs Spaceflight (saline group)

---

## VEN Gene Signature

Curated from Allman et al. (2010), Stimpson et al. (2011), Hodge et al. (2020):

| Category | Genes | Rationale |
|----------|-------|-----------|
| VEN Identity | Fxyd1, Gadd45g, Adra1a, Nos1, Disc1, Bcl11b | Selective VEN markers from literature |
| Fast Signalling | Nefh, Nefm, Nefl, Syt1, Snap25, Vamp2 | Large-axon and fast-release machinery |
| Myelination | Mbp, Mog, Plp1, Mag, Cnp | Fast conduction speed |
| Oxidative Stress | Sod1, Sod2, Cat, Gpx1, Nox4 | VENs are large = high metabolic demand |
| Synaptic | Grin1, Grin2a, Grin2b, Dlg4, Shank3 | Postsynaptic integration |

---

## Repository Structure
```
ven-frontalcortex-spaceflight/
├── data/
│   ├── FCT_DEanalysis.txt <- DE results (GC vs spaceflight)
│   └── FCT_GeneExpression_Q3norm.txt   <- Normalised expression matrix
├── analysis/
│   └── ven_gene_signature_analysis.py
├── figures/
│   ├── fig_ven_volcano.pdf
│   ├── fig_ven_significant.pdf
│   └── fig_ven_categories.pdf
├── results/
│   └── ven_gene_results.json
├── requirements.txt
└── README.md
```

---

## Data Download

1. Go to [GEO GSE239336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239336)
2. Download supplementary files:
   - `GSE239336_FCT_GCvsFLT-SAL_DEanalysis.txt.gz`
   - `GSE239336_FCT_GeneExpression_Q3norm.txt.gz`
3. Decompress and rename to `FCT_DEanalysis.txt` and `FCT_GeneExpression_Q3norm.txt`
4. Place both in `data/`

---

## Usage
```bash
pip install -r requirements.txt
python analysis/ven_gene_signature_analysis.py
```

---

## Relationship to Fast Lane Hypothesis

| Analysis | Data | Question |
|----------|------|----------|
| [Fast Lane Model](https://github.com/esila-keskin/fast-lane-hypothesis) | Computational SNN | Do VENs implement speed-accuracy tradeoff? |
| [Behavioral](https://github.com/esila-keskin/spaceflight-ven-analysis) | NASA OSD-618 | Do spaceflight stressors produce VEN-consistent behavioral signatures? |
| **This repo** | NASA OSD-698 frontal cortex | Are VEN-associated genes dysregulated by real ISS spaceflight? |

---

## Limitations

- Mice do not have VENs - frontal cortex is the closest available proxy to ACC
- n=12 samples limits statistical power
- GeoMx DSP measures RNA at region level, not layer V specifically

---

## Paper

Companion to: **The Fast Lane Hypothesis: Von Economo Neurons Implement a Biological Speed-Accuracy Tradeoff** - Esila Keskin, UWE Bristol (2026)

---

## NASA AWG

Open project under the **Brain AWG**, NASA OSDR Analysis Working Groups.

---

## License

MIT
