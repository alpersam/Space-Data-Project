# GRACE Project Research Outline

## Can Satellites Predict Grain Prices?
**GRACE TWS as a Leading Indicator for Agricultural Commodity Futures**

Space Data FS2026 — Block II | ETH Zürich

---

## 1. Research Question

> *To what extent do GRACE-derived terrestrial water storage anomalies over Brazil's central-southern agricultural belt (Mato Grosso, Goiás, Paraná) anticipate variations in international soybean futures prices, and what is the temporal lag structure between hydrological stress and commodity market response?*

**Why this satisfies the rubric (10% — Research question & topic framing):**
- Specific region, specific signal, specific external dataset
- Grounded in Seo et al. (2025, *Science*) — the regime-shift paper on global soil moisture depletion
- Directly links a GRACE geophysical observable (TWS) to a measurable economic outcome
- Not a yes/no question — asks about *extent* and *lag structure*

---

## 2. Scientific Framing & Literature

### Core references

| Reference | Relevance |
|-----------|-----------|
| **Seo et al. (2025)** — *Science* 387, 1408–1413 | Global soil moisture regime shift; ~1614 Gt loss 2000–2002. South America among most affected. Provides the scientific backbone. |
| **Gonçalves et al. (2020)** — *Sci. Total Environ.* 705, 135845 | GRACE TWS depletion in NE Brazil (Urucuia Aquifer), rate of −6.5 ± 2.6 mm/yr, driven by anthropogenic irrigation for soybean/corn. |
| **Eom et al. (2017)** — *Remote Sens. Environ.* 191, 55–66 | EOF analysis of GRACE over the Amazon/Óbidos sub-basin. EOF1 captures 66.9% variance. Directly relevant methodology for your region. |
| **Landerer & Swenson (2012)** — *Water Resour. Res.* | JPL mascon gain factors for hydrology applications — justifies your scaling approach. |
| **Samaniego (2025)** — *Science* 387, 1348–1350 | Perspective on Seo et al. — calls the observed changes a "permanent" and "irreversible" shift. |

### Supporting context

| Topic | Reference |
|-------|-----------|
| ENSO–TWS linkage in La Plata basin | TWS–TRMM comparison (2003–2017), strong correlation with El Niño/La Niña events |
| Brazil soybean production & water dependency | Brazil is world's largest soybean exporter; growing season Oct–Mar depends on Jul–Sep water availability |
| Satellite alt-data in quant finance | Two Sigma, Citadel, Point72 invest in non-traditional datasets (satellite imagery, shipping data, etc.) |

---

## 3. Analysis Workflow — Mapped to Rubric

### Step-by-step plan

The rubric weights **Methods & Analysis at 30%** — the largest single category. Every step below must be justified in the report.

#### Step 1: GRACE Mascon Preparation
- **Data:** JPL RL06M v03 mascons (NetCDF, already available from HW1)
- **Action:** Mask land pixels over central-southern Brazil (Mato Grosso ~10–18°S, 50–58°W; Goiás; Paraná basin)
- **Gain factors:** Apply JPL hydrology gain factors (valid — no ice sheets in study area)
- **Gap handling:** Linear interpolation for isolated missing months; exclude the 2017–2018 GRACE/GRACE-FO transition gap from cross-correlation (document this choice!)
- **Output:** Area-weighted monthly TWS anomaly time series (2002–2023)
- **Justification for rubric:** Demonstrates understanding of mascon scaling, region selection, and gap treatment

#### Step 2: EOF/PCA Analysis
- **Action:** Run PCA on the masked mascon grid over the study region
- **EOF1:** Should capture the dominant TWS depletion/variability mode
- **PC1 time series → "GRACE drought index"** for the region
- **Check explained variance** — if EOF1 < 50%, consider whether the region is too large or includes mixed signals
- **Justification:** Signal separation technique applied thoughtfully (rubric: "method choices justified, project-specific limitations briefly discussed")

#### Step 3: External Data — ERA5-Land P & ET
- **Source:** Copernicus Climate Data Store (CDS), monthly, free
- **Variables:** Total precipitation (P), evapotranspiration (ET)
- **Action:** Compute P − ET anomaly over the same masked region
- **Purpose:** Physical mechanism linking atmospheric water balance to subsurface TWS changes. This validates *why* TWS changes — not just *that* they change.
- **Justification:** Strengthens the causal chain (precipitation deficit → TWS decline → crop stress → price reaction)

#### Step 4: External Data — Soybean Futures Prices
- **Source options (all free):**
  - Yahoo Finance: `yfinance` Python library, ticker `ZS=F` for CBOT soybean continuous
  - Investing.com: manual CSV download of CBOT soybean monthly close
  - FAO Food Price Index (cereals sub-index): CSV from fao.org/worldfoodsituation
- **Action:** Download monthly close prices, compute log-returns
- **Note:** Use front-month continuous contract (unadjusted monthly close is fine for this analysis)

#### Step 5: Lead-Lag Cross-Correlation Analysis
- **Action:** Compute cross-correlation between GRACE PC1 and soybean price returns at lags τ = 0, 1, 2, ..., 6 months
- **Statistical test:** Granger causality test (bivariate VAR) — does TWS *statistically precede* price movements?
- **Implementation:** `statsmodels.tsa.stattools.grangercausalitytests` in Python
- **Expected result:** Significant correlation at τ = 3–6 months (hypothesis). If weak/insignificant, this is still a valid result — discuss confounding factors (USD strength, Chinese demand, speculation)

#### Step 6: Event Study — Drought Episodes
- **Focus events:** 2005, 2010, 2014–2015, 2021 Brazilian droughts
- **Action:** For each event, plot GRACE TWS decline timeline vs. when soybean futures reacted
- **Question:** How many months in advance does GRACE signal the TWS drop vs. when markets move?
- **This is the "narrative" evidence** to complement the statistical analysis

#### Step 7 (Optional): Backtest Trading Signal
- **Signal:** Go long soybeans when GRACE TWS drops below −1σ of its mean
- **Metric:** Sharpe ratio on simple out-of-sample backtest
- **Caveat:** This is illustrative, not rigorous quant finance. Frame it as a demonstration of signal utility, not a trading strategy.

---

## 4. Data Sources — Complete List

| Dataset | Source | Format | Access | Cost |
|---------|--------|--------|--------|------|
| GRACE/GRACE-FO mascons | JPL RL06M v03 | NetCDF | Already have from HW1 | Free |
| ERA5-Land P & ET | Copernicus CDS | NetCDF | `cdsapi` Python package | Free |
| CBOT soybean futures | Yahoo Finance / Investing.com | CSV | `yfinance` or manual download | Free |
| FAO Food Price Index | fao.org | CSV | Direct download | Free |
| Seo et al. (2025) supplementary data | Science | — | DOI: 10.1126/science.adq6529 | ETH library access |

---

## 5. Expected Figures (4 Key Figures for Report)

The rubric says "key figures instead of a lot of plots" and the report is max 5 pages including figures. Plan exactly 4 figures:

1. **Map: TWS trend (2002–2023)** over central-southern Brazil
   - Spatial pattern of water depletion, with study region boundary overlaid
   - Colorbar in mm/yr (consistent with Seo et al.)

2. **EOF1 spatial pattern + PC1 time series**
   - Two-panel: left = EOF1 map, right = PC1 over time
   - This is your "GRACE drought index"
   - Annotate drought years (2005, 2014, 2021)

3. **Overlay: GRACE TWS vs. soybean futures + ERA5 P−ET**
   - Three signals on the same normalized timeline
   - Dual y-axis: TWS/P−ET on left, price on right
   - Annotate correlation coefficients

4. **Cross-correlogram**
   - Correlation of TWS vs. soybean returns as function of lag (months)
   - Identify the optimal predictive window
   - Include 95% confidence bands

---

## 6. Timeline

| Deadline | Tasks |
|----------|-------|
| **This week (13–19 Mar)** | Region masking + GRACE TWS time series. EOF analysis. Download ERA5 + price data. Start cross-correlation. **Prepare 3-slide kickoff presentation** (due on Moodle by 19.3) |
| **20 Mar** | 3 min presentation in class. Get feedback from Annika. |
| **21 Mar – 2 Apr** | ERA5 P−ET processing. Full cross-correlation & Granger analysis. Event study on drought episodes. Generate all 4 figures. Write report. |
| **3 Apr** | Submit final report (max 5 pages + code) via Moodle |

---

## 7. Risks & Mitigations

| Risk | Mitigation |
|------|-----------|
| **Weak/no correlation** | A null result is valid. Discuss confounding factors: USD index, oil prices, Chinese import policy, speculative positioning. Frame as evidence of market efficiency (prices already incorporate available info). |
| **Short time series (~20 yr monthly)** | Limits statistical power at long lags. Supplement with ERA5 P−ET as a physically motivated proxy. Use bootstrap confidence intervals for correlation estimates. |
| **GRACE/GRACE-FO gap (2017–2018)** | Document gap handling explicitly. Interpolate or exclude. Show sensitivity to this choice. |
| **Multiple confounders in commodity prices** | Control by detrending both series, using returns instead of levels, and acknowledging this as a univariate study in the discussion. |
| **Gain factor choice** | JPL hydrology gain factors are fully appropriate — no ice sheets, no ocean signals in study region. State this in methods. |

---

## 8. Mapping to Grading Rubric

### Presentation (30% of total: 10% + 10% + 10%)

| Category | What to show | Target grade |
|----------|-------------|-------------|
| Research question & topic framing (10%) | RQ grounded in Seo et al. (2025); clear scope (Brazil, soybeans, 3–6 month lag) | 4 — Insightful |
| Analysis strategy & method plan (10%) | EOF justification, gain factors, Granger causality rationale, gap handling | 4 — Clear and appropriate |
| Presentation clarity (10%) | 3 slides max. Lead with the hook ("can satellites predict grain prices?"), show method flowchart, end with expected result | 4 — Very clear & concise |

### Report (70% of total: 30% + 20% + 15% + 5%)

| Category | What to deliver | Target grade |
|----------|----------------|-------------|
| GRACE data handling & analysis (30%) | EOF applied to masked Brazil region, gain factors justified, gap handling documented, PC1 as drought index | 4 — Methods applied appropriately with strong justification |
| External dataset & comparison (20%) | Three external datasets (ERA5 P−ET, soybean futures, FAO index). Cross-correlation with lag analysis. Granger test. Event study. | 4 — Relevant, insightful comparison |
| Discussion & presentation clarity (15%) | Interpret lag structure physically (soil moisture → crop stress → harvest impact → price). Discuss limitations (confounders, short series). "One thing that surprised us." | 4 — Results interpreted critically |
| Code clarity & reproducibility (5%) | Well-documented Jupyter notebook with README section. All data paths relative. Requirements listed. | 2 pts — it works |

---

## 9. Python Package Requirements

```
numpy
scipy
matplotlib
cartopy
netCDF4
xarray
pandas
statsmodels          # Granger causality, VAR
scikit-learn         # PCA
yfinance             # soybean price data
cdsapi               # ERA5 download (if needed)
```

---

## 10. Why This Project Stands Out

Most teams will study ice sheets or straightforward hydrology comparisons. This project:

- Connects GRACE data to a completely different domain (financial markets)
- Satisfies every rubric requirement: clear RQ, EOF signal separation, meaningful external dataset comparison, geophysical interpretation
- Has a direct tie to Seo et al. (2025) — the most recent high-profile GRACE paper
- Demonstrates real-world economic relevance of satellite gravity missions
- The quantitative finance angle (alternative data, lead-lag, Granger causality) adds originality without sacrificing rigor

**In short: two satellites flying 220 km apart might tell you when to buy soybean futures.**
