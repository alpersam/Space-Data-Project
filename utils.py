"""
utils.py -- GRACE Soybean Analysis helpers
==========================================
Data-loading functions, synthetic-data generators, and statistical utilities.

Usage in notebook:
    from utils import (
        load_grace_mascons, make_synthetic_grace,
        load_era5, make_synthetic_era5,
        get_ndvi_gee, make_synthetic_ndvi,
        get_soybean_prices, make_synthetic_soybean,
        to_monthly, adf_test,
        cross_corr_bootstrap, rolling_corr_at_lag,
        partial_corr_at_lag, run_granger,
        drought_event_study, illustrative_backtest,
    )
"""

import numpy as np
import pandas as pd
import xarray as xr
import yfinance as yf
from scipy.stats import pearsonr
from statsmodels.tsa.stattools import grangercausalitytests, adfuller


# GRACE -----------------------------------------------------------------------

def load_grace_mascons(filepath):
    """Load JPL RL06M v04 NetCDF, normalise dims, apply hydrology gain factors."""
    ds = xr.open_dataset(filepath)
    rm = {}
    for d in ds.dims:
        if "lat"  in d.lower(): rm[d] = "lat"
        if "lon"  in d.lower(): rm[d] = "lon"
        if "time" in d.lower(): rm[d] = "time"
    ds = ds.rename(rm)
    tv = next((v for v in ds.data_vars if "lwe" in v.lower()), list(ds.data_vars)[0])
    tws = ds[tv]
    gv = next((v for v in ds.data_vars if "scale" in v.lower() or "gain" in v.lower()), None)
    if gv:
        tws = tws * ds[gv]
    if tws.lon.values.max() > 180:
        tws = tws.assign_coords(lon=(tws.lon % 180) - (tws.lon // 180)*180).sortby("lon")
    tws.attrs["units"] = "cm"
    return tws


def make_synthetic_grace(lat_min, lat_max, lon_min, lon_max,
                          t_start, t_end, gap_start, gap_end, seed=42):
    """Physically plausible synthetic GRACE TWS: trend + seasonal + drought events."""
    rng   = np.random.default_rng(seed)
    times = pd.date_range(t_start, t_end, freq="MS")
    lats  = np.arange(lat_min, lat_max + 0.5, 1.0)
    lons  = np.arange(lon_min, lon_max + 0.5, 1.0)
    nt, nlat, nlon = len(times), len(lats), len(lons)
    t = np.arange(nt)
    trend_map = np.outer(np.linspace(-0.12, -0.06, nlat)[::-1],
                         np.linspace(0.8, 1.2, nlon))
    mi = np.array([d.month for d in times])
    seasonal = 8.0 * np.sin(2 * np.pi * (mi - 1) / 12)
    drought  = np.zeros(nt)
    for yr, amp, dur in [(2005,-18,18),(2010,-22,20),(2015,-30,24),(2021,-25,18)]:
        c_arr = np.where([d.year==yr and d.month==9 for d in times])[0]
        if len(c_arr):
            c = c_arr[0]; w = np.arange(max(0,c-dur//2), min(nt,c+dur//2))
            drought[w] += amp * np.exp(-0.5*((w-c)/(dur/4))**2)
    data = np.zeros((nt, nlat, nlon))
    for i in range(nt):
        data[i] = trend_map*t[i] + seasonal[i] + drought[i] + rng.normal(0, 3., (nlat, nlon))
    gap_mask = (times >= gap_start) & (times <= gap_end)
    data[gap_mask] = np.nan
    return xr.DataArray(data, coords={"time": times, "lat": lats, "lon": lons},
                         dims=["time","lat","lon"],
                         attrs={"units":"cm","note":"SYNTHETIC -- replace with JPL RL06M"})


# ERA5 -----------------------------------------------------------------------

def load_era5(filepath):
    """Load ERA5-Land NetCDF and return area-weighted P-ET time series (mm/month)."""
    ds = xr.open_dataset(filepath)
    pv  = next(v for v in ds.data_vars if "tp" in v.lower() or "prec" in v.lower())
    etv = next(v for v in ds.data_vars if "pev" in v.lower()
               or v.lower() in ("e","et") or "evap" in v.lower())
    P = ds[pv]*1000; ET = ds[etv]*1000
    PET = P - np.abs(ET)
    lat_key = [k for k in ds.dims if "lat" in k.lower()][0]
    lon_key = [k for k in ds.dims if "lon" in k.lower()][0]
    w = np.cos(np.deg2rad(ds[lat_key].values))
    PET_m = (PET * w[np.newaxis,:,np.newaxis]).sum(dim=[lat_key, lon_key]) / \
            (w.sum() * len(ds[lon_key]))
    return PET_m.to_series()


def make_synthetic_era5(times, drought_years):
    """Generate synthetic P-ET series aligned to `times`."""
    rng = np.random.default_rng(123); t = np.arange(len(times))
    mi  = np.array([d.month for d in times])
    seasonal = 50 * np.sin(2*np.pi*(mi-1)/12)
    drought  = np.zeros(len(times))
    for yr, amp, dur in [(2005,-40,16),(2010,-50,18),(2015,-65,22),(2021,-55,16)]:
        c_arr = np.where([d.year==yr and d.month==8 for d in times])[0]
        if len(c_arr):
            c = c_arr[0]; w = np.arange(max(0,c-dur//2), min(len(times),c+dur//2))
            drought[w] += amp*np.exp(-0.5*((w-c)/(dur/4))**2)
    return pd.Series(-0.15*t + seasonal + drought + rng.normal(0, 8, len(times)),
                     index=times, name="P_minus_ET_mm")


# MODIS NDVI (Google Earth Engine) -------------------------------------------

def get_ndvi_gee(lat_min, lat_max, lon_min, lon_max,
                  gee_project, start="2002-01-01", end=None):
    """
    Extract monthly cropland-masked NDVI from MODIS MOD13A3 v061 via GEE.

    Requires: earthengine-api, authenticated GEE project.
    Set USE_GEE = False to load from a local CSV instead.
    """
    import ee, sys
    if end is None: end = pd.Timestamp("today").strftime("%Y-%m-%d")
    try:
        ee.Initialize(project=gee_project)
    except Exception:
        if "google.colab" in sys.modules:
            ee.Authenticate(auth_mode="notebook")
        else:
            ee.Authenticate(auth_mode="localhost")
        ee.Initialize(project=gee_project)
    region = ee.Geometry.Rectangle([lon_min, lat_min, lon_max, lat_max])
    lc = (ee.ImageCollection("MODIS/061/MCD12Q1")
            .filterDate("2020-01-01","2020-12-31").first().select("LC_Type1"))
    crop_mask = lc.eq(12).Or(lc.eq(14))
    col = (ee.ImageCollection("MODIS/061/MOD13A3")
             .filterDate(start, end).select(["NDVI","EVI"]))
    def extract(img):
        masked = img.updateMask(crop_mask)
        s = masked.reduceRegion(reducer=ee.Reducer.mean(), geometry=region,
                                scale=1000, maxPixels=1e9)
        return ee.Feature(None, {"date": img.date().format("YYYY-MM-dd"),
                                  "NDVI": s.get("NDVI"), "EVI": s.get("EVI")})
    feats = col.map(extract).getInfo()["features"]
    df = pd.DataFrame([{"date": f["properties"]["date"],
                        "NDVI": f["properties"]["NDVI"],
                        "EVI":  f["properties"]["EVI"]} for f in feats])
    df["date"] = pd.to_datetime(df["date"])
    df = df.set_index("date").sort_index() / 10000.0
    print(f"GEE: {len(df)} months downloaded ({start} to {end}).")
    return df


def make_synthetic_ndvi(times, drought_years):
    """Synthetic monthly cropland NDVI: seasonal + trend + drought suppression."""
    rng = np.random.default_rng(77); t = np.arange(len(times))
    mi  = np.array([d.month for d in times])
    seasonal = 0.10*np.sin(2*np.pi*(mi-0)/12) + 0.03*np.sin(4*np.pi*(mi-2)/12)
    base = 0.58 - 0.0003*t
    ddrop = np.zeros(len(times))
    for yr, amp, dur in [(2005,-0.07,16),(2010,-0.09,18),(2015,-0.12,22),(2021,-0.10,18)]:
        c_arr = np.where([d.year==yr and d.month==12 for d in times])[0]
        if len(c_arr):
            c = c_arr[0]; w = np.arange(max(0,c-dur//2), min(len(times),c+dur//2))
            ddrop[w] += amp*np.exp(-0.5*((w-c)/(dur/3))**2)
    ndvi = np.clip(base + seasonal + ddrop + rng.normal(0, 0.012, len(times)), 0.1, 0.95)
    df = pd.DataFrame({"NDVI": ndvi, "EVI": ndvi*0.88}, index=times)
    df.index.name = "date"
    return df


# Soybean futures ------------------------------------------------------------

def get_soybean_prices(start="2002-01-01", end="2026-01-01"):
    """Download CBOT soybean futures (ZS=F) monthly prices via yfinance."""
    try:
        s = yf.download("ZS=F", start=start, end=end,
                         interval="1mo", auto_adjust=True, progress=False)
        if s.empty: raise ValueError("empty response")
        c = s["Close"].squeeze()
        c.index = c.index.to_period("M").to_timestamp()
        c.name = "soybean_usdbu"
        print(f"yfinance: {len(c)} months of ZS=F prices.")
        return c
    except Exception as e:
        print(f"yfinance failed: {e}"); return None


def make_synthetic_soybean(times, drought_years):
    """Synthetic soybean price: long-term trend + lagged drought price spikes."""
    rng = np.random.default_rng(99); t = np.arange(len(times))
    base = 500 + 4.5*t; spikes = np.zeros(len(times))
    for yr, amp, lag in [(2005,120,4),(2010,200,5),(2015,160,5),(2021,220,4)]:
        tm = 9+lag; ty = yr+(tm-1)//12; tm = (tm-1)%12+1
        c_arr = np.where([d.year==ty and d.month==tm for d in times])[0]
        if len(c_arr):
            c = c_arr[0]; w = np.arange(max(0,c-6), min(len(times),c+6))
            spikes[w] += amp*np.exp(-0.5*((w-c)/3.5)**2)
    macro = rng.normal(0, 30, len(times)).cumsum()*0.5
    return pd.Series(np.maximum(base+spikes+macro, 200), index=times, name="soybean_centsbu")


# Time-series utilities ------------------------------------------------------

def to_monthly(s, agg="mean"):
    """Resample a pandas Series to monthly start frequency."""
    s = s.copy(); s.index = pd.DatetimeIndex(s.index)
    return s.resample("MS").mean() if agg == "mean" else s.resample("MS").last()


def adf_test(series, name):
    """Augmented Dickey-Fuller test. Prints result and returns True if stationary."""
    clean = series.dropna()
    stat, pval, *_ = adfuller(clean, autolag="AIC")
    status = "STATIONARY" if pval < 0.05 else "non-stationary"
    print(f"  {name:<28s}  ADF={stat:7.3f}  p={pval:.4f}  [{status}]")
    return pval < 0.05


# Statistical analysis -------------------------------------------------------

def cross_corr_bootstrap(x, y, max_lag, n_boot=1000, alpha=0.05, seed=0):
    """
    Cross-correlation at lags 0..max_lag with block-bootstrap 95% CIs.

    Positive lag k means x leads y by k months.
    Returns: lags, cc, lo, hi  (arrays of length max_lag+1).
    """
    rng = np.random.default_rng(seed)
    xy = pd.DataFrame({"x": x, "y": y}).dropna()
    xz = (xy["x"]-xy["x"].mean())/xy["x"].std()
    yz = (xy["y"]-xy["y"].mean())/xy["y"].std()
    n  = len(xz); lags = np.arange(0, max_lag+1)

    def _cc(a, b, lag):
        ne = min(len(a), len(b)) - lag
        if ne < 2: return np.nan
        return np.corrcoef(a[:ne], b[lag:lag+ne])[0,1]

    cc   = np.array([_cc(xz.values, yz.values, lag) for lag in lags])
    blk  = 6; nb = n // blk
    boot = np.zeros((n_boot, len(lags)))
    for b in range(n_boot):
        idx = rng.choice(nb, size=nb, replace=True)
        bi  = np.concatenate([np.arange(i*blk,(i+1)*blk) for i in idx])[:n]
        for j, lag in enumerate(lags):
            boot[b,j] = _cc(xz.values[bi], yz.values[bi], lag)
    lo = np.nanpercentile(boot, 100*alpha/2, axis=0)
    hi = np.nanpercentile(boot, 100*(1-alpha/2), axis=0)
    return lags, cc, lo, hi


def rolling_corr_at_lag(x, y, lag, window=60):
    """Rolling Pearson r between x shifted back by `lag` months and y."""
    df_r = pd.DataFrame({"x": x.shift(lag), "y": y}).dropna()
    return df_r["x"].rolling(window, min_periods=window//2).corr(df_r["y"])


def partial_corr_at_lag(x, y, z, lag):
    """
    Partial correlation r(x[t-lag], y[t] | z[t]).
    Returns (r, p-value).
    """
    df_p = pd.DataFrame({"x": x.shift(lag), "y": y, "z": z}).dropna()
    if len(df_p) < 20: return np.nan, np.nan
    bx = np.polyfit(df_p["z"], df_p["x"], 1)
    by = np.polyfit(df_p["z"], df_p["y"], 1)
    ex = df_p["x"] - np.polyval(bx, df_p["z"])
    ey = df_p["y"] - np.polyval(by, df_p["z"])
    return pearsonr(ex, ey)


def run_granger(df_in, target, predictor, maxlag=6, label=""):
    """
    Granger causality: does `predictor` Granger-cause `target`?
    Gap-interpolated months are excluded when 'Gap_interp' column is present.
    """
    data = df_in[[target, predictor]].dropna()
    if "Gap_interp" in df_in.columns:
        data = data[df_in.loc[data.index, "Gap_interp"] == 0]
    res = grangercausalitytests(data, maxlag=maxlag, verbose=False)
    rows = []
    for lag, r in res.items():
        f, p = r[0]["ssr_ftest"][:2]
        rows.append({"Lag": lag, "F": f, "p": p,
                     "sig": "**" if p<0.05 else ("*" if p<0.10 else "")})
    gc = pd.DataFrame(rows)
    print(f"Granger {label}  (min p={gc.p.min():.4f} at lag {gc.loc[gc.p.idxmin(),'Lag']})")
    for _, row in gc.iterrows():
        print(f"  lag {int(row.Lag)}: F={row.F:.2f}  p={row.p:.4f}  {row.sig}")
    return gc


def drought_event_study(df_in, yr, threshold_sigma=-1.0):
    """
    For a given drought year, identify:
      - TWS onset (first month below threshold_sigma)
      - Soybean price reaction (first month price z > +0.5 after TWS onset)
      - Lag in months between the two.
    """
    tws_z = (df_in["TWS_cm"] - df_in["TWS_cm"].mean()) / df_in["TWS_cm"].std()
    soy_z = (df_in["Soybean"] - df_in["Soybean"].rolling(24).mean()) / \
             df_in["Soybean"].rolling(24).std()
    mask = (df_in.index >= f"{yr-1}-01") & (df_in.index <= f"{yr+2}-12")
    win  = pd.DataFrame({"tws_z": tws_z[mask], "soy_z": soy_z[mask]})
    t_tws = win.index[win.tws_z < threshold_sigma]
    t_tws = t_tws[0] if len(t_tws) else None
    t_soy = None
    if t_tws:
        aft = win.loc[t_tws:]
        t_s = aft.index[aft.soy_z > 0.5]
        t_soy = t_s[0] if len(t_s) else None
    def _lag(t1, t2):
        if t1 and t2: return (t2.year-t1.year)*12 + (t2.month-t1.month)
        return None
    return {"year": yr, "tws_onset": t_tws, "price_reaction": t_soy,
            "lag_tws_soy": _lag(t_tws, t_soy)}


def illustrative_backtest(df_in, signal_col, price_col, best_lag, threshold_sigma=-1.0):
    """
    Long-only backtest: enter when TWS < -1 sigma, hold for `best_lag` months.

    Returns: (strategy_returns, bh_returns, strategy_cumulative,
              bh_cumulative, metrics_dict)

    Disclaimer: academic illustration only -- ignores transaction costs,
    slippage, and information delivery lags.
    """
    sig_z     = (df_in[signal_col] - df_in[signal_col].mean()) / df_in[signal_col].std()
    price_ret = np.log(df_in[price_col]).diff()
    position  = (sig_z < threshold_sigma).astype(float).shift(best_lag).fillna(0)
    if "Gap_interp" in df_in.columns:
        position[df_in["Gap_interp"] == 1] = 0
    strat_ret = position * price_ret
    bh_ret    = price_ret.copy()

    def sharpe(r, ann=12):
        r = r.dropna()
        return r.mean()*ann / (r.std()*np.sqrt(ann)) if r.std() > 0 else 0
    def mdd(c):
        roll = c.cummax(); return ((c - roll) / roll).min()

    sc = (1 + strat_ret.fillna(0)).cumprod()
    bc = (1 + bh_ret.fillna(0)).cumprod()
    return strat_ret, bh_ret, sc, bc, {
        "Strategy Sharpe"  : sharpe(strat_ret),
        "B&H Sharpe"       : sharpe(bh_ret),
        "Strategy Ann. Ret": strat_ret.mean()*12,
        "B&H Ann. Ret"     : bh_ret.mean()*12,
        "Strategy MDD"     : mdd(sc),
        "B&H MDD"          : mdd(bc),
        "Active months"    : int(position.sum()),
        "Entry lag"        : best_lag,
    }
