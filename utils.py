"""
utils.py — GRACE Soybean Analysis helper functions
===================================================
All reusable functions extracted from GRACE_Soybean_Analysis.ipynb.

Usage in notebook / Colab:
    from utils import (
        load_grace_mascons, make_synthetic_grace,
        load_era5, make_synthetic_era5,
        get_ndvi_gee, make_synthetic_ndvi,
        get_soybean_prices, make_synthetic_soybean,
        cross_corr_bootstrap, run_granger,
        drought_event_study, illustrative_backtest,
    )
"""

import numpy as np
import pandas as pd
import xarray as xr
import yfinance as yf
from statsmodels.tsa.stattools import grangercausalitytests

# ---------------------------------------------------------------------------
# GRACE
# ---------------------------------------------------------------------------

def load_grace_mascons(filepath):
    """Load JPL RL06M v03 NetCDF, normalise dims, apply hydrology gain factors."""
    ds = xr.open_dataset(filepath)
    rename_map = {}
    for dim in ds.dims:
        if 'lat'  in dim.lower(): rename_map[dim] = 'lat'
        if 'lon'  in dim.lower(): rename_map[dim] = 'lon'
        if 'time' in dim.lower(): rename_map[dim] = 'time'
    ds = ds.rename(rename_map)
    tws_var = next((v for v in ds.data_vars if 'lwe' in v.lower()), list(ds.data_vars)[0])
    tws = ds[tws_var]
    gf_var = next((v for v in ds.data_vars if 'scale' in v.lower() or 'gain' in v.lower()), None)
    if gf_var:
        tws = tws * ds[gf_var]
        print(f'  Gain factors applied: {gf_var}')
    if tws.lon.values.max() > 180:
        tws = tws.assign_coords(lon=(tws.lon % 180) - (tws.lon // 180)*180).sortby('lon')
    tws.attrs['units'] = 'cm'
    return tws


def make_synthetic_grace(lat_min, lat_max, lon_min, lon_max,
                          t_start, t_end, gap_start, gap_end, seed=42):
    """Physically plausible synthetic GRACE TWS: trend + seasonal + drought events."""
    rng  = np.random.default_rng(seed)
    times = pd.date_range(t_start, t_end, freq='MS')
    lats = np.arange(lat_min, lat_max + 0.5, 1.0)
    lons = np.arange(lon_min, lon_max + 0.5, 1.0)
    nt, nlat, nlon = len(times), len(lats), len(lons)
    t = np.arange(nt)

    trend_map  = np.outer(np.linspace(-0.12, -0.06, nlat)[::-1], np.linspace(0.8, 1.2, nlon))
    month_idx  = np.array([d.month for d in times])
    seasonal   = 8.0 * np.sin(2 * np.pi * (month_idx - 1) / 12)

    drought = np.zeros(nt)
    for yr, amp, dur in [(2005,-18,18),(2010,-22,20),(2015,-30,24),(2021,-25,18)]:
        c_arr = np.where([d.year==yr and d.month==9 for d in times])[0]
        if len(c_arr):
            c = c_arr[0]
            w = np.arange(max(0,c-dur//2), min(nt,c+dur//2))
            drought[w] += amp * np.exp(-0.5*((w-c)/(dur/4))**2)

    data = np.zeros((nt, nlat, nlon))
    for i in range(nt):
        data[i] = trend_map*t[i] + seasonal[i] + drought[i] + rng.normal(0,3.,(nlat,nlon))

    gap_mask = (times >= gap_start) & (times <= gap_end)
    data[gap_mask] = np.nan

    return xr.DataArray(data, coords={'time': times,'lat': lats,'lon': lons},
                         dims=['time','lat','lon'],
                         attrs={'units':'cm','note':'SYNTHETIC – replace with JPL RL06M'})


# ---------------------------------------------------------------------------
# ERA5
# ---------------------------------------------------------------------------

def load_era5(filepath):
    """Load ERA5-Land NetCDF and return area-weighted P−ET time series (mm/month)."""
    ds = xr.open_dataset(filepath)
    p_var  = next(v for v in ds.data_vars if 'tp'  in v.lower() or 'prec' in v.lower())
    et_var = next(v for v in ds.data_vars if 'pev' in v.lower() or v.lower() in ('e', 'et') or 'evap' in v.lower())
    P = ds[p_var]*1000; ET = ds[et_var]*1000
    PET = P - np.abs(ET)
    w = np.cos(np.deg2rad(ds.latitude.values))
    PET_mean = (PET * w[np.newaxis,:,np.newaxis]).sum(dim=['latitude','longitude']) / \
               (w.sum()*len(ds.longitude))
    return PET_mean.to_series()


def make_synthetic_era5(times, drought_years):
    """Generate synthetic P−ET series aligned to `times`."""
    rng = np.random.default_rng(123)
    t = np.arange(len(times))
    month_idx = np.array([d.month for d in times])
    seasonal  = 50 * np.sin(2*np.pi*(month_idx-1)/12)
    drought   = np.zeros(len(times))
    for yr, amp, dur in [(2005,-40,16),(2010,-50,18),(2015,-65,22),(2021,-55,16)]:
        c_arr = np.where([d.year==yr and d.month==8 for d in times])[0]
        if len(c_arr):
            c = c_arr[0]; w = np.arange(max(0,c-dur//2),min(len(times),c+dur//2))
            drought[w] += amp*np.exp(-0.5*((w-c)/(dur/4))**2)
    pet = seasonal + drought - 0.15*t + rng.normal(0,8,len(times))
    return pd.Series(pet, index=times, name='P_minus_ET_mm')


# ---------------------------------------------------------------------------
# MODIS NDVI (Google Earth Engine)
# ---------------------------------------------------------------------------

def get_ndvi_gee(lat_min, lat_max, lon_min, lon_max,
                  gee_project, start='2002-01-01', end=None):
    """
    Extract monthly cropland-masked NDVI from MODIS MOD13A3 v061 via GEE.

    Parameters
    ----------
    gee_project : str
        GEE project ID (e.g. 'your-gee-project').
    """
    import ee

    if end is None:
        end = pd.Timestamp('today').strftime('%Y-%m-%d')

    try:
        ee.Initialize(project=gee_project)
    except Exception:
        ee.Authenticate(auth_mode="localhost")
        ee.Initialize(project=gee_project)

    region = ee.Geometry.Rectangle([lon_min, lat_min, lon_max, lat_max])

    lc = (ee.ImageCollection('MODIS/061/MCD12Q1')
            .filterDate('2020-01-01', '2020-12-31')
            .first()
            .select('LC_Type1'))
    crop_mask = lc.eq(12).Or(lc.eq(14))

    collection = (ee.ImageCollection('MODIS/061/MOD13A3')
                    .filterDate(start, end)
                    .select(['NDVI', 'EVI']))

    def extract(img):
        masked = img.updateMask(crop_mask)
        stats  = masked.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=region,
            scale=1000, maxPixels=1e9)
        return ee.Feature(None, {
            'date': img.date().format('YYYY-MM-dd'),
            'NDVI': stats.get('NDVI'),
            'EVI' : stats.get('EVI'),
        })

    features = collection.map(extract).getInfo()['features']
    rows = [{'date': f['properties']['date'],
             'NDVI': f['properties']['NDVI'],
             'EVI' : f['properties']['EVI']} for f in features]
    ndvi_df = pd.DataFrame(rows)
    ndvi_df['date'] = pd.to_datetime(ndvi_df['date'])
    ndvi_df = ndvi_df.set_index('date').sort_index()
    ndvi_df /= 10000.0
    print(f'GEE: {len(ndvi_df)} months of NDVI downloaded ({start} to {end}).')
    return ndvi_df


def make_synthetic_ndvi(times, drought_years):
    """
    Synthetic monthly NDVI for Brazil agricultural belt.

    Key features:
    - Seasonal peak Jan-Feb (austral summer soybean canopy closure)
    - Drought years: suppressed peak NDVI (~0.05-0.10 units below mean)
    - 1-2 month lag behind TWS signal
    - Slight downward trend consistent with LULCC
    """
    rng = np.random.default_rng(77)
    t = np.arange(len(times))
    month_idx = np.array([d.month for d in times])

    seasonal = (0.10 * np.sin(2*np.pi*(month_idx - 0)/12) +
                0.03 * np.sin(4*np.pi*(month_idx - 2)/12))

    base = 0.58 - 0.0003 * t

    drought_ndvi = np.zeros(len(times))
    for yr, amp, dur in [(2005,-0.07,16),(2010,-0.09,18),(2015,-0.12,22),(2021,-0.10,18)]:
        c_arr = np.where([d.year==yr and d.month==12 for d in times])[0]
        if len(c_arr):
            c = c_arr[0]; w = np.arange(max(0,c-dur//2),min(len(times),c+dur//2))
            drought_ndvi[w] += amp*np.exp(-0.5*((w-c)/(dur/3))**2)

    noise = rng.normal(0, 0.012, len(times))
    ndvi  = base + seasonal + drought_ndvi + noise
    ndvi  = np.clip(ndvi, 0.1, 0.95)

    df = pd.DataFrame({'NDVI': ndvi, 'EVI': ndvi * 0.88}, index=times)
    df.index.name = 'date'
    return df


# ---------------------------------------------------------------------------
# Soybean futures
# ---------------------------------------------------------------------------

def get_soybean_prices(start='2002-01-01', end='2023-12-31'):
    """Download CBOT soybean futures (ZS=F) monthly prices via yfinance."""
    try:
        soy = yf.download('ZS=F', start=start, end=end,
                          interval='1mo', auto_adjust=True, progress=False)
        if soy.empty: raise ValueError('Empty response')
        close = soy['Close'].squeeze()
        close.index = close.index.to_period('M').to_timestamp()
        close.name = 'soybean_usdbu'
        print(f'yfinance: {len(close)} months of ZS=F prices.')
        return close
    except Exception as e:
        print(f'yfinance failed: {e}'); return None


def make_synthetic_soybean(times, drought_years):
    """Generate synthetic soybean price series with trend and drought spikes."""
    rng = np.random.default_rng(99)
    t = np.arange(len(times))
    base = 500 + 4.5*t
    spikes = np.zeros(len(times))
    for yr, amp, lag in [(2005,120,4),(2010,200,5),(2015,160,5),(2021,220,4)]:
        target_month = 9 + lag
        target_year  = yr + (target_month - 1) // 12
        target_month = (target_month - 1) % 12 + 1
        c_arr = np.where([d.year==target_year and d.month==target_month for d in times])[0]
        if len(c_arr):
            c = c_arr[0]; w = np.arange(max(0,c-6),min(len(times),c+6))
            spikes[w] += amp*np.exp(-0.5*((w-c)/3.5)**2)
    macro = rng.normal(0,30,len(times)).cumsum()*0.5
    return pd.Series(np.maximum(base+spikes+macro, 200), index=times, name='soybean_centsbu')


# ---------------------------------------------------------------------------
# Statistical analysis
# ---------------------------------------------------------------------------

def cross_corr_bootstrap(x, y, max_lag, n_boot=1000, alpha=0.05, seed=0):
    """Cross-correlation at lags [-max_lag, max_lag] with block-bootstrap 95% CIs."""
    rng = np.random.default_rng(seed)
    xy = pd.DataFrame({'x': x, 'y': y}).dropna()
    xz = (xy['x']-xy['x'].mean())/xy['x'].std()
    yz = (xy['y']-xy['y'].mean())/xy['y'].std()
    n  = len(xz)
    lags = np.arange(-max_lag, max_lag+1)

    def _cc(a, b, lag):
        nb = n - abs(lag)
        if lag < 0: return np.corrcoef(a[:nb], b[abs(lag):])[0,1]
        else:       return np.corrcoef(a[lag:], b[:nb])[0,1]

    cc = np.array([_cc(xz.values, yz.values, lag) for lag in lags])

    block = 6; nb = n // block
    boot  = np.zeros((n_boot, len(lags)))
    for b in range(n_boot):
        idx = rng.choice(nb, size=nb, replace=True)
        bi  = np.concatenate([np.arange(i*block,(i+1)*block) for i in idx])[:n]
        for j, lag in enumerate(lags):
            boot[b,j] = _cc(xz.values[bi], yz.values[bi], lag)

    lo = np.nanpercentile(boot, 100*alpha/2, axis=0)
    hi = np.nanpercentile(boot, 100*(1-alpha/2), axis=0)
    return lags, cc, lo, hi


def run_granger(df, target, predictor, maxlag=6, label=''):
    """
    Granger causality test: does `predictor` Granger-cause `target`?

    Parameters
    ----------
    df : pd.DataFrame
        Master DataFrame with a 'Gap_interp' column.
    target, predictor : str
        Column names in df.
    """
    data = df[[target, predictor]].dropna()
    data = data[df.loc[data.index,'Gap_interp']==0]
    result = grangercausalitytests(data, maxlag=maxlag, verbose=False)
    rows = []
    for lag, res in result.items():
        f, p = res[0]['ssr_ftest'][:2]
        rows.append({'Lag': lag, 'F': f, 'p': p,
                     'sig': '**' if p<0.05 else ('*' if p<0.10 else '')})
    gc = pd.DataFrame(rows)
    print(f'Granger {label}  (min p={gc.p.min():.4f} at lag {gc.loc[gc.p.idxmin(),"Lag"]})')
    for _, row in gc.iterrows():
        print(f'  lag {row.Lag:d}: F={row.F:.2f}  p={row.p:.4f}  {row.sig}')
    return gc


def drought_event_study(df, yr, threshold_sigma=-1.0):
    """
    Identify TWS onset, NDVI decline, and soybean price reaction for a drought year.

    Parameters
    ----------
    df : pd.DataFrame
        Master DataFrame with columns 'TWS_cm', 'NDVI', 'Soybean'.
    yr : int
        Drought year (e.g. 2005).
    """
    tws_z  = (df['TWS_cm']  - df['TWS_cm'].mean())  / df['TWS_cm'].std()
    ndvi_z = (df['NDVI']    - df['NDVI'].mean())     / df['NDVI'].std()
    soy_z  = (df['Soybean'] - df['Soybean'].rolling(24).mean()) / \
              df['Soybean'].rolling(24).std()

    mask = (df.index >= f'{yr-1}-01') & (df.index <= f'{yr+1}-12')
    win  = pd.DataFrame({'tws_z': tws_z[mask], 'ndvi_z': ndvi_z[mask], 'soy_z': soy_z[mask]})

    t_tws  = win.index[win.tws_z  < threshold_sigma]
    t_tws  = t_tws[0]  if len(t_tws) else None

    t_ndvi = None
    if t_tws is not None:
        after = win.loc[t_tws:]
        t_ndvi_arr = after.index[after.ndvi_z < -0.5]
        t_ndvi = t_ndvi_arr[0] if len(t_ndvi_arr) else None

    t_soy = None
    if t_tws is not None:
        after = win.loc[t_tws:]
        t_soy_arr = after.index[after.soy_z > 0.5]
        t_soy = t_soy_arr[0] if len(t_soy_arr) else None

    def _lag(t1, t2):
        if t1 and t2: return (t2.year-t1.year)*12 + (t2.month-t1.month)
        return None

    return {'year': yr, 'tws_onset': t_tws, 'ndvi_onset': t_ndvi, 'price_reaction': t_soy,
            'lag_tws_ndvi': _lag(t_tws, t_ndvi), 'lag_tws_soy': _lag(t_tws, t_soy)}


def illustrative_backtest(df, signal_col, price_col, best_lag,
                           threshold_sigma=-1.0, entry_lag=None):
    """
    Simple long-only backtest: go long when TWS < -1σ, shifted by best_lag months.

    Parameters
    ----------
    df : pd.DataFrame
        Master DataFrame with 'Gap_interp' column.
    signal_col : str
        Column to use as the drought signal (e.g. 'TWS_cm').
    price_col : str
        Column with soybean price (e.g. 'Soybean').
    best_lag : int
        Lag (months) from cross-correlation analysis.
    """
    if entry_lag is None: entry_lag = int(best_lag)
    sig_z    = (df[signal_col]-df[signal_col].mean())/df[signal_col].std()
    price_ret = np.log(df[price_col]).diff()
    position = (sig_z < threshold_sigma).astype(float).shift(entry_lag).fillna(0)
    position[df['Gap_interp']==1] = 0
    strat_ret = position * price_ret
    bh_ret    = price_ret.copy()

    def sharpe(r, ann=12): r=r.dropna(); return r.mean()*ann/(r.std()*np.sqrt(ann)) if r.std()>0 else 0
    def mdd(c): roll=c.cummax(); return ((c-roll)/roll).min()

    s_cum  = (1+strat_ret.fillna(0)).cumprod()
    bh_cum = (1+bh_ret.fillna(0)).cumprod()
    return strat_ret, bh_ret, s_cum, bh_cum, {
        'Strategy Sharpe': sharpe(strat_ret), 'B&H Sharpe': sharpe(bh_ret),
        'Strategy Ann. Ret': strat_ret.mean()*12, 'B&H Ann. Ret': bh_ret.mean()*12,
        'Strategy MDD': mdd(s_cum), 'B&H MDD': mdd(bh_cum),
        'Active months': int(position.sum()), 'Entry lag': entry_lag
    }
