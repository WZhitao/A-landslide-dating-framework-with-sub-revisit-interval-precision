"""A landslide dating framework integrating remote sensing imagery and natural triggers

Landslide dating from optical time-series data and natural triggers.

Workflow:

  1. periodic adaptive equation  - BIC selection over candidate periods
  2. signature anomaly score     - DeltaSSR segmentation fitting (long term)
                                   times trajectory prediction residuals (short term)
  3. trigger indicators          - accumulated rainfall, earthquake forcing,
                                   accumulated snowmelt
  4. trigger anomaly score       - high-value-biased Isolation Forest
  5. composite anomaly score     - signature score plus the trigger score,
                                   weighted by the signature anomaly window

Input: one CSV per landslide with the columns
    series, date, value, distance_m, magType, satellite
where `series` is one of NDVI / EVI / NBR / NDMI / NIRV for the vegetation
index, Tmask_GREEN / Tmask_SWIR for the Landsat greenness and SWIR masks,
and Rain / EQ / Snow / ERA5 for the trigger series. Vegetation index, Tmask
and trigger rows may be mixed in the same file.

Reported settings:
    TMASK_PROBABILITY = 0.99995, n_trees = 100, psi = number of observations,
    seed = 11
"""

import os
import itertools
import random

import numpy as np
import pandas as pd
import statsmodels.api as sm
from dateutil.relativedelta import relativedelta
from scipy.stats import norm

TMASK_PROBABILITY = 0.99995


def find_anomaly_window(dates_list, scores_list):
    """
    Identify the main anomaly period by locating the most prominent peak
    and marking a +/- 3-month window around that peak.
    Returns a daily binary list (1 = within window, 0 = outside).
    """
    start_date = min(dates_list)
    end_date = max(dates_list)
    total_days = (end_date - start_date).days + 1

    # Find index of the maximum score
    max_idx = np.nanargmax(scores_list)
    peak_date = dates_list[max_idx]

    # Build +/- 3-month window around the peak
    left_date = max(peak_date - relativedelta(months=3), start_date)
    right_date = min(peak_date + relativedelta(months=3), end_date)

    # Convert to day indices and mark the window
    start_idx = (left_date - start_date).days
    end_idx = (right_date - start_date).days

    daily_mask = np.zeros(total_days, dtype=int)
    daily_mask[start_idx:end_idx + 1] = 1

    return daily_mask.tolist()


def normalize_list(lst):
    lst = np.array(lst)
    if lst.size == 0:
        return np.array([])
    min_val = min(lst)
    max_val = max(lst)
    if min_val - max_val != 0:
        lst = np.array([(x - min_val) / (max_val - min_val) for x in lst])
        return lst
    else:
        return [0] * len(lst)


def calcEarthquake(date, distance, mag, mag_type, full_dates):
    """
    Compute daily maximum FEQ (earthquake factor)

    Parameters:
    - date: list of earthquake dates
    - distance: list of distances
    - mag: list of magnitudes
    - mag_type: list of magnitude types

    Returns:
    - list of dates
    - corresponding FEQ values (0.0 filled for missing)
    """
    df = pd.DataFrame({
        "date": pd.to_datetime(date),
        "distance": distance,
        "mag": mag,
        "mag_type": mag_type
    })

    df["mag_type"] = df["mag_type"].astype(str).str.lower()
    df = df[df["mag_type"] != "nan"]

    # ML = Mb + 0.4
    mb_mask = df["mag_type"] == "mb"
    df.loc[mb_mask, "mag"] += 0.4
    df.loc[mb_mask, "mag_type"] = "ml"

    # logM0
    def calculate_logM0(row):
        mag_type = row["mag_type"]
        mag = row["mag"]
        if mag_type.startswith(("mw", "ms")):
            return 1.5 * mag + 16.1
        elif mag_type.startswith("ml"):
            return 1.5 * mag + 16.0
        else:
            return np.nan

    df["logM0"] = df.apply(calculate_logM0, axis=1)
    df = df.dropna(subset=["logM0"])

    df["FEQ"] = df["logM0"] / df["distance"]
    df["D_values"] = 10 ** (0.75 * df["mag"] - 2.89) * 1000
    df = df[df["D_values"] > df["distance"]]

    df_max = df.groupby("date")["FEQ"].max().reset_index() if not df.empty else pd.DataFrame(columns=["date", "FEQ"])
    df_full = pd.DataFrame({"date": full_dates})
    df_result = df_full.merge(df_max, on="date", how="left").fillna(0.0)

    return df_result["date"].tolist(), df_result["FEQ"].tolist()


def calculateSSR(y_true, y_predictions):
    # compute sum of squared residuals (SSR)
    residuals = y_true - y_predictions
    SSR = np.sum(residuals ** 2)
    return SSR


def calcAccPrecipitation(precipitation_list, days=7, a=0.84):
    accprecipitation_list = []
    for idx in range(len(precipitation_list)):
        current_day = precipitation_list[idx]
        for day_idx in range(1, days + 1):
            if idx - day_idx >= 0:
                current_day += precipitation_list[idx - day_idx] * a ** day_idx
        accprecipitation_list.append(current_day)
    return accprecipitation_list


def read_landslide_csv(file_path):
    """
    Read the CSV at the given path and organize into a unified structure:
    VI_time / VI_value / VI_satellite / VI_date_numeric
    Rain_time / Rain_value
    EQ_time / EQ_value / EQ_distance
    Snow_time / Snow_value
    ERA5_time / ERA5_value
    Tmask_time / Tmask_GREEN / Tmask_SWIR
    """
    def most_common_satellite(series):
        if series.empty:
            return None
        return series.mode()[0] if not series.mode().empty else series.iloc[0]

    df = pd.read_csv(file_path)

    base_cols = {'series', 'date', 'value'}
    if not base_cols.issubset(df.columns):
        raise ValueError(f"Missing necessary columns. Please check the column names: {df.columns}")

    df['date'] = pd.to_datetime(df['date'], format='mixed', errors='coerce')

    result = {
        'VI_time': None, 'VI_value': None, 'VI_satellite': None, 'VI_date_numeric': None,
        'Rain_time': None, 'Rain_value': None,
        'EQ_time': None, 'EQ_value': None, 'EQ_distance': None, 'EQ_magType': None,
        'Snow_time': None, 'Snow_value': None,
        'ERA5_time': None, 'ERA5_value': None,
        'Tmask_time': None, 'Tmask_GREEN': None, 'Tmask_SWIR': None
    }

    # ========== Vegetation index (VI) ==========
    vi_types = ['NDVI', 'EVI', 'NBR', 'NDMI', 'NIRV']
    vi_df = df[df['series'].isin(vi_types)].copy()

    if vi_df.empty:
        raise ValueError("No vegetation index (VI) data was found")

    vi_df = (
        vi_df.groupby('date', as_index=False)
        .agg({'value': 'mean', 'satellite': most_common_satellite})
        .sort_values('date')
    )

    base_date = vi_df['date'].iloc[0]
    vi_df['data_numeric'] = (vi_df['date'] - base_date).dt.days

    result['VI_time'] = vi_df['date'].values
    result['VI_value'] = vi_df['value'].values
    result['VI_satellite'] = vi_df['satellite'].values
    result['VI_date_numeric'] = vi_df['data_numeric'].values

    # ========== triggers ==========
    rain_df = df[df['series'] == 'Rain'].sort_values('date')
    if not rain_df.empty:
        result['Rain_time'] = rain_df['date'].values
        result['Rain_value'] = rain_df['value'].values

    eq_df = df[df['series'] == 'EQ'].sort_values('date')
    if not eq_df.empty:
        result['EQ_time'] = eq_df['date'].values
        result['EQ_value'] = eq_df['value'].values
        result['EQ_distance'] = eq_df['distance_m'].values
        result['EQ_magType'] = eq_df['magType'].values

    snow_df = df[df['series'] == 'Snow'].sort_values('date')
    if not snow_df.empty:
        result['Snow_time'] = snow_df['date'].values
        result['Snow_value'] = snow_df['value'].values

    era_df = df[df['series'] == 'ERA5'].sort_values('date')
    if not era_df.empty:
        result['ERA5_time'] = era_df['date'].values
        result['ERA5_value'] = era_df['value'].values

    # ========== Tmask ==========
    # GREEN and SWIR arrive as their own rows of the long-format table,
    # not as columns, and are joined on date.
    tmask_green_df = df[df['series'] == 'Tmask_GREEN'][['date', 'value']].rename(
        columns={'value': 'GREEN'}
    )
    tmask_swir_df = df[df['series'] == 'Tmask_SWIR'][['date', 'value']].rename(
        columns={'value': 'SWIR'}
    )

    if (not tmask_green_df.empty) or (not tmask_swir_df.empty):
        tmask_df = pd.merge(
            tmask_green_df,
            tmask_swir_df,
            on='date',
            how='outer'
        ).sort_values('date')

        tmask_df = tmask_df.dropna(how='all', subset=['GREEN', 'SWIR'])

        if not tmask_df.empty:
            result['Tmask_time'] = tmask_df['date'].values
            result['Tmask_GREEN'] = tmask_df['GREEN'].values
            result['Tmask_SWIR'] = tmask_df['SWIR'].values

    # ========== array & datetime normalization ==========
    for key in result:
        if result[key] is not None:
            result[key] = np.array(result[key])

    date_keys = ['VI_time', 'Rain_time', 'EQ_time', 'Snow_time', 'ERA5_time', 'Tmask_time']
    for key in date_keys:
        if result.get(key) is not None:
            result[key] = pd.to_datetime(result[key], errors='coerce')

    # ========== Compute trigger time range ==========
    trigger_keys = ['Rain_time', 'EQ_time', 'Snow_time', 'ERA5_time']
    all_trigger_dates = []

    for k in trigger_keys:
        arr = result.get(k)
        if arr is None:
            continue
        dts = pd.to_datetime(arr, errors='coerce')
        dts = dts[~pd.isna(dts)]
        if len(dts) > 0:
            all_trigger_dates.extend(list(dts))

    if len(all_trigger_dates) == 0:
        result['trigger_time_range'] = pd.DatetimeIndex([])
    else:
        min_date = pd.to_datetime(min(all_trigger_dates)).normalize()
        max_date = pd.to_datetime(max(all_trigger_dates)).normalize()
        result['trigger_time_range'] = pd.date_range(start=min_date, end=max_date, freq='D')

    return result


def remove_tmask_outliers(data, probability=0.99999, keep_last_n=8, eps=1e-6):
    """
    Remove Tmask outliers (GREEN / SWIR) using harmonic fit residuals and
    normal-tail thresholding with robust sigma (MAD -> sigma).

    - Directionality preserved:
        green_resid = obs - fit  (positive => GREEN higher than model)
        swir_resid  = fit - obs  (positive => SWIR lower than model)
    - Deletion rule: remove if either band signals anomaly (OR).
    """

    # minimal completeness check
    if (data.get('Tmask_time') is None or
        data.get('Tmask_GREEN') is None or
        data.get('Tmask_SWIR') is None):
        return data

    df = pd.DataFrame({
        'DATE': pd.to_datetime(data['Tmask_time']),
        'GREEN': data['Tmask_GREEN'],
        'SWIR': data['Tmask_SWIR']
    }).dropna(subset=['GREEN', 'SWIR'], how='any')

    if len(df) < 10:
        return data

    # time design for harmonic fit
    date_all = df['DATE']
    T = 365.25
    N_years = ((date_all.max() - date_all.min()).days) / T
    if N_years <= 0:
        return data

    omega1 = 2 * np.pi / T
    omega2 = 2 * np.pi / (N_years * T)
    x = (date_all - date_all.min()).dt.days.values
    X = np.column_stack([
        np.ones_like(x),
        np.cos(omega1 * x),
        np.sin(omega1 * x),
        np.cos(omega2 * x),
        np.sin(omega2 * x)
    ])

    # robust sigma from MAD; fallback to std; final fallback eps to avoid zero
    def robust_sigma(arr, eps_local=eps):
        arr = np.asarray(arr)
        if arr.size == 0:
            return float(eps_local)
        mad = np.median(np.abs(arr - np.median(arr)))
        if mad > 0:
            sigma = mad / 0.6745
            return float(sigma)
        std = float(np.std(arr, ddof=1))
        return std if std > 0 else float(eps_local)

    z_p = float(norm.ppf(probability))

    # GREEN (direction: obs - fit)
    green_vals = df['GREEN'].values
    green_fitted = sm.RLM(green_vals, X, M=sm.robust.norms.HuberT()).fit().predict(X)
    green_resid = green_vals - green_fitted
    sigma_g = robust_sigma(green_resid)
    thr_green = z_p * sigma_g
    # anomaly when residual exceeds thr in positive direction
    anom_green = green_resid > thr_green
    # ensure keep last points
    if keep_last_n > 0:
        anom_green[-keep_last_n:] = False

    # SWIR (direction: fit - obs)
    swir_vals = df['SWIR'].values
    swir_fitted = sm.RLM(swir_vals, X, M=sm.robust.norms.HuberT()).fit().predict(X)
    swir_resid = swir_fitted - swir_vals
    sigma_s = robust_sigma(swir_resid)
    thr_swir = z_p * sigma_s
    anom_swir = swir_resid > thr_swir
    if keep_last_n > 0:
        anom_swir[-keep_last_n:] = False

    # deletion rule: remove if either band anomalous (OR)
    anom_any = anom_green | anom_swir
    # keep = not anomalous
    keep_mask = ~anom_any

    df_filt = df.loc[keep_mask].copy()
    data['Tmask_time'] = df_filt['DATE'].values
    data['Tmask_GREEN'] = df_filt['GREEN'].values
    data['Tmask_SWIR'] = df_filt['SWIR'].values

    # masked dates are removed from the vegetation index series as well
    if data.get('VI_time') is not None:
        anomalous_dates = df.loc[anom_any, 'DATE'].values
        keep_vi_mask = ~np.isin(data['VI_time'], anomalous_dates)
        data['VI_time'] = data['VI_time'][keep_vi_mask]
        data['VI_value'] = data['VI_value'][keep_vi_mask]
        if data.get('VI_satellite') is not None:
            data['VI_satellite'] = data['VI_satellite'][keep_vi_mask]
        if data.get('VI_date_numeric') is not None:
            data['VI_date_numeric'] = data['VI_date_numeric'][keep_vi_mask]

    # return modified data; other triggers unchanged
    return data


def calc_triggers(data):
    """
     Automatically compute secondary indicators for rain, earthquake, and snow,
     and add the results to `data`. Skip missing items.
     """
    data = data.copy()
    if data["Rain_value"] is not None:
        df_r = pd.DataFrame({
            "date": pd.to_datetime(data["Rain_time"]),
            "precipitation": np.array(data["Rain_value"], dtype=float)
        })
        df_daily_r = df_r.groupby("date", as_index=False)["precipitation"].sum()
        accPrecipitation = normalize_list(calcAccPrecipitation(df_daily_r["precipitation"].tolist()))
        data["Acc_Precipitation"] = accPrecipitation

    if "EQ_value" in data and data["EQ_value"] is not None:
        FEQ_date, FEQ = calcEarthquake(
            pd.to_datetime(data["EQ_time"]),
            np.array(data["EQ_distance"], dtype=float),
            np.array(data["EQ_value"], dtype=float),
            data["EQ_magType"], data['trigger_time_range']
        )
        data["FEQ"] = np.array(normalize_list(FEQ))

    if "Snow_value" in data and data["Snow_value"] is not None and "ERA5_value" in data:
        df_s = pd.DataFrame({
            "date": pd.to_datetime(data["Snow_time"]),
            "mean_snow": np.array(data["Snow_value"], dtype=float)
        })
        df_lst = pd.DataFrame({
            "date": pd.to_datetime(data["ERA5_time"]),
            "skin_temperature_C": np.array(data["ERA5_value"], dtype=float)
        })

        df_s = df_s.resample("D", on="date")["mean_snow"].mean().reindex(data['trigger_time_range']).ffill().bfill()
        df_lst = (
            df_lst[df_lst["skin_temperature_C"] > 0]
            .resample("D", on="date")["skin_temperature_C"]
            .sum().reindex(data['trigger_time_range'], fill_value=0)
        )

        snowmelt = df_s.values.flatten() * df_lst.values.flatten()
        accSnowmelt = normalize_list(calcAccPrecipitation(snowmelt))
        data["Acc_Snowmelt"] = accSnowmelt

    return data


def select_model_structure(x, y, candidate_periods=[365.25, 182.625, 120], max_seasonal_terms=3):
        """
        Choose the best combination of seasonal periods.

        Parameters:
          x: time (1D array)
          y: observations (1D array)
          candidate_periods: list of candidate periods
          max_seasonal_terms: max number of seasonal terms
        Returns:
          structure_params: dict with model structure {'periods'}
                            empty list means no seasonal terms
        """
        x = np.asarray(x).flatten()
        y = np.asarray(y).flatten()

        best_bic = np.inf
        best_periods = []

        for n_terms in range(0, max_seasonal_terms + 1):
            if n_terms == 0:
                X = sm.add_constant(x)
                model = sm.OLS(y, X).fit()
                bic = model.bic
                if bic < best_bic:
                    best_bic = bic
                    best_periods = []
                continue

            for periods in itertools.product(candidate_periods, repeat=n_terms):
                X = sm.add_constant(x)

                for T in periods:
                    theta = 2 * np.pi * x / T
                    X = np.column_stack([X, np.cos(theta), np.sin(theta)])

                model = sm.OLS(y, X).fit()
                bic = model.bic

                if bic < best_bic:
                    best_bic = bic
                    best_periods = list(periods)
        return {'periods': best_periods}


def fit_with_structure(x_new, y_new, structure_params, x_predict, remove_trend=False):
    """
       Fit data according to structure parameters.

       Parameters:
         x_new: time (1D array)
         y_new: observations (1D array)
         structure_params: dict {'periods'}
         x_predict: array of x values to predict

       Returns:
         dict {'params', 'fitted_values', 'predict_values'}
       """
    x_new = np.asarray(x_new).flatten()
    y_new = np.asarray(y_new).flatten()
    periods = structure_params['periods']
    x_predict = np.asarray(x_predict).flatten()

    X_train = sm.add_constant(np.zeros_like(x_new))
    if not remove_trend:
        X_train = np.column_stack([X_train, x_new])
    for T in periods:
        theta = 2 * np.pi * x_new / T
        X_train = np.column_stack([X_train, np.cos(theta), np.sin(theta)])

    model = sm.OLS(y_new, X_train).fit()
    params = model.params
    X_predict = sm.add_constant(np.zeros_like(x_predict))
    if not remove_trend:
        X_predict = np.column_stack([X_predict, x_predict])
    for T in periods:
        theta = 2 * np.pi * x_predict / T
        X_predict = np.column_stack([X_predict, np.cos(theta), np.sin(theta)])

    return {
        'params': params,
        'fitted_values': model.predict(X_train),  #
        'predict_values': np.dot(X_predict, params)
    }


def calculate_initialization_length(structure_params, date_list):
    periods_length = len(structure_params['periods'])
    base_length = 2 + 2 * periods_length
    start = 3 * base_length
    end = base_length
    while start < len(date_list) and (date_list[start - 1] - date_list[0]).days < 365.25:
        start += 1
    return start, end


def trajectory_prediction(vi_list, data_list_numeric, structure_params, start, end, T=3):
    #Compute prediction residuals for VI time series.
    delta_afterT_T_list = []

    for idx in range(start, len(vi_list) - end):
        data_list_numeric_pre = data_list_numeric[0:idx]
        vi_list_pre = vi_list[0:idx]
        fitted_list = fit_with_structure(data_list_numeric_pre, vi_list_pre, structure_params, data_list_numeric)['predict_values']        # --- mean error over the preceding window ---
        delta_T = fitted_list[idx] - vi_list[idx]
        window_len = min(T, idx-start)
        if window_len > 0:
            start_w = idx - window_len
            deltas_pre = fitted_list[start_w:idx] - np.array(vi_list[start_w:idx])
            delta_T = np.mean(deltas_pre)
        else:
            delta_T = 0

        delta_afterT = (fitted_list[idx: idx + T] - vi_list[idx:idx + T])
        delta_afterT_T = np.mean(delta_afterT) - delta_T

        delta_afterT_T_list.append(delta_afterT_T)

    return normalize_list(delta_afterT_T_list)


def segmentation_fitting(vi_list, data_list_numeric, structure_params, start, end):
    #Compute the change in sum of squared residuals (dSSR) for VI time series.

    SSR_list_pre = []
    SSR_list_post = []
    SSR_list_all = []
    for idx in range(start, len(vi_list) - end):
        vi_list_pre = vi_list[0:idx]
        data_list_numeric_pre = data_list_numeric[0:idx]

        fitted_pre = fit_with_structure(
            data_list_numeric_pre, vi_list_pre, structure_params, data_list_numeric_pre, remove_trend=True
        )["fitted_values"]
        SSR_pre = calculateSSR(vi_list_pre, fitted_pre)
        SSR_list_pre.append(SSR_pre)

        vi_list_post = vi_list[idx:]
        data_list_numeric_post = data_list_numeric[idx:]

        fitted_post = fit_with_structure(
            data_list_numeric_post, vi_list_post, structure_params, data_list_numeric_post, remove_trend=False
        )["fitted_values"]
        SSR_post = calculateSSR(vi_list_post, fitted_post)
        SSR_list_post.append(SSR_post)

        fitted_all = fit_with_structure(
            data_list_numeric, vi_list, structure_params, data_list_numeric
        )["fitted_values"]
        SSR_all = calculateSSR(vi_list, fitted_all)
        SSR_list_all.append(SSR_all)

    deltaSSR_list = normalize_list(
        np.array(SSR_list_all) - np.array(SSR_list_post) - np.array(SSR_list_pre)
    )

    return np.array(deltaSSR_list)


def sequential_backfill(score_dates, signature_score, target_scope):
    score_series = pd.Series(signature_score, index=pd.to_datetime(score_dates))
    full_series = pd.Series(index=pd.to_datetime(target_scope))

    sequential_score = score_series.reindex(full_series.index, method='bfill')  # backward fill
    if pd.isna(sequential_score.iloc[0]):
        sequential_score.iloc[0] = score_series.iloc[0]

    return sequential_score.values


def isolation_forest_scores(X, n_trees=100, psi=None, seed=None, normalize=True, return_forest=False):
    # Compute anomaly scores Isolation Forest .

    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    n, d = X.shape

    if psi is None:
        psi = n
    psi = int(min(max(1, psi), n))

    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)

    def _c(m):
        if m <= 1:
            return 0.0
        return 2.0 * (np.log(m - 1) + 0.5772156649) - 2.0 * (m - 1) / m

    def _iTree(X_sub, depth, max_depth):
        n_sub, d_sub = X_sub.shape
        if depth >= max_depth or n_sub <= 1:
            return [0, n_sub]
        q = np.random.randint(d_sub)
        lo, hi = X_sub[:, q].min(), X_sub[:, q].max()
        beta = random.betavariate(5, 1)
        p = beta * (hi - lo) + lo
        left = _iTree(X_sub[X_sub[:, q] < p], depth + 1, max_depth)
        right = _iTree(X_sub[X_sub[:, q] >= p], depth + 1, max_depth)
        return [1, left, right, q, p]

    def _PathLength(x, T, depth=0):
        if T[0] == 0:
            return depth + _c(T[1])
        _, left, right, q, p = T
        if x[q] < p:
            return _PathLength(x, left, depth + 1)
        else:
            return _PathLength(x, right, depth + 1)

    max_depth = int(np.ceil(np.log2(psi))) if psi > 0 else 0
    forest = []
    for _ in range(n_trees):
        idx = np.random.choice(n, psi, replace=False)
        T = _iTree(X[idx], 0, max_depth)
        forest.append(T)

    c_psi = _c(psi)
    if c_psi == 0.0:
        raw_scores = np.zeros(n, dtype=float)
    else:
        raw_scores = np.empty(n, dtype=float)
        for i in range(n):
            x = X[i]
            avg_path = np.mean([_PathLength(x, T, 0) for T in forest])
            raw_scores[i] = 2.0 ** (-avg_path / c_psi)

    if normalize:
        mn = np.nanmin(raw_scores)
        mx = np.nanmax(raw_scores)
        if np.isclose(mx, mn):
            scores = np.zeros_like(raw_scores)
        else:
            scores = (raw_scores - mn) / (mx - mn)
    else:
        scores = raw_scores

    if return_forest:
        return scores, forest
    return scores


def calc_comprehensive_anomaly(filled_signature_anomaly, trigger_anomaly, signature_anomaly_window, target_scope, vi_time, start, end, ):
    masked_trigger_anomaly =  trigger_anomaly[np.array((target_scope >= vi_time[start])
                                                       & (target_scope <= vi_time[len(vi_time) - end - 1]), dtype=bool)]
    weighted_trigger_anomaly = masked_trigger_anomaly * filled_signature_anomaly * np.array(signature_anomaly_window)
    comprehensive_anomaly = filled_signature_anomaly + weighted_trigger_anomaly
    return comprehensive_anomaly


if __name__ == "__main__":
    # read file: one folder of per-landslide CSVs, see the module docstring
    folder_path = r"PATH"
    data = read_landslide_csv(os.path.join(folder_path, 'slide1.csv'))
    data = remove_tmask_outliers(data, probability=TMASK_PROBABILITY)

    data = calc_triggers(data)
    '''Step 1 signature anomaly detection'''
    #1.1 Construction of the periodic adaptive equation
    structure_params = select_model_structure(data['VI_date_numeric'], data['VI_value'])
    print(structure_params)
    #1.2 segmentation fitting & trajectory prediction
    start, end = calculate_initialization_length(structure_params, data['VI_time'])
    print(start,end)
    long_term_anomaly = segmentation_fitting(data['VI_value'], data['VI_date_numeric'], structure_params, start, end)
    short_term_anomaly = trajectory_prediction(data['VI_value'], data['VI_date_numeric'], structure_params, start, end)
    signature_anomaly = normalize_list(long_term_anomaly * short_term_anomaly)

    startTime_signatureAnomaly = start
    endTime_signatureAnomaly = len(data['VI_time']) - end
    penultTime_signatureAnomaly = len(data['VI_time']) - end - 1

    signature_anomaly_window = np.array(find_anomaly_window(data['VI_time'][startTime_signatureAnomaly:endTime_signatureAnomaly], signature_anomaly), dtype=bool)
    filled_signature_anomaly = normalize_list(
        sequential_backfill(data['VI_time'][startTime_signatureAnomaly:endTime_signatureAnomaly], signature_anomaly,
                                                   pd.date_range(data['VI_time'][startTime_signatureAnomaly], data['VI_time'][penultTime_signatureAnomaly]))
    )
    signature_result = data['VI_time'][startTime_signatureAnomaly:endTime_signatureAnomaly][np.argmax(signature_anomaly)]
    print(signature_result)

    '''Step 2 trigger anomaly detection'''
    # 2.1 calculation of landslide trigger indicators
    trigger_vars = [v for v in ['FEQ', 'Acc_Precipitation', 'Acc_Snowmelt'] if v in data and len(data[v]) > 0]
    print(trigger_vars)
    if trigger_vars:  # only computed when at least one trigger series is present
        trigger_stack = np.column_stack([data[v] for v in trigger_vars])
        scores = isolation_forest_scores(trigger_stack, n_trees=100, psi=len(trigger_stack), seed=11, normalize=True)

        # 2.2 calculation of comprehensive anomaly
        comprehensive_anomaly = calc_comprehensive_anomaly(
            filled_signature_anomaly,
            trigger_anomaly=scores,
            signature_anomaly_window=signature_anomaly_window,
            target_scope=data['trigger_time_range'],
            vi_time=data['VI_time'],
            start=start,
            end=end
        )
        comprehensive_result = pd.date_range(
            data['VI_time'][startTime_signatureAnomaly],
            data['VI_time'][penultTime_signatureAnomaly]
        )[np.argmax(comprehensive_anomaly)]
        print(comprehensive_result)
    else:
        print("no trigger factor, skipping the comprehensive anomaly")
