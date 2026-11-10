
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TRÓPICO | Global Atmospheric Monitor
NOAA CPC CORe daily atmospheric anomalies and original wind vectors.

Author: Alejandro Jaramillo-Moreno
Instituto de Ciencias de la Atmósfera y Cambio Climático, UNAM

Outputs:
    data/latest/{11 atmospheric variables}.json  (anomalies)
    data/latest/wind850.json                      (original U/V)
    data/latest/wind250.json                      (original U/V)
    data/latest/manifest.json
"""

from calendar import monthrange
from datetime import datetime, timedelta, timezone
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import os
import re

import numpy as np
import requests
import xarray as xr

from eccodes import (
    codes_get, codes_get_array, codes_get_values,
    codes_grib_new_from_file, codes_release,
)
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# ==========================================================
# 1. CONFIGURATION
# ==========================================================

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "data" / "latest"

FOLDER_ID = "1iW30G3cAaIZceu2geZ0E9kg26Wlx3zE5"
DRIVE_URL = "https://www.googleapis.com/drive/v3/files"
NOAA_URL = "https://storage.googleapis.com/noaa-nws-ncep-core/grib/day"

TOKEN = os.environ.get("GOOGLE_ACCESS_TOKEN", "")
LOOKBACK_DAYS = 10

if not TOKEN:
    raise RuntimeError("GOOGLE_ACCESS_TOKEN is missing")

# Each entry: NOAA product group, GRIB parameter, level, time descriptor.
SPECS = {
    "t2m":      ("flx", "TMP",   r"^2\s*m\s*above\s*(ground|gnd)$", r"\d+@3\s*hour\s*ave\(anl\)"),
    "prate":    ("flx", "PRATE", r"^(surface|sfc)$", None),
    "pwat":     ("flx", "PWAT",  r"(atmos\s+col|entire\s+atmosphere)", None),
    "u850":     ("pgb", "UGRD",  r"^850\s*(mb|hPa)$", None),
    "v850":     ("pgb", "VGRD",  r"^850\s*(mb|hPa)$", None),
    "u250":     ("pgb", "UGRD",  r"^250\s*(mb|hPa)$", None),
    "v250":     ("pgb", "VGRD",  r"^250\s*(mb|hPa)$", None),
    "z500":     ("pgb", "HGT",   r"^500\s*(mb|hPa)$", None),
    "omega500": ("pgb", "VVEL",  r"^500\s*(mb|hPa)$", None),
}

DERIVED = {
    "wspd850": ("u850", "v850"),
    "wspd250": ("u250", "v250"),
}

# (display units, conversion from NOAA units, JSON decimals)
DISPLAY = {
    "t2m":      ("degC",      1.0,     2),
    "prate":    ("mm day-1",  86400.0, 2),
    "pwat":     ("mm",        1.0,     2),
    "u850":     ("m s-1",     1.0,     2),
    "v850":     ("m s-1",     1.0,     2),
    "u250":     ("m s-1",     1.0,     2),
    "v250":     ("m s-1",     1.0,     2),
    "wspd850":  ("m s-1",     1.0,     2),
    "wspd250":  ("m s-1",     1.0,     2),
    "z500":     ("gpm",       1.0,     2),
    "omega500": ("Pa s-1",    1.0,     5),
}

WIND_LEVELS = (850, 250)
WIND_DECIMALS = 2  # Visualization precision; fields remain float32 in memory.


# ==========================================================
# 2. HTTP SESSION
# ==========================================================

session = requests.Session()
retry = Retry(
    total=4,
    backoff_factor=1,
    status_forcelist=[429, 500, 502, 503, 504],
)
session.mount("https://", HTTPAdapter(max_retries=retry))
session.headers.update({
    "User-Agent": "TROPICO-Global-Atmospheric-Monitor/1.0"
})


# ==========================================================
# 3. NOAA DAILY INVENTORIES
# ==========================================================

def inventory_for_day(group, day):
    stem = (
        f"{NOAA_URL}/{group}/"
        f"{day.year}/{day.month:02d}/{group}.{day:%Y%m%d}"
    )
    response = session.get(stem + ".idx", timeout=(20, 90))

    if response.status_code == 404:
        return None

    response.raise_for_status()
    records = []

    for line in response.text.splitlines():
        parts = line.split(":")
        if len(parts) < 5:
            continue
        records.append({
            "offset": int(parts[1]),
            "parameter": parts[3].strip(),
            "level": parts[4].strip(),
            "timing": parts[5].strip() if len(parts) > 5 else "",
        })

    if not records or any(
        nxt["offset"] <= prev["offset"]
        for prev, nxt in zip(records, records[1:])
    ):
        raise ValueError(f"Invalid NOAA {group} inventory for {day}")

    return stem + ".grb", records


def match_record(records, name):
    _, parameter, level_pattern, time_pattern = SPECS[name]
    matches = [
        i for i, record in enumerate(records)
        if record["parameter"] == parameter
        and re.search(level_pattern, record["level"], re.IGNORECASE)
        and (
            time_pattern is None
            or re.search(time_pattern, record["timing"], re.IGNORECASE)
        )
    ]

    if len(matches) != 1:
        raise ValueError(f"Expected one GRIB message for {name}; got {len(matches)}")

    i = matches[0]
    start = records[i]["offset"]
    end = records[i + 1]["offset"] - 1 if i + 1 < len(records) else None
    return start, end


def latest_complete_day():
    today = datetime.now(timezone.utc).date()

    for lag in range(1, LOOKBACK_DAYS + 1):
        day = today - timedelta(days=lag)
        inventories = {
            group: inventory_for_day(group, day)
            for group in ("flx", "pgb")
        }

        if any(item is None for item in inventories.values()):
            print(f"{day}: FLX/PGB incomplete")
            continue

        try:
            selections = {
                name: match_record(inventories[spec[0]][1], name)
                for name, spec in SPECS.items()
            }
        except ValueError as exc:
            print(f"{day}: {exc}")
            continue

        print(f"Latest common FLX/PGB date: {day}")
        return day, inventories, selections

    raise RuntimeError(
        f"No common complete FLX/PGB date in the last {LOOKBACK_DAYS} days"
    )


# ==========================================================
# 4. GRIB DOWNLOAD AND DECODING
# ==========================================================

def download_grib(url, start, end, path):
    byte_range = f"bytes={start}-{end}" if end is not None else f"bytes={start}-"

    with session.get(
        url,
        headers={"Range": byte_range},
        stream=True,
        timeout=(20, 120),
    ) as response:
        if response.status_code != 206:
            raise RuntimeError(f"Expected HTTP 206; got {response.status_code}: {url}")
        with path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)

    if path.stat().st_size < 20:
        raise ValueError(f"Empty GRIB message: {url}")


def decode_grib(path):
    with path.open("rb") as handle:
        gid = codes_grib_new_from_file(handle)
        if gid is None:
            raise ValueError("Could not decode NOAA GRIB message")
        try:
            ni = int(codes_get(gid, "Ni"))
            nj = int(codes_get(gid, "Nj"))
            units = str(codes_get(gid, "units"))
            lat2d = codes_get_array(gid, "latitudes").reshape(nj, ni)
            lon2d = codes_get_array(gid, "longitudes").reshape(nj, ni) % 360
            values = codes_get_values(gid).reshape(nj, ni)
        finally:
            codes_release(gid)

    if (nj, ni) != (256, 512):
        raise ValueError(f"Unexpected NOAA grid: {(nj, ni)}")
    if not np.allclose(lat2d, lat2d[:, :1], atol=1e-4):
        raise ValueError("NOAA latitude grid is not separable")
    if not np.allclose(lon2d, lon2d[:1, :], atol=1e-4):
        raise ValueError("NOAA longitude grid is not separable")

    values = np.where(
        np.isfinite(values) & (np.abs(values) < 1e15),
        values, np.nan,
    ).astype("float32")

    field = xr.DataArray(
        values,
        dims=("lat", "lon"),
        coords={"lat": lat2d[:, 0], "lon": lon2d[0, :]},
    ).sortby("lat").sortby("lon")

    if not np.isfinite(field.values).any():
        raise ValueError("NOAA field contains no finite values")

    return field.values, field.lat.values, field.lon.values, units


# ==========================================================
# 5. GOOGLE DRIVE CLIMATOLOGIES
# ==========================================================

def drive_child(parent, name):
    response = session.get(
        DRIVE_URL,
        params={
            "q": f"'{parent}' in parents and name = '{name}' and trashed = false",
            "fields": "nextPageToken,files(id,name,mimeType)",
            "pageSize": 100,
        },
        headers={"Authorization": f"Bearer {TOKEN}"},
        timeout=(20, 90),
    )
    response.raise_for_status()
    result = response.json()
    files = result.get("files", [])

    if result.get("nextPageToken") or len(files) != 1:
        raise ValueError(f"Expected one Google Drive child {name!r}; got {len(files)}")
    return files[0]["id"]


def download_drive_file(file_id, path):
    with session.get(
        f"{DRIVE_URL}/{file_id}",
        params={"alt": "media"},
        headers={"Authorization": f"Bearer {TOKEN}"},
        stream=True,
        timeout=(20, 120),
    ) as response:
        response.raise_for_status()
        with path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)

    if path.stat().st_size == 0:
        raise ValueError(f"Empty Drive file: {file_id}")


def read_reference(path, name, month, day):
    with xr.open_dataset(path, decode_times=False) as ds:
        if name not in ds or "lat" not in ds.coords or "lon" not in ds.coords:
            raise ValueError(f"Invalid climatology structure for {name}")
        da = ds[name]
        day_dims = [dim for dim in da.dims if dim not in ("lat", "lon")]
        expected_days = monthrange(2001, month)[1]

        if len(day_dims) != 1 or da.sizes[day_dims[0]] != expected_days:
            raise ValueError(f"Unexpected daily calendar for {name}/{month:02d}.nc")

        field = da.isel({day_dims[0]: day - 1}).transpose("lat", "lon")
        values = np.asarray(field.values, dtype="float32")

        if values.shape != (256, 512) or not np.isfinite(values).any():
            raise ValueError(f"Invalid climatology field for {name}")

        return values, ds.lat.values, ds.lon.values, str(da.attrs.get("units", ""))


# ==========================================================
# 6. FILE VALIDATION / JSON
# ==========================================================

def validate_grid(lat, lon, reference_lat, reference_lon):
    if not (
        np.allclose(lat, reference_lat, rtol=0, atol=1e-5)
        and np.allclose(lon, reference_lon, rtol=0, atol=1e-5)
    ):
        raise ValueError("NOAA and reference grids differ")


def json_values(values, digits):
    return [
        [round(float(x), digits) if np.isfinite(x) else None for x in row]
        for row in values
    ]


def write_json(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")

    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(content, handle, separators=(",", ":"), allow_nan=False)

    temporary.replace(path)


# ==========================================================
# 7. PROCESS ALL DATA FOR THE SAME DAY
# ==========================================================

print("=" * 62)
print("TRÓPICO — GLOBAL ATMOSPHERIC MONITOR")
print("NOAA CPC CORe — DAILY ANOMALIES + ORIGINAL WIND VECTORS")
print("=" * 62)

date, inventories, selections = latest_complete_day()

with TemporaryDirectory() as temporary_directory:
    work = Path(temporary_directory)
    observed = {}
    source_units = {}
    lat = None
    lon = None

    # Download the same nine original fields as before.
    for name, spec in SPECS.items():
        group = spec[0]
        url, _ = inventories[group]
        start, end = selections[name]
        grib_path = work / f"{name}.grb"

        download_grib(url, start, end, grib_path)
        values, field_lat, field_lon, units = decode_grib(grib_path)

        if lat is None:
            lat, lon = field_lat.copy(), field_lon.copy()
        else:
            validate_grid(lat, lon, field_lat, field_lon)

        observed[name] = values
        source_units[name] = units
        print(f"NOAA {name:9s}: {grib_path.stat().st_size:>9,} bytes | {units}")
        grib_path.unlink()

    # Same definition as in our existing pipeline:
    # magnitude of the daily-mean horizontal wind vector.
    for name, (u_name, v_name) in DERIVED.items():
        if source_units[u_name] != source_units[v_name]:
            raise ValueError(f"Inconsistent U/V units for {name}")
        observed[name] = np.hypot(
            observed[u_name], observed[v_name]
        ).astype("float32")
        source_units[name] = source_units[u_name]

    latitude = [round(float(value), 6) for value in lat]
    longitude = [round(float(value), 6) for value in lon]

    # ======================================================
    # 8. NEW: ORIGINAL WIND VECTOR PRODUCTS
    # ======================================================

    wind_outputs = []

    for level in WIND_LEVELS:
        u = observed[f"u{level}"]
        v = observed[f"v{level}"]
        name = f"wind{level}"

        if u.shape != (256, 512) or v.shape != (256, 512):
            raise ValueError(f"Invalid wind grid for {name}")

        good = np.isfinite(u) & np.isfinite(v)
        if not good.any():
            raise ValueError(f"No valid wind vectors for {name}")

        # Ensure the new vectors represent exactly the same
        # original daily-mean winds used by the speed product.
        if not np.allclose(
            np.hypot(u, v), observed[f"wspd{level}"],
            equal_nan=True, rtol=1e-5, atol=1e-5,
        ):
            raise ValueError(f"Wind speed consistency check failed for {name}")

        write_json(
            OUTPUT_DIR / f"{name}.json",
            {
                "variable": name,
                "date": date.isoformat(),
                "level_hpa": level,
                "units": "m/s",
                "source": "NOAA CPC CORe daily ensemble mean",
                "field_type": "original_daily_mean_wind",
                "description": "Original daily-mean horizontal wind components, not anomalies",
                "components": {
                    "u": "eastward",
                    "v": "northward",
                },
                "grid": "native CORe Gaussian grid",
                "shape": [256, 512],
                "precision_decimals": WIND_DECIMALS,
                "lat": latitude,
                "lon": longitude,
                "u": json_values(u, WIND_DECIMALS),
                "v": json_values(v, WIND_DECIMALS),
            },
        )
        wind_outputs.append(name)
        print(f"SUCCESS {name}: {int(good.sum()):,} original U/V vectors")

    # ======================================================
    # 9. EXISTING 11 CLIMATOLOGICAL ANOMALIES
    # ======================================================

    reference_dates = (
        [(2, 28), (3, 1)]
        if (date.month, date.day) == (2, 29)
        else [(date.month, date.day)]
    )

    outputs = []

    for name, (units, scale, digits) in DISPLAY.items():
        drive_folder = drive_child(FOLDER_ID, name)
        climatologies = []
        reference_units = None

        for month, day in reference_dates:
            nc_path = work / f"{name}_{month:02d}.nc"
            drive_file = drive_child(drive_folder, f"{month:02d}.nc")
            download_drive_file(drive_file, nc_path)

            reference, ref_lat, ref_lon, nc_units = read_reference(
                nc_path, name, month, day
            )
            validate_grid(lat, lon, ref_lat, ref_lon)

            if reference_units is None:
                reference_units = nc_units
            elif nc_units != reference_units:
                raise ValueError(f"Inconsistent climatology units for {name}")

            climatologies.append(reference)
            print(f"Drive {name:9s}/{month:02d}.nc: {nc_path.stat().st_size:>9,} bytes")
            nc_path.unlink()

        if (
            name in SPECS
            and reference_units
            and source_units[name] != reference_units
        ):
            raise ValueError(
                f"Units mismatch for {name}: "
                f"NOAA {source_units[name]!r}, climatology {reference_units!r}"
            )

        if len(climatologies) == 1:
            reference = climatologies[0]
        else:
            reference = (
                (
                    climatologies[0].astype("float64")
                    + climatologies[1].astype("float64")
                ) / 2.0
            ).astype("float32")

        anomalies = (observed[name] - reference).astype("float32")
        if scale != 1.0:
            anomalies *= np.float32(scale)

        good = np.isfinite(anomalies)
        if not good.any():
            raise ValueError(f"No valid anomaly values for {name}")

        write_json(
            OUTPUT_DIR / f"{name}.json",
            {
                "variable": name,
                "date": date.isoformat(),
                "units": units,
                "source": "NOAA CPC CORe daily ensemble mean",
                "reference_period": "1991-2020",
                "reference_smoothing": "31-day centered circular smoothing",
                "grid": "native CORe Gaussian grid",
                "shape": [256, 512],
                "precision_decimals": digits,
                "lat": latitude,
                "lon": longitude,
                "anomaly": json_values(anomalies, digits),
            },
        )
        outputs.append(name)
        print(
            f"SUCCESS {name:9s}: {int(good.sum()):,} cells | "
            f"min {float(np.nanmin(anomalies)):+.5f} | "
            f"max {float(np.nanmax(anomalies)):+.5f} {units}"
        )

    # Manifest is written last, after all products succeed.
    write_json(
        OUTPUT_DIR / "manifest.json",
        {
            "date": date.isoformat(),
            "reference_period": "1991-2020",
            "variables": outputs,
            "wind_layers": wind_outputs,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
        },
    )

print("=" * 62)
print(f"SUCCESS: {len(DISPLAY)} anomaly products + {len(WIND_LEVELS)} wind layers")
print("SUCCESS: All products share the same analysis date")
print("SUCCESS: manifest.json")
