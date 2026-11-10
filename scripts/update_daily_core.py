
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TRÓPICO | Global Atmospheric Monitor
Daily NOAA CPC CORe anomalies (1991–2020 climatology).

Author: Alejandro Jaramillo-Moreno
Instituto de Ciencias de la Atmósfera y Cambio Climático, UNAM
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
    codes_get,
    codes_get_array,
    codes_get_values,
    codes_grib_new_from_file,
    codes_release,
)

from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# ==========================================================
# 1. SETTINGS
# ==========================================================

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "data" / "latest"

FOLDER_ID = "1iW30G3cAaIZceu2geZ0E9kg26Wlx3zE5"

DRIVE_URL = (
    "https://www.googleapis.com/drive/v3/files"
)

NOAA_URL = (
    "https://storage.googleapis.com/"
    "noaa-nws-ncep-core/grib/day"
)

TOKEN = os.environ.get(
    "GOOGLE_ACCESS_TOKEN", ""
)

LOOKBACK_DAYS = 10

if not TOKEN:
    raise RuntimeError(
        "GOOGLE_ACCESS_TOKEN is missing"
    )


# ==========================================================
# 2. VARIABLE CATALOG
# ==========================================================

# Original fields:
# group, GRIB parameter, level, timing

SPECS = {

    "t2m": (
        "flx", "TMP",
        r"^2\s*m\s*above\s*(ground|gnd)$",
        r"\d+@3\s*hour\s*ave\(anl\)"
    ),

    "prate": (
        "flx", "PRATE",
        r"^(surface|sfc)$",
        None
    ),

    "pwat": (
        "flx", "PWAT",
        r"(atmos\s+col|entire\s+atmosphere)",
        None
    ),

    "u850": (
        "pgb", "UGRD",
        r"^850\s*(mb|hPa)$",
        None
    ),

    "v850": (
        "pgb", "VGRD",
        r"^850\s*(mb|hPa)$",
        None
    ),

    "u250": (
        "pgb", "UGRD",
        r"^250\s*(mb|hPa)$",
        None
    ),

    "v250": (
        "pgb", "VGRD",
        r"^250\s*(mb|hPa)$",
        None
    ),

    "z500": (
        "pgb", "HGT",
        r"^500\s*(mb|hPa)$",
        None
    ),

    "omega500": (
        "pgb", "VVEL",
        r"^500\s*(mb|hPa)$",
        None
    ),
}


DERIVED = {

    "wspd850": ("u850", "v850"),

    "wspd250": ("u250", "v250"),

}


# Visualization units, scaling and decimal precision.
# Rounding is used only when exporting JSON.

DISPLAY = {

    "t2m": ("degC", 1.0, 2),

    "prate": ("mm day-1", 86400.0, 2),

    "pwat": ("mm", 1.0, 2),

    "u850": ("m s-1", 1.0, 2),
    "v850": ("m s-1", 1.0, 2),

    "u250": ("m s-1", 1.0, 2),
    "v250": ("m s-1", 1.0, 2),

    "wspd850": ("m s-1", 1.0, 2),
    "wspd250": ("m s-1", 1.0, 2),

    "z500": ("gpm", 1.0, 2),

    "omega500": ("Pa s-1", 1.0, 5),

}


# ==========================================================
# 3. HTTP SESSION
# ==========================================================

session = requests.Session()

retry = Retry(
    total=4,
    backoff_factor=1,
    status_forcelist=[
        429, 500, 502, 503, 504
    ],
)

session.mount(
    "https://",
    HTTPAdapter(max_retries=retry)
)

session.headers.update({
    "User-Agent":
        "TROPICO-Global-Atmospheric-Monitor/1.0"
})


# ==========================================================
# 4. NOAA DAILY GRIB INVENTORY
# ==========================================================

def inventory_for_day(group, date):

    stem = (
        f"{NOAA_URL}/{group}/"
        f"{date.year}/{date.month:02d}/"
        f"{group}.{date:%Y%m%d}"
    )

    response = session.get(
        stem + ".idx",
        timeout=(20, 90)
    )

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
            "timing": (
                parts[5].strip()
                if len(parts) > 5
                else ""
            ),
        })

    if (
        not records
        or any(
            b["offset"] <= a["offset"]
            for a, b in zip(
                records, records[1:]
            )
        )
    ):
        raise ValueError(
            f"Invalid {group} inventory for {date}"
        )

    return stem + ".grb", records


def match_record(records, name):

    _, parameter, level, timing = SPECS[name]

    matches = [

        i for i, record in enumerate(records)

        if record["parameter"] == parameter

        and re.search(
            level,
            record["level"],
            re.IGNORECASE
        )

        and (
            timing is None
            or re.search(
                timing,
                record["timing"],
                re.IGNORECASE
            )
        )
    ]

    if len(matches) != 1:
        raise ValueError(
            f"Expected one {name} GRIB record, "
            f"got {len(matches)}"
        )

    i = matches[0]

    end = (
        records[i + 1]["offset"] - 1
        if i + 1 < len(records)
        else None
    )

    return records[i]["offset"], end


def latest_complete_day():

    today = datetime.now(timezone.utc).date()

    for lag in range(1, LOOKBACK_DAYS + 1):

        day = today - timedelta(days=lag)

        inventories = {

            group: inventory_for_day(group, day)

            for group in ("flx", "pgb")
        }

        if any(
            item is None
            for item in inventories.values()
        ):

            print(
                f"{day}: flx/pgb incomplete "
                "(missing inventory)"
            )

            continue

        try:

            selections = {

                name: match_record(
                    inventories[spec[0]][1],
                    name
                )

                for name, spec in SPECS.items()
            }

        except ValueError as exc:

            print(f"{day}: {exc}")
            continue

        print(
            f"Selected common NOAA CORe date: {day}"
        )

        return day, inventories, selections

    raise RuntimeError(
        "No common complete date in the last "
        f"{LOOKBACK_DAYS} days"
    )


# ==========================================================
# 5. GRIB DOWNLOAD AND DECODING
# ==========================================================

def download_grib(url, start, end, path):

    limit = (
        f"bytes={start}-{end}"
        if end is not None
        else f"bytes={start}-"
    )

    with session.get(
        url,
        headers={"Range": limit},
        stream=True,
        timeout=(20, 120)
    ) as response:

        if response.status_code != 206:
            raise RuntimeError(
                "Expected HTTP 206; got "
                f"{response.status_code}: {url}"
            )

        with path.open("wb") as handle:

            for chunk in response.iter_content(
                chunk_size=1024 * 1024
            ):

                if chunk:
                    handle.write(chunk)

    if path.stat().st_size < 20:
        raise ValueError(
            f"Empty GRIB message: {url}"
        )


def decode_grib(path):

    with path.open("rb") as handle:

        gid = codes_grib_new_from_file(
            handle
        )

        if gid is None:
            raise ValueError(
                "Could not decode NOAA GRIB message"
            )

        try:

            ni = int(codes_get(gid, "Ni"))
            nj = int(codes_get(gid, "Nj"))

            units = str(
                codes_get(gid, "units")
            )

            lat2d = codes_get_array(
                gid, "latitudes"
            ).reshape(nj, ni)

            lon2d = (
                codes_get_array(
                    gid, "longitudes"
                ).reshape(nj, ni) % 360
            )

            values = codes_get_values(
                gid
            ).reshape(nj, ni)

        finally:
            codes_release(gid)

    if (nj, ni) != (256, 512):
        raise ValueError(
            f"Unexpected NOAA grid {(nj, ni)}"
        )

    if not np.allclose(
        lat2d,
        lat2d[:, :1],
        atol=1e-4
    ):
        raise ValueError(
            "NOAA latitude grid is not separable"
        )

    if not np.allclose(
        lon2d,
        lon2d[:1, :],
        atol=1e-4
    ):
        raise ValueError(
            "NOAA longitude grid is not separable"
        )

    values = np.where(
        np.isfinite(values)
        & (np.abs(values) < 1e15),
        values,
        np.nan
    ).astype("float32")

    field = xr.DataArray(
        values,
        dims=("lat", "lon"),
        coords={
            "lat": lat2d[:, 0],
            "lon": lon2d[0, :]
        }
    )

    field = (
        field.sortby("lat")
        .sortby("lon")
    )

    if not np.isfinite(field.values).any():
        raise ValueError(
            "NOAA field contains no finite values"
        )

    return (
        field.values,
        field.lat.values,
        field.lon.values,
        units
    )


# ==========================================================
# 6. GOOGLE DRIVE CLIMATOLOGIES
# ==========================================================

def drive_child(parent, name):

    response = session.get(
        DRIVE_URL,
        params={
            "q": (
                f"'{parent}' in parents "
                f"and name = '{name}' "
                "and trashed = false"
            ),
            "fields": (
                "nextPageToken,"
                "files(id,name,mimeType)"
            ),
            "pageSize": 100,
        },
        headers={
            "Authorization": f"Bearer {TOKEN}"
        },
        timeout=(20, 90),
    )

    response.raise_for_status()

    result = response.json()

    files = result.get("files", [])

    if (
        result.get("nextPageToken")
        or len(files) != 1
    ):
        raise ValueError(
            "Expected one Google Drive child "
            f"{name!r}, got {len(files)}"
        )

    return files[0]["id"]


def download_drive_file(file_id, path):

    with session.get(
        f"{DRIVE_URL}/{file_id}",
        params={"alt": "media"},
        headers={
            "Authorization": f"Bearer {TOKEN}"
        },
        stream=True,
        timeout=(20, 120),
    ) as response:

        response.raise_for_status()

        with path.open("wb") as handle:

            for chunk in response.iter_content(
                chunk_size=1024 * 1024
            ):

                if chunk:
                    handle.write(chunk)

    if path.stat().st_size == 0:
        raise ValueError(
            f"Empty Drive file: {file_id}"
        )


def read_reference(path, name, month, day):

    with xr.open_dataset(
        path,
        decode_times=False
    ) as ds:

        if (
            name not in ds
            or "lat" not in ds.coords
            or "lon" not in ds.coords
        ):
            raise ValueError(
                "Invalid climatology structure "
                f"for {name}"
            )

        da = ds[name]

        day_dims = [
            dim for dim in da.dims
            if dim not in ("lat", "lon")
        ]

        expected_days = monthrange(
            2001, month
        )[1]

        if (
            len(day_dims) != 1
            or da.sizes[day_dims[0]]
            != expected_days
        ):
            raise ValueError(
                "Unexpected daily calendar for "
                f"{name}/{month:02d}.nc"
            )

        field = da.isel({
            day_dims[0]: day - 1
        }).transpose("lat", "lon")

        values = np.asarray(
            field.values,
            dtype="float32"
        )

        if (
            values.shape != (256, 512)
            or not np.isfinite(values).any()
        ):
            raise ValueError(
                f"Invalid climatology field for {name}"
            )

        return (
            values,
            ds.lat.values,
            ds.lon.values,
            str(da.attrs.get("units", ""))
        )


# ==========================================================
# 7. VALIDATION AND JSON OUTPUT
# ==========================================================

def validate_grid(
    lat, lon,
    reference_lat, reference_lon
):

    if not (
        np.allclose(
            lat, reference_lat,
            rtol=0, atol=1e-5
        )
        and np.allclose(
            lon, reference_lon,
            rtol=0, atol=1e-5
        )
    ):
        raise ValueError(
            "NOAA and climatology grids differ"
        )


def write_json(path, content):

    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    temporary = path.with_suffix(
        path.suffix + ".tmp"
    )

    with temporary.open(
        "w",
        encoding="utf-8"
    ) as handle:

        json.dump(
            content,
            handle,
            separators=(",", ":"),
            allow_nan=False
        )

    temporary.replace(path)


# ==========================================================
# 8. MAIN PROCESSING
# ==========================================================

print("=" * 62)
print("TRÓPICO — GLOBAL ATMOSPHERIC MONITOR")
print("NOAA CPC CORe — DAILY ATMOSPHERIC ANOMALIES")
print("=" * 62)

date, inventories, selections = (
    latest_complete_day()
)

with TemporaryDirectory() as tmpdir:

    work = Path(tmpdir)

    observed = {}
    source_units = {}

    lat = None
    lon = None

    # ------------------------------------------------------
    # Download 9 original NOAA fields
    # ------------------------------------------------------

    for name, spec in SPECS.items():

        group = spec[0]

        url, _ = inventories[group]

        start, end = selections[name]

        grib_path = work / f"{name}.grb"

        download_grib(
            url,
            start,
            end,
            grib_path
        )

        values, field_lat, field_lon, units = (
            decode_grib(grib_path)
        )

        if lat is None:

            lat = field_lat.copy()
            lon = field_lon.copy()

        else:

            validate_grid(
                lat, lon,
                field_lat, field_lon
            )

        observed[name] = values
        source_units[name] = units

        print(
            f"NOAA {name:9s}: "
            f"{grib_path.stat().st_size:>9,} bytes "
            f"| {units}"
        )

        grib_path.unlink()

    # ------------------------------------------------------
    # Derived daily-mean horizontal wind speeds
    # ------------------------------------------------------

    for name, (u, v) in DERIVED.items():

        observed[name] = np.hypot(
            observed[u],
            observed[v]
        ).astype("float32")

        source_units[name] = source_units[u]

    # ------------------------------------------------------
    # Climatological day selection
    # ------------------------------------------------------

    reference_dates = (

        [(2, 28), (3, 1)]

        if (date.month, date.day) == (2, 29)

        else [(date.month, date.day)]
    )

    outputs = []

    # ------------------------------------------------------
    # Generate all 11 anomaly products
    # ------------------------------------------------------

    for name in DISPLAY:

        drive_folder = drive_child(
            FOLDER_ID,
            name
        )

        climatologies = []

        reference_units = None

        for month, day in reference_dates:

            nc_path = (
                work / f"{name}_{month:02d}.nc"
            )

            drive_file = drive_child(
                drive_folder,
                f"{month:02d}.nc"
            )

            download_drive_file(
                drive_file,
                nc_path
            )

            (
                reference,
                ref_lat,
                ref_lon,
                units
            ) = read_reference(
                nc_path,
                name,
                month,
                day
            )

            validate_grid(
                lat, lon,
                ref_lat, ref_lon
            )

            if reference_units is None:

                reference_units = units

            elif units != reference_units:

                raise ValueError(
                    "Inconsistent climatology units "
                    f"for {name}"
                )

            climatologies.append(
                reference
            )

            print(
                f"Drive {name:9s}/{month:02d}.nc: "
                f"{nc_path.stat().st_size:>9,} bytes"
            )

            nc_path.unlink()

        # --------------------------------------------------
        # Check original-variable units
        # --------------------------------------------------

        if (
            name in SPECS
            and reference_units
            and source_units[name] != reference_units
        ):
            raise ValueError(
                f"Units mismatch for {name}: "
                f"NOAA {source_units[name]!r}, "
                f"climatology {reference_units!r}"
            )

        # --------------------------------------------------
        # Leap-day reference interpolation
        # --------------------------------------------------

        if len(climatologies) == 1:

            reference = climatologies[0]

        else:

            reference = (
                (
                    climatologies[0].astype("float64")
                    + climatologies[1].astype("float64")
                ) / 2
            ).astype("float32")

        # --------------------------------------------------
        # Calculate anomalies at native resolution
        # --------------------------------------------------

        units, scale, digits = DISPLAY[name]

        anomalies = (
            observed[name] - reference
        ).astype("float32")

        if scale != 1.0:

            anomalies *= np.float32(
                scale
            )

        good = np.isfinite(
            anomalies
        )

        if not good.any():

            raise ValueError(
                "No valid anomaly values "
                f"for {name}"
            )

        # --------------------------------------------------
        # JSON for web visualization
        # --------------------------------------------------

        output = (
            OUTPUT_DIR / f"{name}.json"
        )

        write_json(
            output,
            {
                "variable": name,

                "date": date.isoformat(),

                "units": units,

                "source": (
                    "NOAA CPC CORe "
                    "daily ensemble mean"
                ),

                "reference_period": "1991-2020",

                "reference_smoothing": (
                    "31-day centered "
                    "circular smoothing"
                ),

                "grid": (
                    "native CORe Gaussian grid"
                ),

                "shape": [256, 512],

                "precision_decimals": digits,

                "lat": [
                    round(float(x), 6)
                    for x in lat
                ],

                "lon": [
                    round(float(x), 6)
                    for x in lon
                ],

                "anomaly": [
                    [
                        round(float(x), digits)
                        if np.isfinite(x)
                        else None
                        for x in row
                    ]
                    for row in anomalies
                ],
            }
        )

        outputs.append(name)

        print(
            f"SUCCESS {name:9s}: "
            f"{int(good.sum()):,} cells | "
            f"min {float(np.nanmin(anomalies)):+.5f} | "
            f"max {float(np.nanmax(anomalies)):+.5f} "
            f"{units}"
        )

    # ------------------------------------------------------
    # Manifest for the future variable selector
    # ------------------------------------------------------

    write_json(
        OUTPUT_DIR / "manifest.json",
        {
            "date": date.isoformat(),

            "reference_period": "1991-2020",

            "variables": outputs,

            "generated_utc": (
                datetime.now(
                    timezone.utc
                ).isoformat()
            ),
        }
    )

print(
    "SUCCESS: generated "
    f"{len(DISPLAY)} variable files "
    "and manifest.json"
)
