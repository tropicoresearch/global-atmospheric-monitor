
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TRÓPICO — Global Atmospheric Monitor
Daily NOAA CPC CORe 2-m temperature anomalies

Author: Alejandro Jaramillo-Moreno
Instituto de Ciencias de la Atmósfera y
Cambio Climático, UNAM

Reference climatology: 1991–2020
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


# ==========================================================
# 1. SETTINGS
# ==========================================================

ROOT = Path(__file__).resolve().parent.parent

OUTPUT = ROOT / "data" / "latest" / "t2m.json"

FOLDER_ID = "1iW30G3cAaIZceu2geZ0E9kg26Wlx3zE5"

BASE_URL = (
    "https://storage.googleapis.com/"
    "noaa-nws-ncep-core/grib/day/flx"
)

DRIVE_API = "https://www.googleapis.com/drive/v3/files"

TOKEN = os.environ.get("GOOGLE_ACCESS_TOKEN", "")

if not TOKEN:
    raise RuntimeError(
        "Missing GOOGLE_ACCESS_TOKEN. "
        "Run through GitHub Actions."
    )

session = requests.Session()

session.headers.update({
    "User-Agent": "TROPICO-Global-Atmospheric-Monitor/1.0"
})


# ==========================================================
# 2. GOOGLE DRIVE ACCESS
# ==========================================================

def drive_child(parent_id, filename):

    response = session.get(
        DRIVE_API,
        params={
            "q": (
                f"'{parent_id}' in parents "
                f"and name = '{filename}' "
                "and trashed = false"
            ),
            "fields": "nextPageToken,files(id,name,mimeType)",
            "pageSize": 100,
        },
        headers={
            "Authorization": f"Bearer {TOKEN}"
        },
        timeout=60,
    )

    response.raise_for_status()

    result = response.json()
    items = result.get("files", [])

    if result.get("nextPageToken") or len(items) != 1:
        raise ValueError(
            f"Expected one Drive item named {filename!r}; "
            f"found {len(items)}"
        )

    return items[0]["id"]


def download_climatology(month, path):

    folder_id = drive_child(FOLDER_ID, "t2m")
    file_id = drive_child(folder_id, f"{month:02d}.nc")

    with session.get(
        f"{DRIVE_API}/{file_id}",
        params={"alt": "media"},
        headers={
            "Authorization": f"Bearer {TOKEN}"
        },
        stream=True,
        timeout=120,
    ) as response:

        response.raise_for_status()

        with path.open("wb") as output:
            for chunk in response.iter_content(
                chunk_size=1024 * 1024
            ):
                if chunk:
                    output.write(chunk)

    print(
        f"Climatology: t2m/{month:02d}.nc "
        f"({path.stat().st_size:,} bytes)"
    )


def climatology_for_day(path, month, day):

    with xr.open_dataset(
        path, decode_times=False
    ) as ds:

        if (
            "t2m" not in ds
            or "lat" not in ds.coords
            or "lon" not in ds.coords
        ):
            raise ValueError(
                "Unexpected climatology structure"
            )

        da = ds["t2m"]

        day_dims = [
            dim for dim in da.dims
            if dim not in ("lat", "lon")
        ]

        if len(day_dims) != 1:
            raise ValueError(
                f"Expected daily climatology: {da.dims}"
            )

        day_dim = day_dims[0]

        expected_days = monthrange(2001, month)[1]

        if da.sizes[day_dim] != expected_days:
            raise ValueError(
                f"Climatology month {month:02d}: "
                f"expected {expected_days} days, "
                f"got {da.sizes[day_dim]}"
            )

        field = da.isel({
            day_dim: day - 1
        }).transpose("lat", "lon")

        values = np.asarray(
            field.values, dtype="float32"
        )

        lat = np.asarray(ds["lat"].values)
        lon = np.asarray(ds["lon"].values)

        if (
            values.shape != (256, 512)
            or not np.isfinite(values).any()
        ):
            raise ValueError(
                "Invalid climatology grid or data"
            )

        return values, lat, lon


# ==========================================================
# 3. FIND LATEST NOAA DAILY FIELD
# ==========================================================

def find_latest_t2m():

    today = datetime.now(timezone.utc).date()

    for lag in range(1, 11):

        date = today - timedelta(days=lag)

        stem = (
            f"{BASE_URL}/"
            f"{date.year}/{date.month:02d}/"
            f"flx.{date:%Y%m%d}"
        )

        response = session.get(
            stem + ".idx",
            timeout=60
        )

        if response.status_code == 404:
            print(f"{date}: not yet available")
            continue

        response.raise_for_status()

        records = []

        for line in response.text.splitlines():

            parts = line.split(":")

            if len(parts) < 6:
                continue

            records.append({
                "offset": int(parts[1]),
                "parameter": parts[3].strip(),
                "level": parts[4].strip(),
                "timing": parts[5].strip(),
            })

        if not records:
            raise ValueError(
                f"Empty NOAA inventory for {date}"
            )

        if any(
            b["offset"] <= a["offset"]
            for a, b in zip(records, records[1:])
        ):
            raise ValueError(
                f"Invalid NOAA inventory for {date}"
            )

        candidates = [
            i for i, record in enumerate(records)
            if record["parameter"] == "TMP"
            and re.fullmatch(
                r"2\s*m\s*above\s*(ground|gnd)",
                record["level"]
            )
            and re.search(
                r"\d+@3\s*hour\s*ave\(anl\)",
                record["timing"]
            )
        ]

        if len(candidates) != 1:
            raise ValueError(
                f"Expected one daily-mean t2m record "
                f"for {date}; found {len(candidates)}"
            )

        index = candidates[0]

        start = records[index]["offset"]

        end = (
            records[index + 1]["offset"] - 1
            if index + 1 < len(records)
            else None
        )

        print(f"Latest CORe t2m date: {date}")

        return date, stem + ".grb", start, end

    raise RuntimeError(
        "No daily CORe data found in the last 10 days"
    )


def download_t2m_field(url, start, end, path):

    byte_range = (
        f"bytes={start}-{end}"
        if end is not None
        else f"bytes={start}-"
    )

    response = session.get(
        url,
        headers={"Range": byte_range},
        timeout=120
    )

    if response.status_code != 206:
        raise RuntimeError(
            "Expected partial NOAA download "
            f"(HTTP 206), got {response.status_code}"
        )

    response.raise_for_status()

    path.write_bytes(response.content)

    if path.stat().st_size == 0:
        raise ValueError("Empty GRIB message")

    print(
        f"NOAA t2m message: "
        f"{path.stat().st_size:,} bytes"
    )


# ==========================================================
# 4. DECODE NOAA GRIB
# ==========================================================

def decode_t2m(path):

    with path.open("rb") as stream:

        gid = codes_grib_new_from_file(stream)

        if gid is None:
            raise ValueError(
                "Could not decode NOAA GRIB"
            )

        try:
            ni = int(codes_get(gid, "Ni"))
            nj = int(codes_get(gid, "Nj"))

            if (nj, ni) != (256, 512):
                raise ValueError(
                    f"Unexpected NOAA grid: {(nj, ni)}"
                )

            units = str(codes_get(gid, "units"))

            lat2d = codes_get_array(
                gid, "latitudes"
            ).reshape(nj, ni)

            lon2d = (
                codes_get_array(
                    gid, "longitudes"
                ).reshape(nj, ni) % 360.0
            )

            values = codes_get_values(
                gid
            ).reshape(nj, ni)

        finally:
            codes_release(gid)

    if units != "K":
        raise ValueError(
            f"Expected kelvin, got {units!r}"
        )

    if not np.allclose(
        lat2d, lat2d[:, :1],
        rtol=0, atol=1e-4
    ):
        raise ValueError(
            "Latitude grid is not separable"
        )

    if not np.allclose(
        lon2d, lon2d[:1, :],
        rtol=0, atol=1e-4
    ):
        raise ValueError(
            "Longitude grid is not separable"
        )

    values = np.where(
        np.isfinite(values)
        & (np.abs(values) < 1e15),
        values,
        np.nan
    ).astype("float32")

    da = xr.DataArray(
        values,
        dims=("lat", "lon"),
        coords={
            "lat": lat2d[:, 0],
            "lon": lon2d[0, :]
        }
    )

    da = da.sortby("lat").sortby("lon")

    return (
        da.values,
        da.lat.values,
        da.lon.values
    )


# ==========================================================
# 5. PROCESS DAILY TEMPERATURE ANOMALY
# ==========================================================

print("=" * 60)
print("TRÓPICO - GLOBAL ATMOSPHERIC MONITOR")
print("NOAA CPC CORe - DAILY T2M ANOMALY")
print("=" * 60)

date, grib_url, byte_start, byte_end = (
    find_latest_t2m()
)

with TemporaryDirectory() as directory:

    work = Path(directory)

    grib_file = work / "t2m.grb"

    download_t2m_field(
        grib_url,
        byte_start,
        byte_end,
        grib_file
    )

    observed, lat, lon = decode_t2m(grib_file)

    reference_dates = [
        (date.month, date.day)
    ]

    # Leap-day climatology
    if (date.month, date.day) == (2, 29):
        reference_dates = [
            (2, 28),
            (3, 1)
        ]

    references = []

    for month, day in reference_dates:

        nc_file = work / f"{month:02d}.nc"

        download_climatology(
            month, nc_file
        )

        reference, clat, clon = (
            climatology_for_day(
                nc_file, month, day
            )
        )

        if not (
            np.allclose(
                lat, clat, rtol=0, atol=1e-5
            )
            and np.allclose(
                lon, clon, rtol=0, atol=1e-5
            )
        ):
            raise ValueError(
                "NOAA and climatology grids differ"
            )

        references.append(reference)

    if len(references) == 1:
        climatology = references[0]

    else:
        climatology = (
            (
                references[0].astype("float64")
                + references[1].astype("float64")
            ) / 2.0
        ).astype("float32")

    # K - K = degC (temperature anomaly)
    anomaly = (
        observed.astype("float32")
        - climatology
    )

    valid = np.isfinite(anomaly)

    if not valid.any():
        raise ValueError(
            "No finite temperature anomalies"
        )

    print(
        f"Valid grid points: "
        f"{valid.sum():,}/{anomaly.size:,}"
    )

    print(
        f"Anomaly min/max: "
        f"{np.nanmin(anomaly):+.2f}/"
        f"{np.nanmax(anomaly):+.2f} degC"
    )

    # Visualization-ready JSON.
    # Rounding is applied only to this web product.
    result = {
        "variable": "t2m",
        "date": date.isoformat(),
        "units": "degC",
        "source": "NOAA CPC CORe daily ensemble mean",
        "reference_period": "1991-2020",
        "reference_smoothing": (
            "31-day centered circular smoothing"
        ),
        "grid": "native CORe Gaussian grid",
        "shape": list(anomaly.shape),
        "lat": [
            round(float(x), 6) for x in lat
        ],
        "lon": [
            round(float(x), 6) for x in lon
        ],
        "anomaly": [
            [
                round(float(v), 2)
                if np.isfinite(v)
                else None
                for v in row
            ]
            for row in anomaly
        ]
    }

    OUTPUT.parent.mkdir(
        parents=True, exist_ok=True
    )

    temporary = OUTPUT.with_suffix(".json.tmp")

    with temporary.open(
        "w", encoding="utf-8"
    ) as output:
        json.dump(
            result,
            output,
            separators=(",", ":"),
            allow_nan=False
        )

    temporary.replace(OUTPUT)

    print(f"SUCCESS: {OUTPUT}")
