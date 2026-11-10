<p align="center">
  <img src="assets/tropico_logo.png" alt="TRÓPICO Research Group" width="650">
</p>

# Global Atmospheric Monitor

**Daily Global Atmospheric Conditions and Anomalies**

Global Atmospheric Monitor is an interactive scientific visualization project developed within **TRÓPICO — Research Group on Tropical Climatology, Hydrology, and Meteorology**, at the Universidad Nacional Autónoma de México (UNAM).

The project aims to provide an accessible and regularly updated representation of the global atmospheric state using daily analyses and anomalies.

## Overview

Global Atmospheric Monitor is designed to visualize recent atmospheric conditions through an interactive three-dimensional globe.

The project uses atmospheric analyses from **NOAA Climate Prediction Center (CPC) CORe**, together with daily climatological reference fields calculated for the period **1991–2020**.

Atmospheric anomalies are determined by comparing daily analyses with their corresponding calendar-day climatologies.

The monitoring system is designed to update automatically as new daily analyses become available.

## Atmospheric Variables

Initial development includes the following atmospheric variables:

| Variable | Code |
|---|---|
| 2-m air temperature | `t2m` |
| Precipitation rate | `prate` |
| Total precipitable water | `pwat` |
| 850-hPa zonal wind | `u850` |
| 850-hPa meridional wind | `v850` |
| 250-hPa zonal wind | `u250` |
| 250-hPa meridional wind | `v250` |
| 500-hPa geopotential height | `z500` |
| 500-hPa vertical velocity | `omega500` |
| 850-hPa horizontal wind speed | `wspd850` |
| 250-hPa horizontal wind speed | `wspd250` |

Additional atmospheric variables may be incorporated as the project develops.

## Scientific Methodology

The atmospheric reference climatologies are calculated from daily NOAA CPC CORe analyses for 1991–2020.

The methodology includes:

- A fixed 1991–2020 reference period.
- A 365-day climatological calendar.
- A centered 31-day smoothing of the climatological annual cycle.
- Interpolation of the February 29 reference from February 28 and March 1.
- Preservation of the original CORe Gaussian grid during scientific processing.
- Consistent treatment of derived horizontal wind-speed climatologies.

The climatologies are stored as losslessly compressed NetCDF files, preserving their original `float32` values.

The original scientific datasets and climatological processing scripts are maintained separately.

## Project Architecture

Global Atmospheric Monitor consists of three components:

**Scientific processing**

The private `global-atmospheric-monitor-core` repository contains the historical data-processing scripts and climatological calculations.

**Reference data storage**

Daily reference climatologies are maintained in Google Drive and accessed by authorized automated processing tools.

**Public visualization**

This repository contains the interactive web application, daily processing scripts, and visualization-ready atmospheric products.

The public website does not require visitors to authenticate with Google or access the scientific climatology archive.

## Planned Features

- Interactive 3D global atmospheric visualization.
- Daily atmospheric conditions and anomalies.
- Variable-specific color scales.
- Interactive geographic inspection.
- Automatic updates based on NOAA CPC CORe availability.
- Scientific metadata and analysis-date information.
- Responsive interface for desktop and mobile devices.
- Lightweight deployment using GitHub Pages.

## Development Status

**Status: Under active development.**

The 1991–2020 reference climatologies have been calculated and validated.

The public visualization and automated daily processing system are currently under development.

## Data Sources and Disclaimer

Atmospheric data are provided by **NOAA Climate Prediction Center (CPC) CORe**.

Global Atmospheric Monitor is an independent scientific and educational initiative. It is not an official NOAA operational product or weather forecasting service.

Atmospheric anomalies describe departures from climatological reference conditions and do not, by themselves, establish attribution to a particular climate phenomenon.

## Privacy

Google Drive access is used exclusively by the project's authorized processing infrastructure to retrieve scientific reference climatologies.

Visitors do not need Google accounts or Google Drive permissions to use the public visualization.

The project's privacy policy will be available on the public website.

---

## Creator and Maintainer

<table>
  <tr>
    <td width="160" align="center" valign="middle">
      <a href="https://ajaramillomoreno.github.io/">
        <img src="assets/alejandro_jaramillo_moreno.png" alt="Alejandro Jaramillo-Moreno" width="140">
      </a>
    </td>
    <td valign="middle">
      <strong>Alejandro Jaramillo-Moreno</strong><br>
      Research Scientist<br>
      Instituto de Ciencias de la Atmósfera y Cambio Climático (ICAyCC)<br>
      Universidad Nacional Autónoma de México (UNAM)<br>
      Founder of TRÓPICO<br><br>
      <a href="https://ajaramillomoreno.github.io/">Personal Website</a>
    </td>
  </tr>
</table>

This project was created and is maintained by **Alejandro Jaramillo-Moreno**.