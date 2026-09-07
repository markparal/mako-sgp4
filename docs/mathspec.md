# mako-sgp4: Math Specification
**By Mark Paral**


<p align="center">
  <img src="../assets/logo_dark.png" alt="mako-sgp4 logo" width="240">
</p>

## Table of Contents
1. [Introduction](#introduction)
2. [SGP4 Historical Background](#sgp4-historical-background)
3. [GP Element Sets](#gp-element-sets)
    - [Two-Line Element (TLE) Format](#two-line-element-tle-format)
    - [Orbit Mean-Elements Message (OMM) Format](#orbit-mean-elements-message-omm-format)
4. [SGP4 Algorithm](#sgp4-algorithm)
    - [Initialization](#initialization)
    - [Propagation](#propagation)
5. [Thanks](#thanks)
6. [References](#references)

## Introduction
On October 4th, 1957, Sputnik was launched into orbit by the Soviet Union. With the opening of space for artificial satellites, the urgent need to catalog and track all objects around the Earth became apparent. Thanks to the efforts of numerous members of the US military and scientific community, new theories were developed and infrastructure put in place to realize these goals. Today, these tools act as the backbone of modern space situational awareness (SSA). 

Critical to this infrastructure is the Simplified General Perturbations 4 model (SGP4), which, using provided general perturbation element sets (GPs), generates reliably accurate (sub-1 km position error at epoch) position and velocity state estimates for each cataloged orbiting object. Originally implemented in 1970, SGP4 strikes a difficult balance between computational efficiency and accuracy as a semianalytical orbital propagator. The impressiveness of this achievement cannot be overstated, as it remains one of, if not the most, important propagator in the industry nearly 56 years later (as of 2026).

As a member of the space community who makes active use of these SSA resources, [mako-sgp4](https://github.com/markparal/mako-sgp4) ([10]) was developed to gain a deeper understanding of the underlying theory behind this immensely important model. The rest of this document will dive into this theory to facilitate understanding of the codebase.

## SGP4 Historical Background
For a complete historical rundown of the development of SGP4, it is recommended to read *History of Analytical Orbit Modeling in the U.S. Space Surveillance System* by Hoots et al ([9]), which details the creation and evolution of the U.S. Space Surveillance system. This paper also discusses the various theories and works that contributed to the modern SGP4 algorithm. A (non-exhaustive) list includes:
- The effects of the J2, J3, and J4 Earth zonal harmonics on the orbit of a satellite ([1], [2])
- The effects of atmospheric drag on the orbits of satellites ([3], [4], [6])
- The avoidance of small divisors of eccentricity or sine of inclination in propagation ([5])
- The inclusion of lunar and solar gravitational effects as well as Earth tesseral harmonics ([7], [8])

## GP Element Sets

| Field | Units | Description |
| --- | --- | --- |
| Common Name | - | The commonly used name for the satellite |
| Satellite Catalog Number | - | NORAD satellite catalog number (NORAD ID) |
| Classification | - | Security classification (`U` = Unclassified, `C` = Classified, `S` = Secret) |
| International Designator | - | International designator in `Y-NP` form, where `Y` is launch year (4+ digits), `N` is launch number of that year (3+ digits), and `P` is piece of launch (1+ characters) |
| Epoch Datetime | UTC | UTC epoch datetime of the GP elements |
| First Derivative of Mean Motion | revs/day^2 | Brouwer First time derivative of mean motion |
| Second Derivative of Mean Motion | revs/day^3 | Brouwer Second time derivative of mean motion |
| B* | 1/Earth radii | Atmospheric drag coefficient |
| Ephemeris Type | - | Ephemeris type (always zero) |
| Element Set Number | - | Element set number |
| Inclination | degrees | Orbital inclination |
| Right Ascension of Ascending Node | degrees | Orbital right ascension of the ascending node (RAAN) |
| Eccentricity | - | Orbital eccentricity |
| Argument of Perigee | degrees | Orbital argument of perigee |
| Mean Anomaly | degrees | Orbital mean anomaly |
| Mean Motion | revs/day | Kozai Mean motion |
| Revolution Number at Epoch | revs | Revolution number at epoch |

<p align="center"><strong>Table 1.</strong> Standard GP element set fields</p>

The contents of GP element sets are described in Table 1 above. These elements characterize the orbit of a satellite and are what SGP4 uses to propagate the position and velocity over time. There are two primary formats for the distribution and ingestion of GP element sets: the two-line element (TLE) format and the orbit mean-elements message (OMM) format.

### Two-Line Element (TLE) Format
The TLE formats all the necessary GP elements into two lines. An example of this is given below.

```
ISS (ZARYA)
1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537
```

A detailed breakdown of the TLE format and its associated fields can be found at [11]. The standard format follows the convention depicted below, where `N` is a number 0-9 (or a space in certain cases), `C` is the classification character, `P` is any character A-Z (or a space in certain cases), `+` indicates either a '+', '-', or space must be present, `-` indicates a '+' or '-' must be present, `A` is the Alpha-5 encoded digit which can be a character (excluding I and O), digit, or space, and `D` is any printable character or space. Note that line 0 is optional. Table 2 describes what character correspond to what element field.

```
DDDDDDDDDDDDDDDDDDDDDDDD
1 ANNNNC NNNNNPPP NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN
2 ANNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN
```

| Line | Characters | Field | Units | Notes |
| --- | --- | --- | --- | --- |
| 0 | 0-23 | Common Name | - | Optional |
| 1 | 0 | Line Number | - | Always `1` |
| 1 | 2-6 | Satellite Catalog Number | - | NORAD ID, Alpha-5 encoding allowed ([12]) |
| 1 | 7 | Classification | - | `U` = Unclassified, `C` = Classified, `S` = Secret |
| 1 | 9-10 | International Designator | - | Last two digits of launch year |
| 1 | 11-13 | International Designator | - | Launch number of the year |
| 1 | 14-16 | International Designator | - | Piece of the launch |
| 1 | 18-19 | Epoch Datetime | - | Last two digits of year |
| 1 | 20-31 | Epoch Datetime | days | Day of year and fractional portion of the day |
| 1 | 33-42 | First Derivative of Mean Motion | revs/day^2 | TLE stores value divided by two |
| 1 | 44-51 | Second Derivative of Mean Motion | revs/day^3 | TLE stores value divided by six, decimal point assumed |
| 1 | 53-60 | B* | 1/Earth radii | Decimal point assumed |
| 1 | 62 | Ephemeris Type | - | Always `0` |
| 1 | 64-67 | Element Set Number | - | |
| 1 | 68 | Checksum | - | Modulo 10 |
| 2 | 0 | Line Number | - | Always `2` |
| 2 | 2-6 | Satellite Catalog Number | - | Must match line 1 |
| 2 | 8-15 | Inclination | degrees | |
| 2 | 17-24 | Right Ascension of Ascending Node | degrees | |
| 2 | 26-32 | Eccentricity | - | Decimal point assumed |
| 2 | 34-41 | Argument of Perigee | degrees | |
| 2 | 43-50 | Mean Anomaly | degrees | |
| 2 | 52-62 | Mean Motion | revs/day | |
| 2 | 63-67 | Revolution Number at Epoch | revs | |
| 2 | 68 | Checksum | - | Modulo 10 |

<p align="center"><strong>Table 2.</strong> TLE format (zero-based indexing)</p>

While the TLE has traditionally been the GP format of choice, it has several deficiencies which make switching to a new format (OMM) necessary in the coming years. Some examples include:
1. The limiting of epochs and international designators to years between 1957 and 2056
2. The limiting of NORAD IDs to 339999 (even with Alpha-5 encoding)
3. International designators can only represent launches up to 999 in a single year
4. International designators can only represent pieces of a launch up to a 3 character code

All of these issues are alleviated with the OMM format.

### Orbit Mean-Elements Message (OMM) Format
The OMM format assigns values to keywords. These keywords correspond to the necessary GP elements. An example of this is given below.

```
CCSDS_OMM_VERS = 2.0
CREATION_DATE  = 
ORIGINATOR     = 

OBJECT_NAME    = 2026-106A
OBJECT_ID      = 2026-106A
CENTER_NAME    = EARTH
REF_FRAME      = TEME
TIME_SYSTEM    = UTC
MEAN_ELEMENT_THEORY = SGP/SGP4

EPOCH          = 2026-06-14T15:07:48.259488
MEAN_MOTION    = 15.11169557
ECCENTRICITY   = .00147468
INCLINATION    = 97.5103
RA_OF_ASC_NODE = 247.7605
ARG_OF_PERICENTER = 169.6213
MEAN_ANOMALY   = 190.5325

EPHEMERIS_TYPE = 0
CLASSIFICATION_TYPE = U
NORAD_CAT_ID   = 69097
ELEMENT_SET_NO = 999
REV_AT_EPOCH   = 459
BSTAR          = .39221734E-3
MEAN_MOTION_DOT = .6535E-4
MEAN_MOTION_DDOT = 0
```

The standards for the OMM format and its associated fields can be found at [13]. Table 3 describes what keyword correspond to what element field.

| Keyword | Field | Units | Notes |
| --- | --- | --- | --- |
| CCSDS_OMM_VERS | - | - | Typically `2.0` |
| CREATION_DATE | - | - | Optional |
| ORIGINATOR | - | - | Optional |
| OBJECT_NAME | Common Name | - | Optional |
| OBJECT_ID | International Designator | - | Full `Y-NP` form |
| CENTER_NAME | - | - | Always `EARTH` |
| REF_FRAME | - | - | Always `TEME` |
| TIME_SYSTEM | - | - | Always `UTC` |
| MEAN_ELEMENT_THEORY | - | - | Always `SGP/SGP4` |
| EPOCH | Epoch Datetime | UTC | ISO-8601 UTC |
| MEAN_MOTION | Mean Motion | revs/day | |
| ECCENTRICITY | Eccentricity | - | |
| INCLINATION | Inclination | degrees | |
| RA_OF_ASC_NODE | Right Ascension of Ascending Node | degrees | |
| ARG_OF_PERICENTER | Argument of Perigee | degrees | |
| MEAN_ANOMALY | Mean Anomaly | degrees | |
| EPHEMERIS_TYPE | Ephemeris Type | - | Always `0` |
| CLASSIFICATION_TYPE | Classification | - | `U` = Unclassified, `C` = Classified, `S` = Secret |
| NORAD_CAT_ID | Satellite Catalog Number | - | Integer NORAD ID |
| ELEMENT_SET_NO | Element Set Number | - | |
| REV_AT_EPOCH | Revolution Number at Epoch | revs | |
| BSTAR | B* | 1/Earth radii | |
| MEAN_MOTION_DOT | First Derivative of Mean Motion | revs/day^2 | OMM stores value divided by two (Space-Track standard, not CCSDS) |
| MEAN_MOTION_DDOT | Second Derivative of Mean Motion | revs/day^3 | OMM stores value divided by six (Space-Track standard, not CCSDS) |

<p align="center"><strong>Table 3.</strong> OMM format</p>

Because the OMM format follows a keyword-value pattern, it exists in multiple general-purpose data formats as well, including KVN (as the example above shows), XML, JSON, and CSV.

## SGP4 Algorithm
As stated previously, the goal of the SGP4 propagator is to find a balance between accuracy and efficiency. Given the need for efficient propagation, only important orbital perturbations are calculated in SGP4. These perturbations include:
- The J2, J3, and J4 Earth zonal harmonic effects
- The atmospheric drag effects
- The 3rd body effects of the sun and moon
- The J22, J31, J32, J33, J44, J52, and J54 Earth tesseral effects

The SGP4 algorithm can be broken into two primary phases:
1. Initialization - Calculating the time-independent propagation terms
2. Propagation - Calculating the satellite state at a given time

### Initialization

### Propagation

## Thanks
If you've found this section, odds are you've read a lot of what has been written here and made use of mako-sgp4. This has been a passion project of mine for a while, and a lot of hard work went into building this up to the state you see today. I'd like to extend my appreciation to you, the reader, for your attention to this work. I've learned a lot in the process of crafting this codebase, and I hope you've found it useful as well. 

Thank you!

## References
- [1] [Solution of the Problem of Artificial Satellite Theory Without Drag by Brouwer](https://ui.adsabs.harvard.edu/abs/1959AJ.....64..378B/abstract)
- [2] [The Motion of a Close Earth Satellite by Kozai](https://ui.adsabs.harvard.edu/abs/1959AJ.....64..367K/abstract)
- [3] [Theoretical Evaluation of Atmospheric Drag Effects in the Motion of an Artificial Satellite by Brouwer et al](https://scixplorer.org/abs/1961AJ.....66..193B/abstract)
- [4] [An Improved Analytical Drag Theory for the Artificial Satellite Problem by Lane et al](https://arc.aiaa.org/doi/abs/10.2514/6.1969-925)
- [5] [Small Eccentricities or Inclinations in the Brouwer Theory of the Artificial Satellite by Lyddane](https://ui.adsabs.harvard.edu/abs/1963AJ.....68..555L/abstract)
- [6] [General Perturbations Theories Derived from the 1965 Lane Drag Theory by Lane et al]()
- [7] [A First Order Semi-Analytic Perturbation Theory for Highly Eccentric 12 Hour Resonating Satellite Orbits by Bowman]()
- [8] [A Restricted Four Body Solution for Resonating Satellites Without Drag by Hujsak](https://ui.adsabs.harvard.edu/abs/1979aiaa.confV....H/abstract)
- [9] [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
- [10] [The mako-sgp4 Github Repository by Paral](https://github.com/markparal/mako-sgp4)
- [11] [CelesTrak Frequently Asked Questions: Two-Line Element Set Format by Kelso](https://celestrak.org/columns/v04n03/#FAQ04)
- [12] [Space Force Alpha-5 Standard for TLEs](https://www.space-track.org/documentation#/tle-alpha5)
- [13] [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)