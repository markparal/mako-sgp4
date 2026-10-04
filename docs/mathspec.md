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

As a member of the space community who makes active use of these SSA resources, [mako-sgp4](https://github.com/markparal/mako-sgp4) ([10]) was developed to gain a deeper understanding of the underlying theory behind this immensely important model. It combines both the theory detailed in *History of Analytical Orbit Modeling in the U.S. Space Surveillance System* by Hoots et al ([9]) and the practical implementation fixes perscribed in *Revisiting Spacetrack Report #3: Rev 3* by Vallado et al ([14]). The rest of this document will dive into this theory to facilitate understanding of the codebase.

## SGP4 Historical Background
For a complete historical rundown of the development of SGP4, it is recommended to read *History of Analytical Orbit Modeling in the U.S. Space Surveillance System* by Hoots et al ([9]), which details the creation and evolution of the U.S. Space Surveillance system. This paper also discusses the various theories and works that contributed to the modern SGP4 algorithm. A (non-exhaustive) list includes:
- The effects of the J2, J3, and J4 Earth zonal harmonics on the orbit of a satellite ([1], [2])
- The effects of atmospheric drag on the orbits of satellites ([3], [4], [6])
- The avoidance of small divisors of eccentricity or sine of inclination in propagation ([5])
- The inclusion of lunar and solar gravitational effects as well as Earth tesseral harmonics ([7], [8])

## GP Element Sets

| Field | Symbol | Units | Description |
| --- | --- | --- | --- |
| Common Name | - | - | The commonly used name for the satellite |
| Satellite Catalog Number | - | - | NORAD satellite catalog number (NORAD ID) |
| Classification | - | - | Security classification (`U` = Unclassified, `C` = Classified, `S` = Secret) |
| International Designator | - | - | International designator in `Y-NP` form, where `Y` is launch year (4+ digits), `N` is launch number of that year (3+ digits), and `P` is piece of launch (1+ characters) |
| Epoch Datetime | $t_{0}$ | UTC | UTC epoch datetime of the GP elements |
| First Derivative of Mean Motion | $\dot{n}_K$ | revolutions/day^2 | Kozai first time derivative of mean motion (unused in SGP4) |
| Second Derivative of Mean Motion | $\ddot{n}_K$ | revolutions/day^3 | Kozai second time derivative of mean motion (unused in SGP4) |
| B* | $B^{*}$ | 1/Earth radii | Atmospheric drag coefficient |
| Ephemeris Type | - | - | Ephemeris type (always zero) |
| Element Set Number | - | - | Element set number |
| Inclination | $i_B$ | degrees | Orbital inclination |
| Right Ascension of Ascending Node | $\Omega_B$ | degrees | Orbital right ascension of the ascending node (RAAN) |
| Eccentricity | $e_B$ | - | Orbital eccentricity |
| Argument of Perigee | $\omega_B$ | degrees | Orbital argument of perigee |
| Mean Anomaly | $M_B$ | degrees | Orbital mean anomaly |
| Mean Motion | $n_{K}$ | revolutions/day | Kozai Mean motion |
| Revolution Number at Epoch | - | revs | Revolution number at epoch |

<p align="center"><strong>Table 1.</strong> Standard GP element set fields</p>

The units in Table 1 are the units the fields are distributed in. In all equations that follow, the angles $i_B$, $\Omega_B$, $\omega_B$, and $M_B$ are converted to radians and the mean motion $n_K$ is converted to radians/min.

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
| 1 | 33-42 | First Derivative of Mean Motion | revolutions/day^2 | TLE stores value divided by two |
| 1 | 44-51 | Second Derivative of Mean Motion | revolutions/day^3 | TLE stores value divided by six, decimal point assumed |
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
| 2 | 52-62 | Mean Motion | revolutions/day | |
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
| MEAN_ELEMENT_THEORY | - | - | `SGP/SGP4` for the public catalog |
| EPOCH | Epoch Datetime | UTC | ISO-8601 UTC |
| MEAN_MOTION | Mean Motion | revolutions/day | |
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
| MEAN_MOTION_DOT | First Derivative of Mean Motion | revolutions/day^2 | OMM stores value divided by two (Space-Track standard, not CCSDS) |
| MEAN_MOTION_DDOT | Second Derivative of Mean Motion | revolutions/day^3 | OMM stores value divided by six (Space-Track standard, not CCSDS) |

<p align="center"><strong>Table 3.</strong> OMM format</p>

Because the OMM format follows a keyword-value pattern, it exists in multiple general-purpose data formats as well, including KVN (as the example above shows), XML, JSON, and CSV.

## SGP4 Algorithm
As stated previously, the goal of the SGP4 propagator is to find a balance between accuracy and efficiency. Given the need for efficient propagation, only important orbital perturbations are calculated in SGP4. These perturbations include:
- The J2, J3, and J4 Earth zonal harmonic effects
- The atmospheric drag effects
- The 3rd body effects of the sun and moon
- The J22, J31, J32, J33, J44, J52, and J54 Earth tesseral effects

The standard Earth model used by SGP4 is WGS-72 (see Table A1). The standard reference frame used for position and velocity is True Equator Mean Equinox (TEME), an Earth-Centered Inertial (ECI) coordinate frame.

The SGP4 algorithm can be broken into two primary phases:
1. Initialization - Calculating the time-independent propagation terms
2. Propagation - Calculating the satellite state at a given time

### Initialization
Initialization starts from a GP element set (see Table 1). The goal in the initialization process is to calculate the values required for propagation that are independent of time. In mako-sgp4, these values are stored in `Sgp4` structs. 

At a high level, the initialization process can be broken into a series of steps that will be covered individually. These steps are
1. Recover Brouwer mean elements
2. Initialize atmospheric drag parameters
3. Initialize Earth zonal harmonics parameters
4. Initialize lunar and solar third body parameters
5. Initialize Earth half and whole day resonance effects

#### Step 1. Recover Brouwer Mean Elements
We will define the Brouwer mean element set with Table 4 below.

| Element | Symbol | Units | Definition |
| --- | --- | --- | --- |
| Inclination | $i_{B}$ | radians | Orbital inclination |
| Theta | $\theta_{B}$ | - | Cosine of inclination, $\theta_{B} = \cos i_{B}$ |
| Right Ascension of Ascending Node | $\Omega_{B}$ | radians | Right ascension of the ascending node (RAAN) |
| Eccentricity | $e_{B}$ | - | Orbital eccentricity |
| Beta | $\beta_{B}$ | - | $\beta_{B} = \sqrt{1 - e_{B}^{2}}$ |
| Argument of Perigee | $\omega_{B}$ | radians | Argument of perigee |
| Mean Anomaly | $M_{B}$ | radians | Mean anomaly |
| Mean Motion | $n_{B}$ | radians/min | Brouwer mean motion |
| Semi-major Axis | $a_{B}$ | Earth radii | Semi-major axis |
| Period | $T_{B}$ | min | Orbital period |

<p align="center"><strong>Table 4.</strong> Brouwer mean element set</p> 

Many of these values are consistent with what is provided via the GP element set (with unit conversions). The main lift in this step is to convert the Kozai mean motion ($n_{K}$) into the Brouwer mean motion ($n_{B}$). This process is detailed in eqs 1-5.

$$
a_1 = \left(\frac{k_e}{n_{K}}\right)^{2/3} \tag{1}
$$

$$
\delta_1 = \frac{3}{2} \frac{k_2}{a_1^2} \frac{(3 \theta_B^2 - 1)}{(1 - e_B^2)^{3/2}} \tag{2}
$$

$$
a_2 = a_1 (1 - \frac{1}{3} \delta_1 - \delta_1^2 - \frac{134}{81} \delta_1^3) \tag{3}
$$

$$
\delta_0 = \frac{3}{2} \frac{k_2}{a_2^2} \frac{(3 \theta_B^2 - 1)}{(1 - e_B^2)^{3/2}} \tag{4}
$$

$$
n_B = \frac{n_K}{1 + \delta_0} \tag{5}
$$

With the Brouwer mean motion, extract the semi-major axis and orbital period as well with eqs 6-7.

$$
a_B = \left(\frac{k_e}{n_B}\right)^{2/3} \tag{6}
$$

$$
T_B = 2 \pi \frac{\sqrt{(a_B R_e)^3 / \mu_e}}{60} \tag{7}
$$

#### Step 2. Initialize Atmospheric Drag Parameters
We will define the atmospheric drag parameters with Table 5 below.

| Parameter | Symbol | Units | Definition |
| --- | --- | --- | --- |
| Perigee Height | $h_{p}$ | km | Perigee height |
| Reference Distance | $q_{0}$ | Earth radii | Reference distance from the center of the Earth used in the power-law density function (equal to 120km in altitude) |
| s | $s$ | Earth radii | Parameter of the power-law density function |
| Zeta | $\zeta$ | 1 / Earth radii | $\zeta = 1 / (a_{B} - s)$ |
| Eta | $\eta$ | - | $\eta = a_{B} e_{B} \zeta$ |
| C1 | $C_{1}$ | 1 / min | Drag coefficient |
| C3 | $C_{3}$ | Earth radii / min | Drag coefficient |
| C4 | $C_{4}$ | Earth radii / min | Drag coefficient |
| C5 | $C_{5}$ | Earth radii | Drag coefficient |
| D2 | $D_{2}$ | 1 / min^2 | Higher-order drag coefficient |
| D3 | $D_{3}$ | 1 / min^3 | Higher-order drag coefficient |
| D4 | $D_{4}$ | 1 / min^4 | Higher-order drag coefficient |

<p align="center"><strong>Table 5.</strong> Atmospheric drag parameters</p>

Atmospheric drag modeling in SGP4 is based on the power-law density function given in eq 8, where $r$ is the radial distance between the satellite and the center of the Earth, $\rho$ is the atmospheric density, and $\rho_0$ is the reference atmospheric density at a distance $q_0$ from the center of the Earth.

$$
\rho = \rho_0 (q_0 - s)^4 / (r - s)^4 \tag{8}
$$

Because $q_0$ is a reference, it is always the same constant given by eq 9.

$$
q_0 = (120 + R_e) / R_e \tag{9}
$$

The perigee height is given by eq 10.

$$
h_p = R_e \left(a_B (1 - e_B) - 1 \right) \tag{10}
$$

The parameter $s$ is defined as a piecewise function of $h_p$ in eq 11.

$$
s = \begin{cases}
(78 + R_e) / R_e & h_p \ge 156 \\
(h_p - 78 + R_e) / R_e & 98 \le h_p < 156 \\
(20 + R_e) / R_e & h_p < 98
\end{cases} \tag{11}
$$

Additional constants are defined with eqs 12-15 ($A_{3,0}$ is in units of Earth Radii^3)

$$
A_{3,0} = -J_3 R_e^3 / R_e^3 \tag{12}
$$

$$ 
\zeta = \frac{1}{a_B - s} \tag{13}
$$

$$
\eta = a_B e_B \zeta \tag{14}
$$

$$
\psi^2 = |1 - \eta^2| \tag{15}
$$

Drag coefficients are defined with eqs 16-20. $C_3$ is set to zero for small eccentricity to avoid division by $e_B$.

$$
\begin{aligned}
C_2 &= (q_0 - s)^4 \zeta^4 n_B \left(\psi^2\right)^{-7/2} \\
&\qquad \times \Biggl[ a_B \left(1 + \frac{3}{2} \eta^2 + 4 e_B \eta + e_B \eta^3\right) \\
&\qquad + \frac{3}{2} \frac{k_2 \zeta}{\psi^2} \left(-\frac{1}{2} + \frac{3}{2} \theta_B^2\right) \left(8 + 24 \eta^2 + 3 \eta^4\right) \Biggr]
\end{aligned}
\tag{16}
$$

$$
C_1 = B^{*} C_2 \tag{17}
$$

$$
C_3 = \begin{cases}
\dfrac{(q_0 - s)^4 \zeta^5 A_{3,0} n_B \sin i_B}{k_2 e_B} & e_B > 10^{-4} \\
0 & e_B \le 10^{-4}
\end{cases} \tag{18}
$$

$$
\begin{aligned}
C_4 &= 2 n_B (q_0 - s)^4 \zeta^4 a_B \beta_B^2 \left(\psi^2\right)^{-7/2} \\
&\qquad \times \Biggl[ 2\eta \left(1 + e_B \eta\right) + \frac{1}{2} e_B + \frac{1}{2} \eta^3 \\
&\qquad - \frac{2 k_2 \zeta}{a_B \psi^2} \Biggl( \\
&\qquad\qquad 3\left(1 - 3\theta_B^2\right)\left(1 + \frac{3}{2}\eta^2 - 2 e_B \eta - \frac{1}{2} e_B \eta^3\right) \\
&\qquad\qquad + \frac{3}{4}\left(1 - \theta_B^2\right)\left(2\eta^2 - e_B \eta - e_B \eta^3\right) \cos\left(2\omega_B\right) \Biggr) \Biggr]
\end{aligned}
\tag{19}
$$

$$
\begin{aligned}
C_5 &= 2 (q_0 - s)^4 \zeta^4 a_B \beta_B^2 \left(\psi^2\right)^{-7/2} \\
&\qquad \times \left(1 + \frac{11}{4} \eta\left(\eta + e_B\right) + e_B \eta^3\right)
\end{aligned}
\tag{20}
$$

Higher-order drag coefficients are defined with eqs 21-23.

$$
D_2 = 4 a_B \zeta C_1^2 \tag{21}
$$

$$
D_3 = \frac{4}{3} a_B \zeta^2 \left(17 a_B + s\right) C_1^3 \tag{22}
$$

$$
D_4 = \frac{2}{3} a_B^2 \zeta^3 \left(221 a_B + 31 s\right) C_1^4 \tag{23}
$$

#### Step 3. Initialize Earth Zonal Harmonics Parameters
We will define the Earth zonal harmonics parameters with Table 6 below.

| Parameter | Symbol | Units | Definition |
| --- | --- | --- | --- |
| Mean Anomaly Rate | $\dot{M}_{B}$ | radians/min | Secular rate of change of mean anomaly due to the zonal harmonics, excluding the mean motion $n_B$ |
| Argument of Perigee Rate | $\dot{\omega}_{B}$ | radians/min | Secular rate of change of argument of perigee |
| Right Ascension of Ascending Node Rate | $\dot{\Omega}_{B}$ | radians/min | Secular rate of change of the right ascension of the ascending node (RAAN) |

<p align="center"><strong>Table 6.</strong> Earth zonal harmonics parameters</p>

The Earth zonal harmonics in SGP4 consider the impacts of $J_2$ and $J_4$. The resulting secular rates of the Brouwer mean elements are given in eqs 24-26. Note that $\dot{M}_B$ excludes the unperturbed motion $n_B$, so the total secular rate of the mean anomaly is $n_B + \dot{M}_B$.

$$
\begin{aligned}
\dot{M}_B &= n_B \Biggl[
\frac{3 k_2 \left(-1 + 3 \theta_B^{2}\right)}{2 a_B^{2} \beta_B^{3}} \\
&\qquad + \frac{3 k_2^{2} \left(13 - 78 \theta_B^{2} + 137 \theta_B^{4}\right)}{16 a_B^{4} \beta_B^{7}}
\Biggr]
\end{aligned}
\tag{24}
$$

$$
\begin{aligned}
\dot{\omega}_B &= n_B \Biggl[
-\frac{3 k_2 \left(1 - 5 \theta_B^{2}\right)}{2 a_B^{2} \beta_B^{4}} \\
&\qquad + \frac{3 k_2^{2} \left(7 - 114 \theta_B^{2} + 395 \theta_B^{4}\right)}{16 a_B^{4} \beta_B^{8}} \\
&\qquad + \frac{5 k_4 \left(3 - 36 \theta_B^{2} + 49 \theta_B^{4}\right)}{4 a_B^{4} \beta_B^{8}}
\Biggr]
\end{aligned}
\tag{25}
$$

$$
\begin{aligned}
\dot{\Omega}_B &= n_B \Biggl[
-\frac{3 k_2 \theta_B}{a_B^{2} \beta_B^{4}} \\
&\qquad + \frac{3 k_2^{2} \left(4 \theta_B - 19 \theta_B^{3}\right)}{2 a_B^{4} \beta_B^{8}} \\
&\qquad + \frac{5 k_4 \theta_B \left(3 - 7 \theta_B^{2}\right)}{2 a_B^{4} \beta_B^{8}}
\Biggr]
\end{aligned}
\tag{26}
$$

#### Step 4. Initialize Lunar and Solar Third Body Parameters
We will define the lunar and solar third body parameters with Table 7 below.

| Parameter | Symbol | Units | Definition |
| --- | --- | --- | --- |
| Inclination Cosine | $\cos i_{X}$ | - | Cosine of third body orbital inclination |
| Inclination Sine | $\sin i_{X}$ | - | Sine of third body orbital inclination |
| Eccentricity | $e_{X}$ | - | Third body orbital eccentricity |
| Mean Motion | $n_{X}$ | radians/min | Third body mean motion |
| Argument of Perigee Cosine | $\cos \omega_{X}$ | - | Cosine of third body argument of perigee |
| Argument of Perigee Sine | $\sin \omega_{X}$ | - | Sine of third body argument of perigee |
| Right Ascension of Ascending Node | $\Omega_{X}$ | radians | Third body right ascension of the ascending node (RAAN) |
| Mean Anomaly | $M_{X}$ | radians | Third body mean anomaly |
| Beta | $\beta_{X}$ | - | $\beta_{X} = \sqrt{1 - e_{X}^{2}}$ |
| Perturbation Coefficient | $C_{X}$ | radians/min | Third body perturbation coefficient |
| X1 | $x_{1,X}$ | - | Third body geometric coefficient |
| X2 | $x_{2,X}$ | - | Third body geometric coefficient |
| X3 | $x_{3,X}$ | - | Third body geometric coefficient |
| X4 | $x_{4,X}$ | - | Third body geometric coefficient |
| X5 | $x_{5,X}$ | - | Third body geometric coefficient |
| X6 | $x_{6,X}$ | - | Third body geometric coefficient |
| X7 | $x_{7,X}$ | - | Third body geometric coefficient |
| X8 | $x_{8,X}$ | - | Third body geometric coefficient |
| Z1 | $z_{1,X}$ | - | Third body geometric coefficient |
| Z2 | $z_{2,X}$ | - | Third body geometric coefficient |
| Z3 | $z_{3,X}$ | - | Third body geometric coefficient |
| Z11 | $z_{11,X}$ | - | Third body geometric coefficient |
| Z13 | $z_{13,X}$ | - | Third body geometric coefficient |
| Z21 | $z_{21,X}$ | - | Third body geometric coefficient |
| Z23 | $z_{23,X}$ | - | Third body geometric coefficient |
| Z22 | $z_{22,X}$ | - | Third body geometric coefficient |
| Z12 | $z_{12,X}$ | - | Third body geometric coefficient |
| Z31 | $z_{31,X}$ | - | Third body geometric coefficient |
| Z32 | $z_{32,X}$ | - | Third body geometric coefficient |
| Z33 | $z_{33,X}$ | - | Third body geometric coefficient |
| Eccentricity Rate | $\dot{e}_{X}$ | 1 / min | Secular rate of change of satellite eccentricity |
| Inclination Rate | $\dot{i}_{X}$ | radians/min | Secular rate of change of satellite inclination |
| Mean Anomaly Rate | $\dot{M}_{X}$ | radians/min | Secular rate of change of satellite mean anomaly |
| Argument of Perigee Rate | $\dot{\omega}_{X}$ | radians/min | Secular rate of change of satellite argument of perigee |
| Right Ascension of Ascending Node Rate | $\dot{\Omega}_{X}$ | radians/min | Secular rate of change of satellite right ascension of the ascending node (RAAN) |

<p align="center"><strong>Table 7.</strong> Lunar and solar third body parameters (subscript <em>X</em> = <em>M</em> for the Moon, <em>X</em> = <em>S</em> for the Sun)</p>

The lunar and solar third body effects are only considered if the spacecraft has a period greater than or equal to 225 minutes (eq 27). If this is the case, this spacecraft is classified as a "deep space" satellite.

$$
T_B \ge 225 \text{ min} \quad \text{(deep space, } n_B \le 2\pi / 225 \approx 0.0279253 \text{ rad/min)} \tag{27}
$$

The constants used for modeling the orbits and gravitational effects of the Sun and Moon are given in Tables B1 and B2. The time difference between the solar/lunar epoch and the GP element set epoch is defined as $\Delta t = JD_0 - t_{SM}$ in days, where $JD_0$ is the Julian date of the GP element set epoch. We calculate orbital parameters with eqs 28-37.

$$
\Omega_{Me} = \left(\Omega_{Me0} + \dot{\Omega}_{Me0} \Delta t\right) \bmod 2\pi \tag{28}
$$

$$
\cos i_M = 0.91375164 - 0.03568096 \cos \Omega_{Me} \tag{29}
$$

$$
\gamma_M = u_{Me0} + \dot{u}_{Me0} \Delta t \tag{30}
$$

$$
\sin \Omega_M = 0.089683511 \frac{\sin \Omega_{Me}}{\sin i_M} \tag{31}
$$

$$
\Omega_M = \mathrm{atan2}\left(\sin \Omega_M, \cos \Omega_M\right) \tag{32}
$$

$$
z_x = \sin i_S \frac{\sin \Omega_{Me}}{\sin i_M} \tag{33}
$$

$$
z_y = \cos \Omega_M \cos \Omega_{Me} + \cos i_S \sin \Omega_M \sin \Omega_{Me} \tag{34}
$$

$$
\omega_M = \gamma_M + \mathrm{atan2}\left(z_x, z_y\right) - \Omega_{Me} \tag{35}
$$

$$
M_M = \left(M_{M0} + \dot{M}_{M0} \Delta t - \gamma_M\right) \bmod 2\pi \tag{36}
$$

$$
M_S = \left(M_{S0} + \dot{M}_{S0} \Delta t\right) \bmod 2\pi \tag{37}
$$

#### Step 5. Initialize Earth Half and Whole Day Resonance Effects
We will define the Earth resonance parameters with Table 8 below.

| Parameter | Symbol | Units | Definition |
| --- | --- | --- | --- |
| Greenwich Sidereal Time | $\theta_{g}$ | radians | Greenwich mean sidereal time at the GP element set epoch |
| Auxiliary Longitude | $\lambda_{0}$ | radians | Resonance angle at epoch |
| Auxiliary Longitude Rate | $\dot{\lambda}_{0}$ | radians/min | Secular rate of change of the resonance angle, excluding the mean motion $n_B$ |
| D2201 | $D_{2201}$ | radians/min^2 | Half day resonance coefficient |
| D2211 | $D_{2211}$ | radians/min^2 | Half day resonance coefficient |
| D3210 | $D_{3210}$ | radians/min^2 | Half day resonance coefficient |
| D3222 | $D_{3222}$ | radians/min^2 | Half day resonance coefficient |
| D4410 | $D_{4410}$ | radians/min^2 | Half day resonance coefficient |
| D4422 | $D_{4422}$ | radians/min^2 | Half day resonance coefficient |
| D5220 | $D_{5220}$ | radians/min^2 | Half day resonance coefficient |
| D5232 | $D_{5232}$ | radians/min^2 | Half day resonance coefficient |
| D5421 | $D_{5421}$ | radians/min^2 | Half day resonance coefficient |
| D5433 | $D_{5433}$ | radians/min^2 | Half day resonance coefficient |
| Delta1 | $\delta_{1}$ | radians/min^2 | Whole day resonance coefficient |
| Delta2 | $\delta_{2}$ | radians/min^2 | Whole day resonance coefficient |
| Delta3 | $\delta_{3}$ | radians/min^2 | Whole day resonance coefficient |

<p align="center"><strong>Table 8.</strong> Earth resonance parameters</p>

Satellites whose orbital periods are commensurate with the Earth's rotation experience resonant perturbations from the Earth's tesseral harmonics that do not average out over an orbit. SGP4 models two cases: whole day resonance (geosynchronous orbits with periods near 24 hours) and half day resonance (highly eccentric 12 hour orbits, such as Molniya orbits). mako-sgp4 uses the selection criteria of Vallado et al ([14]) given in eqs 38-39. Both ranges correspond to periods above 225 minutes, so resonance effects only apply to deep space satellites.

$$
0.0034906585 < n_B < 0.0052359877 \quad \text{(whole day, } 1200 < T_B < 1800 \text{ min)} \tag{38}
$$

$$
8.26 \times 10^{-3} \le n_B \le 9.24 \times 10^{-3} \text{ and } e_B \ge 0.5 \quad \text{(half day, } 680 \lesssim T_B \lesssim 760.7 \text{ min)} \tag{39}
$$

Both resonances are referenced to the Greenwich mean sidereal time at epoch. This is computed with the IAU-82 model in eqs 40-41, where $T_{UT1}$ is the number of Julian centuries since J2000.0. mako-sgp4 uses the UTC epoch in place of UT1.

$$
T_{UT1} = \frac{JD_0 - 2451545.0}{36525} \tag{40}
$$

$$
\theta_g = \frac{\pi}{180} \cdot \frac{67310.54841 + \left(876600 \cdot 3600 + 8640184.812866\right) T_{UT1} + 0.093104 T_{UT1}^2 - 6.2 \times 10^{-6} T_{UT1}^3}{240} \bmod 2\pi \tag{41}
$$

The Earth rotation rate $\dot{\theta}_E$ and the tesseral resonance constants $Q_{lm}$ and $\lambda_{lm}$ are given in Table C1.

For half day resonance, the functions of inclination are given in eqs 42-51.

$$
F_{220} = \frac{3}{4} \left(1 + \theta_B\right)^2 \tag{42}
$$

$$
F_{221} = \frac{3}{2} \sin^2 i_B \tag{43}
$$

$$
F_{321} = \frac{15}{8} \sin i_B \left(1 - 2\theta_B - 3\theta_B^2\right) \tag{44}
$$

$$
F_{322} = -\frac{15}{8} \sin i_B \left(1 + 2\theta_B - 3\theta_B^2\right) \tag{45}
$$

$$
F_{441} = \frac{105}{4} \sin^2 i_B \left(1 + \theta_B\right)^2 \tag{46}
$$

$$
F_{442} = \frac{315}{8} \sin^4 i_B \tag{47}
$$

$$
F_{522} = \frac{315}{32} \sin i_B \left[\sin^2 i_B \left(1 - 2\theta_B - 5\theta_B^2\right) - \frac{2}{3} + \frac{4}{3}\theta_B + 2\theta_B^2\right] \tag{48}
$$

$$
F_{523} = \frac{105}{16} \sin i_B \left[1 + 2\theta_B - 3\theta_B^2 - \frac{3}{2} \sin^2 i_B \left(1 + 2\theta_B - 5\theta_B^2\right)\right] \tag{49}
$$

$$
F_{542} = \frac{945}{32} \sin i_B \left[2 - 8\theta_B + \theta_B^2 \left(-12 + 8\theta_B + 10\theta_B^2\right)\right] \tag{50}
$$

$$
F_{543} = \frac{945}{32} \sin i_B \left[\theta_B^2 \left(12 + 8\theta_B - 10\theta_B^2\right) - 2 - 8\theta_B\right] \tag{51}
$$

The functions of eccentricity are given below. $G_{201}$ is given as an example by eq 52, and the remaining functions are cubic polynomials in $e_B$ (eq 53) with the coefficients in Table 9.

$$
G_{201} = -0.306 - 0.44 \left(e_B - 0.64\right) \tag{52}
$$

$$
G_{lpq} = g_0 + g_1 e_B + g_2 e_B^2 + g_3 e_B^3 \tag{53}
$$

| Function | Eccentricity Range | $g_0$ | $g_1$ | $g_2$ | $g_3$ |
| --- | --- | --- | --- | --- | --- |
| $G_{211}$ | $e_B \le 0.65$ | 3.616 | -13.247 | 16.29 | 0 |
| $G_{211}$ | $e_B > 0.65$ | -72.099 | 331.819 | -508.738 | 266.724 |
| $G_{310}$ | $e_B \le 0.65$ | -19.302 | 117.39 | -228.419 | 156.591 |
| $G_{310}$ | $e_B > 0.65$ | -346.844 | 1582.851 | -2415.925 | 1246.113 |
| $G_{322}$ | $e_B \le 0.65$ | -18.9068 | 109.7927 | -214.6334 | 146.5816 |
| $G_{322}$ | $e_B > 0.65$ | -342.585 | 1554.908 | -2366.899 | 1215.972 |
| $G_{410}$ | $e_B \le 0.65$ | -41.122 | 242.694 | -471.094 | 313.953 |
| $G_{410}$ | $e_B > 0.65$ | -1052.797 | 4758.686 | -7193.992 | 3651.957 |
| $G_{422}$ | $e_B \le 0.65$ | -146.407 | 841.88 | -1629.014 | 1083.435 |
| $G_{422}$ | $e_B > 0.65$ | -3581.69 | 16178.11 | -24462.77 | 12422.52 |
| $G_{520}$ | $e_B \le 0.65$ | -532.114 | 3017.977 | -5740.032 | 3708.276 |
| $G_{520}$ | $0.65 < e_B < 0.715$ | 1464.74 | -4664.75 | 3763.64 | 0 |
| $G_{520}$ | $e_B \ge 0.715$ | -5149.66 | 29936.92 | -54087.36 | 31324.56 |
| $G_{521}$ | $e_B < 0.7$ | -822.71072 | 4568.6173 | -8491.4146 | 5337.524 |
| $G_{521}$ | $e_B \ge 0.7$ | -51752.104 | 218913.95 | -309468.16 | 146349.42 |
| $G_{532}$ | $e_B < 0.7$ | -853.666 | 4690.25 | -8624.77 | 5341.4 |
| $G_{532}$ | $e_B \ge 0.7$ | -40023.88 | 170470.89 | -242699.48 | 115605.82 |
| $G_{533}$ | $e_B < 0.7$ | -919.2277 | 4988.61 | -9064.77 | 5542.21 |
| $G_{533}$ | $e_B \ge 0.7$ | -37995.78 | 161616.52 | -229838.2 | 109377.94 |

<p align="center"><strong>Table 9.</strong> Half day resonance eccentricity function coefficients</p>

The half day resonance coefficients are then given by eqs 54-63. Note that several of the $D_{44pq}$ and $D_{54pq}$ expressions printed in Hoots et al ([9]) contain typos.

$$
D_{2201} = \frac{3 n_B^2}{a_B^2} Q_{22} F_{220} G_{201} \tag{54}
$$

$$
D_{2211} = \frac{3 n_B^2}{a_B^2} Q_{22} F_{221} G_{211} \tag{55}
$$

$$
D_{3210} = \frac{3 n_B^2}{a_B^3} Q_{32} F_{321} G_{310} \tag{56}
$$

$$
D_{3222} = \frac{3 n_B^2}{a_B^3} Q_{32} F_{322} G_{322} \tag{57}
$$

$$
D_{4410} = \frac{6 n_B^2}{a_B^4} Q_{44} F_{441} G_{410} \tag{58}
$$

$$
D_{4422} = \frac{6 n_B^2}{a_B^4} Q_{44} F_{442} G_{422} \tag{59}
$$

$$
D_{5220} = \frac{3 n_B^2}{a_B^5} Q_{52} F_{522} G_{520} \tag{60}
$$

$$
D_{5232} = \frac{3 n_B^2}{a_B^5} Q_{52} F_{523} G_{532} \tag{61}
$$

$$
D_{5421} = \frac{6 n_B^2}{a_B^5} Q_{54} F_{542} G_{521} \tag{62}
$$

$$
D_{5433} = \frac{6 n_B^2}{a_B^5} Q_{54} F_{543} G_{533} \tag{63}
$$

The half day auxiliary longitude and its secular rate are given by eqs 64-65. The rate combines the zonal (Table 6) and third body (Table 7) secular rates.

$$
\lambda_0 = \left(M_B + 2\Omega_B - 2\theta_g\right) \bmod 2\pi \tag{64}
$$

$$
\dot{\lambda}_0 = \dot{M}_B + \dot{M}_M + \dot{M}_S + 2\left(\dot{\Omega}_B + \dot{\Omega}_M + \dot{\Omega}_S\right) - 2\dot{\theta}_E \tag{65}
$$

For whole day resonance, the functions of inclination and eccentricity are given in eqs 66-71. $F_{220}$ is the same function as in the half day resonance (eq 42), and is repeated here for completeness. However, note that the whole day $G_{310}$ (eq 70) is a different function from the half day $G_{310}$ (Table 9), despite sharing the same name. 

$$
F_{220} = \frac{3}{4} \left(1 + \theta_B\right)^2 \tag{66}
$$

$$
F_{311} = \frac{15}{16} \sin^2 i_B \left(1 + 3\theta_B\right) - \frac{3}{4} \left(1 + \theta_B\right) \tag{67}
$$

$$
F_{330} = \frac{15}{8} \left(1 + \theta_B\right)^3 \tag{68}
$$

$$
G_{200} = 1 - \frac{5}{2} e_B^2 + \frac{13}{16} e_B^4 \tag{69}
$$

$$
G_{310} = 1 + 2 e_B^2 \tag{70}
$$

$$
G_{300} = 1 - 6 e_B^2 + \frac{423}{64} e_B^4 \tag{71}
$$

The whole day resonance coefficients are given by eqs 72-74.

$$
\delta_1 = \frac{3 n_B^2}{a_B^3} F_{311} G_{310} Q_{31} \tag{72}
$$

$$
\delta_2 = \frac{6 n_B^2}{a_B^2} F_{220} G_{200} Q_{22} \tag{73}
$$

$$
\delta_3 = \frac{9 n_B^2}{a_B^3} F_{330} G_{300} Q_{33} \tag{74}
$$

Finally, the whole day auxiliary longitude and its secular rate are given by eqs 75-76.

$$
\lambda_0 = M_B + \Omega_B + \omega_B - \theta_g \tag{75}
$$

$$
\dot{\lambda}_0 = \dot{M}_B + \dot{M}_M + \dot{M}_S + \dot{\Omega}_B + \dot{\Omega}_M + \dot{\Omega}_S + \dot{\omega}_B + \dot{\omega}_M + \dot{\omega}_S - \dot{\theta}_E \tag{76}
$$

During propagation, $\lambda_0$ and $n_B$ are numerically integrated forward in time using these coefficients (see Propagation step 3).

### Propagation
Once a time is provided at which to propagate to, the state of the spacecraft can be calculated using the values found in the initialization process (stored in the `Sgp4` struct).

The propagation process can be broken into a series of steps that will be covered individually. These steps are
1. Account for Earth zonal gravity and partial atmospheric drag effects
2. Account for Lunar and Solar third body effects
3. Account for the whole and half day resonance effects of Earth's gravity
4. Account for remaining atmospheric drag effects
5. Account for long-period periodic effects of lunar and solar gravity
6. Account for long-period periodic effects of Earth's gravity
7. Account for short-period periodic effects of Earth's gravity (solve Kepler's equation)
8. Calculate position and velocity vectors in the TEME frame

## Appendix A: World Geodetic System (WGS) Models

| Variable | Symbol | Units | Value | Description |
| --- | --- | --- | --- | --- |
| `mu` | $\mu_{e}$ | km^3/s^2 | 398600.8 | Standard gravitational parameter |
| `r_earth_eq` | $R_{e}$ | km | 6378.135 | Earth's equatorial radius |
| `j2` | $J_{2}$ | - | 0.001082616 | Second zonal harmonic (Earth's oblateness) |
| `k2` | $k_{2}$ | Earth radii^2 | 0.000541308 | $k_{2} = \frac{1}{2} J_{2}$ |
| `j3` | $J_{3}$ | - | -0.00000253881 | Third zonal harmonic (pear-shaped component) |
| `j4` | $J_{4}$ | - | -0.00000165597 | Fourth zonal harmonic |
| `k4` | $k_{4}$ | Earth radii^4 | 0.00000062098875 | $k_{4} = -\frac{3}{8} J_{4}$ |
| `ke` | $k_{e}$ | Earth radii^1.5 / min | 0.07436691613317 | $k_{e} = 60 \sqrt{\mu_{e} / R_{e}^{3}}$, the square root of $\mu_{e}$ in Earth radii^1.5 / min |

<p align="center"><strong>Table A1.</strong> WGS-72 constants (SGP4 default)</p>

| Variable | Symbol | Units | Value | Description |
| --- | --- | --- | --- | --- |
| `mu` | $\mu_{e}$ | km^3/s^2 | 398600.5 | Standard gravitational parameter |
| `r_earth_eq` | $R_{e}$ | km | 6378.137 | Earth's equatorial radius |
| `j2` | $J_{2}$ | - | 0.00108262998905 | Second zonal harmonic (Earth's oblateness) |
| `k2` | $k_{2}$ | Earth radii^2 | 0.000541314994525 | $k_{2} = \frac{1}{2} J_{2}$ |
| `j3` | $J_{3}$ | - | -0.00000253215306 | Third zonal harmonic (pear-shaped component) |
| `j4` | $J_{4}$ | - | -0.00000161098761 | Fourth zonal harmonic |
| `k4` | $k_{4}$ | Earth radii^4 | 0.0000006041203538 | $k_{4} = -\frac{3}{8} J_{4}$ |
| `ke` | $k_{e}$ | Earth radii^1.5 / min | 0.07436685316871 | $k_{e} = 60 \sqrt{\mu_{e} / R_{e}^{3}}$, the square root of $\mu_{e}$ in Earth radii^1.5 / min |

<p align="center"><strong>Table A2.</strong> WGS-84 constants</p>

## Appendix B: Constants for the Sun and Moon

| Variable | Symbol | Units | Value | Description |
| --- | --- | --- | --- | --- |
| `epoch_sm` | $t_{SM}$ | days | 2415020.0 | Lunar/solar element epoch (12/31/1899 12:00:00 UTC), Julian date |
| `sin_i_s` | $\sin i_{S}$ | - | 0.39785416 | Sine of solar orbital inclination |
| `cos_i_s` | $\cos i_{S}$ | - | 0.91744867 | Cosine of solar orbital inclination |
| `e_s` | $e_{S}$ | - | 0.01675 | Solar orbital eccentricity |
| `n_s` | $n_{S}$ | radians/min | 1.19459e-5 | Solar mean motion |
| `raan_s` | $\Omega_{S}$ | radians | 0.0 | Solar right ascension of the ascending node (RAAN) |
| `sin_omega_s` | $\sin \omega_{S}$ | - | -0.98088458 | Sine of solar argument of perigee |
| `cos_omega_s` | $\cos \omega_{S}$ | - | 0.1945905 | Cosine of solar argument of perigee |
| `c_s` | $C_{S}$ | radians/min | 2.9864797e-6 | Solar perturbation coefficient |
| `m_s0` | $M_{S0}$ | radians | 6.2565837 | Solar mean anomaly at the lunar/solar element epoch |
| `m_s0_dot` | $\dot{M}_{S0}$ | radians/day | 0.017201977 | Solar mean anomaly time rate of change at the lunar/solar element epoch |

<p align="center"><strong>Table B1.</strong> Solar model constants</p>

| Variable | Symbol | Units | Value | Description |
| --- | --- | --- | --- | --- |
| `epoch_sm` | $t_{SM}$ | days | 2415020.0 | Lunar/solar element epoch (12/31/1899 12:00:00 UTC), Julian date |
| `e_m` | $e_{M}$ | - | 0.05490 | Lunar orbital eccentricity |
| `n_m` | $n_{M}$ | radians/min | 1.5835218e-4 | Lunar mean motion |
| `c_m` | $C_{M}$ | radians/min | 4.7968065e-7 | Lunar perturbation coefficient |
| `raan_me0` | $\Omega_{Me0}$ | radians | 4.5236020 | Lunar RAAN with respect to the ecliptic plane at the lunar/solar element epoch |
| `raan_me0_dot` | $\dot{\Omega}_{Me0}$ | radians/day | -9.2422029e-4 | Lunar RAAN with respect to the ecliptic plane time rate of change at the lunar/solar element epoch |
| `u_me0` | $u_{Me0}$ | radians | 5.8351514 | Lunar longitude of perigee with respect to the ecliptic plane at the lunar/solar element epoch |
| `u_me0_dot` | $\dot{u}_{Me0}$ | radians/day | 0.0019443680 | Lunar longitude of perigee with respect to the ecliptic plane time rate of change at the lunar/solar element epoch |
| `m_m0` | $M_{M0}$ | radians | 4.7199672 | Lunar mean anomaly at the lunar/solar element epoch |
| `m_m0_dot` | $\dot{M}_{M0}$ | radians/day | 0.22997150 | Lunar mean anomaly time rate of change at the lunar/solar element epoch |

<p align="center"><strong>Table B2.</strong> Lunar model constants</p>

## Appendix C: Constants for Earth Resonance

| Variable | Symbol | Units | Value | Description |
| --- | --- | --- | --- | --- |
| `RPTIM` | $\dot{\theta}_{E}$ | radians/min | 4.3752690880113e-3 | Earth's rotation rate |
| `c22s22`, `q22` | $Q_{22}$ | - | 1.7891679e-6 | Resonance amplitude of the (2, 2) tesseral harmonic |
| `q31` | $Q_{31}$ | - | 2.1460748e-6 | Resonance amplitude of the (3, 1) tesseral harmonic |
| `c32s32` | $Q_{32}$ | - | 3.7393792e-7 | Resonance amplitude of the (3, 2) tesseral harmonic |
| `q33` | $Q_{33}$ | - | 2.2123015e-7 | Resonance amplitude of the (3, 3) tesseral harmonic |
| `c44s44` | $Q_{44}$ | - | 7.3636953e-9 | Resonance amplitude of the (4, 4) tesseral harmonic |
| `c52s52` | $Q_{52}$ | - | 1.1428639e-7 | Resonance amplitude of the (5, 2) tesseral harmonic |
| `c54s54` | $Q_{54}$ | - | 2.1765803e-9 | Resonance amplitude of the (5, 4) tesseral harmonic |
| `lam31` | $\lambda_{31}$ | radians | 0.13130908 | Phase angle of the (3, 1) tesseral harmonic (whole day resonance) |
| `lam22` | $\lambda_{22}$ | radians | 2.88431980 | Phase angle of the (2, 2) tesseral harmonic (whole day resonance) |
| `lam33` | $\lambda_{33}$ | radians | 0.37448087 | Phase angle of the (3, 3) tesseral harmonic (whole day resonance) |

<p align="center"><strong>Table C1.</strong> Earth resonance constants</p>

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
- [14] [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)