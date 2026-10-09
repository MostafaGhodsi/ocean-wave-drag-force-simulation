# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:42:23 2026

@author: FATEMEH Moradimoradpour
         Mahsa Zareibagher
         Mostafa 
         Shayan Mohammadsharif
"""
"""
Semester Task - Ocean Waves

This file is based on the original forcesSeaState_FillIn.py template.

It:
1. Creates the Bretschneider spectrum.
2. Checks the resulting significant wave height.
3. Calculates the time series for three phase-shift cases.
4. Calculates surface elevation, horizontal velocity, and drag force.
5. Saves all plots as clean PNG images in the same folder as this file.

Required file in the same folder:
    airyLib.py
"""

from pathlib import Path
from math import ceil

import numpy as np
import matplotlib.pyplot as plt

import airyLib as alib


# =============================================================================
# Output folder
# =============================================================================

# Save all figures in the same folder as this Python file.
outputFolder = Path(__file__).resolve().parent


# =============================================================================
# Phase-shift modes
# =============================================================================

RND_MODE_NONE = 0
RND_MODE_RAND = 1
RND_MODE_PERIODS = 2


# =============================================================================
# Configuration
# =============================================================================

# Gravitational acceleration, water density, and water depth
gAcc = 9.81
rhoW = 1020.0
depth = 98.1

# Wave spectrum
TPeak = 10.0
HSig = 0.98
deltaF = 0.005
fMax = 0.5

# Evaluation position
xPosEval = 0.0
zPosEval = -2.0

# Housing parameters
cDHous = 0.45
diamHous = 0.50
ARefHous = np.pi / 4.0 * diamHous**2


# =============================================================================
# Check parameter validity
# =============================================================================

steepness = HSig / (gAcc * TPeak**2)
relDepth = depth / (gAcc * TPeak**2)

print("Steepness:      ", steepness)
print("Relative depth: ", relDepth)
print("Selected wave theory: Linear Airy wave theory")
print("Selected spectrum: Bretschneider spectrum")
print()


# =============================================================================
# Calculate spectrum
# =============================================================================

nWaves, periods, oms, spectrum, amplis, wLengths, kNums = \
    alib.calcSpectrumBretschneider(
        fMax,
        deltaF,
        HSig,
        TPeak,
        depth,
        gAcc
    )

# Backwards check of the discretised spectrum
deltaOm = 2.0 * np.pi * deltaF
m0 = np.sum(spectrum * deltaOm)
HSigResult = 4.0 * np.sqrt(m0)
HSigError = abs(HSigResult - HSig) / HSig * 100.0

print("Number of wave components:", nWaves)
print("Given Hs:                 ", HSig, "m")
print("Resulting Hs:             ", HSigResult, "m")
print("Relative Hs error:        ", HSigError, "%")
print()


# =============================================================================
# Calculate time values
# =============================================================================

# The lowest frequency is deltaF, therefore the longest period is 1/deltaF.
tEnd = 1.0 / deltaF

# The highest frequency must be resolved with 10 steps per cycle.
shortestPeriod = 1.0 / fMax
dT = shortestPeriod / 10.0

nSteps = ceil(tEnd / dT) + 1
tEval = np.linspace(0.0, tEnd, nSteps)

print("Analysis end time:", tEnd, "s")
print("Time step:        ", dT, "s")
print("Number of steps:  ", nSteps)
print()


# =============================================================================
# Function for one phase-shift case
# =============================================================================

def calculateTimeSeries(rndPhi, randomSeed=2026):
    """
    Calculate surface elevation, horizontal velocity, and drag force.

    Parameters
    ----------
    rndPhi : int
        Phase-shift mode.
    randomSeed : int
        Seed used for the random phase-shift case.

    Returns
    -------
    phis : numpy.ndarray
        Phase shifts [rad].
    etaSer : numpy.ndarray
        Surface elevation time series [m].
    uXSer : numpy.ndarray
        Horizontal velocity time series [m/s].
    fXSer : numpy.ndarray
        Horizontal drag-force time series [N].
    """

    # Define phase shifts
    if rndPhi == RND_MODE_NONE:
        phis = np.zeros(nWaves)

    elif rndPhi == RND_MODE_RAND:
        randomGenerator = np.random.default_rng(randomSeed)
        phis = 2.0 * np.pi * randomGenerator.random(nWaves)

    elif rndPhi == RND_MODE_PERIODS:
        phis = 2.0 * np.pi * periods**2

    else:
        raise ValueError("Unknown phase-shift mode.")

    # Allocate output arrays
    uXSer = np.zeros(nSteps)
    etaSer = np.zeros(nSteps)

    # Calculate the time series
    for iSt in range(nSteps):
        uXSt = 0.0
        etaSt = 0.0

        for iW in range(nWaves):
            etaSt += amplis[iW] * np.cos(
                phis[iW] - oms[iW] * tEval[iSt]
            )

            uXSti, uZSti = alib.calcVelAiry(
                amplis[iW],
                kNums[iW],
                oms[iW],
                phis[iW],
                depth,
                xPosEval,
                zPosEval,
                tEval[iSt]
            )

            uXSt += uXSti

        etaSer[iSt] = etaSt
        uXSer[iSt] = uXSt

    # Calculate the horizontal drag force
    fXSer = (
        0.5
        * rhoW
        * cDHous
        * ARefHous
        * uXSer
        * np.abs(uXSer)
    )

    return phis, etaSer, uXSer, fXSer


# =============================================================================
# Calculate all three phase-shift cases
# =============================================================================

phaseCases = {
    "case_a_zero_phases": {
        "mode": RND_MODE_NONE,
        "title": "Case A - Zero phase shifts"
    },
    "case_b_random_phases": {
        "mode": RND_MODE_RAND,
        "title": "Case B - Random phase shifts"
    },
    "case_c_period_based_phases": {
        "mode": RND_MODE_PERIODS,
        "title": "Case C - Period-based phase shifts"
    }
}

results = {}

for caseName, caseInfo in phaseCases.items():
    phis, etaSer, uXSer, fXSer = calculateTimeSeries(caseInfo["mode"])

    results[caseName] = {
        "title": caseInfo["title"],
        "phis": phis,
        "eta": etaSer,
        "velocity": uXSer,
        "force": fXSer
    }

    print(caseInfo["title"])
    print("Maximum absolute eta:  ", np.max(np.abs(etaSer)), "m")
    print("Maximum absolute ux:   ", np.max(np.abs(uXSer)), "m/s")
    print("Maximum absolute force:", np.max(np.abs(fXSer)), "N")
    print()


# =============================================================================
# Plot 1: Spectrum and amplitude
# =============================================================================

fig, axes = plt.subplots(2, 1, figsize=(10, 9))

axes[0].plot(oms, amplis, "-o", markersize=3)
axes[0].set_xlabel("Angular frequency, omega [rad/s]")
axes[0].set_ylabel("Wave amplitude [m]")
axes[0].set_title("Discrete wave-amplitude spectrum")
axes[0].grid(True)

axes[1].plot(oms, spectrum, "-o", markersize=3)
axes[1].set_xlabel("Angular frequency, omega [rad/s]")
axes[1].set_ylabel("Energy density S(omega) [m^2 s]")
axes[1].set_title("Bretschneider energy-density spectrum")
axes[1].grid(True)

fig.tight_layout()

spectrumFile = outputFolder / "01_spectrum_and_amplitudes.png"
fig.savefig(spectrumFile, dpi=300, bbox_inches="tight")
plt.close(fig)


# =============================================================================
# Plot 2: One figure for each phase-shift case
# =============================================================================

for caseName, caseResult in results.items():
    etaSer = caseResult["eta"]
    uXSer = caseResult["velocity"]
    fXSer = caseResult["force"]

    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)

    axes[0].plot(tEval, etaSer)
    axes[0].axhline(0.5 * HSig, linestyle="--", linewidth=1)
    axes[0].axhline(-0.5 * HSig, linestyle="--", linewidth=1)
    axes[0].set_ylabel("Surface elevation [m]")
    axes[0].set_title(caseResult["title"])
    axes[0].grid(True)

    axes[1].plot(tEval, uXSer)
    axes[1].set_ylabel("Horizontal velocity [m/s]")
    axes[1].grid(True)

    axes[2].plot(tEval, fXSer)
    axes[2].set_xlabel("Time [s]")
    axes[2].set_ylabel("Horizontal drag force [N]")
    axes[2].grid(True)

    fig.tight_layout()

    outputFile = outputFolder / f"02_{caseName}.png"
    fig.savefig(outputFile, dpi=300, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Plot 3: Comparison of surface elevations
# =============================================================================

fig, ax = plt.subplots(figsize=(11, 5))

for caseName, caseResult in results.items():
    ax.plot(
        tEval,
        caseResult["eta"],
        label=caseResult["title"]
    )

ax.set_xlabel("Time [s]")
ax.set_ylabel("Surface elevation [m]")
ax.set_title("Comparison of surface elevation")
ax.grid(True)
ax.legend()

fig.tight_layout()

etaComparisonFile = outputFolder / "03_surface_elevation_comparison.png"
fig.savefig(etaComparisonFile, dpi=300, bbox_inches="tight")
plt.close(fig)


# =============================================================================
# Plot 4: Comparison of horizontal velocities
# =============================================================================

fig, ax = plt.subplots(figsize=(11, 5))

for caseName, caseResult in results.items():
    ax.plot(
        tEval,
        caseResult["velocity"],
        label=caseResult["title"]
    )

ax.set_xlabel("Time [s]")
ax.set_ylabel("Horizontal velocity [m/s]")
ax.set_title("Comparison of horizontal fluid velocity")
ax.grid(True)
ax.legend()

fig.tight_layout()

velocityComparisonFile = outputFolder / "04_velocity_comparison.png"
fig.savefig(velocityComparisonFile, dpi=300, bbox_inches="tight")
plt.close(fig)


# =============================================================================
# Plot 5: Comparison of drag forces
# =============================================================================

fig, ax = plt.subplots(figsize=(11, 5))

for caseName, caseResult in results.items():
    ax.plot(
        tEval,
        caseResult["force"],
        label=caseResult["title"]
    )

ax.set_xlabel("Time [s]")
ax.set_ylabel("Horizontal drag force [N]")
ax.set_title("Comparison of horizontal drag force")
ax.grid(True)
ax.legend()

fig.tight_layout()

forceComparisonFile = outputFolder / "05_drag_force_comparison.png"
fig.savefig(forceComparisonFile, dpi=300, bbox_inches="tight")
plt.close(fig)


# =============================================================================
# Final information
# =============================================================================

print("All calculations were completed successfully.")
print("All plots were saved in:")
print(outputFolder)
print()
print("Saved image files:")
print("01_spectrum_and_amplitudes.png")
print("02_case_a_zero_phases.png")
print("02_case_b_random_phases.png")
print("02_case_c_period_based_phases.png")
print("03_surface_elevation_comparison.png")
print("04_velocity_comparison.png")
print("05_drag_force_comparison.png")
