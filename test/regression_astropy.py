"""
Optional Astropy comparison for the checked-in wcslib regression fixtures.

Run from the repository root with:

    python3.10 test/regression_astropy.py

The script uses Astropy's FITS 1-based pixel origin (`origin=1`) to match the
public FITSWCS.jl convention and the values stored in regression_wcslib.jl.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


@dataclass(frozen=True)
class RegressionCase:
    name: str
    header: dict[str, object]
    points: tuple[tuple[tuple[float, ...], tuple[float, ...]], ...]
    world_atol: float = 1e-7
    pixel_atol: float = 1e-6


def fits_header(values: dict[str, object]) -> fits.Header:
    """Build a FITS Header while preserving the literal keyword values."""
    header = fits.Header()

    # Astropy's WCS constructor reads the standard FITS keyword/value mapping.
    for key, value in values.items():
        header[key] = value

    return header


def angular_diff_deg(a: float, b: float) -> float:
    """Return the absolute angular separation modulo 360 degrees."""
    return abs((a - b + 180.0) % 360.0 - 180.0)


def point_differences(actual: np.ndarray, expected: tuple[float, ...], wrap_first: bool) -> list[float]:
    """Compare one Astropy result point with one stored reference point."""
    diffs: list[float] = []

    # Longitude-like first axes must be compared modulo a full turn.
    for index, (actual_value, expected_value) in enumerate(zip(actual, expected)):
        if index == 0 and wrap_first:
            diffs.append(angular_diff_deg(float(actual_value), float(expected_value)))
        else:
            diffs.append(abs(float(actual_value) - float(expected_value)))

    return diffs


def case_wcs(case: RegressionCase) -> WCS:
    """Construct an Astropy WCS object for one regression header."""
    return WCS(fits_header(case.header), relax=True)


def check_case(case: RegressionCase) -> list[str]:
    """Return human-readable failure messages for one regression case."""
    wcs = case_wcs(case)
    failures: list[str] = []
    pixels = np.array([point[0] for point in case.points], dtype=float)
    worlds = np.array([point[1] for point in case.points], dtype=float)
    wrap_first = str(case.header.get("CTYPE1", "")).startswith(("RA", "GLON", "ELON"))

    # Forward comparison: pixel -> world must match stored wcslib values.
    astropy_worlds = np.asarray(wcs.all_pix2world(pixels, 1), dtype=float)
    for pixel, expected_world, actual_world in zip(pixels, worlds, astropy_worlds):
        diffs = point_differences(actual_world, tuple(expected_world), wrap_first)
        if any(diff > case.world_atol for diff in diffs):
            failures.append(
                f"{case.name} forward pixel={pixel.tolist()} "
                f"expected={expected_world.tolist()} actual={actual_world.tolist()} diffs={diffs}"
            )

    # Inverse comparison: stored world values must map back to fixture pixels.
    astropy_pixels = np.asarray(wcs.all_world2pix(worlds, 1), dtype=float)
    for expected_pixel, world, actual_pixel in zip(pixels, worlds, astropy_pixels):
        diffs = [abs(float(a) - float(e)) for a, e in zip(actual_pixel, expected_pixel)]
        if any(diff > case.pixel_atol for diff in diffs):
            failures.append(
                f"{case.name} inverse world={world.tolist()} "
                f"expected={expected_pixel.tolist()} actual={actual_pixel.tolist()} diffs={diffs}"
            )

    return failures


def regression_cases() -> tuple[RegressionCase, ...]:
    """Return the wcslib regression fixtures mirrored from regression_wcslib.jl."""
    cdelt_v = 1.0 / 3600.0
    rho = math.radians(45.0)
    theta = math.radians(0.1)
    scale = 0.05 / 3600.0

    # Each case stores representative pixels and the wcslib world values.
    return (
        RegressionCase(
            "TAN CDELT-form",
            {
                "NAXIS": 2,
                "CTYPE1": "RA---TAN",
                "CTYPE2": "DEC--TAN",
                "CRPIX1": 512.0,
                "CRPIX2": 512.0,
                "CRVAL1": 83.8221,
                "CRVAL2": -5.3911,
                "CDELT1": -2.7778e-4,
                "CDELT2": 2.7778e-4,
            },
            (
                ((512.0, 512.0), (83.8221, -5.3911)),
                ((1.0, 1.0), (83.96470930312101, -5.533028256996351)),
                ((1024.0, 1.0), (83.67921161916392, -5.533028190267628)),
                ((1.0, 1024.0), (83.96464257062250, -5.248860779320227)),
                ((1024.0, 1024.0), (83.67927848225278, -5.248860716038349)),
                ((100.0, 200.0), (83.93707010775421, -5.477756332944454)),
            ),
        ),
        RegressionCase(
            "TAN PC-matrix 45 degree rotation",
            {
                "NAXIS": 2,
                "CTYPE1": "RA---TAN",
                "CTYPE2": "DEC--TAN",
                "CRPIX1": 256.0,
                "CRPIX2": 256.0,
                "CRVAL1": 45.0,
                "CRVAL2": 20.0,
                "CDELT1": -cdelt_v,
                "CDELT2": cdelt_v,
                "PC1_1": math.cos(rho),
                "PC1_2": -math.sin(rho),
                "PC2_1": math.sin(rho),
                "PC2_2": math.cos(rho),
            },
            (
                ((256.0, 256.0), (45.0, 20.0)),
                ((300.0, 256.0), (44.99080242788927, 20.008642178799985)),
                ((256.0, 300.0), (45.00919757211073, 20.008642178799985)),
            ),
        ),
        RegressionCase(
            "TAN CD-matrix HST ACS/WFC-like",
            {
                "NAXIS": 2,
                "CTYPE1": "RA---TAN",
                "CTYPE2": "DEC--TAN",
                "CRPIX1": 2048.0,
                "CRPIX2": 1024.0,
                "CRVAL1": 150.0,
                "CRVAL2": 2.5,
                "CD1_1": -scale * math.cos(theta),
                "CD1_2": -scale * math.sin(theta),
                "CD2_1": -scale * math.sin(theta),
                "CD2_2": scale * math.cos(theta),
            },
            (
                ((2048.0, 1024.0), (150.0, 2.5)),
                ((1.0, 1.0), (150.02848210976424, 2.485841002491482)),
                ((4096.0, 2048.0), (149.9715033488136, 2.514172244812679)),
            ),
        ),
        RegressionCase(
            "AIT galactic all-sky",
            {
                "NAXIS": 2,
                "CTYPE1": "GLON-AIT",
                "CTYPE2": "GLAT-AIT",
                "CRPIX1": 360.5,
                "CRPIX2": 180.5,
                "CRVAL1": 0.0,
                "CRVAL2": 0.0,
                "CDELT1": -0.5,
                "CDELT2": 0.5,
            },
            (
                ((360.5, 180.5), (0.0, 0.0)),
                ((300.0, 150.0), (31.167480197179597, -15.155848828445626)),
                ((400.0, 190.0), (340.1744667748081, 4.733615120590103)),
                ((100.0, 180.5), (138.5334189628714, 0.0)),
            ),
        ),
        RegressionCase(
            "SIN projection",
            {
                "NAXIS": 2,
                "CTYPE1": "RA---SIN",
                "CTYPE2": "DEC--SIN",
                "CRPIX1": 100.0,
                "CRPIX2": 100.0,
                "CRVAL1": 180.0,
                "CRVAL2": 30.0,
                "CDELT1": -0.01,
                "CDELT2": 0.01,
            },
            (
                ((100.0, 100.0), (180.0, 30.0)),
                ((150.0, 80.0), (179.42380496785952, 29.79874251298449)),
                ((60.0, 120.0), (180.46281699820287, 30.199192628819894)),
            ),
        ),
        RegressionCase(
            "CAR projection",
            {
                "NAXIS": 2,
                "CTYPE1": "RA---CAR",
                "CTYPE2": "DEC--CAR",
                "CRPIX1": 181.0,
                "CRPIX2": 91.0,
                "CRVAL1": 0.0,
                "CRVAL2": 0.0,
                "CDELT1": -1.0,
                "CDELT2": 1.0,
            },
            (
                ((181.0, 91.0), (0.0, 0.0)),
                ((100.0, 50.0), (81.0, -41.0)),
                ((50.0, 20.0), (131.0, -71.0)),
            ),
            world_atol=1e-10,
            pixel_atol=1e-8,
        ),
        RegressionCase(
            "STG projection",
            {
                "NAXIS": 2,
                "CTYPE1": "RA---STG",
                "CTYPE2": "DEC--STG",
                "CRPIX1": 64.0,
                "CRPIX2": 64.0,
                "CRVAL1": 270.0,
                "CRVAL2": -45.0,
                "CDELT1": -0.1,
                "CDELT2": 0.1,
            },
            (
                ((64.0, 64.0), (270.0, -45.0)),
                ((70.0, 64.0), (269.1515106303453, -44.996858579590125)),
                ((64.0, 70.0), (270.0, -44.400005483023364)),
            ),
        ),
        RegressionCase(
            "3D cube RA+DEC+FREQ",
            {
                "NAXIS": 3,
                "CTYPE1": "RA---TAN",
                "CTYPE2": "DEC--TAN",
                "CTYPE3": "FREQ",
                "CRPIX1": 50.0,
                "CRPIX2": 50.0,
                "CRPIX3": 1.0,
                "CRVAL1": 10.0,
                "CRVAL2": 25.0,
                "CRVAL3": 1.42e9,
                "CDELT1": -0.01,
                "CDELT2": 0.01,
                "CDELT3": 1.0e6,
            },
            (
                ((50.0, 50.0, 1.0), (10.0, 25.0, 1.42e9)),
                ((60.0, 40.0, 5.0), (9.8897520707026, 24.89995959401809, 1.424e9)),
            ),
            pixel_atol=1e-5,
        ),
        RegressionCase(
            "TAN-SIP distortion",
            {
                "NAXIS": 2,
                "CTYPE1": "RA---TAN-SIP",
                "CTYPE2": "DEC--TAN-SIP",
                "CRPIX1": 512.0,
                "CRPIX2": 512.0,
                "CRVAL1": 150.0,
                "CRVAL2": 2.5,
                "CDELT1": -2.7778e-4,
                "CDELT2": 2.7778e-4,
                "A_ORDER": 2,
                "A_2_0": 5.0e-6,
                "A_0_2": 2.0e-6,
                "A_1_1": 0.0,
                "B_ORDER": 2,
                "B_2_0": 1.0e-6,
                "B_0_2": 0.0,
                "B_1_1": 3.0e-6,
            },
            (
                ((512.0, 512.0), (150.0, 2.5)),
                ((600.0, 500.0), (149.97552128963324, 2.4966676835562676)),
                ((400.0, 600.0), (150.03111983052793, 2.524439537711708)),
                ((300.0, 300.0), (150.05885532834, 2.4411593124884754)),
            ),
            pixel_atol=1e-4,
        ),
    )


def main() -> int:
    """Run all optional Astropy comparisons."""
    failures: list[str] = []
    cases = regression_cases()

    # Accumulate every mismatch so a single run reports the complete drift set.
    for case in cases:
        failures.extend(check_case(case))

    if failures:
        print("Astropy regression mismatches:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    npoints = sum(len(case.points) for case in cases)
    print(f"Astropy regression check passed: {len(cases)} cases, {npoints} points")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
