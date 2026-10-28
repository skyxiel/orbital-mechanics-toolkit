"""
Orbital mechanics toolkit — Keplerian orbits, Hohmann transfers,
ground tracks, and delta-v budgets.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from dataclasses import dataclass

# ── Constants ────────────────────────────────────────────────────────────────

MU_EARTH = 3.986004418e14      # m^3/s^2
R_EARTH = 6.371e6              # m
G0 = 9.80665                   # m/s^2
OMEGA_EARTH = 7.2921159e-5     # rad/s  (Earth rotation rate)

# ── Data classes ─────────────────────────────────────────────────────────────

@dataclass
class Orbit:
    a: float          # semi-major axis (m)
    e: float = 0.0    # eccentricity
    i: float = 0.0    # inclination (deg)

    @property
    def period(self):
        return 2 * np.pi * np.sqrt(self.a**3 / MU_EARTH)

    @property
    def v_circular(self):
        return np.sqrt(MU_EARTH / self.a)

    @property
    def altitude(self):
        return self.a - R_EARTH

    @property
    def apoapsis(self):
        return self.a * (1 + self.e)

    @property
    def periapsis(self):
        return self.a * (1 - self.e)


@dataclass
class HohmannResult:
    dv1: float
    dv2: float
    transfer_time: float
    transfer_orbit: Orbit

    @property
    def total_dv(self):
        return abs(self.dv1) + abs(self.dv2)

# ── Core functions ───────────────────────────────────────────────────────────

def hohmann_transfer(r1: float, r2: float) -> HohmannResult:
    a_transfer = (r1 + r2) / 2
    v1_circular = np.sqrt(MU_EARTH / r1)
    v2_circular = np.sqrt(MU_EARTH / r2)
    v1_transfer = np.sqrt(MU_EARTH * (2 / r1 - 1 / a_transfer))
    v2_transfer = np.sqrt(MU_EARTH * (2 / r2 - 1 / a_transfer))

    dv1 = v1_transfer - v1_circular
    dv2 = v2_circular - v2_transfer
    t_transfer = np.pi * np.sqrt(a_transfer**3 / MU_EARTH)
    e_transfer = abs(r2 - r1) / (r2 + r1)

    return HohmannResult(dv1, dv2, t_transfer, Orbit(a_transfer, e_transfer))


def orbital_velocity(r: float, a: float) -> float:
    return np.sqrt(MU_EARTH * (2 / r - 1 / a))


def vis_viva(r: float, a: float) -> float:
    return MU_EARTH * (2 / r - 1 / a)


def escape_velocity(r: float) -> float:
    return np.sqrt(2 * MU_EARTH / r)


def orbit_positions(orbit: Orbit, n_points: int = 360) -> tuple:
    theta = np.linspace(0, 2 * np.pi, n_points)
    r = orbit.a * (1 - orbit.e**2) / (1 + orbit.e * np.cos(theta))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y


def ground_track(orbit: Orbit, n_orbits: int = 3, n_points: int = 1000):
    t_total = n_orbits * orbit.period
    t = np.linspace(0, t_total, n_points)
    inc = np.radians(orbit.i)

    lats, lons = [], []
    for ti in t:
        M = 2 * np.pi * ti / orbit.period
        # solve Kepler for circular/low-e orbits (Newton iteration)
        E = M
        for _ in range(10):
            E = E - (E - orbit.e * np.sin(E) - M) / (1 - orbit.e * np.cos(E))
        nu = 2 * np.arctan2(
            np.sqrt(1 + orbit.e) * np.sin(E / 2),
            np.sqrt(1 - orbit.e) * np.cos(E / 2),
        )
        lat = np.degrees(np.arcsin(np.sin(inc) * np.sin(nu)))
        lon_inertial = np.degrees(np.arctan2(
            np.cos(inc) * np.sin(nu), np.cos(nu)
        ))
        lon = (lon_inertial - np.degrees(OMEGA_EARTH * ti) + 180) % 360 - 180
        lats.append(lat)
        lons.append(lon)

    return np.array(lats), np.array(lons)


def tsiolkovsky(isp: float, m0: float, mf: float) -> float:
    return isp * G0 * np.log(m0 / mf)


def required_mass_ratio(dv: float, isp: float) -> float:
    return np.exp(dv / (isp * G0))

# ── Plotting ─────────────────────────────────────────────────────────────────

def plot_hohmann(r1: float, r2: float):
    result = hohmann_transfer(r1, r2)
    fig, ax = plt.subplots(figsize=(8, 8))

    earth = Circle((0, 0), R_EARTH / 1e6, color="#4a90d9", alpha=0.8, label="Earth")
    ax.add_patch(earth)

    orbit1 = Orbit(r1)
    orbit2 = Orbit(r2)
    x1, y1 = orbit_positions(orbit1)
    x2, y2 = orbit_positions(orbit2)
    ax.plot(x1 / 1e6, y1 / 1e6, "w-", linewidth=1.5, label=f"LEO ({orbit1.altitude/1e3:.0f} km)")
    ax.plot(x2 / 1e6, y2 / 1e6, "w-", linewidth=1.5, alpha=0.6, label=f"Target ({orbit2.altitude/1e3:.0f} km)")

    xt, yt = orbit_positions(result.transfer_orbit)
    half = len(xt) // 2
    ax.plot(xt[:half] / 1e6, yt[:half] / 1e6, "#ff6b35", linewidth=2.5, linestyle="--", label="Transfer orbit")

    ax.plot(r1 / 1e6, 0, "o", color="#2ecc71", markersize=10, zorder=5)
    ax.plot(-r2 / 1e6, 0, "s", color="#e74c3c", markersize=10, zorder=5)

    ax.set_facecolor("#1a1a2e")
    fig.patch.set_facecolor("#1a1a2e")
    ax.set_aspect("equal")
    ax.set_xlabel("x (thousands of km)", color="white")
    ax.set_ylabel("y (thousands of km)", color="white")
    ax.set_title("Hohmann Transfer Orbit", color="white", fontsize=14, fontweight="bold")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_color("white")
    ax.legend(facecolor="#2a2a3e", edgecolor="white", labelcolor="white", fontsize=9)

    textbox = (
        f"Δv₁ = {result.dv1:.1f} m/s\n"
        f"Δv₂ = {result.dv2:.1f} m/s\n"
        f"Total Δv = {result.total_dv:.1f} m/s\n"
        f"Transfer time = {result.transfer_time/3600:.1f} h"
    )
    ax.text(0.02, 0.98, textbox, transform=ax.transAxes, fontsize=10,
            verticalalignment="top", color="white",
            bbox=dict(boxstyle="round", facecolor="#2a2a3e", alpha=0.9))

    plt.tight_layout()
    return fig


def plot_ground_track(orbit: Orbit, n_orbits: int = 3):
    lats, lons = ground_track(orbit, n_orbits)
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.set_facecolor("#1a1a2e")
    fig.patch.set_facecolor("#1a1a2e")

    # split track at longitude wraps
    split_indices = np.where(np.abs(np.diff(lons)) > 180)[0] + 1
    lon_segments = np.split(lons, split_indices)
    lat_segments = np.split(lats, split_indices)

    for lo, la in zip(lon_segments, lat_segments):
        ax.plot(lo, la, "#ff6b35", linewidth=1.5)

    ax.plot(lons[0], lats[0], "o", color="#2ecc71", markersize=8, zorder=5, label="Start")

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude (°)", color="white")
    ax.set_ylabel("Latitude (°)", color="white")
    ax.set_title(
        f"Ground Track — {orbit.altitude/1e3:.0f} km, {orbit.i}° inclination, {n_orbits} orbits",
        color="white", fontsize=13, fontweight="bold",
    )
    ax.tick_params(colors="white")
    ax.grid(True, alpha=0.2, color="white")
    for spine in ax.spines.values():
        spine.set_color("white")
    ax.legend(facecolor="#2a2a3e", edgecolor="white", labelcolor="white")

    plt.tight_layout()
    return fig

# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  ORBITAL MECHANICS TOOLKIT")
    print("=" * 60)

    # ── LEO to GEO Hohmann transfer ──
    leo_alt = 400e3
    geo_alt = 35_786e3
    r_leo = R_EARTH + leo_alt
    r_geo = R_EARTH + geo_alt

    result = hohmann_transfer(r_leo, r_geo)
    print(f"\n--- Hohmann Transfer: LEO ({leo_alt/1e3:.0f} km) → GEO ({geo_alt/1e3:.0f} km) ---")
    print(f"  Δv₁ (departure burn)  : {result.dv1:,.1f} m/s")
    print(f"  Δv₂ (arrival burn)    : {result.dv2:,.1f} m/s")
    print(f"  Total Δv              : {result.total_dv:,.1f} m/s")
    print(f"  Transfer time         : {result.transfer_time/3600:.2f} hours")

    # ── Orbit comparisons ──
    print(f"\n--- Orbit Properties ---")
    for name, alt in [("ISS (408 km)", 408e3), ("GPS (20200 km)", 20_200e3), ("GEO (35786 km)", 35_786e3)]:
        orb = Orbit(R_EARTH + alt)
        print(f"  {name:20s}  v = {orb.v_circular:,.0f} m/s   T = {orb.period/3600:.2f} h")

    # ── Escape velocity ──
    v_esc = escape_velocity(r_leo)
    print(f"\n--- Escape velocity from {leo_alt/1e3:.0f} km ---")
    print(f"  v_escape = {v_esc:,.0f} m/s")

    # ── Tsiolkovsky rocket equation ──
    isp = 311         # s (Merlin 1D vacuum)
    m0 = 549_054      # kg (Falcon 9 full)
    mf = 75_000       # kg (dry + payload)
    dv = tsiolkovsky(isp, m0, mf)
    print(f"\n--- Tsiolkovsky (Falcon 9-class) ---")
    print(f"  Isp = {isp} s,  m₀ = {m0:,} kg,  mf = {mf:,} kg")
    print(f"  Δv = {dv:,.0f} m/s")

    # ── Required mass ratio for LEO → GEO ──
    mr = required_mass_ratio(result.total_dv, isp)
    print(f"\n  Mass ratio needed for LEO→GEO Δv ({result.total_dv:.0f} m/s): {mr:.2f}")

    # ── Plots ──
    print("\nGenerating plots...")
    fig1 = plot_hohmann(r_leo, r_geo)
    fig1.savefig("hohmann_transfer.png", dpi=150, bbox_inches="tight")
    print("  → hohmann_transfer.png")

    iss = Orbit(R_EARTH + 408e3, e=0.0002, i=51.6)
    fig2 = plot_ground_track(iss)
    fig2.savefig("ground_track_iss.png", dpi=150, bbox_inches="tight")
    print("  → ground_track_iss.png")

    molniya = Orbit(R_EARTH + 24_730e3, e=0.74, i=63.4)
    fig3 = plot_ground_track(molniya, n_orbits=2)
    fig3.savefig("ground_track_molniya.png", dpi=150, bbox_inches="tight")
    print("  → ground_track_molniya.png")

    plt.show()
    print("\nDone.")


if __name__ == "__main__":
    main()
