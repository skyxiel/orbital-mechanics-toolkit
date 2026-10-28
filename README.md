# Orbital Mechanics Toolkit

Keplerian orbit calculator with Hohmann transfer planning, ground track projection and delta-v budgeting — built to get a feel for the maths behind mission design.

![Hohmann transfer](hohmann_transfer.png)
![ISS ground track](ground_track_iss.png)

## What it does

| Feature | Detail |
|---|---|
| **Hohmann transfers** | Two-impulse transfer between circular orbits — departure/arrival burns, transfer time, orbit geometry |
| **Orbit properties** | Period, circular velocity, apoapsis/periapsis for any Keplerian orbit |
| **Ground tracks** | Kepler-equation solver → latitude/longitude projection with Earth rotation, works for circular and eccentric orbits |
| **Rocket equation** | Tsiolkovsky delta-v from Isp and mass ratio, or required mass ratio for a given delta-v |
| **Escape velocity** | From any altitude |
| **Visualisation** | Dark-theme Matplotlib plots for transfer orbits and ground tracks |

## Example output

```
--- Hohmann Transfer: LEO (400 km) → GEO (35786 km) ---
  Δv₁ (departure burn)  : 2,427.7 m/s
  Δv₂ (arrival burn)    : 1,467.0 m/s
  Total Δv              : 3,894.7 m/s
  Transfer time         : 5.26 hours

--- Orbit Properties ---
  ISS (408 km)            v = 7,667 m/s   T = 1.54 h
  GPS (20200 km)          v = 3,874 m/s   T = 11.97 h
  GEO (35786 km)          v = 3,075 m/s   T = 23.93 h
```

## Usage

```bash
pip install -r requirements.txt
python orbital_mechanics.py
```

Generates text output plus PNG plots for:
- LEO → GEO Hohmann transfer orbit
- ISS ground track (51.6° inclination, 3 orbits)
- Molniya ground track (63.4° inclination, e = 0.74)

## Built with

- Python 3.10+
- NumPy
- Matplotlib

## What I learned

Kepler's equation and why you need Newton iteration to solve it. The vis-viva equation tying speed to position and orbit shape. How Earth's rotation under a satellite creates the sinusoidal ground track pattern, and why Molniya orbits look completely different. The Tsiolkovsky equation and why staging exists.
