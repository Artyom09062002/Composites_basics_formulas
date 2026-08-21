"""
src/materials_db.py — Phase 3: Materials database loader and CLT screening.

Sources: data/materials/materials_db.csv
         SAND2013-2569 (MAT-002), DOE/MSU database (MAT-001, MAT-003),
         manufacturer data (MAT-005, MAT-006).

Usage:
    from src.materials_db import load_materials, get_ply, clt_screen, specific_props
"""

import os
import pandas as pd
import numpy as np

# Path relative to repo root
_DB_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "materials_db.csv")


def load_materials(db_path: str = None) -> pd.DataFrame:
    """
    Load the materials database CSV.
    Returns a DataFrame indexed by material name.
    Skips comment lines (starting with '#').
    """
    path = db_path or _DB_PATH
    df = pd.read_csv(path, comment="#")
    df = df.set_index("name")
    return df


def get_ply(name: str, db: pd.DataFrame = None) -> dict:
    """
    Return a dict of ply properties for a named material.
    Keys match src/clt_engine.py conventions: E1, E2, G12, nu12,
    Xt, Xc, Yt, Yc, S, rho, t_ply.

    Parameters
    ----------
    name : str
        Material name as in materials_db.csv (e.g. 'ELT5500_UD').
    db   : DataFrame, optional
        Pre-loaded DB (avoids re-reading file on every call).
    """
    if db is None:
        db = load_materials()
    if name not in db.index:
        raise KeyError(f"Material '{name}' not found in database. "
                       f"Available: {list(db.index)}")
    row = db.loc[name]
    return {
        "name"     : name,
        "type"     : row["type"],
        "E1"       : row["E1_GPa"] * 1e9,        # [Pa]
        "E2"       : row["E2_GPa"] * 1e9,
        "G12"      : row["G12_GPa"] * 1e9,
        "nu12"     : row["nu12"],
        "Xt"       : row["Xt_MPa"] * 1e6,         # [Pa]
        "Xc"       : row["Xc_MPa"] * 1e6,
        "Yt"       : row["Yt_MPa"] * 1e6,
        "Yc"       : row["Yc_MPa"] * 1e6,
        "S"        : row["S_MPa"] * 1e6,
        "rho"      : row["rho_kgm3"],              # [kg/m³]
        "t_ply"    : row["t_ply_mm"] * 1e-3,       # [m]
        "varim_ok" : row.get("varim_ok", "yes") == "yes",
    }


def specific_props(db: pd.DataFrame = None) -> pd.DataFrame:
    """
    Compute specific stiffness and specific strength for each material.
    Returns DataFrame with columns:
        E1_rho   [MN·m/kg] = (E1/rho) * 1e-3   (conventional unit for composites)
        Xt_rho   [kN·m/kg] = Xt/rho
        Xc_rho   [kN·m/kg] = Xc/rho
    These are the primary CLT screening metrics.
    """
    if db is None:
        db = load_materials()
    out = pd.DataFrame(index=db.index)
    out["E1_GPa"]    = db["E1_GPa"]
    out["rho_kgm3"]  = db["rho_kgm3"]
    out["Xt_MPa"]    = db["Xt_MPa"]
    out["Xc_MPa"]    = db["Xc_MPa"]
    # Specific stiffness: E/ρ in GPa·m³/kg = 10³ kN·m/kg; display as MN·m/kg
    out["E1_specific_MNm_kg"]  = (db["E1_GPa"] * 1e9) / db["rho_kgm3"] / 1e6
    out["Xt_specific_kNm_kg"]  = (db["Xt_MPa"] * 1e6) / db["rho_kgm3"] / 1e3
    out["Xc_specific_kNm_kg"]  = (db["Xc_MPa"] * 1e6) / db["rho_kgm3"] / 1e3
    return out


def clt_screen(db: pd.DataFrame = None, verbose: bool = True) -> pd.DataFrame:
    """
    CLT screening: rank materials by specific stiffness and strength.
    Prints a summary table and returns the DataFrame.

    Screening metrics (higher = better for structural efficiency):
        1. E1/rho — controls tip deflection per unit mass
        2. Xt/rho — controls tensile failure load per unit mass
        3. Xc/rho — controls spar cap compressive failure per unit mass
    """
    if db is None:
        db = load_materials()
    sp = specific_props(db)

    # VARIM processability flag
    sp["varim_ok"] = db.get("varim_ok", pd.Series("yes", index=db.index))

    if verbose:
        print("\n" + "=" * 72)
        print("  CLT MATERIAL SCREENING — Phase 3")
        print("=" * 72)
        print(f"  {'Material':<22} {'E1/rho':>10} {'Xt/rho':>10} {'Xc/rho':>10} {'VARIM':>6}")
        print(f"  {'':22} {'[MN·m/kg]':>10} {'[kN·m/kg]':>10} {'[kN·m/kg]':>10} {'':>6}")
        print("  " + "-" * 62)
        for nm, row in sp.iterrows():
            varim = "yes" if row["varim_ok"] == "yes" else "no"
            print(f"  {nm:<22} {row['E1_specific_MNm_kg']:>10.2f} "
                  f"{row['Xt_specific_kNm_kg']:>10.1f} "
                  f"{row['Xc_specific_kNm_kg']:>10.1f} {varim:>6}")
        print("=" * 72)
        carbon = sp.loc["Newport307_CarbonUD"]
        glass  = sp.loc["ELT5500_UD"]
        print(f"\n  Carbon vs Glass UD:")
        print(f"    E1/rho ratio: {carbon['E1_specific_MNm_kg'] / glass['E1_specific_MNm_kg']:.2f}×")
        print(f"    Xt/rho ratio: {carbon['Xt_specific_kNm_kg'] / glass['Xt_specific_kNm_kg']:.2f}×")
        print(f"    Xc/rho ratio: {carbon['Xc_specific_kNm_kg'] / glass['Xc_specific_kNm_kg']:.2f}×")
        print()
    return sp


if __name__ == "__main__":
    db = load_materials()
    print(f"Loaded {len(db)} materials: {list(db.index)}")
    clt_screen(db)

    print("\nSample ply dict for ELT5500_UD:")
    ply = get_ply("ELT5500_UD", db)
    for k, v in ply.items():
        print(f"  {k:10s} = {v}")
