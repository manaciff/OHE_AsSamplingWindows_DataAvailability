# -*- coding: utf-8 -*-
"""
Export the area of every sampling unit at full precision
========================================================

WHY THIS SCRIPT EXISTS

The column of unit areas in Data_Raw_FINAL.csv (UA_Area_ha) is stored as whole
numbers. Of the 67 units, 34 are exact multiples of 10, 28 of 100 and five of
1,000 (27,000; 36,000; 51,000; 73,000 and 150,000 ha). Areas computed by ArcGIS
do not produce that pattern: the values were rounded somewhere between ArcGIS
and the table.

This matters for one argument of the manuscript, an algebraic identity between
the response and patch density:

    log(OHE) = log(NP) + log(100) - log(PD)

Part of the discrepancy in this identity is unavoidable, because the landscape
metrics are counted on 30 m pixels while the response is the area of the
polygon, so the two do not describe exactly the same surface. The other part
came from the rounding of the response, and it is removed by using the areas at
full precision. With them, the largest discrepancy across the 67 units is 0.028
on the logarithmic scale (Table S28e).

The expected effect on the results is negligible. Rounding to the nearest
100 ha is a relative error of at most about 2.6% in the smallest unit
(1,910 ha), and much less in all the others.

--------------------------------------------------------------------------------
HOW TO USE

1. Open ArcGIS Pro 3.1 and go to Analysis > Python > Python Window, or use the
   project notebook.
2. Adjust the paths of the CONFIGURATION block below, if needed. By default
   they point to the folders of this project.
3. Paste the whole script and run it.
4. It writes Dados/Processados/AOH_area_full_precision.csv and prints a
   unit-by-unit comparison with the stored values.
5. Then run script 01 again (and run_all.R). Script 01 detects the file and
   uses the full-precision areas; if the file is absent, it continues with the
   stored values and says so.

The script does not modify any shapefile; it works on copies in memory.
--------------------------------------------------------------------------------
"""

import os
import csv
import arcpy

# ==============================================================================
# CONFIGURATION
# ==============================================================================
PROJECT_ROOT = r"D:\Duda_Nacif_TCC"

# Folder with the nine shapefiles of sampling units, one per taxon, already
# split into single parts (singlepart). They are the same files from which
# script 01 extracts the centroids.
SHP_DIR = os.path.join(PROJECT_ROOT, "08_Dados_Especies",
                       "Dados_geo_especies", "Sp_data_singlepart")

OUT_CSV = os.path.join(PROJECT_ROOT, "Dados", "Processados",
                       "AOH_area_full_precision.csv")

# Current table, used only for the comparison printed at the end. If it does
# not exist, the script still runs and skips the comparison.
CURRENT_CSV = os.path.join(PROJECT_ROOT, "Dados", "Processados",
                           "Data_Raw_FINAL.csv")

# South America Albers Equal Area Conic, the system stated in Section 2.4 of
# the manuscript. In an equal-area projection, the planar area is the correct
# measure.
ALBERS_WKID = 102033

arcpy.env.overwriteOutput = True

# ==============================================================================
# 1. FIND THE SHAPEFILES
# ==============================================================================
if not os.path.isdir(SHP_DIR):
    raise RuntimeError("Shapefile folder not found: {}".format(SHP_DIR))

shapefiles = sorted(
    os.path.join(SHP_DIR, f) for f in os.listdir(SHP_DIR)
    if f.lower().endswith(".shp")
)
if not shapefiles:
    raise RuntimeError("No .shp file in: {}".format(SHP_DIR))

print("Shapefiles found: {}".format(len(shapefiles)))
for s in shapefiles:
    print("   {}".format(os.path.basename(s)))

albers = arcpy.SpatialReference(ALBERS_WKID)
print("\nTarget coordinate system: {}".format(albers.name))


# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================
def find_field(fields, patterns):
    """Return the first field whose name matches one of the patterns, ignoring
    case. Return None if none matches."""
    lower = {f.name.lower(): f.name for f in fields}
    for p in patterns:
        if p in lower:
            return lower[p]
    for p in patterns:
        for name_l, name in lower.items():
            if p in name_l:
                return name
    return None


def species_from_filename(path):
    """Taxon name from the file name, without the projection suffix. Used only
    when the shapefile has no species field."""
    base = os.path.splitext(os.path.basename(path))[0]
    for suffix in ("_utm", "_albers", "_wgs", "_UTM", "_ALBERS", "_WGS"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    return base


# ==============================================================================
# 3. COMPUTE THE AREAS
# ==============================================================================
rows = []
for shp in shapefiles:
    name = os.path.basename(shp)
    fields = arcpy.ListFields(shp)

    sp_field = find_field(fields, ["sp_id", "species"])
    ua_field = find_field(fields, ["orig_fid", "origfid", "orig_id"])
    if ua_field is None:
        raise RuntimeError(
            "Field ORIG_FID not found in {}. Available fields: {}"
            .format(name, ", ".join(f.name for f in fields))
        )

    # Reprojected copy in memory. The original shapefile is not touched.
    tmp = "memory\\aoh_{}".format(
        os.path.splitext(name)[0].replace("-", "_").replace(".", "_")
    )
    src_sr = arcpy.Describe(shp).spatialReference
    if src_sr.factoryCode != ALBERS_WKID:
        arcpy.management.Project(shp, tmp, albers)
    else:
        arcpy.management.CopyFeatures(shp, tmp)

    # Two measures. In an equal-area projection they almost coincide; the
    # difference between them is a useful check of the coordinate system.
    arcpy.management.AddField(tmp, "AREA_PLN", "DOUBLE")
    arcpy.management.AddField(tmp, "AREA_GEO", "DOUBLE")
    arcpy.management.CalculateGeometryAttributes(
        tmp,
        [["AREA_PLN", "AREA"], ["AREA_GEO", "AREA_GEODESIC"]],
        area_unit="HECTARES",
        coordinate_system=albers,
    )

    read_fields = [ua_field, "AREA_PLN", "AREA_GEO"]
    if sp_field:
        read_fields.insert(0, sp_field)

    n = 0
    with arcpy.da.SearchCursor(tmp, read_fields) as cur:
        for rec in cur:
            if sp_field:
                sp, ua, a_pln, a_geo = rec
            else:
                sp = species_from_filename(shp)
                ua, a_pln, a_geo = rec
            rows.append({
                "SPECIES": str(sp).strip(),
                "UA_ID": int(ua),
                "area_planar_ha": float(a_pln),
                "area_geodesic_ha": float(a_geo),
                "source_shapefile": name,
            })
            n += 1
    arcpy.management.Delete(tmp)
    print("   {}: {} unit(s)".format(name, n))

print("\nTotal units: {}".format(len(rows)))
if len(rows) != 67:
    print("WARNING: the design has 67 units and this script read {}. "
          "Check that the folder holds exactly the nine singlepart "
          "shapefiles used in the analysis.".format(len(rows)))


# ==============================================================================
# 4. WRITE THE CSV AT FULL PRECISION
# ==============================================================================
# repr() of a float writes every significant digit. Do not round here: the
# purpose of the file is to keep the full precision.
out_dir = os.path.dirname(OUT_CSV)
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

rows.sort(key=lambda r: (r["SPECIES"], r["UA_ID"]))
with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["SPECIES", "UA_ID", "area_planar_ha",
                "area_geodesic_ha", "source_shapefile"])
    for r in rows:
        w.writerow([r["SPECIES"], r["UA_ID"],
                    repr(r["area_planar_ha"]), repr(r["area_geodesic_ha"]),
                    r["source_shapefile"]])
print("\nWritten: {}".format(OUT_CSV))


# ==============================================================================
# 5. COMPARISON WITH THE STORED VALUES
# ==============================================================================
if not os.path.isfile(CURRENT_CSV):
    print("\n{} not found; comparison skipped.".format(CURRENT_CSV))
else:
    current = {}
    with open(CURRENT_CSV, "r", encoding="utf-8-sig") as fh:
        rd = csv.DictReader(fh)
        area_col = None
        for cand in ("AOH_unit_area_ha", "AOH_ha", "UA_Area_ha", "EHA_ha"):
            if rd.fieldnames and cand in rd.fieldnames:
                area_col = cand
                break
        if area_col is None:
            print("\nArea column not found in {}. Fields: {}"
                  .format(os.path.basename(CURRENT_CSV), rd.fieldnames))
        else:
            for rec in rd:
                try:
                    key = (str(rec["SPECIES"]).strip(), int(float(rec["UA_ID"])))
                    current[key] = float(rec[area_col])
                except (KeyError, ValueError, TypeError):
                    continue

    if current:
        print("\n{:<16} {:>4} {:>14} {:>16} {:>10}".format(
            "SPECIES", "UA", "stored (ha)", "full precision", "diff %"))
        print("-" * 66)
        diffs, unmatched = [], 0
        for r in rows:
            key = (r["SPECIES"], r["UA_ID"])
            old = current.get(key)
            if old is None:
                unmatched += 1
                continue
            new = r["area_planar_ha"]
            pct = 100.0 * (new - old) / old if old else float("nan")
            diffs.append(abs(pct))
            print("{:<16} {:>4d} {:>14.1f} {:>16.4f} {:>9.3f}%".format(
                r["SPECIES"][:16], r["UA_ID"], old, new, pct))
        if diffs:
            print("-" * 66)
            print("absolute difference: median {:.3f}%  maximum {:.3f}%".format(
                sorted(diffs)[len(diffs) // 2], max(diffs)))
        if unmatched:
            print("{} unit(s) without a match on (SPECIES, UA_ID). This is "
                  "usually a difference in the spelling of the taxon; script 01 "
                  "normalizes the names when it reads the file.".format(unmatched))

print("\nDone. Run script 01 (and run_all.R) next.")
