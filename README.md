# RoadRAT — Road Risk Assessment Tool

RoadRAT is a physics-based, reduced-complexity modeling framework for assessing the probability of inundation, wave runup, and erosion impacts on coastal roads. It is designed as a regional-scale screening tool to identify vulnerable road segments under current and future climate conditions.

This framework is detailed in:
> Hallin, C., Adell, A., Almström, B., Kroon, A., Larson, M. (2025).  
> RoadRAT – A new framework to assess the probability of inundation, wave runup, and erosion impacting coastal roads.  
> *Coastal Engineering*, 199, 104741. DOI: [10.1016/j.coastaleng.2025.104741](https://doi.org/10.1016/j.coastaleng.2025.104741)

---

## 🔧 Features

- Transect generation in computation points along roads
- Extraction of morphological parameters from input files (DEM, coastline/vegetation line, shoreline)
- Extreme value analysis (GEV) for water levels and runup
- Shoreline change calculation based on historical trends and Bruun rule response to sea-level rise (SLR)
- Protective sediment volume and erosion probability calculations
- User-friendly dictionary-based configuration (`input.txt`)
- Possibility to run timeseries for calibration of erosion coefficient
- Visualizations of outputs (maps, cdf and pdf distributions)
- NetCDF and text export of results

---

## 📦 Requirements

- Python 3.7+
- Required libraries:
  - numpy, pandas, matplotlib, scipy, rasterio, shapefile (`pyshp`), geopandas, shapely, numba, xarray, pyproj

Install all with:

```bash
pip install -r requirements.txt
```

---

## 📁 Input Data

- `input.txt` in the `in/` folder specifies file paths and parameters.
- Required shapefiles:
  - Road polyline (`road_file`)
  - Coastal schematization lines: vegetation (`CL_file`), shoreline (`SH_file`), seaward foreshore limit (`FS_file`)
- Required tif-file:
  - Digital Elevation Model (`DEM_file`)
- Required text files:
  - Water level file(s): `SWL_file_1` and optionally `SWL_file_2`
- Optional: bathyline, SWAN wave data files, line representing non-erodible features, historical coastline

---

## 🚀 Usage

1. Edit `input.txt` in the `in/` folder to set parameters and file paths.
2. Run the model:

```bash
python main_program.py
```

3. Logs are written to:
   - Console (INFO level)
   - File: `logs/road_rat.log`

4. Outputs go to the `out/` folder:
   - Text file
   - NetCDF file
   - Plots (e.g.maps, distributions)

---

## 📊 Outputs

- Maps of road segments with flood/runup/erosion return levels
- NetCDF and `.txt` summaries per computation point
- Distributions and model fit SWL, runup, erosion
---

## 📖 Reference

If you use RoadRAT in your work, please cite:

> Hallin et al. (2025), RoadRAT — A new framework to assess the probability of inundation, wave runup, and erosion impacting coastal roads. *Coastal Engineering*.

---

## 👩‍💻 Contributing

Contributions are welcome. Fork the repository and submit a pull request.

---

## 📄 License

MIT License — see the `LICENSE` file.