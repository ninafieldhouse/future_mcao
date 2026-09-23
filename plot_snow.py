
#%%
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import xarray as xr


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Arial"

map_extent = {
	"JS": [60, 260, 25, 90],
	"NA": [130, 330, 25, 90],
}

central_longitudes = {
	"JS": 120,
	"NA": 270,
}

present_file = "/net/litho/atmosdyn2/nfieldhouse/data/cao/clim/snow_djf_mean_present_1990_2000.nc"
future_file = "/net/litho/atmosdyn2/nfieldhouse/data/cao/clim/snow_djf_mean_future_2090_2100.nc"
fig_path = "/home/nfieldhouse/ma/plots/statistics"


def load_fsno_data(file_path):
	ds = xr.open_dataset(file_path, decode_times=True)
	return ds["FSNO"].isel(time=0)


def create_map_projection(region):
	return ccrs.LambertAzimuthalEqualArea(
		central_longitude=central_longitudes[region],
		central_latitude=60,
	)


def setup_map_axes(ax, region):
    ax.add_feature(cfeature.LAND, facecolor=(217/256, 217/256, 217/256))
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.set_extent(map_extent[region], crs=ccrs.PlateCarree())

    gl = ax.gridlines(draw_labels=True, dms=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = True
    gl.bottom_labels = True
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 30))
    gl.ylocator = mticker.FixedLocator(np.arange(0, 100, 20))
    gl.xlabel_style = {"size": 20, "color": "k"}
    gl.ylabel_style = {"size": 20, "color": "k"}


def get_figsize_from_map_extent(extent, base_height=10):
	lon_min, lon_max, lat_min, lat_max = extent
	if lon_max < lon_min:
		lon_span = (360 + lon_max) - lon_min
	else:
		lon_span = lon_max - lon_min
	lat_span = lat_max - lat_min
	aspect_ratio = lon_span / lat_span
	width = round(aspect_ratio * base_height)
	return width, base_height


def get_colorbar_shrink(extent, base_width=12.0, base_shrink=0.6):
	fig_width = get_figsize_from_map_extent(extent)[0]
	shrink = base_shrink * (base_width / fig_width)
	return min(1.0, round(shrink, 2))


def subset_region(data_array, region):
	lon_min, lon_max, lat_min, lat_max = map_extent[region]
	return data_array.sel(lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max))


def plot_single_map(region, data, title, cmap, vmin, vmax, cbar_label, output_base):
	fig = plt.figure(figsize=(10, 8), layout="constrained")
	fig = plt.figure(
		figsize=get_figsize_from_map_extent(map_extent[region], base_height=6),
		layout="constrained",
	)
	ax = fig.add_subplot(1, 1, 1, projection=create_map_projection(region))
	setup_map_axes(ax, region)
	lon_2d, lat_2d = np.meshgrid(data.lon, data.lat)
	lon_flat = lon_2d.flatten()
	lat_flat = lat_2d.flatten()
	data_flat = data.values.flatten()
	valid = ~np.isnan(data_flat)

	mesh = ax.scatter(
		lon_flat[valid],
		lat_flat[valid],
		c=data_flat[valid],
		cmap=cmap,
		s=6,
		vmin=vmin,
		vmax=vmax,
		transform=ccrs.PlateCarree(),
	)

	cbar = plt.colorbar(
		mesh,
		ax=ax,
		orientation="horizontal",
		pad=0.05,
		aspect=50,
		extend="both",
		shrink=get_colorbar_shrink(map_extent[region], base_width=10.0, base_shrink=1.0),
	)
	cbar.set_label(cbar_label, size=20)
	cbar.ax.tick_params(labelsize=16)
	# ax.set_title(f"{region} {title}", fontsize=20)

	plt.savefig(f"{output_base}.png", dpi=300, bbox_inches="tight")
	plt.savefig(f"{output_base}.pdf", bbox_inches="tight")
	plt.show()


def plot_region_fsno(region, present_fsno, future_fsno):
	present_region = subset_region(present_fsno, region)
	present_region = present_fsno
	future_region = subset_region(future_fsno, region)
	future_region = future_fsno
	diff_region = future_region - present_region
	diff_region = diff_region

	vmax_diff = float(np.nanmax(np.abs(diff_region.values)))
	vmax_diff = max(0.05, np.ceil(vmax_diff / 0.05) * 0.05)

	os.makedirs(fig_path, exist_ok=True)

	plot_single_map(
		region,
		present_region,
		"Present",
		plt.get_cmap("Blues"),
		0.0,
		1.0,
		"fraction of ground covered by snow [-]",
		os.path.join(fig_path, f"snow_fsno_{region}_present"),
	)
	plot_single_map(
		region,
		future_region,
		"Future",
		plt.get_cmap("Blues"),
		0.0,
		1.0,
		"fraction of ground covered by snow[-]",
		os.path.join(fig_path, f"snow_fsno_{region}_future"),
	)
	plot_single_map(
		region,
		diff_region,
		"Future - Present",
		plt.get_cmap("RdBu_r"),
		-vmax_diff,
		vmax_diff,
		"Δ fraction of ground covered by snow [-]",
		os.path.join(fig_path, f"snow_fsno_{region}_diff"),
	)


def main():
	present_fsno = load_fsno_data(present_file)
	future_fsno = load_fsno_data(future_file)

	for region in ["JS", "NA"]:
		plot_region_fsno(region, present_fsno, future_fsno)


if __name__ == "__main__":
	main()

#%%
present_fsno = load_fsno_data(present_file)
future_fsno = load_fsno_data(future_file)
plot_region_fsno("JS", present_fsno, future_fsno)
plot_region_fsno("NA", present_fsno, future_fsno)
# %%
