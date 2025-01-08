
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio as rio
from shapely.geometry import box

class Visualisation:
    def __init__(self, center_line_path: str):
        self.center_line_path = center_line_path
        self.dataset_rio: rio.DatasetReader
        self.gdf_buffer: gpd.GeoDataFrame
        self._create_buffer_river_center_line()

    def intersect_tif_and_shp(self) -> pd.DataFrame:
        """_summary_
        Args:
            dataset_rio (rio.DatasetReader): Rasterio dataset reader object containing the raster data
            gdf_buffer (gpd.GeoDataFrame): GeoDataFrame containing the shapefile buffer geometry

        Returns:
            pd.DataFrame: DataFrame containing the intersection bounds and pixel indexes
                Columns:
                - CRS: Coordinates in the raster CRS (minx, miny, maxx, maxy)
                - Index: Corresponding pixel indexes (row_min, col_min, row_max, col_max)
        
        Finds the intersection between a raster dataset and shapefile buffer geometry.
        Projects the shapefile to match the raster CRS if needed.
        Returns the bounds of the intersection in both coordinate space and pixel indexes.
        """
        raster_bounds = self.dataset_rio.bounds
        raster_crs = self.dataset_rio.crs
        shapefile_crs = self.gdf_buffer.crs  # CRS of the shapefile
        # Check if CRS are different
        if shapefile_crs != raster_crs:
            self.gdf_buffer = self.gdf_buffer.to_crs(raster_crs)
            # shapefile_bounds = shapefile.total_bounds
            # shapefile_box = box(*shapefile_bounds)
        shapefile_bounds = self.gdf_buffer.total_bounds  # [minx, miny, maxx, maxy]
        shape_geom :gpd.GeoSeries = self.gdf_buffer.geometry
        # Create a geometry for the raster bounds
        raster_box = box(*raster_bounds)
        # Calculate the intersection
        intersection: gpd.GeoSeries = shape_geom.intersection(raster_box)
        intersection = intersection[~intersection.is_empty]
        bounds = intersection.geometry.bounds.values[0]
        indexes = self.get_indexes_intersection(bounds, self.dataset_rio)
        df = pd.DataFrame(data=np.c_[bounds,indexes],columns=["CRS","Index"])
        df = df.astype({"CRS": float, "Index": int})
        return df

    @classmethod
    def get_indexes_intersection(cls, bounds: list, dataset_rio: rio.DatasetReader) -> list:
        """Gets the pixel indexes corresponding to geographic bounds in a raster.

        Args:
            bounds (list): Geographic bounds [minx, miny, maxx, maxy] in the raster's CRS
            dataset_rio (rio.DatasetReader): Rasterio dataset reader object containing the raster

        Returns:
            list: Pixel indexes [row_min, col_min, row_max, col_max] corresponding to the bounds
                The indexes define a rectangular window containing the bounds.
                row_min/max are the minimum/maximum row indexes
                col_min/max are the minimum/maximum column indexes
        """
        row_min, col_min = dataset_rio.index(np.amin([bounds[0], bounds[2]]),
                                        np.amin([bounds[1], bounds[3]]))
        row_max, col_max = dataset_rio.index(np.amax([bounds[0], bounds[2]]),
                                        np.amax([bounds[1], bounds[3]]))
        return [np.amin([row_min,row_max]), np.amin([col_min,col_max]),
                np.amax([row_min,row_max]), np.amax([col_min,col_max])]
    
    @classmethod
    def get_RGB(cls, src: rio.DatasetReader, extent_clip: pd.DataFrame) -> np.ndarray:
        """Gets the RGB bands from a raster dataset and clips them to the specified extent.

        Args:
            src (rio.DatasetReader): Rasterio dataset reader object containing the raster
            extent_clip (pd.DataFrame): DataFrame with the clipping extent, containing:
                - Index: Row/column indexes defining the clip window [row_min, col_min, row_max, col_max]

        Returns:
            np.ndarray: 3D array containing the clipped RGB bands stacked along the third axis
                       Shape is (height, width, 3) where the last dimension contains (red, green, blue)
        """
        red = src.read(1)
        green = src.read(2)
        blue = src.read(3)
        red = red[extent_clip.loc[0, "Index"]:extent_clip.loc[2, "Index"],
                  extent_clip.loc[1, "Index"]:extent_clip.loc[3, "Index"]]
        green = green[extent_clip.loc[0, "Index"]:extent_clip.loc[2, "Index"],
                      extent_clip.loc[1, "Index"]:extent_clip.loc[3, "Index"]]
        blue = blue[extent_clip.loc[0, "Index"]:extent_clip.loc[2, "Index"],
                    extent_clip.loc[1, "Index"]:extent_clip.loc[3, "Index"]]
        
        return np.dstack((red, green, blue))
    
    def _create_buffer_river_center_line(self):
        """Creates a buffer around the river centerline.
        
        Reads the centerline shapefile specified in self.center_line_path,
        creates a 100m buffer around it, and stores the buffered geometry
        in self.gdf_buffer as a GeoDataFrame.
        """
        gdf_centerline: gpd.GeoDataFrame = gpd.read_file(self.center_line_path)
        geom: gpd.GeoSeries = gdf_centerline.geometry
        buffer: gpd.GeoSeries = geom.buffer(100)
        gdf_centerline.geometry = buffer
        self.gdf_buffer = gdf_centerline
    
    @classmethod
    def get_mask(cls, thereshold: float, array: str, log_val="l")->np.ndarray:
        """Creates a binary mask from an array based on a threshold value.

        Args:
            thereshold (float): Threshold value to compare array values against
            array (str): Input array to create mask from
            log_val (str, optional): Logic for threshold comparison. Defaults to "l".
                "l" means mask where array < threshold
                Any other value means mask where array > threshold

        Returns:
            np.ndarray: Binary mask array where:
                1 indicates pixels meeting the threshold condition
                0 indicates pixels not meeting the threshold condition
        """
        # raster_name = os.path.split(tif_path)
        # with rio.open(tif_path) as src:
        # array = src.read(1)
        if log_val == "l":
            return np.where(array < thereshold, 1, 0)
        else:
            return np.where(array > thereshold, 1, 0)
    
    @classmethod
    def update_plot(cls, nir_threshold: int, 
                    ndvi_threshold: float, 
                    nir: np.ndarray, ndvi: np.ndarray, 
                    im, fig: plt.Figure):
        """Updates the plot display with new threshold masks.

        Updates the displayed mask overlay based on new NIR and NDVI threshold values.
        Creates binary masks from the NIR and NDVI arrays using the threshold values,
        combines them, and updates the plot.

        Args:
            nir_threshold (int): Threshold value for NIR band masking
            ndvi_threshold (float): Threshold value for NDVI masking  
            nir (np.ndarray): Array of NIR band values
            ndvi (np.ndarray): Array of NDVI values
            im: The image object to update
            fig (plt.Figure): The figure containing the plot

        Returns:
            None
        """
        mask_nir = cls.get_mask(nir_threshold, nir)
        mask_ndvi = cls.get_mask(ndvi_threshold, ndvi)
        mask_nir = np.where((mask_ndvi == 1) & (mask_nir == 1), 1, 0)
        mask_ = mask_nir + mask_ndvi
        mask_ = mask_.astype(float)
        mask_[mask_ == 0] = np.nan

        im.set_data(mask_)
        fig.canvas.draw_idle()

    @classmethod
    def calibrate_threshold(cls, rgb: np.ndarray, mask: np.ndarray,
                            ndvi: np.ndarray, nir: np.ndarray,
                            nir_threshold_init: int,
                            ndvi_threshold_init: float):
        fig, ax = plt.subplots(nrows=1,
                            ncols=1,
                            layout="constrained",
                            figsize=(15, 10))
        """Calibrates threshold values for NIR and NDVI masking through an interactive plot.

        Creates an interactive matplotlib figure with sliders to adjust NIR and NDVI threshold values.
        The plot shows the RGB image with a mask overlay that updates in real-time as the thresholds
        are adjusted using the sliders.

        Args:
            rgb (np.ndarray): RGB image array to display as background
            mask (np.ndarray): Initial mask array to overlay
            ndvi (np.ndarray): NDVI values array used for threshold calculations
            nir (np.ndarray): NIR band values array used for threshold calculations  
            nir_threshold_init (int): Initial NIR threshold value
            ndvi_threshold_init (float): Initial NDVI threshold value

        Returns:
            tuple: Final NDVI and NIR threshold values selected via the sliders
        """
        # f = plt.figure()
        ax.imshow(rgb)
        im = ax.imshow(mask,alpha=0.4)
        # im = ax.contour(mask, colors=["white", 'blue', 'red'],
        #               linewidths=[0, 0.5, 0.5])
        # define the values to use for snapping
        allowed_NIR = np.arange(10,200,1)
        # Create sliders
        ax_nir = plt.axes([0.25, 0.1, 0.65, 0.03])
        nir_slider = Slider(
            ax=ax_nir,
            label='NIR Threshold',
            valmin=np.amin(allowed_NIR),
            valmax=np.amax(allowed_NIR),
            valinit=nir_threshold_init,
            valstep=allowed_NIR,
            color="blue"
        )

        ax_ndvi = plt.axes([0.25, 0.05, 0.65, 0.03])
        ndvi_slider = Slider(
            ax=ax_ndvi,
            label='NDVI Threshold',
            valmin=-1,
            valmax=1,
            valinit=ndvi_threshold_init,
            color="green"
        )

        def update(val):
            cls.update_plot(nir_slider.val, ndvi_slider.val, nir, ndvi, im, fig)

        nir_slider.on_changed(update)
        ndvi_slider.on_changed(update)

        cls.update_plot(nir_threshold_init, ndvi_threshold_init, nir, ndvi, im, fig)
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        plt.show()
        return ndvi_slider.val, nir_slider.val
