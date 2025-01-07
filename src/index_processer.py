import os
import geopandas as gpd
import rasterio as rio
from shapely.geometry import box

class IndexProcesser:
    def __init__(self, index_path: str, center_line_path: str, photos_path: str):
        self.index_path: str = index_path
        self.center_line_path: str = center_line_path
        self.photos_path: str = photos_path
        self.epsg: str
        self.index_gdf: gpd.GeoDataFrame = gpd.read_file(self.index_path)
        self.trimmed_index: gpd.GeoDataFrame
        self.raster_box_index: gpd.GeoDataFrame

    def run(self):
        self._trim_index()
        self.create_box_from_extent()

    def _trim_index(self):
        # Open index gdf
        # self.index_gdf = gpd.read_file(self.index_path)
        # Check that the name exist
        try:
            self.trimmed_index = self.index_gdf.sort_values(by=["DATE_PHOTO"])
        except Exception as exc:
            raise ValueError("The name DATE_PHOTO does not exist in your dataset.") from exc
        # Get unique values
        names_compararison = self._get_image_name()
        self.trimmed_index["names_short"] = names_compararison
        self.trimmed_index = self.trimmed_index.drop_duplicates(subset=['names_short'])
        self.trimmed_index["ind"] = [i for i in range(len(self.trimmed_index))]
        self.trimmed_index.index = [i for i in range(len(self.trimmed_index))]
        self.trimmed_index = self._get_neighbours()

    def create_box_from_extent(self):
        modified_gdf = self.trimmed_index.copy()
        # with rio.open(raster_path) as src:
        for i, name in enumerate(["NOM_IMAGE_","NOM_IMAGE"]):
            try:
                self.trimmed_index.sort_values(by=[name])
                break
            except AttributeError as exc:
                if i > 2:
                    raise ValueError("The name NOM_IMAGE_ or NOM_IMAGE_ do not exist in your dataset.") from exc
                # repeat the loop on failure
                continue

        for idx, feature in self.trimmed_index.iterrows():

            raster_name = feature["NOM_IMAGE_"]
            raster_name = raster_name.split("_")
            raster_name.pop()
            raster_name = "_".join(raster_name)
            raster_path = self.complete_filename(self.photos_path,raster_name)
            if raster_path is None:
                raise ValueError(f"The file {raster_path} does not exist")
            # Open the associated raster
            # raster_path = os.path.join(self.photos_path, raster_name)
            with rio.open(raster_path) as src:
                bounds = src.bounds
                geom = box(*bounds)
                modified_gdf.loc[idx, "geometry"] = geom
        shp_path = os.path.split(raster_path)
        modified_gdf.to_file(os.path.join(shp_path[0],"Index/RasterBox_Index.shp"))
        self.raster_box_index = modified_gdf

    @classmethod
    def complete_filename(cls,folder_path: str, filename_start: str)-> str:
        """
        Finds and returns the full name of the first file in the given folder 
        that starts with the specified filename_start.

        Args:
            folder_path: The path to the folder to search.
            filename_start: The starting part of the filename to search for.

        Returns:
            The full name of the file if found, otherwise None.
        """
        for filename in os.listdir(folder_path):
            if filename.startswith(filename_start) and filename.endswith(".tif"):
                return os.path.join(folder_path, filename)
        return None
    
    def _get_neighbours(self) -> gpd.GeoDataFrame:
        # https://gis.stackexchange.com/questions/281652/finding-all-neighbors-using-geopandas
        gdf = self.trimmed_index
        # Drop the column if it exist
        if 'NEIGHBORS' in gdf.columns:
            gdf = gdf.drop(columns=["NEIGHBORS", "KEEP"])
        # add NEIGHBORS column
        gdf = gdf.reindex(columns=gdf.columns.tolist() + ['NEIGH', 'NEIGH_INDEX'])
        gdf["NEIGH"] = ''
        gdf["NEIGH_INDEX"] = ''
        for index, feature in gdf.iterrows():
            # get 'not disjoint' countries
            neighbors = gdf[~gdf.geometry.disjoint(feature.geometry)]["NOM_IMAGE_"].tolist()
            neighbors_index = gdf[~gdf.geometry.disjoint(feature.geometry)].index.tolist()
            neighbors = [name for name in neighbors if feature["NOM_IMAGE_"] != name]
            neighbors_index = [ind for ind in neighbors_index if feature["ind"] != ind]
            # except ValueError:
            if isinstance(neighbors, list):
                gdf.at[index, "NEIGH"] = str(neighbors)
                gdf.at[index, "NEIGH_INDEX"] = str(neighbors_index)
            else:
                gdf.at[index, "NEIGHBORS"] = str(list([int(neighbors)]))
                gdf.at[index, "NEIGH_INDEX"] = str(list([int(neighbors_index)]))
        return gdf
    
    def _get_image_name(self) -> list:
        self.trimmed_index.index = [i for i in range(len(self.trimmed_index))]
        # Check which name is correct
        for i, name in enumerate(["NOM_IMAGE_","NOM_IMAGE"]):
            try:
                name_list = self.trimmed_index[name].values
                break
            except AttributeError as exc:
                if i > 2:
                    raise ValueError("The name NOM_IMAGE_ or NOM_IMAGE_ do not exist in your dataset.") from exc
                # repeat the loop on failure
                continue
            
        names_compararison = []
        for name in name_list:
            name_split = name.split("_")
            name_split.pop()
            name_numbered = "_".join(name_split)
            names_compararison.append(name_numbered)
        # name_list = []
        return names_compararison