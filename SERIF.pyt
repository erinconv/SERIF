# -*- coding: utf-8 -*-

import os
import shutil
import math
import arcpy
from arcpy.sa import *
import numpy as np
import pandas as pd
from tqdm import tqdm


def print_and_message(arcpy_message, text_to_print, warning = False):
    if warning:
        arcpy_message.addWarningMessage(text_to_print)
    else:
        arcpy_message.addMEssage(text_to_print)
    print(text_to_print)



class Toolbox:
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "SERIF"
        self.alias = "SERIF"

        # List of tool classes associated with this toolbox
        self.tools = [Segmentation, Hierarchy]


class Segmentation:
    def __init__(self):
        # Predifined Arcpy api parameters
        """Define the tool (tool name is the name of the class)."""
        self.label = "Segmentation"
        self.description = "Automatically perform the river segmentation from an input Ortho-photo"
        
        # Custom parameters
        self._photos_path : str
        self._center_line_path : str
        self._index_path : str
        self._NIR : str
        self._RED : str
        self._NDVI : str
        self._streets : str
        self._GRHQ : str
        self._ndvi_threshold = 0.2
        self._nir_threshold = 77
        self._files_list : list

    def getParameterInfo(self):
        """Define the tool parameters."""
        param0 = arcpy.Parameter(displayName="Path to photos",
                                 name="photos", datatype="DEFolder",
                                 parameterType="Required", direction="Input")
        param1 = arcpy.Parameter(displayName="Centerline path",
                                 name="center_line", datatype="DEFeatureClass", 
                                 parameterType="Required", direction="Input")
        param2 = arcpy.Parameter(displayName="GRHQ lines",
                                 name="grhq", datatype="DEFeatureClass", 
                                 parameterType="Required", direction="Input")
        param3 = arcpy.Parameter(displayName="Street network",
                                 name="street", datatype="DEFeatureClass", 
                                 parameterType="Required", direction="Input")
        param4 = arcpy.Parameter(displayName="NIR thereshold",
                                 name="nir", datatype="GPLong", 
                                 parameterType="Optional", direction="Input")
        param5 = arcpy.Parameter(displayName="NDVI thereshold",
                                 name="ndvi", datatype="GPDouble", 
                                 parameterType="Optional", direction="Input")
        params = [param0, param1, param2, param3, param4, param5]
        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        self._photos_path = parameters[0].valueAsText
        self._center_line_path = parameters[1].valueAsText
        self._GRHQ = parameters[2].valueAsText
        self._streets = parameters[3].valueAsText
        if len(parameters[4].valueAsText) != 0:
            self._nir_threshold = int(parameters[5].valueAsText)
        if len(parameters[4].valueAsText) != 0:
            self._ndvi_threshold = int(parameters[4].valueAsText)
        self._delete_temps()
        self._manage_files()
        self._segmentation(messages)
        print_and_message(messages, "Starting the preparation for manual corrections")
        self._manual_cor_prep()
        print_and_message(messages, "Finishing the preparation for manual corrections")
        self._delete_temps()
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

    def _segmentation(self, messages):
        out_feature_class = os.path.split(self._center_line_path)[0]
        buffer_name = os.path.join(out_feature_class,"buffer_centerline.shp")
        arcpy.analysis.Buffer(self._center_line_path, buffer_name,
                              "100 meters", dissolve_option="ALL")
        # Get raster bands
        self._NIR = os.path.join(self._photos_path,"temp/NIR.tif")
        self._RED = os.path.join(self._photos_path,"temp/RED.tif")
        self._NDVI = os.path.join(self._photos_path,"temp/NDVI.tif")
        mask_path = os.path.join(self._photos_path,"temp/MASK.tif")

        for i in tqdm(range(len(self._files_list))):
            file  = self._files_list[i]
            print_and_message(messages, f"\nDoing segmentation in file: {file} {i+1}/{len(self._files_list)}")
            # Create intersection area of buffered river and the photo to optimize the further calculations
            extent_clip = self.intersect_tif_and_shp(file)

            red = self.get_raster_band(os.path.join(self._photos_path,file),0)
            nir = self.get_raster_band(os.path.join(self._photos_path,file),3)

            # Slice arrays
            red = red[extent_clip.loc[0, "Index"]:extent_clip.loc[2, "Index"],
                      extent_clip.loc[1, "Index"]:extent_clip.loc[3, "Index"]]
            nir = nir[extent_clip.loc[0, "Index"]:extent_clip.loc[2, "Index"],
                      extent_clip.loc[1, "Index"]:extent_clip.loc[3, "Index"]]

            ndvi = (nir.astype(float)-red.astype(float))/(nir.astype(float)+red.astype(float))
            ndvi[np.isnan(ndvi)] = np.nan
            # Save the clipped tif
            self.save_geotif(extent_clip, file, red, self._RED)
            self.save_geotif(extent_clip, file, nir, self._NIR)
            self.save_geotif(extent_clip, file, ndvi, self._NDVI,np.nan)
            # Free the memory
            nir = None
            red = None
            ndvi = None
            # Check if there is files with theresholds in the folder.
            # Otherwise, use the same input values for all the images
            file_csv = os.path.join(self._photos_path,file.replace(".tif",".csv"))
            if os.path.exists(file_csv):
                self._ndvi_threshold, self._nir_threshold = self.read_thresholds(file_csv)
            # Get and save mask
            mask_nir = self.get_mask(self._NIR, self._nir_threshold)
            mask_ndvi = self.get_mask(self._NDVI, self._ndvi_threshold)
            
            mask_nir = np.where((mask_ndvi== 1) &
                             (mask_nir == 1),1,0)
            # Get total mask
            mask_ = mask_nir + mask_ndvi
            # Create buffer mask to filter far pixels from water course
            self.save_geotif(extent_clip, file, mask_, mask_path)
            self.clip_mask_by_buffer(mask_path)
            # Get paths for filteres files
            majority_filtered_path = os.path.join(self._photos_path,f"temp/MASK_majority_{file}")
            clean_bounds_path = os.path.join(self._photos_path,f"temp/MASK_clean_{file}")
            # arcpy.sa.MajorityFilter(mask_path, "FOUR", "MAJORITY")
            majority_filter = arcpy.sa.MajorityFilter(mask_path, "EIGHT", "MAJORITY")
            majority_filter.save(majority_filtered_path)
            clean_bounds = arcpy.sa.BoundaryClean(majority_filtered_path, "NO_SORT", "TWO_WAY")
            clean_bounds.noDataValue = 0
            clean_bounds.save(clean_bounds_path)
            # Save segments
            file_name_shp = file.replace(".tif",".shp")
            poly_seg_path = os.path.join(self._photos_path,f"shp/SEG_{file_name_shp}")
            arcpy.conversion.RasterToPolygon(in_raster=clean_bounds, out_polygon_features=poly_seg_path)
            # Dissolve all features to create a single feature shp
            chanel_feature = os.path.join(self._photos_path,f"shp/MK_{file_name_shp}")
            arcpy.management.Dissolve(poly_seg_path,chanel_feature)
        return 0

    def _manual_cor_prep(self):
        # Prepare the data for manual correction
        # List all the SEG files
        shp_path = os.path.join(self._photos_path, "shp")
        SEG_list = self.list_files_starting_with(shp_path, "SEG")
        MK_list = self.list_files_starting_with(shp_path, "MK")
        merged_features = os.path.join(shp_path,"Merged_SEG.shp")
        # Dissolve features
        if os.path.isfile(merged_features):
            self.delete_files_starting_with(shp_path,"Merged_SEG")
        arcpy.management.Merge(SEG_list, merged_features)
        dissolved_features = os.path.join(shp_path,"Dissolved_SEG.shp")
        if os.path.isfile(dissolved_features):
            self.delete_files_starting_with(shp_path,"Dissolved_SEG")
        arcpy.Dissolve_management(merged_features, dissolved_features, "gridcode")
        # Compute field
        # field_DN = arcpy.management.AddField(merged_features, 
        field_DN = arcpy.management.AddField(dissolved_features, 
                                             field_name="DN", 
                                             field_type="SHORT", 
                                             field_is_nullable="NULLABLE", 
                                             field_is_required="REQUIRED")[0]
        vertex_number = arcpy.management.CalculateGeometryAttributes(in_features=field_DN, 
                                                                     geometry_property=[["area", "AREA"], ["NVER", "POINT_COUNT"], ["length", "PERIMETER_LENGTH"]], 
                                                                     area_unit="SQUARE_METERS")[0]
        complexity = arcpy.management.CalculateField(in_table=vertex_number,
                                                     field="COMPLEX", 
                                                     expression="!LENGTH! / !NVER!")[0]
        
        # Create folders within the shp folder to keep order when creating files
        if not os.path.exists(os.path.join(self._photos_path,"shp/water")):
            os.mkdir(os.path.join(self._photos_path, "shp/water"))
        
        if not os.path.exists(os.path.join(self._photos_path,"shp/bank")):
            os.mkdir(os.path.join(self._photos_path, "shp/bank"))
        
        # start creating the BANKS and WATER features
        bank_layer, _ = arcpy.management.SelectLayerByAttribute(in_layer_or_view=complexity,
                                                                    where_clause="gridcode = 1")
        bank_raw_name = os.path.join(os.path.join(self._photos_path, "shp/bank/bank_raw.shp"))
        arcpy.management.CopyFeatures(bank_layer, bank_raw_name)
        # Get all the polygons that intersects the river center line
        polygon_intersects_centerline, _, _ = arcpy.management.SelectLayerByLocation(in_layer=complexity,
                                                                                     select_features=self._center_line_path)
        
        # Get the valid polygons that are at leats 10000m2
        valid_poly, _ = arcpy.management.SelectLayerByAttribute(in_layer_or_view=polygon_intersects_centerline,
                                                                      selection_type="SUBSET_SELECTION",
                                                                      where_clause="gridcode = 2 And area >= 10000")
        # Detele polygons that are not valid
        water_polygon = os.path.join(self._photos_path,"shp/water/WATER_L1000.shp")
        arcpy.management.EliminatePolygonPart(in_features=valid_poly,
                                              out_feature_class=water_polygon,
                                              part_area="1000 SquareMeters",
                                              part_option="ANY")
        # Get water polygons
        water_DN_1 = arcpy.management.CalculateField(in_table=water_polygon, field="DN", expression="1")[0]
        water_filtered = os.path.join(self._photos_path,"shp/water/water_filtered.shp")
        arcpy.management.CopyFeatures(in_features=water_DN_1, out_feature_class=water_filtered)
        clipped_watter_GRHQ = os.path.join(self._photos_path,"shp/water/water_GRHQ.shp")
        arcpy.analysis.Clip(in_features=water_filtered,
                            clip_features=self._GRHQ,
                            out_feature_class=clipped_watter_GRHQ)
        
        # Intersect banks with water
        bank_water_intersect, _, _ = arcpy.management.SelectLayerByLocation(in_layer=bank_layer,
                                                                            overlap_type="INTERSECT",
                                                                            select_features=clipped_watter_GRHQ,
                                                                            selection_type="SUBSET_SELECTION")
        clipped_watter_GRHQ_uniques = os.path.join(self._photos_path,"shp/water/water.shp")
        arcpy.management.MultipartToSinglepart(in_features=clipped_watter_GRHQ,
                                               out_feature_class=clipped_watter_GRHQ_uniques)
        # Intersect banks with streets
        bank_street_intersect, _, _ = arcpy.management.SelectLayerByLocation(in_layer=bank_water_intersect,
                                                                             overlap_type="INTERSECT", 
                                                                             select_features=self._streets,
                                                                             selection_type="REMOVE_FROM_SELECTION")
        banks_clean = os.path.join(self._photos_path,"shp/bank/banks_clean.shp")
        arcpy.management.CopyFeatures(in_features=bank_street_intersect,
                                      out_feature_class=banks_clean)
        
        banks_clean_unique = os.path.join(self._photos_path,"shp/bank/banks.shp")
        arcpy.management.MultipartToSinglepart(in_features=banks_clean,
                                               out_feature_class=banks_clean_unique)
        # Get the final layers
        water_final, banks_final = arcpy.cartography.SmoothSharedEdges(in_features=[clipped_watter_GRHQ_uniques],
                                                                       algorithm="PAEK",
                                                                       tolerance="5 Meters",
                                                                       shared_edge_features=[banks_clean_unique])
        # Fusion water and banks
        field_mappings=f"gridcode \"gridcode\" true true false 10 Long 0 10,First,#,{clipped_watter_GRHQ_uniques},gridcode,-1,-1,{banks_clean_unique},gridcode,-1,-1;DN \"DN\" true true false 5 Long 0 5,First,#,{clipped_watter_GRHQ_uniques},DN,-1,-1,{banks_clean_unique},DN,-1,-1;area \"area\" true true false 19 Double 0 0,First,#,{clipped_watter_GRHQ_uniques},area,-1,-1,{banks_clean_unique},area,-1,-1;NVER \"NVER\" true true false 10 Long 0 10,First,#,{clipped_watter_GRHQ_uniques},NVER,-1,-1,{banks_clean_unique},NVER,-1,-1;COMPLEX \"COMPLEX\" true true false 254 Text 0 0,First,#,{clipped_watter_GRHQ_uniques},COMPLEX,0,253,{banks_clean_unique},COMPLEX,0,253;ORIG_FID \"ORIG_FID\" true true false 10 Long 0 10,First,#,{clipped_watter_GRHQ_uniques},ORIG_FID,-1,-1,{banks_clean_unique},ORIG_FID,-1,-1;length \"length\" true true false 19 Double 0 0,First,#,{clipped_watter_GRHQ_uniques},length,-1,-1,{banks_clean_unique},length,-1,-1;area \"area\" true true false 19 Double 0 0,First,#,{clipped_watter_GRHQ_uniques},area,-1,-1,{banks_clean_unique},area,-1,-1"
        fusioned_features = os.path.join(self._photos_path,"shp/fusion_active.shp")
        arcpy.management.Merge(inputs=[banks_clean_unique, clipped_watter_GRHQ_uniques],
                               output=fusioned_features,
                               field_mappings=field_mappings)
        active_shp = os.path.join(self._photos_path,"shp/Active.shp")
        arcpy.management.Dissolve(in_features=fusioned_features,
                                  out_feature_class=active_shp,
                                  multi_part="SINGLE_PART")
        return 0
    

    def _fusion_features(self):
        
        # inputs=[BANKS, WATER],
        # output=MERGED,
        return 0

    def _manage_files(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        # List all the .tif files in the folder
        files_list = os.listdir(self._photos_path)
        self._files_list = [file for file in files_list if file.endswith(".tif")]
        # check if folder shp exist. If it does delete it and create it from scratch
        if os.path.exists(os.path.join(self._photos_path, "shp")):
            shutil.rmtree(os.path.join(self._photos_path, "shp"),ignore_errors=True)
            os.mkdir(os.path.join(self._photos_path, "shp"))
        else:
            os.mkdir(os.path.join(self._photos_path, "shp"))

        # Create temporal folder to process files. This will be removed after finishing all the process
        if not os.path.exists(os.path.join(self._photos_path, "temp")):
            os.mkdir(os.path.join(self._photos_path, "temp"))
        return 0

    def _delete_temps(self):
        if os.path.exists(os.path.join(self._photos_path, "temp")):
            shutil.rmtree(os.path.join(self._photos_path, "temp"),ignore_errors=True)


    def save_geotif(self, extent: pd.DataFrame, raster_name: str,
                    array_values: np.ndarray, path: str, noDataValue=0):
        
        src_raster = arcpy.Raster(os.path.join(self._photos_path,raster_name))
        new_origin_x = np.amin([extent.loc[0,"CRS"],extent.loc[2,"CRS"]])
        new_origin_y = np.amin([extent.loc[1,"CRS"],extent.loc[3,"CRS"]])
        myRasterBlock = arcpy.NumPyArrayToRaster(array_values, arcpy.Point(new_origin_x, new_origin_y),
                                                 src_raster.meanCellWidth,
                                                 src_raster.meanCellHeight)
        myRasterBlock.noDataValue = noDataValue
        myRasterBlock.save(path)
        # whereClause = "Value = " + str(noDataValue)
        # SetNull(path,path,whereClause) 
        crs = src_raster.spatialReference
        arcpy.DefineProjection_management(path, crs)
        return 0

    def intersect_tif_and_shp(self, file: str):
        # Read the reference tif file
        src_raster = arcpy.Raster(os.path.join(self._photos_path,file))
        # Get the extent of the raster to clip after the buffered river
        clip_shp  = self.create_raster_extent_shp(src_raster, self._photos_path)
        out_feature_class = os.path.split(self._center_line_path)[0]
        buffer_name = os.path.join(out_feature_class,"buffer_centerline.shp")
        clipped_buffer = os.path.join(self._photos_path,"temp/clipped_buffer.shp")
        arcpy.analysis.Clip(buffer_name, clip_shp, clipped_buffer)
        
        # Create a polygon with the extent of the clipped buffer
        extent_clipped_name = self.create_shp_extent_poly(clipped_buffer,self._photos_path)
        
        # Get the indexes
        indexes = self.get_raster_extent_indexes(src_raster, extent_clipped_name)
        return indexes

    def clip_mask_by_buffer(self, mask_path: str):
        # Read the reference tif file
        # src_raster = arcpy.Raster(mask_path)
        buffer_name = os.path.join(self._photos_path,"temp/clipped_buffer.shp")
        extracted_raster = ExtractByMask(mask_path, buffer_name)
        extracted_raster.save(mask_path)

    @classmethod
    def read_thresholds(cls, file_csv: str)->tuple:
        df = pd.read_csv(file_csv)
        a = (df["ndvi"].values[0], df["nir"].values[0])
        return df["ndvi"].values[0], df["nir"].values[0]


    @classmethod
    def get_mask(cls,raster_name: str, thereshold: float, log_val="l") -> np.ndarray:
        band = arcpy.RasterToNumPyArray(raster_name)
        # band = rgb[:,:]
        if log_val == "l":
            return np.where(band < thereshold, 1, 0)
        else:
            return np.where(band > thereshold, 1, 0)

    @classmethod
    def create_shp_extent_poly(cls, in_shp: str, output_path: str):
        desc = arcpy.Describe(in_shp)
        extent = desc.extent
        # Create a polygon from the extent coordinates
        polygon = arcpy.Polygon(
            arcpy.Array([
                arcpy.Point(extent.XMin, extent.YMin),
                arcpy.Point(extent.XMin, extent.YMax),
                arcpy.Point(extent.XMax, extent.YMax),
                arcpy.Point(extent.XMax, extent.YMin),
                arcpy.Point(extent.XMin, extent.YMin)  # Close the polygon
            ]),
            desc.spatialReference  # Use the spatial reference of the shapefile
        )
        # Save the polygon to a shapefile or feature class
        output_name = os.path.join(output_path,"temp/clipped_extent.shp")
        arcpy.CopyFeatures_management(polygon, output_name)
        return output_name

    @classmethod
    def create_raster_extent_shp(cls, raster: arcpy.Raster, output_path: str):
        extent = raster.extent
        # Create a polygon from the extent
        polygon = arcpy.Polygon(
            arcpy.Array([
                arcpy.Point(extent.XMin, extent.YMin),
                arcpy.Point(extent.XMin, extent.YMax),
                arcpy.Point(extent.XMax, extent.YMax),
                arcpy.Point(extent.XMax, extent.YMin),
                arcpy.Point(extent.XMin, extent.YMin)  # Close the ring
            ]),
            raster.spatialReference  # Match the raster's spatial reference
        )
        # Save the polygon to a shapefile or feature class
        output_name = os.path.join(output_path,"temp/extent_poly.shp")
        arcpy.CopyFeatures_management(polygon, output_name)
        return output_name

    @classmethod
    def get_raster_extent_indexes(cls, raster_src: arcpy.Raster, shp_name: str) -> list:
        # Load the shp
        desc = arcpy.Describe(shp_name)
        extent = desc.extent
        # Get the extent indexes
        col_min = np.amax([math.floor((extent.XMin - raster_src.extent.XMin) / raster_src.meanCellWidth),0])
        row_min = np.amax([math.ceil((raster_src.extent.YMax - extent.YMin) / raster_src.meanCellHeight),0])
        col_max = np.amax([math.floor((extent.XMax - raster_src.extent.XMin) / raster_src.meanCellWidth),0])
        row_max = np.amax([math.ceil((raster_src.extent.YMax - extent.YMax) / raster_src.meanCellHeight),0])
        
        indexes = [np.amin([row_min,row_max]), np.amin([col_min,col_max]),
                   np.amax([row_min,row_max]), np.amax([col_min,col_max])]
        bounds = [extent.XMin,extent.YMin,
                  extent.XMax, extent.YMax]
        
        df = pd.DataFrame(data=np.c_[bounds,indexes],columns=["CRS","Index"])
        df = df.astype({"CRS": float, "Index": int})
        
        # col, row = raster_src.mapToPixel(x_coord, y_coord)
        return df

    @classmethod
    def get_raster_band(cls, raster_name: str, band_index: int):
        rgb = arcpy.RasterToNumPyArray(raster_name)
        band = rgb[band_index,:,:]
        return band

    @ classmethod
    def list_files_starting_with(cls, directory: str, prefix: str):
        """
        Lists all files in a directory that start with the given prefix.

        Args:
            directory: The path to the directory to search.
            prefix: The prefix to match.

        Returns:
            A list of file paths.
        """

        files = [os.path.join(directory,filename) for filename in os.listdir(directory) if filename.startswith(prefix)]
        files = [filename for filename in files if filename.endswith(".shp")]
        # for filename in os.listdir(directory):
        #     if filename.startswith(prefix): files.append()
        return files

    @classmethod
    def delete_files_starting_with(cls, directory: str, prefix: str):
        files = [os.path.join(directory,filename) for filename in os.listdir(directory) if filename.startswith(prefix)]
        # files = [filename for filename in files if filename.endswith(".shp")]
        for fname in files:
            if os.path.isfile(fname): # this makes the code more robust
                os.remove(fname)
                

class Hierarchy:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Hierarchy"
        self.description = "Hierarchy"

    def getParameterInfo(self):
        """Define the tool parameters."""
        params = None
        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

# Class template for future developments
class Tool:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tool"
        self.description = ""

    def getParameterInfo(self):
        """Define the tool parameters."""
        params = None
        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return