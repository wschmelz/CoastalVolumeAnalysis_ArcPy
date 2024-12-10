import sys
import os
import arcpy
import math
backslash = '\\'
wkspc = str(os.getcwd()).replace(backslash,"/")
comppar = "Compartment_params.txt"
comppars = "Compartment_params.py"
os.rename(comppar, comppars)
from Compartment_params import *
from arcpy.sa import *
from arcpy import env
import time
os.rename(comppars, comppar)
arcpy.CheckOutExtension("Spatial")
arcpy.CheckOutExtension("3D")
env.overwriteOutput = True
arcpy.overwriteOutput = True
backslash = '\\'
sr = str("PROJCS['NAD_1983_UTM_Zone_18N',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-75.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]")

wkspc = str(os.getcwd()).replace(backslash,"/") 

locbmk = wkspc + '/01_Benchmarks/'
locrwd = wkspc + '/02_RawData/'
locshp = wkspc + '/03_Shapefiles/'
loctin = wkspc + '/04_TINs/'
locmap = wkspc + '/05_ArcMap/'
locras = wkspc + '/06_Rasters/'
locdif = wkspc + '/07_Differences/'
loctab = wkspc + '/08_Tables/'
loccomp = wkspc + '/09_Compartments/'
locvec = wkspc + '/10_Vectors/'
loclyr = wkspc + '/11_Layers/'
locndif = wkspc + '/12_TotalDifferences/'
locntab = wkspc + '/13_TotalTables/'
locnvec = wkspc + '/14_TotalVectors/'
locimg = wkspc + '/15_Images/'

print" "
print "Directory:"
print str(str(os.getcwd()).replace(backslash,"/") + "/")
print" "

baselineradian = math.radians(baselineangle)
baselinelength = compartments*distance
baselinex2 = easting + baselinelength*math.cos(baselineradian)
baseliney2 = northing + baselinelength*math.sin(baselineradian)

baselname = "Baseline.shp"
baseline1 = locvec + baselname
feature_info = [[[easting, northing], [baselinex2, baseliney2]],]

features = []

for feature in feature_info:
	features.append(
		arcpy.Polyline(
			arcpy.Array([arcpy.Point(*coords) for coords in feature])))

arcpy.CopyFeatures_management(features, baseline1)

arcpy.AddField_management(baseline1, 'StartX', 'DOUBLE')
arcpy.AddField_management(baseline1, 'StartY', 'DOUBLE')
arcpy.AddField_management(baseline1, 'EndX', 'DOUBLE')
arcpy.AddField_management(baseline1, 'EndY', 'DOUBLE')

arcpy.CalculateField_management(baseline1, 'StartX', "!SHAPE.firstPoint.X!", "PYTHON_9.3")
arcpy.CalculateField_management(baseline1, 'StartY', "!SHAPE.firstPoint.Y!", "PYTHON_9.3")
arcpy.CalculateField_management(baseline1, 'EndX', "!SHAPE.lastPoint.X!", "PYTHON_9.3")
arcpy.CalculateField_management(baseline1, 'EndY', "!SHAPE.lastPoint.Y!", "PYTHON_9.3")

if compdirection == "left":
	compradian = baselineradian + (math.pi)/2
elif compdirection == "right":
	compradian = baselineradian - (math.pi)/2
else:
	compdirection = raw_input("Shoreline must be specified as located to the left or right of the offshore boundary for compartments, re-enter: LEFT or RIGHT? ").lower()
	if compdirection == "left":
		compradian = baselineradian + (math.pi)/2
	elif compdirection == "right":
		compradian = baselineradian - (math.pi)/2
	else:
		print "Invalid answer, direction must be left or right."

print "Compartment angle is %f." % math.degrees(compradian)

print "Draw compartments."


q = easting

w = northing

e = q + length * (math.cos(compradian))

r = w + length * (math.sin(compradian))

t = e + distance * (math.cos(baselineradian))

y = r + distance * (math.sin(baselineradian))

u = q + distance * (math.cos(baselineradian))

i = w + distance * (math.sin(baselineradian))

n = 2

compname = "Compartments.shp"

comp = loccomp + compname

arcpy.CreateFeatureclass_management(loccomp, compname, "POLYGON")

feature_info = [[[q, w], [e, r], [t, y], [u, i]]]

for n in range (2, compartments + 1):
	iterable = [[[q + distance * (math.cos(baselineradian)), w + distance * (math.sin(baselineradian))], [e + distance * (math.cos(baselineradian)), r + distance * (math.sin(baselineradian))], [t + distance * (math.cos(baselineradian)), y + distance * (math.sin(baselineradian))], [u + distance * (math.cos(baselineradian)), i + distance * (math.sin(baselineradian))]],]
	feature_info.extend(iterable)

	q = q + distance * (math.cos(baselineradian))

	w = w + distance * (math.sin(baselineradian))

	e = e + distance * (math.cos(baselineradian))

	r = r + distance * (math.sin(baselineradian))

	t = t + distance * (math.cos(baselineradian))

	y = y + distance * (math.sin(baselineradian))

	u = u + distance * (math.cos(baselineradian))

	i = i + distance * (math.sin(baselineradian))

	sys.stdout.write("\rProcessed %i of %i" % (n, compartments))

	n = n + 1

print "."

features = []

for feature in feature_info:

	features.append(
		arcpy.Polygon(
			arcpy.Array([arcpy.Point(*coords) for coords in feature])))


arcpy.CopyFeatures_management(features, comp)
arcpy.AddField_management(comp, 'COMP', 'DOUBLE')
arcpy.CalculateField_management(comp, 'COMP', '!FID! +1', "PYTHON_9.3")

easting2 = easting + (-5)*math.cos(baselineradian)
northing2 = northing + (-5)*math.sin(baselineradian)
baselinex22 = easting + (baselinelength+5)*math.cos(baselineradian)
baseliney22 = northing + (baselinelength+5)*math.sin(baselineradian)
compartmentsx1 = baselinex2 + (length+5) * (math.cos(compradian))
compartmentsy1 = baseliney2 + (length+5) * (math.sin(compradian))
compartmentsx2 = easting + (length+5) * (math.cos(compradian))
compartmentsy2 = northing + (length+5) * (math.sin(compradian))

shotrend_spacing = [[easting2,northing2],]
shotrend_spacing.append([baselinex22,baseliney22],)
shotrend_spacing.append([compartmentsx1,compartmentsy1],)
shotrend_spacing.append([compartmentsx2,compartmentsy2],)

point = arcpy.Point()

pointGeometryList = []

for pt in shotrend_spacing:
	point.X = pt[0]
	point.Y = pt[1]

	pointGeometry = arcpy.PointGeometry(point)
	pointGeometryList.append(pointGeometry)

tpts = wkspc + 'Temp.shp'

arcpy.CopyFeatures_management(pointGeometryList, tpts)
arcpy.DefineProjection_management (tpts, sr)	

arcpy.AddField_management(tpts, 'z', 'DOUBLE')
arcpy.CalculateField_management(tpts, 'z', elev_threshold, "PYTHON_9.3")

TIN = wkspc + "tempTIN"
RAS = locras + "base_ras"
arcpy.CreateTin_3d(TIN, sr, [[tpts, "z", "masspoints"]],"CONSTRAINED_DELAUNAY")
arcpy.TinRaster_3d(TIN, RAS, "FLOAT", "LINEAR", "CELLSIZE 5.0") #1 m cell size 

n1 = -1
cursor2 = arcpy.SearchCursor(tpts)
for row in cursor2:
	
	arcpy.CalculateField_management(tpts, 'z', n, "PYTHON_9.3")
	
for n in xrange (0,4):
	if n == 0:
		w = -3.5
	else:
		w = 3.5
	compid = str('"FID" = ' + str(n))
	cursor3 = arcpy.UpdateCursor(tpts, compid)
	for row in cursor3:
		row.setValue('z', w)
		cursor3.updateRow(row)

TIN = wkspc + "tempTIN"
RAS2 = locdif + "Differences"
arcpy.CreateTin_3d(TIN, sr, [[tpts, "z", "masspoints"]],"CONSTRAINED_DELAUNAY")
arcpy.TinRaster_3d(TIN, RAS2, "FLOAT", "LINEAR", "CELLSIZE 5.0") #1 m cell size 

arcpy.CreateFileGDB_management(locvec, 'Vectors')
wkspcgdb = locvec + 'Vectors.gdb'

vecdlname = "Vectors_D"
vecelname = "Vectors_E"
mrklname = "Markers"
baselname = "Baseline"

arcpy.CreateFeatureclass_management(wkspcgdb, vecdlname, "POLYLINE")
arcpy.CreateFeatureclass_management(wkspcgdb, vecelname, "POLYLINE")
arcpy.CreateFeatureclass_management(wkspcgdb, mrklname, "POLYLINE")
arcpy.CreateFeatureclass_management(wkspcgdb, baselname, "POLYLINE")


arcpy.Delete_management(TIN)
arcpy.Delete_management(tpts)

compparsc = "Compartment_params.pyc"
os.remove(compparsc)

print " "
print raw_input("Compartments created, hit Enter to exit. ;)")