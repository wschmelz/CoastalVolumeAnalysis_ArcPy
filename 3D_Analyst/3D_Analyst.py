import os
import arcpy
from arcpy.sa import *
import datetime			
import glob
import random
import numpy
import math
import shutil

comppar = "Compartment_params.txt"
comppars = "Compartment_params.py"
an3D = "3D_Analyst_params.txt"
ans3D = "Analyst_params.py"

os.rename(comppar, comppars)
os.rename(an3D, ans3D)
from Compartment_params import *
from Analyst_params import *

from numpy import matrix
from numpy import linalg
from numpy import genfromtxt

os.rename(comppars, comppar)
os.rename(ans3D, an3D)

print (time.strftime("%H:%M:%S"))

arcpy.CheckOutExtension("Spatial")
arcpy.CheckOutExtension("3D")
arcpy.env.overwriteOutput = True
workspace = os.getcwd()
arcpy.env.workspace = workspace
arcpy.env.overwriteOutput = True
backslash = '\\'

sr = arcpy.SpatialReference(26918) #UTM18N
wkspc = str(os.getcwd()).replace(backslash,"/") + "/"

locbmk = wkspc + '01_Benchmarks/'
locrwd = wkspc + '02_RawData/'
locshp = wkspc + '03_Shapefiles/'
loctin = wkspc + '04_TINs/'
locmap = wkspc + '05_ArcMap/'
locras = wkspc + '06_Rasters/'
locdif = wkspc + '07_Differences/'
loctab = wkspc + '08_Tables/'
loccomp = wkspc + '09_Compartments/'
locvec = wkspc + '10_Vectors/'
loclyr = wkspc + '11_Layers/'
locndif = wkspc + '12_TotalDifferences/'
locntab = wkspc + '13_TotalTables/'
locnvec = wkspc + '14_TotalVectors/'
locimg = wkspc + '15_Images/'

locvecgdb = locvec + 'Vectors.gdb'
locnvecgdb = locnvec + 'Vectors.gdb'

baselinedistance = distance 

vectordirection = compdirection

vecdlname = "Vectors_D"
vecelname = "Vectors_E"
diflname = "Differences"
imglname = "Background"
mrklname = "Markers"
baselname = "Baseline"

compname = "Compartments.shp"
integ=str(random.randint(1,12000000))
in_dem = locras + "base_ras"
baseline1 = locvec + baselname + ".shp"

mapname = "Template/Template.aprx"

compfile = loccomp + compname

print ("Directory:")
print (str(str(os.getcwd()).replace(backslash,"/") + "/"))
print (" ")

compcount = int(arcpy.GetCount_management(loccomp + compname).getOutput(0))

gpsdata = glob.glob(locshp + "*.shp")

n=1
data = []

for n in range (0,len(gpsdata)):
	data.append(gpsdata[n][-12:-4])
print ("Shapefiles:")
print (data)
print (" ")

data2 = []

for n in range (0,len(gpsdata)):
	data2.append(gpsdata[n][-19:-4])
print ("Shapefiles:")
print (data2)
print (" ")

months = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]

titles = []
for n in range (1, len(data)):
	shoname = data[n-1]
	shoname2 = data[n]
	shomm = shoname[-4:-2]
	shomm2 = shoname2[-4:-2]
	shoyyyy = shoname[-8:-4]
	shoyyyy2 = shoname2[-8:-4]
	title = str(months[int(shomm)-1] + " " + shoyyyy + " to " + months[int(shomm2)-1] + " " + shoyyyy2)
	titles.append(title)

labels = []
for n in range (0, len(data)):
	shoname = data[n]
	shomm = shoname[-4:-2]
	shoyyyy = shoname[-8:-4]
	shoyyyy2 = shoname2[-8:-4]
	label = str(months[int(shomm)-1] + " " + shoyyyy)

	labels.append(label)

print ("Surveys:")
print (list(labels))
print (" ")
print ("Analyses:")
print (list(titles))
print (" ")
		
#set environment
arcpy.MakeRasterLayer_management(in_dem, "extlyr")
arcpy.env.extent = arcpy.Describe("extlyr").extent
arcpy.env.overwriteOutput = True
#+++#
#Create TINs

n = 1

for n in range (1, len(data) + 1):
	shp = locshp + data2[n-1] + ".shp"
	TIN = loctin + data[n-1]
	RAS = locras + data[n-1]
	arcpy.CreateTin_3d(TIN, sr, [[shp, "Z", "masspoints"]], "DELAUNAY")
	arcpy.DelineateTinDataArea_3d (TIN, 4000)
	print (data[n-1] + ' TIN created')

#Tins to rasters

n = 1

for n in range (1, len(data) + 1):
	TIN = loctin + data[n-1]
	RAS = locras + data[n-1]# + ".tif"
	arcpy.TinRaster_3d(TIN, RAS, "FLOAT", "LINEAR", "CELLSIZE 1.0") #1 m cell size 
	print (data[n-1] + ' raster created')
 
#Merge rasters


n=1

for n in range (1, len(data) + 1):
	RAS = locras + data[n-1]# + ".tif"
	mRAS = locras + 'm_' + data[n-1]
	temp = "in_memory/t1" + integ
	tempn = "t1" + integ
	mRASname = 'm_' + data[n-1]
	arcpy.MosaicToNewRaster_management([RAS,'extlyr'], "in_memory", tempn, sr, "32_BIT_FLOAT", "1.0", "1", "MAXIMUM","FIRST") #1 m cell size
	arcpy.gp.ExtractByMask_sa(temp,compfile,mRAS)
	arcpy.Delete_management("in_memory")
	
	print (data[n-1] + ' data raster merged with base')
		  
 

#Subtract rasters

shutil.rmtree(locdif)
os.makedirs(locdif)

n = 2

for n in range (2, len(data) + 1):

	mRAS = locras + 'm_' + data[n-1]
	mRASname = 'm_' + data[n-1]
	mRAS2 = locras + 'm_' + data[n-2]
	mRAS2name = 'm_' + data[n-2]
	RASname = data[n-1]
	RAS2name = data[n-2]
	RASdif = locdif + RASname[:6] + '_' + RAS2name[:6]
	RASdifname = RASname[:6] + '_' + RAS2name[:6]
	outMinus = arcpy.sa.Minus(mRAS, mRAS2)
	outMinus.save(RASdif)
	
	print (RASdifname + ' difference raster created')
	n = n + 1

if arcpy.Exists(locvec + 'Vectors.gdb'):
	arcpy.Delete_management(locvec + 'Vectors.gdb')

arcpy.CreateFileGDB_management(locvec, 'Vectors')
	
arcpy.env.workspace = locvec + 'Vectors.gdb'
arcpy.env.overwriteOutput = True
wkspcgdb = locvec + 'Vectors.gdb'
baseline = 'Baseline'
arcpy.CopyFeatures_management(baseline1, baseline)
arcpy.DefineProjection_management (baseline, sr)

baselineradian = math.radians(baselineangle)
baselinelength = compartments*distance
baselinex2 = easting + baselinelength*math.cos(baselineradian)
baseliney2 = northing + baselinelength*math.sin(baselineradian)
	
mark_dists = numpy.array([(marker_labels)])
invarray = -1*mark_dists
markdists = numpy.append(mark_dists,invarray)

print ("")

print ("Creating marker lines for maps.")

n=0
markers = 'Markers'
arcpy.CreateFeatureclass_management(wkspcgdb, markers, "POLYLINE", "", "","", sr)
arcpy.AddField_management(markers, "Label", "Double")
markertemp = "tmp"
for n in range (0, len(markdists)):
	distval =  int(markdists[n] * -1)
	distmarker = vectormultiplier * distval * -1
	
	polylineX = []

	z=1
	for z in range(1, 2):
		tx1 = easting
		ty1 = northing
		tx2 = baselinex2
		ty2 = baseliney2	
		x2_x1_1 = tx2 - tx1
		y2_y1_1 = ty2 - ty1
		
		id = z
		x2_x1_1 = tx2 - tx1
		y2_y1_1 = ty2 - ty1
		if z==1:
			seg1ang = math.atan2(y2_y1_1,x2_x1_1)
			baselineangle = (seg1ang)
			if vectordirection == "left":
				transectangle = baselineangle - math.pi/2
			elif vectordirection == "right":
				transectangle = baselineangle + math.pi/2				
			q = tx1
			w = ty1
			e = q + distmarker * (math.cos(transectangle))
			r = w + distmarker * (math.sin(transectangle))
			polylineX.append([[e,r],])
			sys.stdout.write("\rProcessed %i of %i vertices on %i of %i marker lines" % (1,2,n+1, len(markdists)))
		if 2 == 2:
			seg1ang = math.atan2(y2_y1_1,x2_x1_1)
			baselineangle = (seg1ang)
			if vectordirection == "left":
				transectangle = baselineangle - math.pi/2
			elif vectordirection == "right":
				transectangle = baselineangle + math.pi/2				
			q = tx2
			w = ty2
			e = q + distmarker * (math.cos(transectangle))
			r = w + distmarker * (math.sin(transectangle))
			polylineX.append([[e,r],])
			sys.stdout.write("\rProcessed %i of %i vertices on %i of %i marker lines" % (2,2,n+1, len(markdists)))			
		else:

			seg1ang = math.atan2(y2_y1_1,x2_x1_1)
			baselineangle = (seg1ang)
			
			if vectordirection == "left":
				transectangle = baselineangle - math.pi/2
			elif vectordirection == "right":
				transectangle = baselineangle + math.pi/2
				
			q = tx2
			w = ty2
			e = q + distmarker * (math.cos(transectangle))
			r = w + distmarker * (math.sin(transectangle))

			polylineX.append([[e,r],])
			sys.stdout.write("\rProcessed %i of %i vertices on %i of %i marker lines" % (2,2,n+1, len(markdists)))

	point = arcpy.Point()  
	array = arcpy.Array()  
	features = []  
	for feature in polylineX:
		for coordPair in feature:  
			point.X = coordPair[0]  
			point.Y = coordPair[1]  
			array.add(point)
	
	features.append(
		arcpy.Polyline(array))

	markertmp = 'marktmp'
	arcpy.CopyFeatures_management(features, markertmp)
	arcpy.DefineProjection_management (markertmp, sr)
	arcpy.AddField_management(markertmp, "Label", "Double")
	
	arcpy.CalculateField_management(markertmp, "Label", distval, "PYTHON3")			

	arcpy.Append_management (markertmp, markers, "NO_TEST")		



print (".")
	
#Create data table

difs = [x[0][-13:] for x in os.walk(locdif)]

voltab = loctab + 'Volumes'

arcpy.CreateTable_management(loctab, 'Volumes')
arcpy.AddField_management(voltab, "COMP", "SHORT")
arcpy.AddField_management(voltab, "DURATION", "TEXT")
arcpy.AddField_management(voltab, "DEPOSITION", "Double")
arcpy.AddField_management(voltab, "EROSION", "Double")
arcpy.AddField_management(voltab, "VOLUME", "Double")
arcpy.DeleteField_management (voltab, "FIELD1")
arcpy.DeleteField_management (voltab, "OBJECTID")
n = 2
for n in range (2, len(data) + 1):
	dif = locdif + difs[n-1]
	difname = str(difs[n-1])
	y = 1
	for y in range (1, compcount + 1):
		rowid = str(((n-2)*compcount)+y)
		comp = loccomp + 'Compartments.shp'
		compnum = str(y)
		compid  = str('"COMP" = ' + str(y))
		temp1 = "in_memory/t2" + integ
		arcpy.MakeFeatureLayer_management(comp, "lyr")
		arcpy.SelectLayerByAttribute_management("lyr", "NEW_SELECTION", compid)
		arcpy.gp.ExtractByMask_sa(dif,"lyr",temp1)
		result = arcpy.SurfaceVolume_3d(temp1,"#","ABOVE","0","1","0")
		depvolp = str(arcpy.GetMessages())
		volpos = depvolp.find('Volume=')
		depvol1 = str(depvolp[volpos+7:volpos+17])
		try:
			depvol = float(depvol1)
		except:
			depvol = 0
		result = arcpy.SurfaceVolume_3d(temp1,"#","BELOW","0","1","0")
		arcpy.GetMessages()
		erovolp = str(arcpy.GetMessages())
		volpos = erovolp.find('Volume=')
		erovol1 = str(erovolp[volpos+7:volpos+17])
		try:
			erovol = float(erovol1)
		except:
			erovol = 0
		totvol = depvol - erovol
		fields = ["Rowid", "COMP","DURATION", "DEPOSITION", "EROSION", "VOLUME"]
		cursor = arcpy.da.InsertCursor(voltab, fields)
		cursor.insertRow([rowid, compnum, difname, depvol, erovol, totvol])
		sys.stdout.write("\rProcessed %i of %i for duration %s" % (y, compcount, difname))
       
		arcpy.Delete_management("in_memory")
			
			
	print ('.')	

#Pivot data table
voltab = loctab + 'Volumes'
voltabp = loctab + 'Volumes_Pivot'
voltabpname = 'Volumes_Pivot'
difs = [x[0][-13:] for x in os.walk(locdif)]

arcpy.CreateTable_management(loctab, 'Volumes_Pivot')
arcpy.AddField_management(voltabp, "COMP", "SHORT")
arcpy.DeleteField_management (voltabp, "FIELD1")
arcpy.DeleteField_management (voltabp, "OBJECTID")

for n in range (2, len(data) + 1):
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]
	arcpy.AddField_management(voltabp, dif_F_name, "Double")
	

	

y = 1
for y in range (1, compcount + 1):
	compnum = str(y)
	cursor = arcpy.InsertCursor(voltabp)
	row = cursor.newRow()
	row.setValue("COMP", compnum)
	cursor.insertRow(row)
	del(cursor)
	for n in range (2, len(data) + 1):
		dif = locdif + difs[n-1]
		difname = str(difs[n-1])
		dif_F_name = 'D' + difname[-12:]
		compid = str('"COMP" = ' + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		cursor2 = arcpy.SearchCursor(voltab,difcompid)
		for row in cursor2:
			comptimvol = row.getValue("VOLUME")
		cursor3 = arcpy.UpdateCursor(voltabp, compid)
		for row in cursor3:
			row.setValue(dif_F_name, comptimvol)
			cursor3.updateRow(row)
		del(cursor2, cursor3)
		n=n+1
	sys.stdout.write("\rProcessed %i of %i" % (y, compcount))
	y=y+1	
arcpy.TableToTable_conversion (voltabp, loctab, voltabpname)
arcpy.TableToExcel_conversion (voltabp, voltabp + '.xls') 	
print ('.')
		
#Pivot deposition data table
voltab = loctab + 'Volumes'
deptabp = loctab + 'Deposition_Pivot'
deptabpname = 'Deposition_Pivot'
difs = [x[0][-13:] for x in os.walk(locdif)]

arcpy.CreateTable_management(loctab, 'Deposition_Pivot')
arcpy.AddField_management(deptabp, "COMP", "SHORT")
arcpy.DeleteField_management (deptabp, "FIELD1")
arcpy.DeleteField_management (deptabp, "OBJECTID")

for n in range (2, len(data) + 1):
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]
	arcpy.AddField_management(deptabp, dif_F_name, "Double")

	

y = 1
for y in range (1, compcount + 1):
	compnum = str(y)
	cursor = arcpy.InsertCursor(deptabp)
	row = cursor.newRow()
	row.setValue("COMP", compnum)
	cursor.insertRow(row)
	del(cursor)
	for n in range (2, len(data) + 1):
		dif = locdif + difs[n-1]
		difname = str(difs[n-1])
		dif_F_name = 'D' + difname[-12:]
		compid = str('"COMP" = ' + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		cursor2 = arcpy.SearchCursor(voltab,difcompid)
		for row in cursor2:
			comptimdep = row.getValue("DEPOSITION")
		cursor3 = arcpy.UpdateCursor(deptabp, compid)
		for row in cursor3:
			row.setValue(dif_F_name, comptimdep)
			cursor3.updateRow(row)
		del(cursor2, cursor3)
		
	sys.stdout.write("\rProcessed %i of %i" % (y, compcount))
	
arcpy.TableToTable_conversion (deptabp, loctab, deptabpname)
arcpy.TableToExcel_conversion (deptabp, deptabp + '.xls') 	
print ('.')

#Pivot erosion data table
voltab = loctab + 'Volumes'
erotabp = loctab + 'Erosion_Pivot'
erotabpname = 'Erosion_Pivot'
difs = [x[0][-13:] for x in os.walk(locdif)]

arcpy.CreateTable_management(loctab, 'Erosion_Pivot')
arcpy.AddField_management(erotabp, "COMP", "SHORT")
arcpy.DeleteField_management (erotabp, "FIELD1")
arcpy.DeleteField_management (erotabp, "OBJECTID")

for n in range (2, len(data) + 1):
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]
	arcpy.AddField_management(erotabp, dif_F_name, "Double")
	
	
	

y = 1
for y in range (1, compcount + 1):
	compnum = str(y)
	cursor = arcpy.InsertCursor(erotabp)
	row = cursor.newRow()
	row.setValue("COMP", compnum)
	cursor.insertRow(row)
	del(cursor)
	for n in range (2, len(data) + 1):
		dif = locdif + difs[n-1]
		difname = str(difs[n-1])
		dif_F_name = 'D' + difname[-12:]
		compid = str('"COMP" = ' + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		cursor2 = arcpy.SearchCursor(voltab,difcompid)
		for row in cursor2:
			comptimero = row.getValue("EROSION")
		cursor3 = arcpy.UpdateCursor(erotabp, compid)
		for row in cursor3:
			row.setValue(dif_F_name, comptimero)
			cursor3.updateRow(row)
		del(cursor2, cursor3)
		
	sys.stdout.write("\rProcessed %i of %i" % (y, compcount))
	
arcpy.TableToTable_conversion (erotabp, loctab, erotabpname)
arcpy.TableToExcel_conversion (erotabp, erotabp + '.xls') 	
print ('.')

print (time.strftime("%H:%M:%S"))

#Create Vectors
arcpy.env.workspace = locvec + 'Vectors.gdb'
arcpy.env.overwriteOutput = True
arcpy.env.extent = "MAXOF"

arcpy.AddField_management(baseline, 'StartX', 'DOUBLE')
arcpy.AddField_management(baseline, 'StartY', 'DOUBLE')
arcpy.AddField_management(baseline, 'EndX', 'DOUBLE')
arcpy.AddField_management(baseline, 'EndY', 'DOUBLE')

arcpy.CalculateField_management(baseline, 'StartX', "!SHAPE.firstPoint.X!", "PYTHON3")
arcpy.CalculateField_management(baseline, 'StartY', "!SHAPE.firstPoint.Y!", "PYTHON3")
arcpy.CalculateField_management(baseline, 'EndX', "!SHAPE.lastPoint.X!", "PYTHON3")
arcpy.CalculateField_management(baseline, 'EndY', "!SHAPE.lastPoint.Y!", "PYTHON3")

cursor = arcpy.SearchCursor(baseline)

for row in cursor:
	startx = row.getValue("StartX")
	starty = row.getValue("StartY")
	endx = row.getValue("EndX")
	endy = row.getValue("EndY")

deltax = endx - startx
deltay = endy - starty

baselineradian = math.atan2(deltay,deltax)

print ("Shoreline is located in which direction relative to the baseline: left or right?")

if vectordirection == "left":
	vectorangle1 = baselineradian + math.pi/2
	vectorangle2 = baselineradian - math.pi/2
elif vectordirection == "right":
	vectorangle1 = baselineradian - math.pi/2
	vectorangle2 = baselineradian + math.pi/2
else:
	print ("Invalid answer, direction must be left or right.")

#Draw Vectors

print ("Draw vectors.")

voltab = loctab + 'Volumes'
deptabp = loctab + 'Deposition_Pivot'
difs = [x[0][-13:] for x in os.walk(locdif)]
for n in range (2, len(data) + 1):
	
	dif = locdif + difs[n-1]
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]

	vectors = "Vectors_Deposition_" + difname
	dunedifout = voltab


	y = 1

	compid = str('"COMP" = ' + str(y))
	difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")

	q = startx + (baselinedistance * 0.5) * (math.cos(baselineradian))

	w = starty + (baselinedistance * 0.5) * (math.sin(baselineradian))

	cursor = arcpy.SearchCursor(voltab, difcompid)

	for row in cursor:
		volumedif = float(row.getValue("Deposition"))

	if volumedif >= 0:
		vectordistance = vectormultiplier * volumedif * - 1
		vectorradian = vectorangle2
	else:
		vectordistance = vectormultiplier * volumedif
		vectorradian = vectorangle1
		
	e = q + vectordistance * (math.cos(vectorradian))

	r = w + vectordistance * (math.sin(vectorradian))

	feature_info = [[[q, w], [e, r]],]
	y = 2
	for y in range (2, compcount + 1):
		compid = str('"COMP" = ' + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		cursor = arcpy.SearchCursor(voltab, difcompid)
		for row in cursor:
			volumedif = float(row.getValue("Deposition"))
		if volumedif >= 0:
			vectordistance = vectormultiplier * volumedif * - 1
			vectorradian = vectorangle2
		else:
			vectordistance = vectormultiplier * volumedif
			vectorradian = vectorangle1
		iterable = [[[q + baselinedistance * (math.cos(baselineradian)), w + baselinedistance * (math.sin(baselineradian))], [(q + baselinedistance * (math.cos(baselineradian))) + vectordistance * (math.cos(vectorradian)), (w + baselinedistance * (math.sin(baselineradian))) + vectordistance * (math.sin(vectorradian))]],]
		feature_info.extend(iterable)

		q = q + baselinedistance * (math.cos(baselineradian))
		
		w = w + baselinedistance * (math.sin(baselineradian))

		sys.stdout.write("\rProcessed %i of %i vectors for duration %s" % (y, compcount, difname))
		
	print (".")
	
	features = []

	for feature in feature_info:
		features.append(
			arcpy.Polyline(
				arcpy.Array([arcpy.Point(*coords) for coords in feature]),sr))
			
	arcpy.CopyFeatures_management(features, vectors)
	arcpy.AddField_management(vectors, "Comp", "Integer")
	arcpy.AddField_management(vectors, "Deposition", "Double")

	n = 1

	for y in range (1, compcount + 1):
		zone = str(y)
		compid = str('"COMP" = ' + str(y))
		objid = str("OBJECTID = " + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		compid = str('"COMP" = ' + str(y))
		cursor = arcpy.SearchCursor(voltab, difcompid)
		for row in cursor:
			volumedif = float(row.getValue("Deposition"))
		arcpy.MakeFeatureLayer_management(vectors, "lyr1")
		arcpy.SelectLayerByAttribute_management("lyr1", "NEW_SELECTION", objid)
		arcpy.CalculateField_management("lyr1", "Comp", zone, "PYTHON3")
		arcpy.CalculateField_management("lyr1", "Deposition", volumedif, "PYTHON3")
		sys.stdout.write("\rCaluculated %i of %i vector attributes for duration %s" % (y, compcount, difname))
		
	print (".")

print ("Vectors drawn.")

#Draw Vectors 2

print ("Draw vectors.")

voltab = loctab + 'Volumes'
erotabp = loctab + 'Erosion_Pivot'
difs = [x[0][-13:] for x in os.walk(locdif)]
for n in range (2, len(data) + 1):
	
	dif = locdif + difs[n-1]
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]

	vectors = "Vectors_Erosion_" + difname
	dunedifout = voltab

	y = 1

	compid = str('"COMP" = ' + str(y))
	difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")

	q = startx + (baselinedistance * 0.5) * (math.cos(baselineradian))

	w = starty + (baselinedistance * 0.5) * (math.sin(baselineradian))

	cursor = arcpy.SearchCursor(voltab, difcompid)

	for row in cursor:
		volumedif = float(row.getValue("Erosion")*-1)

	if volumedif >= 0:
		vectordistance = vectormultiplier * volumedif * - 1
		vectorradian = vectorangle2
	else:
		vectordistance = vectormultiplier * volumedif
		vectorradian = vectorangle1
		
	e = q + vectordistance * (math.cos(vectorradian))

	r = w + vectordistance * (math.sin(vectorradian))

	feature_info = [[[q, w], [e, r]],]
	y = 2

	for y in range (2, compcount + 1):
		compid = str('"COMP" = ' + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		cursor = arcpy.SearchCursor(voltab, difcompid)
		for row in cursor:
			volumedif = float(row.getValue("Erosion")*-1)
		if volumedif >= 0:
			vectordistance = vectormultiplier * volumedif * - 1
			vectorradian = vectorangle2
		else:
			vectordistance = vectormultiplier * volumedif
			vectorradian = vectorangle1
		iterable = [[[q + baselinedistance * (math.cos(baselineradian)), w + baselinedistance * (math.sin(baselineradian))], [(q + baselinedistance * (math.cos(baselineradian))) + vectordistance * (math.cos(vectorradian)), (w + baselinedistance * (math.sin(baselineradian))) + vectordistance * (math.sin(vectorradian))]],]
		feature_info.extend(iterable)

		q = q + baselinedistance * (math.cos(baselineradian))
		
		w = w + baselinedistance * (math.sin(baselineradian))

		sys.stdout.write("\rProcessed %i of %i vectors for duration %s" % (y, compcount, difname))
		y = y + 1 
	print (".")
	
	features = []

	for feature in feature_info:
		features.append(
			arcpy.Polyline(
				arcpy.Array([arcpy.Point(*coords) for coords in feature]),sr))
			
	arcpy.CopyFeatures_management(features, vectors)
	arcpy.AddField_management(vectors, "Comp", "Integer")
	arcpy.AddField_management(vectors, "Erosion", "Double")

	n = 1

	for y in range (1, compcount + 1):
		zone = str(y)
		compid = str('"COMP" = ' + str(y))
		objid = str("OBJECTID = " + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		compid = str('"COMP" = ' + str(y))
		cursor = arcpy.SearchCursor(voltab, difcompid)
		for row in cursor:
			volumedif = float(row.getValue("Erosion")*-1)
		arcpy.MakeFeatureLayer_management(vectors, "lyr1")
		arcpy.SelectLayerByAttribute_management("lyr1", "NEW_SELECTION", objid)
		arcpy.CalculateField_management("lyr1", "Comp", zone, "PYTHON3")
		arcpy.CalculateField_management("lyr1", "Erosion", volumedif, "PYTHON3")
		sys.stdout.write("\rCaluculated %i of %i vector attributes for duration %s" % (y, compcount, difname))
		
	print (".")

print ("Vectors drawn.")
#+++#
locvecgdb = locvec + 'Vectors.gdb'

#Create Maps

#file locations for datalayers

baslyr_path = loclyr + "Baseline.lyrx"
veclyre_path = loclyr + "Vectors_E.lyrx"
veclyrd_path = loclyr + "Vectors_D.lyrx"
diflyr_path = loclyr + "Differences.lyrx"
mrklyr_path = loclyr + "Markers.lyrx"

baslyr_lf = arcpy.mp.LayerFile(baslyr_path)
veclyre_lf = arcpy.mp.LayerFile(veclyre_path)
veclyrd_lf = arcpy.mp.LayerFile(veclyrd_path)
diflyr_lf = arcpy.mp.LayerFile(diflyr_path)
mrklyr_lf = arcpy.mp.LayerFile(mrklyr_path)

difs = [x[0][-13:] for x in os.walk(locdif)]

n=2
voltabp = loctab + 'Volumes_Pivot'

for n in range (2, len(data) + 1):
	dif = locdif + difs[n-1]
	difname = str(difs[n-1])
	comp = loccomp + 'Compartments.shp'
	compid  = str('"COMP" < ' + str(1000))
	temp1 = "in_memory/t2" + integ
	arcpy.MakeFeatureLayer_management(comp, "lyr")
	arcpy.SelectLayerByAttribute_management("lyr", "NEW_SELECTION", compid)
	outRaster = ExtractByMask(dif,"lyr")
	outRaster.save(temp1)
	arcpy.CopyRaster_management(temp1, dif) 
	sys.stdout.write("\rClipped difference raster %s" % (difname))
	arcpy.Delete_management("in_memory")
	
	
print (".")

site = gpsdata[0][-19:-12]	

for n in range (2, len(data) + 1):
	dif = locdif + difs[n-1]
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]

	searchid = "MEAN_" + str(dif_F_name)

	compid1 = str('"COMP" <= ' + str(compbreak1))
	temp_table = "in_memory/tempt" + integ
	arcpy.TableSelect_analysis (voltabp, temp_table, compid1)
	stats_table = "in_memory/stats" + integ
	arcpy.Statistics_analysis(temp_table,stats_table,[[dif_F_name,"MEAN"]])
	cursor = arcpy.SearchCursor(stats_table)
	for row in cursor:
		vol1 = str(round(row.getValue(searchid),2))

	compid2 = str('"COMP" > ' + str(compbreak1) + 'AND "COMP" <= ' + str(compbreak2))
	temp_table = "in_memory/tempt" + integ
	arcpy.TableSelect_analysis (voltabp, temp_table, compid2)
	stats_table = "in_memory/stats" + integ
	arcpy.Statistics_analysis(temp_table,stats_table,[[dif_F_name,"MEAN"]])
	cursor = arcpy.SearchCursor(stats_table)
	for row in cursor:
		vol2 = str(round(row.getValue(searchid),2))

	compid3 = str('"COMP" > ' + str(compbreak2))
	temp_table = "in_memory/tempt" + integ
	arcpy.TableSelect_analysis (voltabp, temp_table, compid3)
	stats_table = "in_memory/stats" + integ
	arcpy.Statistics_analysis(temp_table,stats_table,[[dif_F_name,"MEAN"]])
	cursor = arcpy.SearchCursor(stats_table)
	for row in cursor:
		vol3 = str(round(row.getValue(searchid),2))

	stats_table = "in_memory/stats" + integ
	arcpy.Statistics_analysis(voltabp,stats_table,[[dif_F_name,"MEAN"]])
	cursor = arcpy.SearchCursor(stats_table)
	for row in cursor:
		vol4 = str(round(row.getValue(searchid),2))
	totalm3_1 = float(float(vol4) * float(compcount))
	totalm3 = str(round(totalm3_1,2))

	Comp_1 = "Mean Volumetric Change of" + "\n" + "Compartments 1 - " + str(compbreak1) + ": " + str(vol1 + " m<SUP>3</SUP>")
	Comp_2 = "Mean Volumetric Change of" + "\n" + "Compartments " + str(compbreak1+1) + " - " + str(compbreak2) + ": " + str(vol2 + " m<SUP>3</SUP>")
	Comp_3 = "Mean Volumetric Change of" + "\n" + "Compartments " + str(compbreak2+1) + " - " + str(compcount) + ": " + str(vol3 + " m<SUP>3</SUP>")
	Comp_All = "Mean Volumetric Change" + "\n" + "per Compartment: " + str(vol4 + " m<SUP>3</SUP>")
	Total_Change = "Total Volumetric Change: " + str(totalm3 + " m<SUP>3</SUP>")

	title = site[0:3] + " Erosion and Deposition" + "\n" + titles[n-2]


	dif = locdif + difs[n-1]
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]

	vectorse = "Vectors_Erosion_" + difname
	vectorsd = "Vectors_Deposition_" + difname

	mapname = "Layers"  # Make sure a map named "Template" exists in the APRX
	aprx_path = locmap + "Template/Template.aprx"
	aprx = arcpy.mp.ArcGISProject(aprx_path)
	maps = aprx.listMaps("Layers")
	if not maps:
		raise RuntimeError("No map named 'Layers' found in the project.")

	map_obj = maps[0]

	veclyre_name = veclyre_lf.listLayers()[0].name
	veclyrd_name = veclyrd_lf.listLayers()[0].name
	diflyr_name = diflyr_lf.listLayers()[0].name
	mrklyr_name = mrklyr_lf.listLayers()[0].name
	baslyr_name = baslyr_lf.listLayers()[0].name
	
	layers_to_remove = [vecdlname, vecelname, diflname, mrklname, baselname]
	for layer in map_obj.listLayers():
		if layer.name in layers_to_remove:
			map_obj.removeLayer(layer)

	map_obj.addDataFromPath(diflyr_path)
	map_obj.addDataFromPath(veclyre_path)
	map_obj.addDataFromPath(veclyrd_path)
	map_obj.addDataFromPath(mrklyr_path)
	map_obj.addDataFromPath(baslyr_path)
	
	for lyr in map_obj.listLayers():
		if lyr.name == veclyre_name:
			
			print ("Layer:",veclyre_name)
			old_conn = lyr.connectionProperties
			print("Old connection:", lyr.connectionProperties)
			new_conn = old_conn.copy()			
			new_conn['connection_info']['database'] = locvecgdb
			new_conn['dataset'] = vectorse
			lyr.updateConnectionProperties(old_conn, new_conn)
			print("New connection:", lyr.connectionProperties)
						
		elif lyr.name == veclyrd_name:
			
			print ("Layer:",veclyrd_name)
			old_conn = lyr.connectionProperties
			print("Old connection:", lyr.connectionProperties)
			new_conn = old_conn.copy()			
			new_conn['connection_info']['database'] = locvecgdb
			new_conn['dataset'] = vectorsd			
			lyr.updateConnectionProperties(old_conn, new_conn)
			print("New connection:", lyr.connectionProperties)	
			
		elif lyr.name == diflyr_name:
		
			print ("Layer:",veclyrd_name)
			old_conn = lyr.connectionProperties
			print("Old connection:", lyr.connectionProperties)
			new_conn = old_conn.copy()			
			new_conn['connection_info']['database'] = locdif
			new_conn['dataset'] = difname			
			lyr.updateConnectionProperties(old_conn, new_conn)
			print("New connection:", lyr.connectionProperties)			

	layouts = aprx.listLayouts()
	if layouts:
		layout = layouts[0] 
		for element in layout.listElements("TEXT_ELEMENT"):
			if element.name == "Title":
				element.text = title
			elif element.name == "Comp_1":
				element.text = Comp_1
			elif element.name == "Comp_2":
				element.text = Comp_2
			elif element.name == "Comp_3":
				element.text = Comp_3
			elif element.name == "Comp_All":
				element.text = Comp_All
			elif element.name == "Total_Change":
				element.text = Total_Change

		print("Updated layout elements.")

	aprx.saveACopy(locmap + difname + ".aprx")
	print (titles[n-2] + " map created")
	del aprx

	
################################################################################################################

###############################################
###############################################



print ("Directory:")
print (str(str(os.getcwd()).replace(backslash,"/") + "/"))
print (" ")

compcount = int(arcpy.GetCount_management(loccomp + compname).getOutput(0))

gpsdata = glob.glob(locshp + "*.shp")

n=1

data3 = [startdate, enddate]

print ("Shapefiles:")
print (data3)
print (" ")

site = gpsdata[0][-19:-12]

startdatefile = site + startdate
enddatefile = site + enddate
data4 = [startdatefile, enddatefile]

print ("Shapefiles:")
print (data4)
print (" ")

months = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]

titles2 = []
for n in range (1, len(data3)):
	shoname = data3[n-1]
	shoname2 = data3[n]
	shomm = shoname[-4:-2]
	shomm2 = shoname2[-4:-2]
	shoyyyy = shoname[-8:-4]
	shoyyyy2 = shoname2[-8:-4]
	title = str(months[int(shomm)-1] + " " + shoyyyy + " to " + months[int(shomm2)-1] + " " + shoyyyy2)
	titles2.append(title)

labels = []
for n in range (0, len(data3)):
	shoname = data3[n]
	shomm = shoname[-4:-2]
	shoyyyy = shoname[-8:-4]
	shoyyyy2 = shoname2[-8:-4]
	label = str(months[int(shomm)-1] + " " + shoyyyy)

	labels.append(label)

print ("Surveys:")
print (list(labels))
print (" ")
print ("Analyses:")
print (list(titles2))
print (" ")
	
	
	
#set environment
arcpy.MakeRasterLayer_management(in_dem, "extlyr")
arcpy.env.extent = arcpy.Describe("extlyr").extent
arcpy.env.overwriteOutput = True


#Subtract rasters

shutil.rmtree(locndif)
os.makedirs(locndif)

n = 2

for n in range (2, len(data3) + 1):

	mRAS = locras + 'm_' + data3[n-1]
	mRASname = 'm_' + data3[n-1]
	mRAS2 = locras + 'm_' + data3[n-2]
	mRAS2name = 'm_' + data3[n-2]
	RASname = data3[n-1]
	RAS2name = data3[n-2]
	RASdif = locndif + RASname[:6] + '_' + RAS2name[:6]
	RASdifname = RASname[:6] + '_' + RAS2name[:6]
	outMinus = Minus(mRAS,mRAS2)
	outMinus.save(RASdif)
	
	
	print (RASdifname + ' difference raster created')
	
#Create data3 table

difs = [x[0][-13:] for x in os.walk(locndif)]

voltab = locntab + 'Volumes'

arcpy.CreateTable_management(locntab, 'Volumes')
arcpy.AddField_management(voltab, "COMP", "SHORT")
arcpy.AddField_management(voltab, "DURATION", "TEXT")
arcpy.AddField_management(voltab, "DEPOSITION", "Double")
arcpy.AddField_management(voltab, "EROSION", "Double")
arcpy.AddField_management(voltab, "VOLUME", "Double")
arcpy.DeleteField_management (voltab, "FIELD1")
arcpy.DeleteField_management (voltab, "OBJECTID")
n = 2
for n in range (2, len(data3) + 1):
	dif = locndif + difs[n-1]
	difname = str(difs[n-1])
	y = 1
	for y in range (1, compcount + 1):
		rowid = str(((n-2)*compcount)+y)
		comp = loccomp + 'Compartments.shp'
		compnum = str(y)
		compid  = str('"COMP" = ' + str(y))
		temp1 = "in_memory/t2" + integ
		arcpy.MakeFeatureLayer_management(comp, "lyr")
		arcpy.SelectLayerByAttribute_management("lyr", "NEW_SELECTION", compid)
		arcpy.gp.ExtractByMask_sa(dif,"lyr",temp1)
		result = arcpy.SurfaceVolume_3d(temp1,"#","ABOVE","0","1","0")
		depvolp = str(arcpy.GetMessages())
		volpos = depvolp.find('Volume=')
		depvol1 = str(depvolp[volpos+7:volpos+17])
		try:
			depvol = float(depvol1)
		except:
			depvol = 0
		result = arcpy.SurfaceVolume_3d(temp1,"#","BELOW","0","1","0")
		arcpy.GetMessages()
		erovolp = str(arcpy.GetMessages())
		volpos = erovolp.find('Volume=')
		erovol1 = str(erovolp[volpos+7:volpos+17])
		try:
			erovol = float(erovol1)
		except:
			erovol = 0
		totvol = depvol - erovol
		fields = ["Rowid", "COMP","DURATION", "DEPOSITION", "EROSION", "VOLUME"]
		cursor = arcpy.da.InsertCursor(voltab, fields)
		cursor.insertRow([rowid, compnum, difname, depvol, erovol, totvol])
		del(cursor)
		sys.stdout.write("\rProcessed %i of %i for duration %s" % (y, compcount, difname))
		y = y + 1
		arcpy.Delete_management("in_memory")
	
	print (".")


#Pivot data3 table
voltab = locntab + 'Volumes'
voltabp = locntab + 'Volumes_Pivot'
voltabpname = 'Volumes_Pivot'
difs = [x[0][-13:] for x in os.walk(locndif)]

arcpy.CreateTable_management(locntab, 'Volumes_Pivot')
arcpy.AddField_management(voltabp, "COMP", "SHORT")
arcpy.DeleteField_management (voltabp, "FIELD1")
arcpy.DeleteField_management (voltabp, "OBJECTID")

for n in range (2, len(data3) + 1):
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]
	arcpy.AddField_management(voltabp, dif_F_name, "Double")
	
	n = n + 1
	

y = 1
for y in range (1, compcount + 1):
	compnum = str(y)
	cursor = arcpy.InsertCursor(voltabp)
	row = cursor.newRow()
	row.setValue("COMP", compnum)
	cursor.insertRow(row)
	del(cursor)
	for n in range (2, len(data3) + 1):
		dif = locndif + difs[n-1]
		difname = str(difs[n-1])
		dif_F_name = 'D' + difname[-12:]
		compid = str('"COMP" = ' + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		cursor2 = arcpy.SearchCursor(voltab,difcompid)
		for row in cursor2:
			comptimvol = row.getValue("VOLUME")
		cursor3 = arcpy.UpdateCursor(voltabp, compid)
		for row in cursor3:
			row.setValue(dif_F_name, comptimvol)
			cursor3.updateRow(row)
		del(cursor2, cursor3)
		n=n+1
	sys.stdout.write("\rProcessed %i of %i" % (y, compcount))
	y=y+1	
arcpy.TableToTable_conversion (voltabp, locntab, voltabpname)
arcpy.TableToExcel_conversion (voltabp, voltabp + '.xls') 	
print ('.')
		
#Pivot deposition data3 table
voltab = locntab + 'Volumes'
deptabp = locntab + 'Deposition_Pivot'
deptabpname = 'Deposition_Pivot'
difs = [x[0][-13:] for x in os.walk(locndif)]

arcpy.CreateTable_management(locntab, 'Deposition_Pivot')
arcpy.AddField_management(deptabp, "COMP", "SHORT")
arcpy.DeleteField_management (deptabp, "FIELD1")
arcpy.DeleteField_management (deptabp, "OBJECTID")

for n in range (2, len(data3) + 1):
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]
	arcpy.AddField_management(deptabp, dif_F_name, "Double")
	
	n = n + 1
	

y = 1
for y in range (1, compcount + 1):
	compnum = str(y)
	cursor = arcpy.InsertCursor(deptabp)
	row = cursor.newRow()
	row.setValue("COMP", compnum)
	cursor.insertRow(row)
	del(cursor)
	for n in range (2, len(data3) + 1):
		dif = locndif + difs[n-1]
		difname = str(difs[n-1])
		dif_F_name = 'D' + difname[-12:]
		compid = str('"COMP" = ' + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		cursor2 = arcpy.SearchCursor(voltab,difcompid)
		for row in cursor2:
			comptimdep = row.getValue("DEPOSITION")
		cursor3 = arcpy.UpdateCursor(deptabp, compid)
		for row in cursor3:
			row.setValue(dif_F_name, comptimdep)
			cursor3.updateRow(row)
		del(cursor2, cursor3)
		n=n+1
	sys.stdout.write("\rProcessed %i of %i" % (y, compcount))
	y=y+1	
arcpy.TableToTable_conversion (deptabp, locntab, deptabpname)
arcpy.TableToExcel_conversion (deptabp, deptabp + '.xls') 	
print ('.')

#Pivot erosion data3 table
voltab = locntab + 'Volumes'
erotabp = locntab + 'Erosion_Pivot'
erotabpname = 'Erosion_Pivot'
difs = [x[0][-13:] for x in os.walk(locndif)]

arcpy.CreateTable_management(locntab, 'Erosion_Pivot')
arcpy.AddField_management(erotabp, "COMP", "SHORT")
arcpy.DeleteField_management (erotabp, "FIELD1")
arcpy.DeleteField_management (erotabp, "OBJECTID")

for n in range (2, len(data3) + 1):
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]
	arcpy.AddField_management(erotabp, dif_F_name, "Double")
	
y = 1
for y in range (1, compcount + 1):
	compnum = str(y)
	cursor = arcpy.InsertCursor(erotabp)
	row = cursor.newRow()
	row.setValue("COMP", compnum)
	cursor.insertRow(row)
	del(cursor)
	for n in range (2, len(data3) + 1):
		dif = locndif + difs[n-1]
		difname = str(difs[n-1])
		dif_F_name = 'D' + difname[-12:]
		compid = str('"COMP" = ' + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		cursor2 = arcpy.SearchCursor(voltab,difcompid)
		for row in cursor2:
			comptimero = row.getValue("EROSION")
		cursor3 = arcpy.UpdateCursor(erotabp, compid)
		for row in cursor3:
			row.setValue(dif_F_name, comptimero)
			cursor3.updateRow(row)
		del(cursor2, cursor3)

	sys.stdout.write("\rProcessed %i of %i" % (y, compcount))

arcpy.TableToTable_conversion (erotabp, locntab, erotabpname)
arcpy.TableToExcel_conversion (erotabp, erotabp + '.xls') 	
print (".")

print (time.strftime("%H:%M:%S"))

locvec

#Create Vectors

if arcpy.Exists(locnvec + 'Vectors.gdb'):
	arcpy.Delete_management(locnvec + 'Vectors.gdb')

arcpy.CreateFileGDB_management(locnvec, 'Vectors')
	
arcpy.env.workspace = locnvec + 'Vectors.gdb'
arcpy.env.overwriteOutput = True
baseline = 'Baseline'
arcpy.CopyFeatures_management(baseline1, baseline)
arcpy.DefineProjection_management (baseline, sr)

arcpy.env.extent = "MAXOF"

arcpy.AddField_management(baseline, 'StartX', 'DOUBLE')
arcpy.AddField_management(baseline, 'StartY', 'DOUBLE')
arcpy.AddField_management(baseline, 'EndX', 'DOUBLE')
arcpy.AddField_management(baseline, 'EndY', 'DOUBLE')

arcpy.CalculateField_management(baseline, 'StartX', "!SHAPE.firstPoint.X!", "PYTHON3")
arcpy.CalculateField_management(baseline, 'StartY', "!SHAPE.firstPoint.Y!", "PYTHON3")
arcpy.CalculateField_management(baseline, 'EndX', "!SHAPE.lastPoint.X!", "PYTHON3")
arcpy.CalculateField_management(baseline, 'EndY', "!SHAPE.lastPoint.Y!", "PYTHON3")

cursor = arcpy.SearchCursor(baseline)

for row in cursor:
	startx = row.getValue("StartX")
	starty = row.getValue("StartY")
	endx = row.getValue("EndX")
	endy = row.getValue("EndY")

deltax = endx - startx
deltay = endy - starty

baselineradian = math.atan2(deltay,deltax)

print ("Shoreline is located in which direction relative to the baseline: left or right?")

if vectordirection == "left":
	vectorangle1 = baselineradian + math.pi/2
	vectorangle2 = baselineradian - math.pi/2
elif vectordirection == "right":
	vectorangle1 = baselineradian - math.pi/2
	vectorangle2 = baselineradian + math.pi/2
else:
	print ("Invalid answer, direction must be left or right.")

#Draw Vectors

print ("Draw vectors.")

voltab = locntab + 'Volumes'
deptabp = locntab + 'Deposition_Pivot'
difs = [x[0][-13:] for x in os.walk(locndif)]
for n in range (2, len(data3) + 1):
	
	dif = locndif + difs[n-1]
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]

	vectors = "Vectors_Deposition_" + difname
	dunedifout = voltab

	y = 1

	compid = str('"COMP" = ' + str(y))
	difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")

	q = startx + (baselinedistance * 0.5) * (math.cos(baselineradian))

	w = starty + (baselinedistance * 0.5) * (math.sin(baselineradian))

	cursor = arcpy.SearchCursor(voltab, difcompid)

	for row in cursor:
		volumedif = float(row.getValue("Deposition"))

	if volumedif >= 0:
		vectordistance = vectormultiplier * volumedif * - 1
		vectorradian = vectorangle2
	else:
		vectordistance = vectormultiplier * volumedif
		vectorradian = vectorangle1
		
	e = q + vectordistance * (math.cos(vectorradian))

	r = w + vectordistance * (math.sin(vectorradian))

	feature_info = [[[q, w], [e, r]],]
	y = 2
	for y in range (2, compcount + 1):
		compid = str('"COMP" = ' + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		cursor = arcpy.SearchCursor(voltab, difcompid)
		for row in cursor:
			volumedif = float(row.getValue("Deposition"))
		if volumedif >= 0:
			vectordistance = vectormultiplier * volumedif * - 1
			vectorradian = vectorangle2
		else:
			vectordistance = vectormultiplier * volumedif
			vectorradian = vectorangle1
		iterable = [[[q + baselinedistance * (math.cos(baselineradian)), w + baselinedistance * (math.sin(baselineradian))], [(q + baselinedistance * (math.cos(baselineradian))) + vectordistance * (math.cos(vectorradian)), (w + baselinedistance * (math.sin(baselineradian))) + vectordistance * (math.sin(vectorradian))]],]
		feature_info.extend(iterable)

		q = q + baselinedistance * (math.cos(baselineradian))
		
		w = w + baselinedistance * (math.sin(baselineradian))

		sys.stdout.write("\rProcessed %i of %i vectors for duration %s" % (y, compcount, difname))
		
	print (".")
	
	features = []

	for feature in feature_info:
		features.append(
			arcpy.Polyline(
				arcpy.Array([arcpy.Point(*coords) for coords in feature]),sr))
	
		
	arcpy.CopyFeatures_management(features, vectors)
	arcpy.AddField_management(vectors, "Comp", "Integer")
	arcpy.AddField_management(vectors, "Deposition", "Double")

	n = 1

	for y in range (1, compcount + 1):
		zone = str(y)
		compid = str('"COMP" = ' + str(y))
		objid = str("OBJECTID = " + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		compid = str('"COMP" = ' + str(y))
		cursor = arcpy.SearchCursor(voltab, difcompid)
		for row in cursor:
			volumedif = float(row.getValue("Deposition"))
		arcpy.MakeFeatureLayer_management(vectors, "lyr1")
		arcpy.SelectLayerByAttribute_management("lyr1", "NEW_SELECTION", objid)
		arcpy.CalculateField_management("lyr1", "Comp", zone, "PYTHON3")
		arcpy.CalculateField_management("lyr1", "Deposition", volumedif, "PYTHON3")
		sys.stdout.write("\rCaluculated %i of %i vector attributes for duration %s" % (y, compcount, difname))
		
	print (".")

print ("Vectors drawn.")
#Draw Vectors 2

print ("Draw vectors.")

voltab = locntab + 'Volumes'
erotabp = locntab + 'Erosion_Pivot'
difs = [x[0][-13:] for x in os.walk(locndif)]
for n in range (2, len(data3) + 1):
	
	dif = locndif + difs[n-1]
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]

	vectors = "Vectors_Erosion_" + difname
	dunedifout = voltab

	y = 1

	compid = str('"COMP" = ' + str(y))
	difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")

	q = startx + (baselinedistance * 0.5) * (math.cos(baselineradian))

	w = starty + (baselinedistance * 0.5) * (math.sin(baselineradian))

	cursor = arcpy.SearchCursor(voltab, difcompid)

	for row in cursor:
		volumedif = float(row.getValue("Erosion")*-1)

	if volumedif >= 0:
		vectordistance = vectormultiplier * volumedif * - 1
		vectorradian = vectorangle2
	else:
		vectordistance = vectormultiplier * volumedif
		vectorradian = vectorangle1
		
	e = q + vectordistance * (math.cos(vectorradian))

	r = w + vectordistance * (math.sin(vectorradian))

	feature_info = [[[q, w], [e, r]],]
	y = 2

	for y in range (2, compcount + 1):
		compid = str('"COMP" = ' + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		cursor = arcpy.SearchCursor(voltab, difcompid)
		for row in cursor:
			volumedif = float(row.getValue("Erosion")*-1)
		if volumedif >= 0:
			vectordistance = vectormultiplier * volumedif * - 1
			vectorradian = vectorangle2
		else:
			vectordistance = vectormultiplier * volumedif
			vectorradian = vectorangle1
		iterable = [[[q + baselinedistance * (math.cos(baselineradian)), w + baselinedistance * (math.sin(baselineradian))], [(q + baselinedistance * (math.cos(baselineradian))) + vectordistance * (math.cos(vectorradian)), (w + baselinedistance * (math.sin(baselineradian))) + vectordistance * (math.sin(vectorradian))]],]
		feature_info.extend(iterable)

		q = q + baselinedistance * (math.cos(baselineradian))
		
		w = w + baselinedistance * (math.sin(baselineradian))

		sys.stdout.write("\rProcessed %i of %i vectors for duration %s" % (y, compcount, difname))
		
	print (".")
	
	features = []

	for feature in feature_info:
		features.append(
			arcpy.Polyline(
				arcpy.Array([arcpy.Point(*coords) for coords in feature]),sr))
			
	arcpy.CopyFeatures_management(features, vectors)
	arcpy.AddField_management(vectors, "Comp", "Integer")
	arcpy.AddField_management(vectors, "Erosion", "Double")

	n = 1

	for y in range (1, compcount + 1):
		zone = str(y)
		compid = str('"COMP" = ' + str(y))
		objid = str("OBJECTID = " + str(y))
		difcompid = str('"COMP" = ' + str(y)) + "AND" + str('"DURATION" = ' + "'" + difname + "'")
		compid = str('"COMP" = ' + str(y))
		cursor = arcpy.SearchCursor(voltab, difcompid)
		for row in cursor:
			volumedif = float(row.getValue("Erosion")*-1)
		arcpy.MakeFeatureLayer_management(vectors, "lyr1")
		arcpy.SelectLayerByAttribute_management("lyr1", "NEW_SELECTION", objid)
		arcpy.CalculateField_management("lyr1", "Comp", zone, "PYTHON3")
		arcpy.CalculateField_management("lyr1", "Erosion", volumedif, "PYTHON3")
		sys.stdout.write("\rCaluculated %i of %i vector attributes for duration %s" % (y, compcount, difname))
		
	print (".")

print ("Vectors drawn.")

locnvecgdb = locnvec + 'Vectors.gdb'

#Create Maps

#file locations for datalayers

baslyr_path = loclyr + "Baseline.lyrx"
veclyre_path = loclyr + "Vectors_E.lyrx"
veclyrd_path = loclyr + "Vectors_D.lyrx"
diflyr_path = loclyr + "Differences.lyrx"
mrklyr_path = loclyr + "Markers.lyrx"

baslyr_lf = arcpy.mp.LayerFile(baslyr_path)
veclyre_lf = arcpy.mp.LayerFile(veclyre_path)
veclyrd_lf = arcpy.mp.LayerFile(veclyrd_path)
diflyr_lf = arcpy.mp.LayerFile(diflyr_path)
mrklyr_lf = arcpy.mp.LayerFile(mrklyr_path)

difs = [x[0][-13:] for x in os.walk(locndif)]

n=2
voltabp = locntab + 'Volumes_Pivot'

for n in range (2, len(data3) + 1):
	dif = locndif + difs[n-1]
	difname = str(difs[n-1])
	comp = loccomp + 'Compartments.shp'
	compid  = str('"COMP" < ' + str(1000))
	temp1 = "in_memory/t2" + integ
	arcpy.MakeFeatureLayer_management(comp, "lyr")
	arcpy.SelectLayerByAttribute_management("lyr", "NEW_SELECTION", compid)
	outRaster = ExtractByMask(dif,"lyr")
	outRaster.save(temp1)	
	arcpy.CopyRaster_management(temp1, dif) 
	sys.stdout.write("\rClipped difference raster %s" % (difname))
	arcpy.Delete_management("in_memory")
	
	
print (".")	
	

for n in range (2, len(data3) + 1):
	dif = locndif + difs[n-1]
	difname = str(difs[n-1])
	dif_F_name = 'D' + difname[-12:]
	
	searchid = "MEAN_" + str(dif_F_name)
	
	compid1 = str('"COMP" <= ' + str(compbreak1))
	temp_table = "in_memory/tempt" + integ
	arcpy.TableSelect_analysis (voltabp, temp_table, compid1)
	stats_table = "in_memory/stats" + integ
	arcpy.Statistics_analysis(temp_table,stats_table,[[dif_F_name,"MEAN"]])
	cursor = arcpy.SearchCursor(stats_table)
	for row in cursor:
		vol1 = str(round(row.getValue(searchid),2))

	compid2 = str('"COMP" > ' + str(compbreak1) + 'AND "COMP" <= ' + str(compbreak2))
	temp_table = "in_memory/tempt" + integ
	arcpy.TableSelect_analysis (voltabp, temp_table, compid2)
	stats_table = "in_memory/stats" + integ
	arcpy.Statistics_analysis(temp_table,stats_table,[[dif_F_name,"MEAN"]])
	cursor = arcpy.SearchCursor(stats_table)
	for row in cursor:
		vol2 = str(round(row.getValue(searchid),2))
	
	compid3 = str('"COMP" > ' + str(compbreak2))
	temp_table = "in_memory/tempt" + integ
	arcpy.TableSelect_analysis (voltabp, temp_table, compid3)
	stats_table = "in_memory/stats" + integ
	arcpy.Statistics_analysis(temp_table,stats_table,[[dif_F_name,"MEAN"]])
	cursor = arcpy.SearchCursor(stats_table)
	for row in cursor:
		vol3 = str(round(row.getValue(searchid),2))
	
	stats_table = "in_memory/stats" + integ
	arcpy.Statistics_analysis(voltabp,stats_table,[[dif_F_name,"MEAN"]])
	cursor = arcpy.SearchCursor(stats_table)
	for row in cursor:
		vol4 = str(round(row.getValue(searchid),2))
		
	totalm3_1 = float(float(vol4) * float(compcount))
	totalm3 = str(round(totalm3_1,2))
		
	Comp_1 = "Mean Volumetric Change of" + "\n" + "Compartments 1 - " + str(compbreak1) + ": " + str(vol1 + " m<SUP>3</SUP>")
	Comp_2 = "Mean Volumetric Change of" + "\n" + "Compartments " + str(compbreak1+1) + " - " + str(compbreak2) + ": " + str(vol2 + " m<SUP>3</SUP>")
	Comp_3 = "Mean Volumetric Change of" + "\n" + "Compartments " + str(compbreak2+1) + " - " + str(compcount) + ": " + str(vol3 + " m<SUP>3</SUP>")
	Comp_All = "Mean Volumetric Change" + "\n" + "per Compartment: " + str(vol4 + " m<SUP>3</SUP>")
	Total_Change = "Total Volumetric Change: " + str(totalm3 + " m<SUP>3</SUP>")

	title = site[0:3] + " Erosion and Deposition" + "\n" + titles2[n-2]
	
	vectorse = "Vectors_Erosion_" + difname
	vectorsd = "Vectors_Deposition_" + difname

	mapname = "Layers"  # Make sure a map named "Template" exists in the APRX
	aprx_path = locmap + "Template/Template.aprx"
	aprx = arcpy.mp.ArcGISProject(aprx_path)
	maps = aprx.listMaps("Layers")
	if not maps:
		raise RuntimeError("No map named 'Layers' found in the project.")

	map_obj = maps[0]

	veclyre_name = veclyre_lf.listLayers()[0].name
	veclyrd_name = veclyrd_lf.listLayers()[0].name
	diflyr_name = diflyr_lf.listLayers()[0].name
	mrklyr_name = mrklyr_lf.listLayers()[0].name
	baslyr_name = baslyr_lf.listLayers()[0].name
	
	layers_to_remove = [vecdlname, vecelname, diflname, mrklname, baselname]
	for layer in map_obj.listLayers():
		if layer.name in layers_to_remove:
			map_obj.removeLayer(layer)

	map_obj.addDataFromPath(diflyr_path)
	map_obj.addDataFromPath(veclyre_path)
	map_obj.addDataFromPath(veclyrd_path)
	map_obj.addDataFromPath(mrklyr_path)
	map_obj.addDataFromPath(baslyr_path)
		
	for lyr in map_obj.listLayers():
			
		if lyr.name == veclyre_name:

			old_conn = lyr.connectionProperties
			new_conn = old_conn.copy()			
			new_conn['connection_info']['database'] = locnvecgdb
			new_conn['dataset'] = vectorse
			lyr.updateConnectionProperties(old_conn, new_conn)
			
		elif lyr.name == veclyrd_name:

			old_conn = lyr.connectionProperties
			new_conn = old_conn.copy()			
			new_conn['connection_info']['database'] = locnvecgdb
			new_conn['dataset'] = vectorsd
			lyr.updateConnectionProperties(old_conn, new_conn)
					
		elif lyr.name == diflyr_name:
			
			old_conn = lyr.connectionProperties
			new_conn = old_conn.copy()			
			new_conn['connection_info']['database'] = locndif
			new_conn['dataset'] = difname
			lyr.updateConnectionProperties(old_conn, new_conn)
		
	layouts = aprx.listLayouts()
	if layouts:
		layout = layouts[0] 
		for element in layout.listElements("TEXT_ELEMENT"):
			if element.name == "Title":
				element.text = title
			elif element.name == "Comp_1":
				element.text = Comp_1
			elif element.name == "Comp_2":
				element.text = Comp_2
			elif element.name == "Comp_3":
				element.text = Comp_3
			elif element.name == "Comp_All":
				element.text = Comp_All
			elif element.name == "Total_Change":
				element.text = Total_Change

		print("Updated layout elements.")

	aprx.saveACopy(locmap + difname + ".aprx")
	print (titles2[n-2] + " map created")		
	
	del aprx
	
difs = [x[0][-13:] for x in os.walk(locndif)]

#Export Maps

n = 2

for n in range (2, len(data3) + 1):
	difname = str(difs[n-1])
	mapname = locmap + difname + ".aprx" 
	jpggname = locimg + difname + ".jpg"
	aprx = arcpy.mp.ArcGISProject(mapname)
	maps = aprx.listMaps()
	if maps:
		map_obj = maps[0]
		layout = aprx.listLayouts()[0]
		layout.exportToJPEG(jpggname, resolution=300)
	del aprx
	print (titles2[n-2] + " map exported as jpeg")
	
difs = [x[0][-13:] for x in os.walk(locdif)]

#Export Maps

n = 2
for n in range (2, len(data) + 1):
	difname = str(difs[n-1])
	mapname = locmap + difname + ".aprx" 
	jpggname = locimg + difname + ".jpg"
	aprx = arcpy.mp.ArcGISProject(mapname)
	maps = aprx.listMaps()
	if maps:
		map_obj = maps[0]
		layout = aprx.listLayouts()[0]
		layout.exportToJPEG(jpggname, resolution=300)
	del aprx
	print (titles[n-2] + " map exported as jpeg")