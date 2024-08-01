# Plots the optical efficiency files in this directory using polar coordinate system. 
# Reads declination-hra data
# Also overlays sunpath for a given latitude.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import csv

def ThetaRValues(filename,index,lat): #lat is inputted in radians
	"""Generates values of theta, r and efficiency by reading a file"""
	#f = open(filename+"_"+index+".txt","r")
	f = open(filename+".motab","r")

	data = f.readlines()

	#data = data[2:]
	data = data[index-1:index+9]
	matrix = []


	i = 2
	for line in data:
		matrix.append(line[0:-3].split(" "))

	#print(matrix)
	matrix = np.array(matrix).astype(np.float)
	#print(np.shape(matrix))
	#print(matrix)

	#extract not azimuths but HRAs
	HRAs = (matrix)[0][1:]
	HRAs = HRAs*np.pi/180.0 #convert to rad
	#print(HRAs)
	#azimuths = (matrix)[0][1:]
	#print(HRAs)
	

	#extract not elevations but declinations
	dec = (matrix[:,0])[1:]
	dec = dec*np.pi/180.0 #convert to rad
	#elevations = (matrix[:,0])[1:]
	#print(dec)


	#Now we need to convert dec-hra to Azimuth elevation
	#elevations = np.arcsin(np.sin(dec)*np.sin(lat)+np.cos(dec)*np.cos(lat)*np.cos(HRAs))
	#azimuths = np.arccos((np.sin(dec)*np.cos(lat)-np.cos(dec)*np.sin(lat)*np.cos(HRAs))/(np.cos(elevations)))

	#for row in azimuths: #correct azimuths to fit within a range
	#	for col in row:
	#		if col > 0.0:
	#			col = 2.0*np.pi-col


	values = (matrix[1:,1:]).T #transpose so it is rows of zeniths stacked up


	#zeniths = 90.0 - elevations
	#zeniths = 0.5*np.pi - elevations #radians
	#azimuths = np.radians(azimuths) #already converted
	#r , theta = np.meshgrid(zeniths,azimuths)

	r_dec , theta_HRAs = np.meshgrid(dec,HRAs) #Not yet converted to ele,azi

	r = np.arcsin(np.sin(r_dec)*np.sin(lat)+np.cos(r_dec)*np.cos(lat)*np.cos(theta_HRAs)) # r is now elevation (rad) between negative pi and pi
	
	for row in r:
		for col in row:
			if col < 0.0:
				col = col + 2.0*np.pi #now between 0 and 2pi
	r = r*(180.0/np.pi) #convert to deg
	r = 90.0 - r #convert to zenith

	argument = (np.sin(r_dec)*np.cos(lat)-np.cos(r_dec)*np.sin(lat)*np.cos(theta_HRAs))/(np.cos((90.0-r)*np.pi/180.0))

	#print(argument)
	i = 0
	j = 0
	while i < np.shape(argument)[0]:
		while j < np.shape(argument)[1]:
			if argument[i][j] < -0.999:
				argument[i][j] = -0.999
			if argument[i][j] > 0.999:
				argument[i][j] = 0.999
			j += 1
		j = 0
		i += 1

	
	theta = np.arccos(argument) #theta is now azimuth



	#print(theta)
	#print(np.shape(theta))

	i = 0
	j = 0

	while i < np.shape(theta_HRAs)[0]:
		while j < np.shape(theta_HRAs)[1]:
			if theta_HRAs[i][j] >= 0:
				theta[i][j] = 2.0*np.pi - theta[i][j] #correct for evening
			if theta[i][j] > np.pi:
				theta[i][j] = theta[i][j] - 2.0*np.pi #bring it between neg pi and pi
			j += 1
		j = 0
		i += 1

	#for row in theta:
		#for col in row:
			#if col < 0.0:
				#col = col + 2.0*np.pi #now between 0 and 2pi	
	#theta = theta-(2.0*np.pi)
	print(theta)

	return [theta,r,values]

def Polar_Single(filename,latitude): #latitude is in degrees, lat is in radians
	lat = np.radians(latitude)
	HRA_Summer = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(0.40928)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(0.40928)),100)
	HRA_Spring = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(0.0)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(0.0)),100)
	HRA_Winter = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(-0.40928)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(-0.40928)),100)

	Elev_Summer = np.arcsin(np.sin(0.40928)*np.sin(lat)+np.cos(0.40928)*np.cos(lat)*np.cos(HRA_Summer))
	Elev_Spring = np.arcsin(np.sin(0.0)*np.sin(lat)+np.cos(0.0)*np.cos(lat)*np.cos(HRA_Spring))
	Elev_Winter = np.arcsin(np.sin(-0.40928)*np.sin(lat)+np.cos(-0.40928)*np.cos(lat)*np.cos(HRA_Winter))

	Azi_Summer = np.arccos((np.sin(0.40928)*np.cos(lat)-np.cos(0.40928)*np.sin(lat)*np.cos(HRA_Summer))/(np.cos(Elev_Summer)))
	Azi_Spring = np.arccos((np.sin(0.0)*np.cos(lat)-np.cos(0.0)*np.sin(lat)*np.cos(HRA_Spring))/(np.cos(Elev_Spring)))
	Azi_Winter = np.arccos((np.sin(-0.40928)*np.cos(lat)-np.cos(-0.40928)*np.sin(lat)*np.cos(HRA_Winter))/(np.cos(Elev_Winter)))

	i = 0
	while i < len(HRA_Summer):
		if HRA_Summer[i] > 0.0:
			Azi_Summer[i] = 2.0*np.pi-Azi_Summer[i]
		i += 1
	i = 0
	while i < len(HRA_Spring):
		if HRA_Spring[i] > 0.0:
			Azi_Spring[i] = 2.0*np.pi-Azi_Spring[i]
		i += 1
	i = 0
	while i < len(HRA_Winter):
		if HRA_Winter[i] > 0.0:
			Azi_Winter[i] = 2.0*np.pi-Azi_Winter[i]
		i += 1

	Zenith_Summer = 90.0 - np.degrees(Elev_Summer)
	Zenith_Spring = 90.0 - np.degrees(Elev_Spring)
	Zenith_Winter = 90.0 - np.degrees(Elev_Winter)

	#index of optical efficiency is 7
	opt_list = ThetaRValues(filename,7,lat)
	fig = plt.figure()
	grid = plt.GridSpec(1,12,wspace=0.1,hspace=0.1)

	ax1 = fig.add_subplot(grid[0,0:10],projection="polar")
	ax2 = fig.add_subplot(grid[0,11]) #1st Colorbar
	ax1.set_ylim(0.0,90.0)

	ax1.contourf(opt_list[0],opt_list[1],opt_list[2],np.linspace(0.0,1.00,41),extend='neither',cmap="jet",)

	ax1.plot(Azi_Summer,Zenith_Summer,"k-",linewidth=1)
	ax1.plot(Azi_Spring,Zenith_Spring,"k--",linewidth=1)
	ax1.plot(Azi_Winter,Zenith_Winter,"k-",linewidth=1)

	#print(opt_list[0])
	#print(opt_list[1])
	
	#print(Azi_Spring)`
	#print(Zenith_Spring)

	ax1.set_theta_offset(np.pi/2.0)

	cNorm1 = matplotlib.colors.Normalize(vmin=0.0,vmax=1.0)
	cbar1 = matplotlib.colorbar.ColorbarBase(ax2, norm=cNorm1,cmap="jet",ticks = np.linspace(0.0,1.0,21),label="Optical Efficiency")
	plt.show()
	return


def Polar(filename,field_info,location,latitude):

	lat = np.radians(latitude)
	HRA_Summer = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(0.40928)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(0.40928)),100)
	HRA_Spring = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(0.0)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(0.0)),100)
	HRA_Winter = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(-0.40928)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(-0.40928)),100)

	Elev_Summer = np.arcsin(np.sin(0.40928)*np.sin(lat)+np.cos(0.40928)*np.cos(lat)*np.cos(HRA_Summer))
	Elev_Spring = np.arcsin(np.sin(0.0)*np.sin(lat)+np.cos(0.0)*np.cos(lat)*np.cos(HRA_Spring))
	Elev_Winter = np.arcsin(np.sin(-0.40928)*np.sin(lat)+np.cos(-0.40928)*np.cos(lat)*np.cos(HRA_Winter))

	Azi_Summer = np.arccos((np.sin(0.40928)*np.cos(lat)-np.cos(0.40928)*np.sin(lat)*np.cos(HRA_Summer))/(np.cos(Elev_Summer)))
	Azi_Spring = np.arccos((np.sin(0.0)*np.cos(lat)-np.cos(0.0)*np.sin(lat)*np.cos(HRA_Spring))/(np.cos(Elev_Spring)))
	Azi_Winter = np.arccos((np.sin(-0.40928)*np.cos(lat)-np.cos(-0.40928)*np.sin(lat)*np.cos(HRA_Winter))/(np.cos(Elev_Winter)))

	i = 0
	while i < len(HRA_Summer):
		if HRA_Summer[i] > 0.0:
			Azi_Summer[i] = 2.0*np.pi-Azi_Summer[i]
		i += 1
	i = 0
	while i < len(HRA_Spring):
		if HRA_Spring[i] > 0.0:
			Azi_Spring[i] = 2.0*np.pi-Azi_Spring[i]
		i += 1
	i = 0
	while i < len(HRA_Winter):
		if HRA_Winter[i] > 0.0:
			Azi_Winter[i] = 2.0*np.pi-Azi_Winter[i]
		i += 1

	Zenith_Summer = 90.0 - np.degrees(Elev_Summer)
	Zenith_Spring = 90.0 - np.degrees(Elev_Spring)
	Zenith_Winter = 90.0 - np.degrees(Elev_Winter)
	#Indices are 1.Cos, 2.Shade, 3.Block, 4.Spill, 5.OptEff
	#7.DNI 8.OptEff
	cos_list = ThetaRValues(filename,"CosEff")
	shade_list = ThetaRValues(filename,"ShadeEff")
	block_list = ThetaRValues(filename,"BlockEff")
	spill_list = ThetaRValues(filename,"SpillEff")
	opt_list = ThetaRValues(filename,"OptEff")

	fig = plt.figure()

	grid = plt.GridSpec(20,34,wspace=0.1,hspace=0.1)

	ax1 = fig.add_subplot(grid[0:4,0:7],projection="polar")
	ax2 = fig.add_subplot(grid[0:4,8:15],projection="polar")
	ax3 = fig.add_subplot(grid[5:9,0:7],projection="polar")
	ax4 = fig.add_subplot(grid[5:9,8:15],projection="polar")
	ax5 = fig.add_subplot(grid[0:9,17:32],projection="polar")
	ax6 = fig.add_subplot(grid[0:9,33]) #1st Colorbar
	ax7 = fig.add_subplot(grid[10:19,0:15],projection = "polar")
	ax8 = fig.add_subplot(grid[10:19,17:32],projection="polar")
	ax9 = fig.add_subplot(grid[10:19,33]) #2nd ColorBar

	ax1.contourf(cos_list[0],cos_list[1],cos_list[2],np.linspace(0.0,1.00,41),extend='neither',cmap="jet")
	ax2.contourf(shade_list[0],shade_list[1],shade_list[2],np.linspace(0.0,1.00,41),extend='max',cmap="jet")
	ax3.contourf(block_list[0],block_list[1],block_list[2],np.linspace(0.0,1.00,41),extend='neither',cmap="jet")
	ax4.contourf(spill_list[0],spill_list[1],spill_list[2],np.linspace(0.0,1.00,41),extend='neither',cmap="jet")
	ax5.contourf(opt_list[0],opt_list[1],opt_list[2],np.linspace(0.0,1.00,41),extend='neither',cmap="jet")

	#Polt interpolation dots
	ax5.plot(opt_list[0],opt_list[1],"ko",markersize=1)

	#Plot Sunpath
	ax1.plot(Azi_Summer,Zenith_Summer,"k-",linewidth=1)
	ax1.plot(Azi_Spring,Zenith_Spring,"k--",linewidth=1)
	ax1.plot(Azi_Winter,Zenith_Winter,"k-",linewidth=1)

	ax2.plot(Azi_Summer,Zenith_Summer,"k-",linewidth=1)
	ax2.plot(Azi_Spring,Zenith_Spring,"k--",linewidth=1)
	ax2.plot(Azi_Winter,Zenith_Winter,"k-",linewidth=1)

	ax3.plot(Azi_Summer,Zenith_Summer,"k-",linewidth=1)
	ax3.plot(Azi_Spring,Zenith_Spring,"k--",linewidth=1)
	ax3.plot(Azi_Winter,Zenith_Winter,"k-",linewidth=1)

	ax4.plot(Azi_Summer,Zenith_Summer,"k-",linewidth=1)
	ax4.plot(Azi_Spring,Zenith_Spring,"k--",linewidth=1)
	ax4.plot(Azi_Winter,Zenith_Winter,"k-",linewidth=1)

	ax5.plot(Azi_Summer,Zenith_Summer,"k-",linewidth=1)
	ax5.plot(Azi_Spring,Zenith_Spring,"k--",linewidth=1)
	ax5.plot(Azi_Winter,Zenith_Winter,"k-",linewidth=1)

	#DNI Calculaions
	#Remember that zeniths are in degrees
	AM = (np.cos(np.radians(opt_list[1])))**-1.0
	DNI = 1367*(0.7**(AM**0.678))

	ax7.contourf(opt_list[0],opt_list[1],DNI,np.linspace(0.0,1000.0,41),extend='neither',cmap="jet")
	ax7.plot(Azi_Summer,Zenith_Summer,"k-",linewidth=1)
	ax7.plot(Azi_Spring,Zenith_Spring,"k--",linewidth=1)
	ax7.plot(Azi_Winter,Zenith_Winter,"k-",linewidth=1)

	ax8.contourf(opt_list[0],opt_list[1],opt_list[2]*DNI,np.linspace(0.0,1000.0,41),extend='neither',cmap="jet")
	ax8.plot(Azi_Summer,Zenith_Summer,"k-",linewidth=1)
	ax8.plot(Azi_Spring,Zenith_Spring,"k--",linewidth=1)
	ax8.plot(Azi_Winter,Zenith_Winter,"k-",linewidth=1)
	
	#Labels and hide the zeniths
	ax1.set_yticklabels([])
	ax2.set_yticklabels([])
	ax3.set_yticklabels([])
	ax4.set_yticklabels([])
	ax1.set_xticklabels([])
	ax2.set_xticklabels([])
	ax3.set_xticklabels([])
	ax4.set_xticklabels([])

	ax1.set_title("Cos")
	ax2.set_title("Shade")
	ax3.set_title("Block")
	ax4.set_title("Spill")
	ax5.set_title("OptEff")
	ax7.set_title(r"DNI = $1367 \times 0.7^{AM^{0.678}}$")
	ax8.set_title("DNI*OptEff")
	

	fig.suptitle("Collector_RecvArea = "+field_info+"   "+"Location = "+location+"\n"+"Latitude = "+str(np.degrees(lat)),fontsize="x-large")

	cNorm1 = matplotlib.colors.Normalize(vmin=0.0,vmax=1.0)
	cNorm2 = matplotlib.colors.Normalize(vmin=0.0,vmax=1000.0)
	cbar1 = matplotlib.colorbar.ColorbarBase(ax6, norm=cNorm1,cmap="jet",ticks = np.linspace(0.0,1.0,21),label="Efficiency")
	cbar2 = matplotlib.colorbar.ColorbarBase(ax9, norm=cNorm2,cmap="jet",ticks = np.linspace(0.0,1000.0,21),label="Power (W/m2)")

	fig.set_size_inches(8.27,11.69)
	plt.subplots_adjust(left=0.05, right=0.90, top=0.90, bottom=0.10)
	plt.savefig("PolarContour.png",dpi=100)
	return

def Sun_Path(latitude):
	lat = np.radians(latitude)
	declination = 23.45
	dec = declination*np.pi/180.0
	HRA_Summer = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(dec)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(dec)),100)
	HRA_Spring = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(0.0)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(0.0)),100)
	HRA_Winter = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(-1.0*dec)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(-1.0*dec)),100)

	Elev_Summer = np.arcsin(np.sin(dec)*np.sin(lat)+np.cos(dec)*np.cos(lat)*np.cos(HRA_Summer))
	Elev_Spring = np.arcsin(np.sin(0.0)*np.sin(lat)+np.cos(0.0)*np.cos(lat)*np.cos(HRA_Spring))
	Elev_Winter = np.arcsin(np.sin(-1.0*dec)*np.sin(lat)+np.cos(-1.0*dec)*np.cos(lat)*np.cos(HRA_Winter))

	Azi_Summer = np.arccos((np.sin(dec)*np.cos(lat)-np.cos(dec)*np.sin(lat)*np.cos(HRA_Summer))/(np.cos(Elev_Summer)))
	Azi_Spring = np.arccos((np.sin(0.0)*np.cos(lat)-np.cos(0.0)*np.sin(lat)*np.cos(HRA_Spring))/(np.cos(Elev_Spring)))
	Azi_Winter = np.arccos((np.sin(-1.0*dec)*np.cos(lat)-np.cos(-1.0*dec)*np.sin(lat)*np.cos(HRA_Winter))/(np.cos(Elev_Winter)))

	i = 0
	while i < len(HRA_Summer):
		if HRA_Summer[i] > 0.0:
			Azi_Summer[i] = 2.0*np.pi-Azi_Summer[i]
		i += 1
	i = 0
	while i < len(HRA_Spring):
		if HRA_Spring[i] > 0.0:
			Azi_Spring[i] = 2.0*np.pi-Azi_Spring[i]
		i += 1
	i = 0
	while i < len(HRA_Winter):
		if HRA_Winter[i] > 0.0:
			Azi_Winter[i] = 2.0*np.pi-Azi_Winter[i]
		i += 1

	Zenith_Summer = 90.0 - np.degrees(Elev_Summer)
	Zenith_Spring = 90.0 - np.degrees(Elev_Spring)
	Zenith_Winter = 90.0 - np.degrees(Elev_Winter)

	
	grid = plt.GridSpec(1,12,wspace=0.1,hspace=0.1)
	fig = plt.figure()
	ax1 = fig.add_subplot(grid[0,0:10],projection="polar")

	#ax1 = fig.add_subplot(projection="polar")
	ax1.set_ylim(0.0,90.0)

	#ax1.contourf(opt_list[0],opt_list[1],opt_list[2],np.linspace(0.0,1.00,41),extend='neither',cmap="jet",)

	if latitude > 0:
		ax1.plot(Azi_Summer,Zenith_Summer,"r-",linewidth=1)
		ax1.plot(Azi_Spring,Zenith_Spring,"k--",linewidth=1)
		ax1.plot(Azi_Winter,Zenith_Winter,"b-",linewidth=1)
	else:
		ax1.plot(Azi_Summer,Zenith_Summer,"b-",linewidth=1)
		ax1.plot(Azi_Spring,Zenith_Spring,"k--",linewidth=1)
		ax1.plot(Azi_Winter,Zenith_Winter,"r-",linewidth=1)

	ax1.set_theta_offset(np.pi/2.0)
	plt.show()
	return

def Sun_Path2(latitude,declination=23.45): #Plots the 12 solar months of the year
	lat = np.radians(latitude)
	dec = declination*np.pi/180.0
	HRA_Summer = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(dec)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(dec)),100)
	HRA_Spring = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(0.0)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(0.0)),100)
	HRA_Winter = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(-1.0*dec)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(-1.0*dec)),100)

	dec2 = dec*np.sin(((31.0-1.0)/364.0)*2.0*np.pi) 
	dec3 = dec*np.sin(((61.0-1.0)/364.0)*2.0*np.pi) 
	dec8 = dec*np.sin(((213.0-1.0)/364.0)*2.0*np.pi) 
	dec9 = dec*np.sin(((243.0-1.0)/364.0)*2.0*np.pi) 

	print(dec2)
	print(dec3)
	print(dec8)
	print(dec9)

	Elev_Summer = np.arcsin(np.sin(dec)*np.sin(lat)+np.cos(dec)*np.cos(lat)*np.cos(HRA_Summer))
	Elev_Spring = np.arcsin(np.sin(0.0)*np.sin(lat)+np.cos(0.0)*np.cos(lat)*np.cos(HRA_Spring))
	Elev_Winter = np.arcsin(np.sin(-1.0*dec)*np.sin(lat)+np.cos(-1.0*dec)*np.cos(lat)*np.cos(HRA_Winter))

	Azi_Summer = np.arccos((np.sin(dec)*np.cos(lat)-np.cos(dec)*np.sin(lat)*np.cos(HRA_Summer))/(np.cos(Elev_Summer)))
	Azi_Spring = np.arccos((np.sin(0.0)*np.cos(lat)-np.cos(0.0)*np.sin(lat)*np.cos(HRA_Spring))/(np.cos(Elev_Spring)))
	Azi_Winter = np.arccos((np.sin(-1.0*dec)*np.cos(lat)-np.cos(-1.0*dec)*np.sin(lat)*np.cos(HRA_Winter))/(np.cos(Elev_Winter)))

	HRA_2 = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(dec2)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(dec2)),100)
	HRA_3 = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(dec3)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(dec3)),100)
	HRA_8 = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(dec8)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(dec8)),100)
	HRA_9 = np.linspace(-1.0*np.arccos(-1.0*np.tan(lat)*np.tan(dec9)),+1.0*np.arccos(-1.0*np.tan(lat)*np.tan(dec9)),100)

	print(HRA_2[-1])
	print(HRA_3[-1])
	print(HRA_8[-1])
	print(HRA_9[-1])


	Elev_2 = np.arcsin(np.sin(dec2)*np.sin(lat)+np.cos(dec2)*np.cos(lat)*np.cos(HRA_2))
	Elev_3 = np.arcsin(np.sin(dec3)*np.sin(lat)+np.cos(dec3)*np.cos(lat)*np.cos(HRA_3))
	Elev_8 = np.arcsin(np.sin(dec8)*np.sin(lat)+np.cos(dec8)*np.cos(lat)*np.cos(HRA_8))
	Elev_9 = np.arcsin(np.sin(dec9)*np.sin(lat)+np.cos(dec9)*np.cos(lat)*np.cos(HRA_9))

	Azi_2 = np.arccos((np.sin(dec2)*np.cos(lat)-np.cos(dec2)*np.sin(lat)*np.cos(HRA_2))/(np.cos(Elev_2)))	
	Azi_3 = np.arccos((np.sin(dec3)*np.cos(lat)-np.cos(dec3)*np.sin(lat)*np.cos(HRA_3))/(np.cos(Elev_3)))	
	Azi_8 = np.arccos((np.sin(dec8)*np.cos(lat)-np.cos(dec8)*np.sin(lat)*np.cos(HRA_8))/(np.cos(Elev_8)))	
	Azi_9 = np.arccos((np.sin(dec9)*np.cos(lat)-np.cos(dec9)*np.sin(lat)*np.cos(HRA_9))/(np.cos(Elev_9)))	

	i = 0
	while i < len(HRA_Summer):
		if HRA_Summer[i] > 0.0:
			Azi_Summer[i] = 2.0*np.pi-Azi_Summer[i]
		i += 1
	i = 0
	while i < len(HRA_Spring):
		if HRA_Spring[i] > 0.0:
			Azi_Spring[i] = 2.0*np.pi-Azi_Spring[i]
		i += 1
	i = 0
	while i < len(HRA_Winter):
		if HRA_Winter[i] > 0.0:
			Azi_Winter[i] = 2.0*np.pi-Azi_Winter[i]
		i += 1

	i = 0
	while i < len(HRA_2):
		if HRA_2[i] > 0.0:
			Azi_2[i] = 2.0*np.pi-Azi_2[i]
		i += 1
	i = 0
	while i < len(HRA_3):
		if HRA_3[i] > 0.0:
			Azi_3[i] = 2.0*np.pi-Azi_3[i]
		i += 1
	i = 0
	while i < len(HRA_8):
		if HRA_8[i] > 0.0:
			Azi_8[i] = 2.0*np.pi-Azi_8[i]
		i += 1
	i = 0
	while i < len(HRA_9):
		if HRA_9[i] > 0.0:
			Azi_9[i] = 2.0*np.pi-Azi_9[i]
		i += 1

	Zenith_Summer = 90.0 - np.degrees(Elev_Summer)
	Zenith_Spring = 90.0 - np.degrees(Elev_Spring)
	Zenith_Winter = 90.0 - np.degrees(Elev_Winter)

	Zenith_2 = 90.0 - np.degrees(Elev_2)
	Zenith_3 = 90.0 - np.degrees(Elev_3)
	Zenith_8 = 90.0 - np.degrees(Elev_8)
	Zenith_9 = 90.0 - np.degrees(Elev_9)
	
	grid = plt.GridSpec(1,12,wspace=0.1,hspace=0.1)
	fig = plt.figure()
	ax1 = fig.add_subplot(grid[0,0:10],projection="polar")

	#ax1 = fig.add_subplot(projection="polar")
	ax1.set_ylim(0.0,90.0)

	#ax1.contourf(opt_list[0],opt_list[1],opt_list[2],np.linspace(0.0,1.00,41),extend='neither',cmap="jet",)

	if latitude > 0:
		ax1.plot(Azi_Summer,Zenith_Summer,"r-",linewidth=1)
		ax1.plot(Azi_2,Zenith_2,"r--",linewidth=1)
		ax1.plot(Azi_3,Zenith_3,"r--",linewidth=1)
		ax1.plot(Azi_Spring,Zenith_Spring,"k--",linewidth=1)
		ax1.plot(Azi_Winter,Zenith_Winter,"b-",linewidth=1)
		ax1.plot(Azi_8,Zenith_8,"b--",linewidth=1)
		ax1.plot(Azi_9,Zenith_9,"b--",linewidth=1)
	else:
		ax1.plot(Azi_Summer,Zenith_Summer,"b-",linewidth=1)
		ax1.plot(Azi_2,Zenith_2,"b--",linewidth=1)
		ax1.plot(Azi_3,Zenith_3,"b--",linewidth=1)
		ax1.plot(Azi_Spring,Zenith_Spring,"k--",linewidth=1)
		ax1.plot(Azi_Winter,Zenith_Winter,"r-",linewidth=1)
		ax1.plot(Azi_8,Zenith_8,"r--",linewidth=1)
		ax1.plot(Azi_9,Zenith_9,"r--",linewidth=1)

	ax1.set_theta_offset(np.pi/2.0)
	ax1.set_thetagrids(range(0,360,30))
	plt.show()
	return

filename = "124%phi_100%HT_120%Arecv_optics"
#latitude = 19.27#-35.3
latitude = 37.24
#Polar_Single(filename,latitude)
Sun_Path(latitude)
#Sun_Path2(latitude,declination=46.0)
