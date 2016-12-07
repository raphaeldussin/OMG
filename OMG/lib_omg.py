import OMG.utils as lu
import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.cm as cm
from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                       AltitudeMode, Camera)

class gearth_anim():

	def __init__(self,gridfile,plotdir,model='ROMS'):
		self.gridfile = gridfile
		self.reference = dt.datetime(1900,1,1,0,0)
		self.maxdigit=6
		self.dpi = 600
		self.figure_size = 8 # *100 pixels
		self.plotdir = plotdir
		self.model = model
		if self.model == 'ROMS':
			self.read_coords()
			self.timevar = 'ocean_time'
		elif self.model == 'NEMO':
			self.read_coords(lonvar='nav_lon',latvar='nav_lat')
			self.timevar = 'time_counter'
			self.spval = 0.

	def create_animation(self,listfiles,output,showing='velocity_module'):
		times = [] ; frame_list = [] ; ct=1
		for file in listfiles:
			print 'working on file ', file
			# define png output
			current_png = self.plotdir + 'output_' + str(ct).zfill(self.maxdigit) + '.png'
			frame_list.append(current_png)
			# read time
			try:
				times.append(self.read_time(file))
			except:
				times.append(self.reference)
			# compute data in subroutine
			if showing == 'velocity_module':
				data = self.velocity_module(file)
			elif showing == 'sst':
				data = self.sst(file)
			elif showing == 'bathy_nemo':
				data = self.bathy_nemo(file)
			elif showing == 'bathy_roms':
				data = self.bathy_roms(file)
			else:
				exit('no such kind of animation')

			if self.model == 'NEMO':
				# we need to sort longitude to be monotonic
				lonmono, latmono, datamono = lu.sort_monotonic(self.lon, self.lat, data)
				# remove the northfold
				lonplt = lonmono[:-1,:] ; latplt = latmono[:-1,:] ; dataplt = datamono[:-1,:]
				dataplt = np.ma.masked_values(dataplt,self.spval)	
			else:
				 lonplt = self.lon ; latplt = self.lat ; dataplt = data

			# make the plot
			fig, ax = self.gearth_fig()
			cs = ax.pcolormesh(lonplt, latplt, dataplt, cmap=self.colormap,vmin=self.vmin,vmax=self.vmax)
			ax.set_axis_off()
			fig.savefig(current_png, transparent=True, format='png')
			plt.close(fig)
			ct = ct + 1
		# create final kmz file
		self.make_kml(times,frame_list,output)
		return None

	def velocity_module(self,file):
		var_u = self.read_data(file,'u')
		var_v = self.read_data(file,'v')
		u = 0.5 * (var_u[1:,:] + var_u[:-1,:])
		v = 0.5 * (var_v[:,1:] + var_v[:,:-1])
		data = np.sqrt(u*u+v*v)
		# plot parameters
		self.vmin = 0. ; self.vmax = 2.
		self.colormap = cm.Blues_r
		return data

	def sst(self,file):
		data = self.read_data(file,'temp')
		# plot parameters
		self.vmin = -2. ; self.vmax = 35.
		self.colormap = cm.gist_ncar
		return data

	def bathy_nemo(self,file):
		data = self.read_data(file,'hdept',level=None,time=None)
		# plot parameters
		self.vmin = 0. ; self.vmax = 5800.
		self.colormap = cm.jet
		self.spval = data.min() # ugly but works
		return data

	def bathy_roms(self,file):
		data = self.read_data(file,'h',level=None,time=None)
		mask = self.read_data(self.gridfile,'mask_rho',level=None,time=None)
		logicalmask = np.empty(mask.shape)
		logicalmask[:] = False
		logicalmask[np.where(mask == 0)] = True
		data_mask = ma.masked_array(data,mask=logicalmask)
		# plot parameters
		self.vmin = 0. ; self.vmax = 1000.
		self.colormap = cm.jet
		return data_mask

	def read_coords(self,lonvar='lon_rho',latvar='lat_rho'):
		nc_grd = nc.Dataset(self.gridfile,'r')
		self.lat = nc_grd.variables[latvar][:]
		self.lon = nc_grd.variables[lonvar][:]
		nc_grd.close()	
		self.llcrnrlon=self.lon.min()
		self.llcrnrlat=self.lat.min()
		self.urcrnrlon=self.lon.max()
		self.urcrnrlat=self.lat.max()
		return None

	def read_data(self,datafile,datavar,level=-1,time=0):
		nc_data = nc.Dataset(datafile,'r')
		if time is None:
			if level is None:
				data = nc_data.variables[datavar][:].squeeze()
			else:
				data = nc_data.variables[datavar][level,:].squeeze()
		else:
			if level is None:
				data = nc_data.variables[datavar][time,:,:].squeeze()
			else:
				data = nc_data.variables[datavar][time,level,:,:].squeeze()
		nc_data.close()
		return data

	def read_time(self,datafile,ref=None):
		nc_data = nc.Dataset(datafile,'r')
		seconds_from_ref = nc_data.variables[self.timevar][0]
		nc_data.close()
		if ref is None:
			ref = self.reference
		time = ref + dt.timedelta(seconds=seconds_from_ref)
		return time

	def make_kml(self,times,figs,fileout,colorbar=None,debug=False,**kw): 
		"""TODO: LatLon bbox, list of figs, optional colorbar figure,
		and several simplekml kw..."""

		kml = Kml()
		altitude = kw.pop('altitude', 2e7)
		roll = kw.pop('roll', 0)
		tilt = kw.pop('tilt', 0)
		altitudemode = kw.pop('altitudemode', AltitudeMode.relativetoground)
		camera = Camera(latitude=np.mean([self.urcrnrlat, self.llcrnrlat]),
		                longitude=np.mean([self.urcrnrlon, self.llcrnrlon]),
		                altitude=altitude, roll=roll, tilt=tilt,
		                altitudemode=altitudemode)

		# we need another date to close last interval
		dt = times[1] - times[0]
		next_time = times[-1] + dt
		times.append(next_time)

		kml.document.camera = camera
		draworder = 0
		for fig in figs:  # NOTE: Overlays are limited to the same bbox.
			draworder += 1
			ground = kml.newgroundoverlay(name='GroundOverlay')
			ground.draworder = draworder
			ground.visibility = kw.pop('visibility', 1)
			ground.name = kw.pop('name', 'overlay')
			ground.color = kw.pop('color', '9effffff')
			ground.atomauthor = kw.pop('author', 'esm')
			ground.latlonbox.rotation = kw.pop('rotation', 0)
			ground.description = kw.pop('description', 'Matplotlib figure')
			ground.gxaltitudemode = kw.pop('gxaltitudemode',
			                               'clampToSeaFloor')
			ground.icon.href = fig
			ground.latlonbox.east = self.llcrnrlon
			ground.latlonbox.south = self.llcrnrlat
			ground.latlonbox.north = self.urcrnrlat
			ground.latlonbox.west = self.urcrnrlon
			# date span
			ground.timespan.begin = times[draworder-1].strftime(format="%Y-%m-%d")
			ground.timespan.end = times[draworder].strftime(format="%Y-%m-%d")
	
		kmzfile = kw.pop('kmzfile', self.plotdir + fileout)
		kml.savekmz(kmzfile)
		return None
	
	def gearth_fig(self):
		"""Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image."""
		aspect = np.cos(np.mean([self.llcrnrlat, self.urcrnrlat]) * np.pi/180.0)
		xsize = np.ptp([self.urcrnrlon, self.llcrnrlon]) * aspect
		ysize = np.ptp([self.urcrnrlat, self.llcrnrlat])
		aspect = ysize / xsize
	
		if aspect > 1.0:
			figsize = (self.figure_size / aspect, self.figure_size)
		else:
			figsize = (self.figure_size, self.figure_size * aspect)
	
		fig = plt.figure(figsize=figsize,frameon=False,dpi=self.dpi)
		ax = fig.add_axes([0, 0, 1, 1])
		ax.set_xlim(self.llcrnrlon, self.urcrnrlon)
		ax.set_ylim(self.llcrnrlat, self.urcrnrlat)
		return fig, ax
	
