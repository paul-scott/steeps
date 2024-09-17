from __future__ import division, print_function, unicode_literals
import os
from pathlib import Path
import shutil
import warnings
from pipes import quote as sh_quote
import subprocess as sp
import xml.etree.ElementTree as ET
import multiprocessing as mp
import re
import tempfile
import sysconfig

if os.environ.get('ST_DEBUG'):
	import colorama
	colorama.init()
	def bright(str1):
		return colorama.Style.BRIGHT+colorama.Fore.GREEN  + str1 + colorama.Style.RESET_ALL 
	def sp_run(call):
		print(bright(" ".join(call)))
		sp.check_call(call)
else:
	sp_run = sp.check_call	

# TODO: Add in option for different result file output
# TODO: Need to add in error checking for calls (possibly use in tests)

var_re = re.compile(r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)(\S+)')
# Assuming all linear relations between units (factor, offset)
unit_conv = {
	's': {
		's': (1,0),
	},
	'm': {
		's': (60,0),
	},
	'h': {
		's': (60*60,0),
		'm': (60,0),
	},
	'd': {
		's': (24*60*60,0),
		'm': (24*60,0),
		'h': (24,0),
	},
	'y': {
		's': (365*24*60*60,0),
		'm': (365*24*60,0),
		'h': (365*24,0),
		'd': (365,0),
	},
	'W': {
		'W': (1,0),
	},
	'kW': {
		'W': (1000,0),
	},
	'MW': {
		'W': (1e6,0),
		'kW': (1000,0),
	},
	'GW': {
		'W': (1e9,0),
		'kW': (1e6,0),
		'MW': (1000,0),
	},
	'J': {
		'J': (1,0),
	},
	'kWh': {
		'J': (3.6e6,0),
	},
	'MWh': {
		'J': (3.6e9,0),
		'kWh': (1000,0),
	},
	'GWh': {
		'J': (3.6e12,0),
		'MWh': (1000,0),
		'kWh': (1e6,0),
	},
	'$': {
		'$': (1,0),
	},
	'k$': {
		'$': (1000,0),
	},
	'M$': {
		'$': (1e6,0),
		'k$': (1000,0),
	},
}


def convert_val(v1, u1, u2):
	try:
		fac, off = unit_conv[u1][u2]
		return fac*v1 + off
	except KeyError:
		fac, off = unit_conv[u2][u1]
		return (v1 - off)/fac


def move_overwrite(src,dst):
	"""If `src` exists and is a normal file, move it to `dst`, overwiting a file
	that is there already if necessary. Raise exception if `src` is a directory.
	If `dst` is a directory, then place the `src` file in that directory."""
	if not os.access(src,os.R_OK):
		return
	assert not os.path.isdir(src)
	if os.path.isdir(dst):
		dst = os.path.join(dst,os.path.basename(src))
	if os.path.exists(dst):
		assert os.access(dst,os.W_OK) and not os.path.isdir(dst)
		os.unlink(dst)
	# print("Moving '%s' to '%s'"%(src,dst))
	shutil.move(src,dst)


def parse_var_val(vstr, unit):
	"""Convert variable value from string with unit to target unit.

	The string valstr must be a number which has an optional suffix:

		'<number>[unit]'

	Raises an exception if not in the correct format or if there is no known
	conversion between the unit types.
	"""
	# First we assume it has no suffix, and hence already in the target units
	try:
		return float(vstr)
	except ValueError:
		pass

	# Try to get the unit suffix and convert the first part to a number
	res = var_re.match(vstr)
	if res is None:
		raise ValueError('Cannot parse variable value ' + vstr)

	unit_old = res.group(2)  # original unit
	val_old = float(res.group(1))  # original value

	if unit == unit_old:
		return val_old

	try:
		return convert_val(val_old, unit_old, unit)
	except KeyError:
		raise ValueError('Can\'t convert from unit ' + unit_old + ' to ' + unit)


from pathlib import Path
import os, sys

def in_dir_of(destination: Path):
	"""decorator for use in testing functions"""
	def decorator(func):
		def wrapper(*args, **kwargs):
			dest = Path(destination)
			orig_cwd = Path.cwd()
			os.chdir(dest if dest.is_dir() else dest.parent)
			#print("ENVIRON:",os.environ)
			print("SYS.PATH:",sys.path)			
			print("CHECKING IMPORT NUMPY")
			import numpy
			print("...OK")			
			print("ENTERING",Path.cwd())		    
			try:
				return func(*args, **kwargs)
			finally:
				print("LEAVING BACK TO",orig_cwd)
				os.chdir(orig_cwd)
		return wrapper
	return decorator

UNIONFS = "/usr/bin/unionfs-fuse"
FUSERMOUNT = "/bin/fusermount"


class Simulator(object):
	"""Compilation and simulation of a modelica model."""
	def __init__(self, fn, model=None, suffix=None, fusemount=False, reuse=False):
		"""Constructor. `fn` is a .mo filename, `suffix` is a init file suffix,
		`tempdir` is FIXME the location where temporary files should be stored,
		`True` if a machine-generated temporary location should be generated and
		used, or `None` if the current directory should be used.

		The `suffix` is applied to the init file names written and read by the
		simulator.  This can be changed after construction and thereafter
		changes	calls to `write_init`, `update_pars` and `simulate`.

		`fusemount` instructs Simulator to run its operations inside a
		unionfs-fuse filesystem, which enables us to segregate and remove all
		of the temporary files created by OMC. Set False if you want to disable
		this feature.

		`reuse` instructs Simulator to re-use an earlier-created mount. If
		passes, it should be the tuple returned by Simulator::get_fuse_dirs()

		During cleanup, Simulator will keep the init_out_fn and res_fn (if they
		exist) which are assumed to be the only two important output files from
		the calculation.
		"""
		self.fn = os.path.abspath(fn)
		if not os.path.exists(fn):
			raise RuntimeError("Model file '%s' does not exist"%fn)

		if fusemount and not reuse:
			if not os.access(UNIONFS,os.X_OK) or not os.access(FUSERMOUNT,os.X_OK):
				warnings.warn(
					"'%s' or '%s are not executable (try st_simulate --nofuse, or sudo apt install unionfs-fuse)"
						%(UNIONFS,FUSERMOUNT))
				fusemount = 0

		if model is None:
			self.model = os.path.splitext(os.path.split(self.fn)[1])[0]
		else:
			self.model = model

		self.makefile_fn = self.model + '.makefile'
		self.init_in_fn = self.model + '_init.xml'
		self.fusemount = fusemount
		self.enter_fuse(reuse)
		self.suffix = suffix
		self.init_et = None

	def __del__(self):
		if hasattr(self,'fusemount') and self.fusemount:
			self.leave_fuse()

	@property
	def init_out_fn(self):
		if self.suffix is None:
			return self.init_in_fn
		else:
			return self.model + '_init_' + self.suffix + '.xml'

	@property
	def res_fn(self):
		if self.suffix is None:
			return self.model + '_res.mat'
		else:
			return self.model + '_res_' + self.suffix + '.mat'

	def enter_fuse(self,reuse=False):
		assert not getattr(self,'entered_fuse',0)
		if reuse:
			self.mountdir = reuse[0]
			self.tempdir = reuse[1]
			self.init_cwd = reuse[2]
			self.reuse = True
		elif self.fusemount:
			self.mountdir = tempfile.mkdtemp(prefix="solartherm-mount-")
			self.tempdir = tempfile.mkdtemp(prefix="solartherm-temp-")
			self.init_cwd = os.getcwd()
			print("Temporary directory '%s'" % self.tempdir)
			print("Mount directory '%s'" % self.mountdir)
			self.makefile_fn = os.path.join(self.mountdir,self.makefile_fn)
			self.init_in_fn = os.path.join(self.mountdir,self.init_in_fn)
			sp.check_call([UNIONFS,'-o','cow','%s=RW:%s=RO'%(sh_quote(self.tempdir),sh_quote(self.init_cwd)), self.mountdir])
			os.chdir(self.mountdir)
			self.entered_fuse = True
			self.reuse = False

	def leave_fuse(self):		
		os.chdir(self.init_cwd)
		if self.fusemount:
			if not self.reuse:		
				assert getattr(self,'entered_fuse')
				# print("REMOVING FUSE MOUNTPOINT")
				try:
					sp.check_call([FUSERMOUNT,'-uz',self.mountdir])
				except sp.CalledProcessError as e:
					print("UNABLE TO UNMOUNT: %s", str(e))
					print("cwd = ",os.getcwd())
					print("init_cwd =",self.init_cwd)
					import psutil
					proc = psutil.Process()
					print(proc.open_files())
				self.cleanup_fuse()
			self.entered_fuse = False

	def cleanup_fuse(self):
		assert(self.fusemount)
		if self.reuse:
			print("Copying results files to original directory")
			print("result file: ",self.res_fn)
			shutil.move(os.path.join(self.tempdir,self.res_fn),self.init_cwd)
			print("init_out file: ",self.init_out_fn)
			shutil.move(os.path.join(self.tempdir,self.init_out_fn),self.init_cwd)
		else:			
			#print("Cleaning up tempdir '%s'" % self.tempdir)
			move_overwrite(os.path.join(self.tempdir,self.res_fn),self.init_cwd)
			move_overwrite(os.path.join(self.tempdir,self.init_out_fn),self.init_cwd)
			assert(self.tempdir[0:4]=="/tmp") # just being cautious!
			assert(os.path.isdir(self.tempdir))
			shutil.rmtree(self.tempdir)
			os.rmdir(self.mountdir)

	def get_fuse_dirs(self):
		if self.fusemount:
			return (self.mountdir,self.tempdir,self.init_cwd)
		else:
			return None

	def realpath(self,fn):
		"""Return the fuse-mounted version of a local path, if required"""
		if self.fusemount:
			assert self.entered_fuse
			return os.path.join(self.mountdir,fn)
		else:
			return fn

	def compile_model(self, n_proc=0, libs=['Modelica', 'SolarTherm'], args=['-d=nonewInst']):
		"""Compile modelica model in .mo file."""
		call = ['omc', '-s', '-q', '-n='+str(n_proc)] + args + ['-i='+self.model, self.fn] + libs
		if os.environ.get('ST_DEBUG',0):
			call += ['--debug=gendebugsymbols'] # add debug flags
		sp_run(call)
			
		#TODO solve the issue of linker flags in the latest msys2 (v20210228), ASLR enabled by default
		# Ref: 
		# https://www.msys2.org/news/#2021-01-31-aslr-enabled-by-default 
		# https://github.com/msys2/MINGW-packages/issues/7023
		if sysconfig.get_platform()=='mingw':
			makefile=self.model+'.makefile'	
			extraflags = " -Wl,--disable-dynamicbase,--disable-high-entropy-va,--default-image-base-low\n"

			f=open(makefile, "r")
			s=f.readlines()
			f.close()
			i=0
			for r in s:
				if 'LDFLAGS=' in r[:8]:
					s[i]=r[:-1]+extraflags
				i+=1
			f=open(makefile, "w")
			f.writelines(s)
			f.close()		

	def compile_sim(self, n_jobs=(1 + mp.cpu_count()//2), args=[]):
		"""Compile model source code into a simulation executable."""
		call = ['make', '-j', str(n_jobs), '-f', self.makefile_fn] + args
		sp_run(call)

	def load_init(self):
		"""Load in init XML."""
		self.init_et = ET.parse(self.init_in_fn)

	def write_init(self):
		"""Write new init XML with optional suffix."""
		self.init_et.write(self.init_out_fn)

	def clear_init(self):
		"""Clear init XML."""
		self.init_et = None

	def update_et_pars(self, par_n, par_v):
		"""Update parameter values in internal XML.

		The values of parameters will fail to change in the simulation if they
		are final (protected).  This can happen indirectly if a final parameter
		derives its value from a non-final changed parameter.
		"""
		root = self.init_et.getroot()
		for i, n in enumerate(par_n):
			root.find('*ScalarVariable[@name=\''+n+'\']/*[@start]').attrib['start'] = par_v[i]

	def update_pars(self, par_n, par_v):
		"""Update parameters values in XML and write results to file."""
		if self.init_et is None:
			self.load_init()
		self.update_et_pars(par_n, par_v)
		self.write_init()
	
	def get_unit(self, var_n):
		root = self.init_et.getroot()

		node = root.find('*ScalarVariable[@name=\''+var_n+'\']/*[@unit]')
		return '' if node is None else node.attrib['unit']

	def simulate(self, start='0', stop='86400', step='60', tolerance = '1e-04', initStep=None, maxStep=None, integOrder=None, solver='rungekutta', nls='newton', lv='-LOG_SUCCESS,-stdout', noEventEmit=True, args=[]):
		"""Run simulation.

		If running an optimisation then 'optimization' needs to be used as
		solver type.
		"""
		start = str(parse_var_val(start, 's'))
		stop = str(parse_var_val(stop, 's'))
		step = str(parse_var_val(step, 's'))
		tolerance = str(tolerance)
		
		if initStep!=None:
			initStep = str(parse_var_val(initStep, 's'))
		if maxStep!=None:
			maxStep = str(parse_var_val(maxStep, 's'))

		sim_args = [
			'-override',
			'startTime='+start+',stopTime='+stop+',stepSize='+step+',tolerance='+tolerance,
			'-s', solver,
			'-nls', nls, #Nonlinear solver
			'-initialStepSize', initStep,
			'-maxStepSize', maxStep,
			'-maxIntegrationOrder', integOrder,
			'-lv', lv, # Specifies which logging levels to enable
			'-f', self.init_out_fn,
			'-r', self.res_fn,
			]
		
		if noEventEmit:
			sim_args.append('-noEventEmit',)	

		if initStep==None:
			sim_args = [e for e in sim_args if e not in ('-initialStepSize', initStep)]
		if maxStep==None:
			sim_args = [e for e in sim_args if e not in ('-maxStepSize', maxStep)]
		if integOrder==None:
			sim_args = [e for e in sim_args if e not in ('-maxIntegrationOrder', integOrder)]
		if lv==None:
			sim_args = [e for e in sim_args if e not in ('-lv', lv)]

		call = ['./'+self.model] + sim_args + args
		sp_run(call)
		#sp.call(['./'+self.model] + sim_args + args)
		# assert also that there must be a result file
		assert os.access(self.res_fn,os.R_OK)

# vim: ts=4:sw=4:noet:tw=80
