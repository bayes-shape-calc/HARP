import logging
import sys
import os
import struct
from . import __version__

from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import Qt,QSharedMemory
from PyQt5.QtGui import QIcon

from .interface import main_window


def setup_logging():
	
	logger = logging.getLogger('harpgui')
	logger.setLevel(logging.DEBUG)
	log_fmt = ("[(%(lineno)4s:%(filename)s %(funcName)s)  %(message)s]")
	formatter = logging.Formatter(log_fmt)
	
	logging.getLogger("matplotlib").setLevel(logging.WARNING)
	logging.getLogger("numba").setLevel(logging.WARNING)
	logging.getLogger("h5py").setLevel(logging.WARNING)

	stdout_handler = logging.StreamHandler()
	stdout_handler.setFormatter(formatter)
	stdout_handler.setLevel(logging.DEBUG)
	logger.addHandler(stdout_handler)

	return logger
	
def excepthook(*exc_args):
	"""
	Log exception and exit cleanly.
	"""
	logging.error("Unrecoverable error", exc_info=(exc_args))
	# Very important to release shared memory used to signal an app instance is running
	# as we are going to exit below
	_shared_memory.release()
	if exc_args[0] != KeyboardInterrupt:
		logging.error("Failed to report crash", exc_info=e)
		sys.__excepthook__(*exc_args)
		sys.exit(1)
	else:  # It's harmless, don't sound the alarm.
		sys.exit(0)

def setup_exception_handler():
	"""
	Install global exception handler
	"""
	sys.excepthook = excepthook

class MutexError(BaseException):
	pass

def check_only_running_once():
	"""If the application is already running log the error and exit"""
	try:
		with _shared_memory:
			_shared_memory.acquire()
	except MutexError as exc:
		[message] = exc.args
		logging.error(message)
		sys.exit(2)
		

class SharedMemoryMutex(object):
	"""Simple wrapper around the QSharedMemory object, adding a context
	handler which uses the built in Semaphore as a locking mechanism
	and raises an error if the shared memory object is already in use

	TODO: The *nix implementation doesn't release the shared memory if a
	process attached to it crashes. There is code to attempt to detect this
	but it doesn't seem worth implementing for the moment: we're only talking
	about a page at most.
	"""

	NAME = "harp_gui_mutex"

	def __init__(self):
		sharedAppName = self.NAME
		self._shared_memory = QSharedMemory(sharedAppName)

	def __enter__(self):
		self._shared_memory.lock()
		return self

	def __exit__(self, *args, **kwargs):
		self._shared_memory.unlock()

	def acquire(self):
		#
		# The attach-detach dance is a shim from
		# https://stackoverflow.com/questions/42549904/qsharedmemory-is-not-getting-deleted-on-application-crash
		# If the existing shared memory is not held by any active application
		# (eg because an appimage has hard-crashed) then it will be released
		# If the memory is held by an active application it will have no effect
		#
		self._shared_memory.attach()
		self._shared_memory.detach()

		if self._shared_memory.attach():
			pid = struct.unpack("q", self._shared_memory.data()[:8])
			raise MutexError("HARP_GUI mutex: HARP_GUI is already running with pid %d" % pid)
		else:
			self._shared_memory.create(8)
			self._shared_memory.data()[:8] = struct.pack("q", os.getpid())

	def release(self):
		self._shared_memory.detach()


_shared_memory = SharedMemoryMutex()



def run():
	logger = setup_logging()
	logger.info("Starting HARP GUI ({})".format(__version__))
	
	try:
		import harp
		logger.info('Load HARP from %s'%(str(harp)))
	except:
		logger.info('Failed to load HARP')
		sys.exit(1)
	
	setup_exception_handler()
	check_only_running_once()

	# Remove the VIRTUAL_ENV env var: when set, prevents pgzero mode from working
	# if HARP_GUI is installed from source/PyPI into a virtual environment.
	os.environ.pop("VIRTUAL_ENV", None)

	# Images (such as toolbar icons) aren't scaled nicely on retina/4k displays
	# unless this flag is set
	os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
	if hasattr(Qt, "AA_EnableHighDpiScaling"):
		QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
	QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps)

	# An issue in PyQt5 v5.13.2 to v5.15.1 makes PyQt5 application
	# hang on Mac OS 11 (Big Sur)
	# Setting this environment variable fixes the problem.
	# See issue #1147 for more information
	os.environ["QT_MAC_WANTS_LAYER"] = "1"

	# The app object is the application running on your computer.
	app = QApplication(sys.argv)
	# By default PyQt uses the script name (run.py)
	app.setApplicationName("HARP")
	# Set hint as to the .desktop files name
	app.setDesktopFileName("harp.gui")
	app.setApplicationVersion(__version__)
	app.setAttribute(Qt.AA_DontShowIconsInMenus)

	# Create the "window" we'll be looking at.
	harp_window = main_window()

	# Make sure all windows have the icon 
	ipath = os.path.dirname(os.path.realpath(__file__))
	app.setWindowIcon(QIcon(ipath + os.path.sep + 'HARP_logo.png'))

	# Save the exit code for sys.exit call below.
	exit_status = app.exec_()
	# Clean up the shared memory used to signal an app instance is running
	_shared_memory.release()

	# Stop the program after the application finishes executing.
	sys.exit(exit_status)