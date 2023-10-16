import logging
logger = logging.getLogger('harpgui')
import sys
import os

from PyQt5.QtCore import QSize, Qt
from PyQt5.QtWidgets import QMainWindow,QLabel,QVBoxLayout,QHBoxLayout,QWidget,QSizePolicy,QPushButton,QFileDialog,QLineEdit,QTabWidget,QCheckBox,QDoubleSpinBox,QSpinBox,QFormLayout
from PyQt5.QtGui import QIcon

from . import __version__
from .. import bin,bayes_model_select

class main_window(QMainWindow):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.initialize_widgets()
		
		self.setWindowTitle('HARP (%s)'%(str(__version__)))
		self.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)

		self.show()
	
	def initialize_widgets(self):
		tabs = QTabWidget()
		tabs.setStyleSheet('''QTabWidget::tab-bar {left: 0;}''')
		
		
		## Tab one
		pdblabel = QLabel('PDB ID:')
		self.le_pdbid = QLineEdit()
		self.b_runpdbid = QPushButton('Run HARP')
		self.b_openpdbid = QPushButton('Open on RCSB')
		
		self.le_pdbid.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
		self.b_runpdbid.clicked.connect(self.run_pdbid)
		self.b_openpdbid.clicked.connect(self.open_pdbid)
		
		qw = QWidget()
		hbox = QHBoxLayout()
		hbox.addWidget(pdblabel)
		hbox.addWidget(self.le_pdbid)
		hbox.addWidget(self.b_runpdbid)
		hbox.addWidget(self.b_openpdbid)
		qw.setLayout(hbox)
		
		qwtab = QWidget()
		vbox = QVBoxLayout()
		vbox.addWidget(qw)
		vbox.addStretch()
		qwtab.setLayout(vbox)
		tabs.insertTab(0,qwtab,'PDB ID')




		## Tab two

		l_map = QLabel('Map:')
		self.le_map = QLineEdit()
		self.le_map.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
		self.b_selectmap = QPushButton('Select File')
		qmap = QWidget()
		hbox = QHBoxLayout()
		hbox.addWidget(l_map)
		hbox.addWidget(self.le_map)
		hbox.addWidget(self.b_selectmap)
		qmap.setLayout(hbox)
		
		l_model = QLabel('Model:')
		self.le_model = QLineEdit()
		self.le_model.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
		self.b_selectmodel = QPushButton('Select File')
		qmodel = QWidget()
		hbox = QHBoxLayout()
		hbox.addWidget(l_model)
		hbox.addWidget(self.le_model)
		hbox.addWidget(self.b_selectmodel)
		qmodel.setLayout(hbox)
		
		qwr = QWidget()
		hbox = QHBoxLayout()
		hbox.addStretch()
		self.b_runfile = QPushButton('Run HARP')
		hbox.addWidget(self.b_runfile)
		qwr.setLayout(hbox)
		
		qw = QWidget()
		vbox = QVBoxLayout()
		vbox.addWidget(qmodel)
		vbox.addWidget(qmap)
		vbox.addWidget(qwr)
		vbox.addStretch()
		qw.setLayout(vbox)
		tabs.insertTab(1,qw,'Files')

		self.b_runfile.clicked.connect(self.run_file)
		self.b_selectmap.clicked.connect(self.select_map)
		self.b_selectmodel.clicked.connect(self.select_model)
		
		
		
		## Tab three
		
		l_wdir = QLabel('Working Directory:')
		self.le_wdir = QLineEdit()
		self.le_wdir.setText(os.getcwd())
		self.le_wdir.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
		self.b_selectwdir = QPushButton('Select Directory')
		qwdir = QWidget()
		hbox = QHBoxLayout()
		hbox.addWidget(l_wdir)
		hbox.addWidget(self.le_wdir)
		hbox.addWidget(self.b_selectwdir)
		qwdir.setLayout(hbox)
		
		qfl = QFormLayout()
		
		self.sb_voxels = QDoubleSpinBox()
		self.sb_voxels.setMinimum(0)
		self.sb_voxels.setMaximum(1.)
		self.sb_voxels.setValue(0.5)
		self.sb_voxels.setDecimals(2)

		self.sb_nsig0 = QSpinBox()
		self.sb_nsig0.setMinimum(1)
		self.sb_nsig0.setSingleStep(1)
		self.sb_nsig0.setValue(10)

		self.sb_sig0min = QDoubleSpinBox()
		self.sb_sig0min.setMinimum(0.000001)
		self.sb_sig0min.setMaximum(100.)
		self.sb_sig0min.setValue(0.25)
		self.sb_sig0min.setDecimals(4)

		self.sb_sig0max = QDoubleSpinBox()
		self.sb_sig0max.setMinimum(0.000001)
		self.sb_sig0max.setMaximum(100.)
		self.sb_sig0max.setValue(1.00)
		self.sb_sig0max.setDecimals(4)

		self.sb_nsig1 = QSpinBox()
		self.sb_nsig1.setMinimum(1)
		self.sb_nsig1.setSingleStep(1)
		self.sb_nsig1.setValue(20)

		self.sb_sig1min = QDoubleSpinBox()
		self.sb_sig1min.setMinimum(0.000001)
		self.sb_sig1min.setMaximum(100.)
		self.sb_sig1min.setValue(0.25)
		self.sb_sig1min.setDecimals(4)
		
		self.sb_sig1max = QDoubleSpinBox()
		self.sb_sig1max.setMinimum(0.000001)
		self.sb_sig1max.setMaximum(100.)
		self.sb_sig1max.setValue(2.8)
		self.sb_sig1max.setDecimals(4)
		

		
		self.cb_authid = QCheckBox('Use Auth ID instead of Label ID?')
		self.cb_onlypolymers = QCheckBox('Only calculate chains that come from entities that are polymers?')
		self.cb_removemetals = QCheckBox('Remove common metal/salt ions from molecule before calculation?')
		self.cb_removewater = QCheckBox('Remove waters from molecule before calculation?')
		self.cb_removemetals.setChecked(True)
		self.cb_removewater.setChecked(True)
		
		qfl.addRow(self.cb_authid)
		qfl.addRow(self.cb_onlypolymers)
		qfl.addRow(self.cb_removemetals)
		qfl.addRow(self.cb_removewater)
		qfl.addRow('Voxel Offset',self.sb_voxels)
		qfl.addRow('Number of sigma_0 points',self.sb_nsig0)
		qfl.addRow('Minimum sigma_0 value',self.sb_sig0min)
		qfl.addRow('Maximum sigma_0 value',self.sb_sig0max)
		qfl.addRow('Number of sigma_1 points',self.sb_nsig1)
		qfl.addRow('Minimum sigma_1 value',self.sb_sig1min)
		qfl.addRow('Maximum sigma_1 value',self.sb_sig1max)
		
		qfw = QWidget()
		qfw.setLayout(qfl)
		
		
		qw = QWidget()
		vbox = QVBoxLayout()
		vbox.addWidget(qwdir)

		vbox.addWidget(qfw)
		
		qw.setLayout(vbox)
		tabs.insertTab(2,qw,'Options')
		

		self.b_selectwdir.clicked.connect(self.select_wdir)

		tabs.setCurrentIndex(0)
		self.setCentralWidget(tabs)
		
		
	def select_map(self):
		logger.info('select map')

		fname = QFileDialog.getOpenFileName(self,"Select map file (.mrc)")[0]
		if not fname == "":
			self.le_map.setText(fname)

	def select_model(self):
		logger.info('select model')
		
		fname = QFileDialog.getOpenFileName(self,"Select model file (.cif,.cif.gz)")[0]
		if not fname == "":
			self.le_model.setText(fname)
			
	def select_wdir(self):
		logger.info('select working dir')
		
		fname = QFileDialog.getExistingDirectory(self,"Select working directory")
		if not fname == "":
			self.le_wdir.setText(fname)


	def run_file(self):
		logger.info('Run File')
		
		fname_map = self.le_map.text()
		fname_model = self.le_model.text()
		
		if os.path.isfile(fname_map) and os.path.isfile(fname_model):
			basedir = self.le_wdir.text()
			if not os.path.isdir(basedir):
				basedir='./'

			adfs,blobs = bayes_model_select.gen_adfsblobs(
				self.sb_sig0min.value(),self.sb_sig0max.value(),self.sb_nsig0.value(),
				self.sb_sig1min.value(),self.sb_sig1max.value(),self.sb_nsig1.value()
			)

			bin.harpcalc(
				None,
				basedir,
				input_files=[fname_model,fname_map],
				adfs = adfs,
				blobs = blobs,
				offset=self.sb_voxels.value(),
				authid=self.cb_authid.isChecked(),
				only_polymers=self.cb_onlypolymers.isChecked(),
				remove_metals=self.cb_removemetals.isChecked(),
				remove_water=self.cb_removewater.isChecked(),
			)
		
	
	def run_pdbid(self):
		pdbid = self.le_pdbid.text()
		logger.info('Running HARP on %s'%(pdbid))
		
		basedir = self.le_wdir.text()
		if not os.path.isdir(basedir):
			basedir='./'
			
		adfs,blobs = bayes_model_select.gen_adfsblobs(
			self.sb_sig0min.value(),self.sb_sig0max.value(),self.sb_nsig0.value(),
			self.sb_sig1min.value(),self.sb_sig1max.value(),self.sb_nsig1.value()
		)
		
		bin.harpcalc(
			pdbid,
			basedir,
			adfs = adfs,
			blobs = blobs,
			offset=self.sb_voxels.value(),
			authid=self.cb_authid.isChecked(),
			only_polymers=self.cb_onlypolymers.isChecked(),
			remove_metals=self.cb_removemetals.isChecked(),
			remove_water=self.cb_removewater.isChecked(),
		)
		
		
	def open_pdbid(self):
		pdbid = self.le_pdbid.text()
		logger.info('Opening %s (PDB ID) on RCSB'%(pdbid))
		bin.openrcsb(pdbid)
		
		
		
		
		