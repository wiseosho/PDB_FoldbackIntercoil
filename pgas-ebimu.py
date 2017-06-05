# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/home/june/pgas-ebimu-draft.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
import serial, glob
from time import sleep
import numpy as np
import pandas as pd
from datetime import datetime

IMU_raw = np.empty([1,24]) # Declare as a Global Const
Ser = serial.Serial()
FILESTAMP = ''
RECORD = False
TIMESTAMP = ''

class EBIMU_thread(QtCore.QThread):
    update = QtCore.pyqtSignal(list)
    def __init__(self):
        QtCore.QThread.__init__(self)
    
    def run(self):        
        while Ser.isOpen():
            packet = []
            for i in range(8):
                packet.append(Ser.readline().decode('utf-8'))
            self.update.emit(packet)
            
class Ui_PGAS_demo(object):
    def setupUi(self, PGAS_demo):
        PGAS_demo.setObjectName("PGAS_demo")
        PGAS_demo.resize(672, 499)
        self.centralwidget = QtWidgets.QWidget(PGAS_demo)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(50, 100, 551, 291))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setSizeConstraint(QtWidgets.QLayout.SetMaximumSize)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_15 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_15.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_15.sizePolicy().hasHeightForWidth())
        self.label_15.setSizePolicy(sizePolicy)
        self.label_15.setObjectName("label_15")
        self.verticalLayout.addWidget(self.label_15)
        self.label_5 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_5.setObjectName("label_5")
        self.verticalLayout.addWidget(self.label_5)
        self.label_6 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_6.setObjectName("label_6")
        self.verticalLayout.addWidget(self.label_6)
        self.label_7 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_7.setObjectName("label_7")
        self.verticalLayout.addWidget(self.label_7)
        self.label_8 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_8.setObjectName("label_8")
        self.verticalLayout.addWidget(self.label_8)
        self.label_9 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_9.setObjectName("label_9")
        self.verticalLayout.addWidget(self.label_9)
        self.label_10 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_10.setObjectName("label_10")
        self.verticalLayout.addWidget(self.label_10)
        self.label_11 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_11.setObjectName("label_11")
        self.verticalLayout.addWidget(self.label_11)
        self.label_12 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_12.setObjectName("label_12")
        self.verticalLayout.addWidget(self.label_12)
        self.horizontalLayout_2.addLayout(self.verticalLayout)
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.imu2y = QtWidgets.QTextBrowser(self.gridLayoutWidget)
        self.imu2y.setObjectName("imu2y")
        self.gridLayout_2.addWidget(self.imu2y, 4, 0, 1, 1)
        self.label_13 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_13.setObjectName("label_13")
        self.gridLayout_2.addWidget(self.label_13, 0, 0, 1, 1)
        self.imu5y = QtWidgets.QTextBrowser(self.gridLayoutWidget)
        self.imu5y.setObjectName("imu5y")
        self.gridLayout_2.addWidget(self.imu5y, 7, 0, 1, 1)
        
        self.imu1 = []
        self.imu = []
        for j in range(8):
            self.imu.append([])    
            for i in range (3):
                self.imu[j].append(QtWidgets.QTextBrowser(self.gridLayoutWidget))
                self.imu[j][i].setEnabled(True)
                self.imu[j][i].setObjectName("imu"+str(j)+str(i))
                self.gridLayout_2.addWidget(self.imu[j][i], 3+j, i, 1, 1)
        

        self.label_16 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_16.setObjectName("label_16")
        self.gridLayout_2.addWidget(self.label_16, 0, 2, 1, 1)

        self.label_14 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_14.setObjectName("label_14")
        self.gridLayout_2.addWidget(self.label_14, 0, 1, 1, 1)
        self.horizontalLayout_2.addLayout(self.gridLayout_2)
        self.gridLayout.addLayout(self.horizontalLayout_2, 0, 0, 1, 1)
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(70, 30, 96, 27))
        self.pushButton.setObjectName("pushButton")
        self.comboBox = QtWidgets.QComboBox(self.centralwidget)
        self.comboBox.setGeometry(QtCore.QRect(180, 30, 85, 27))
        self.comboBox.setObjectName("comboBox")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(280, 30, 151, 31))
        self.label.setObjectName("label")
        PGAS_demo.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(PGAS_demo)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 672, 25))
        self.menubar.setObjectName("menubar")
        PGAS_demo.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(PGAS_demo)
        self.statusbar.setObjectName("statusbar")
        PGAS_demo.setStatusBar(self.statusbar)
        self.Record = QtWidgets.QPushButton(self.centralwidget)
        self.Record.setGeometry(QtCore.QRect(70, 60, 96, 27))
        self.Record.setObjectName("Record")
        self.Record.setCheckable(True)
        self.Analyse = QtWidgets.QPushButton(self.centralwidget)
        self.Analyse.setGeometry(QtCore.QRect(330, 60, 96, 27))
        self.Analyse.setObjectName("Analyse")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(180, 63, 141, 21))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(440, 63, 141, 21))
        self.label_3.setObjectName("label_3")
        
        self.retranslateUi(PGAS_demo)
        
        self.pushButton.clicked.connect(self.Initialize_Port)
        self.Record.clicked.connect(self.recordbuttonclicked)
        #self.pushButton.clicked.connect(self.UpdateImus)
        self.thread = EBIMU_thread()
        self.thread.start()
        
        self.thread.update.connect(self.testdisplay)
        #self.UpdateImus()

        QtCore.QMetaObject.connectSlotsByName(PGAS_demo)
        #self.imu7p.setText('Null')
    def retranslateUi(self, PGAS_demo):
        _translate = QtCore.QCoreApplication.translate
        PGAS_demo.setWindowTitle(_translate("PGAS_demo", "MainWindow"))
        self.label_15.setText(_translate("PGAS_demo", "NoName"))
        self.label_5.setText(_translate("PGAS_demo", "IMU1"))
        self.label_6.setText(_translate("PGAS_demo", "IMU2"))
        self.label_7.setText(_translate("PGAS_demo", "IMU3"))
        self.label_8.setText(_translate("PGAS_demo", "IMU4"))
        self.label_9.setText(_translate("PGAS_demo", "IMU5"))
        self.label_10.setText(_translate("PGAS_demo", "IMU6"))
        self.label_11.setText(_translate("PGAS_demo", "IMU7"))
        self.label_12.setText(_translate("PGAS_demo", "IMU8"))
        self.label_13.setText(_translate("PGAS_demo", "YAW"))
        self.label_16.setText(_translate("PGAS_demo", "Roll"))
        self.label_14.setText(_translate("PGAS_demo", "Pitch"))
        self.pushButton.setText(_translate("PGAS_demo", "Connect"))
        self.Record.setText(_translate("PGAS_demo", "Record"))
        self.Analyse.setText(_translate("PGAS_demo", "Analyse"))
        self.label_2.setText(_translate("PGAS_demo", "RecName"))
        self.label_3.setText(_translate("PGAS_demo", "TextLabel"))

        ports = self.Detect_serial_ports()
        for idx,port in enumerate(ports):
            self.comboBox.addItem("")
            self.comboBox.setItemText(idx, _translate("PGAS_demo", port))
        self.label.setText(_translate("PGAS_demo", "Closed"))
    
    def recordbuttonclicked(self):
        global IMU_raw
        global FILESTAMP
        global TIMESTAMP
        if self.Record.isChecked() == True:
            self.Record.setText("Stop")
            TIMESTAMP = datetime.now().isoformat()
            FILESTAMP = TIMESTAMP.split('.')[0].replace(':', '-')
            
            IMU_raw = np.empty([1, 24])
            with open("./cache"+FILESTAMP+".csv", 'w') as cache:
                cache.write(FILESTAMP + '\n')
                cache.close()
        elif self.Record.isChecked() == False:
            self.Record.setText("Record")
            self.savetofile()
            
    def testdisplay(self, packet):
        global IMU_raw
        Imus = np.empty([8,3])
        for line in packet:
            
            elements = line.split(',')
            Nimu = int(elements[0][-1])-1
            Imus[Nimu] = np.array(elements[1:4])
                        
            for i in range(3):
                try:
                    self.imu[Nimu][i].setText(str(elements[i+1]))
                except Exception as ex:
                    print (ex.args)
                    continue
        IMU_raw = np.vstack((IMU_raw, np.append(Imus[0], Imus[1:])))
        
    def savetofile(self):
        global FILESTAMP
        global TIMESTAMP
        global IMU_raw
        #This operation occurs once Save/ Analysis is called.
        #Convert np.array into the pandas DataFrame and save into the CSV file.
        arrays = [['imu1','imu2', 'imu3', 'imu4', 'imu5', 'imu6', 'imu7', 'imu8' ],['roll','pitch','yaw'] ]
        columns = pd.MultiIndex.from_product(arrays, names=['IMUs','Coord'])
        ENDSTAMP = datetime.now().isoformat()
        Timedelta = (datetime.strptime(ENDSTAMP, "%Y-%m-%dT%H:%M:%S.%f")
        - datetime.strptime(TIMESTAMP, "%Y-%m-%dT%H:%M:%S.%f")).total_seconds()
        with open('./'+FILESTAMP+'.csv', 'wb') as arr:
            np.savetxt(arr, IMU_raw, delimiter=',', header = TIMESTAMP, footer=ENDSTAMP)
            
        #IMU_raw = np.loadtxt('./'+FILESTAMP + '.csv', delimiter=',', comments = '#')
        IMUdf = pd.DataFrame(IMU_raw, columns=columns)
        print (Timedelta)
        #IMUdf.to_csv('./'+FILESTAMP+'.csv')

        #IMUdf.columns = pd.MultiIndex.from_tuples(IMUdf.columns)
        #IMUdf = pd.read_csv('./imu_raw.csv', header=[0,1], skipinitialspace=False, tupleize_cols= True)
        #IMUdf = pd.read_csv('./'+FILESTAMP+'.csv', index_col=False, names=columns, skiprows=[0,1])
        #IMUdf.columns = pd.MultiIndex.from_tuples(IMUdf.columns)
  
    def Detect_serial_ports(self):
        """ Lists serial port names
    
            :raises EnvironmentError:
                On unsupported or unknown platforms
            :returns:
                A list of the serial ports available on the system
        """
        if sys.platform.startswith('win'):
            ports = ['COM%s' % (i + 1) for i in range(256)]
        elif sys.platform.startswith('linux') or sys.platform.startswith('cygwin'):
            # this excludes your current terminal "/dev/tty"
            ports = glob.glob('/dev/tty[A-Za-z]*')
        elif sys.platform.startswith('darwin'):
            ports = glob.glob('/dev/tty.*')
        else:
            raise EnvironmentError('Unsupported platform')
        
        result = []
        for port in ports:
            try:
                s = serial.Serial(port)
                s.close()
                result.append(port)
            except (OSError, serial.SerialException):
                pass
        return result
    
    def Initialize_Port(self):
        ser = Ser#ser = serial.Serial()
        ser.baudrate = 115200
        port = self.comboBox.currentText()
        ser.port = port
        if self.pushButton.text() == "Disconnect":
            ser.close()
            if ser.isOpen() == False:
                self.label.setText("Closed")
                self.comboBox.setEnabled(True)
                self.pushButton.setText("Connect")
                #EBIMU_thread.stop()
                return ser
            else:
                self.label.setText("Not closed, Check port.")
                return ser
        
        if (self.pushButton.text() == "Connect" and ser.isOpen()) :
            self.label.setText("port already Opened, closing..")
            ser.close()
        else:
            ser.open()
            if ser.isOpen() == True:
                self.label.setText(port + " Opened")
                self.comboBox.setEnabled(False)
                self.pushButton.setText("Disconnect")
                self.thread.start()
            
            else:
                self.label.setText(port + " Cannot be open")
        return ser

    def UpdateImus(self):
        for i in range(100):
            sleep(1)
            self.imu1p.setText(str(i))
        return 0
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    PGAS_demo = QtWidgets.QMainWindow()
    ui = Ui_PGAS_demo()
    ui.setupUi(PGAS_demo)
    PGAS_demo.show()
    sys.exit(app.exec_())
    

