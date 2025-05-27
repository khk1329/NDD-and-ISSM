# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'ISSMrCnMmD.ui'
##
## Created by: Qt User Interface Compiler version 6.9.0
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QCheckBox, QDoubleSpinBox,
    QFrame, QLabel, QLineEdit, QListWidget,
    QListWidgetItem, QMainWindow, QPushButton, QRadioButton,
    QSizePolicy, QTextEdit, QWidget)
from ui import resource_rc

class Ui_In_silico_sequence_mining(object):
    def setupUi(self, In_silico_sequence_mining):
        if not In_silico_sequence_mining.objectName():
            In_silico_sequence_mining.setObjectName(u"In_silico_sequence_mining")
        In_silico_sequence_mining.resize(1200, 601)
        In_silico_sequence_mining.setCursor(QCursor(Qt.CursorShape.ArrowCursor))
        In_silico_sequence_mining.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        In_silico_sequence_mining.setStyleSheet(u"")
        self.centralwidget = QWidget(In_silico_sequence_mining)
        self.centralwidget.setObjectName(u"centralwidget")
        self.centralwidget.setStyleSheet(u"QWidget{background-color: #F9FAFA; color: black;}\n"
"     QPushButton {\n"
"                background-color: #1D83D5;\n"
"                color: white;\n"
"                font-weight: bold;\n"
"                border: 1px solid #4183BC;\n"
"                border-radius: 6px;\n"
"            }\n"
"            QPushButton:hover {\n"
"                background-color: #306999;\n"
"            }")
        self.radioButton = QRadioButton(self.centralwidget)
        self.radioButton.setObjectName(u"radioButton")
        self.radioButton.setGeometry(QRect(150, 20, 61, 21))
        font = QFont()
        font.setBold(True)
        self.radioButton.setFont(font)
        self.radioButton.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.radioButton.setToolTipDuration(-1)
        self.radioButton.setAutoFillBackground(False)
        self.radioButton.setStyleSheet(u"QRadioButton::indicator {\n"
"    width: 16px;\n"
"    height: 16px;\n"
"}\n"
"QRadioButton::indicator:checked{image:url(:/icons/icons/circle-check.png)}\n"
"QRadioButton::indicator:unchecked{image:url(:/icons/icons/circle-dashed.png)}")
        self.radioButton.setIconSize(QSize(15, 15))
        self.radioButton_2 = QRadioButton(self.centralwidget)
        self.radioButton_2.setObjectName(u"radioButton_2")
        self.radioButton_2.setGeometry(QRect(230, 20, 61, 21))
        self.radioButton_2.setFont(font)
        self.radioButton_2.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.radioButton_2.setStyleSheet(u"QRadioButton::indicator {\n"
"    width: 16px;\n"
"    height: 16px;\n"
"}\n"
"QRadioButton::indicator:checked{image:url(:/icons/icons/circle-check.png)}\n"
"QRadioButton::indicator:unchecked{image:url(:/icons/icons/circle-dashed.png)}")
        self.label = QLabel(self.centralwidget)
        self.label.setObjectName(u"label")
        self.label.setGeometry(QRect(30, 20, 101, 20))
        self.selectFolderButton = QPushButton(self.centralwidget)
        self.selectFolderButton.setObjectName(u"selectFolderButton")
        self.selectFolderButton.setGeometry(QRect(30, 50, 111, 24))
        self.selectFolderButton.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.inputPathLineEdit = QLineEdit(self.centralwidget)
        self.inputPathLineEdit.setObjectName(u"inputPathLineEdit")
        self.inputPathLineEdit.setGeometry(QRect(150, 50, 401, 22))
        self.inputPathLineEdit.setStyleSheet(u"QLineEdit{background-color: white;}")
        self.allFilesList = QListWidget(self.centralwidget)
        self.allFilesList.setObjectName(u"allFilesList")
        self.allFilesList.setGeometry(QRect(30, 110, 231, 441))
        self.allFilesList.setStyleSheet(u"QListWidget {\n"
"                background-color: white;\n"
"                color: #2c3e50;\n"
"                border: 1px solid #ccc;\n"
"            }\n"
"QListWidget::item:selected {\n"
"                background-color: #1663AE;\n"
"                color: white;\n"
"            }\n"
"QListWidget::item:hover {\n"
"                background-color: #EBF5FD;\n"
"                color : black;\n"
"            }\n"
"QScrollBar:vertical {\n"
"        background: transparent;\n"
"        width: 6px;\n"
"        margin: 2px 0 2px 0;\n"
"    }\n"
"\n"
"QScrollBar::handle:vertical {\n"
"        background: #1D83D5;\n"
"        border-radius: 3px;\n"
"        min-height: 20px;\n"
"    }\n"
"\n"
"QScrollBar::add-line:vertical,\n"
"QScrollBar::sub-line:vertical {\n"
"        height: 0px;\n"
"        background: none;\n"
"    }\n"
"\n"
"QScrollBar::add-page:vertical,\n"
"QScrollBar::sub-page:vertical {\n"
"        background: none;\n"
"    }\n"
"\n"
"QScrollBar:horizontal {\n"
"    background: transparent;\n"
"    h"
                        "eight: 6px;\n"
"    margin: 0 2px 0 2px;\n"
"}\n"
"\n"
"QScrollBar::handle:horizontal {\n"
"    background: #1D83D5;\n"
"    border-radius: 3px;\n"
"    min-width: 20px;\n"
"}\n"
"\n"
"QScrollBar::add-line:horizontal,\n"
"QScrollBar::sub-line:horizontal {\n"
"    width: 0px;\n"
"    background: none;\n"
"}\n"
"\n"
"QScrollBar::add-page:horizontal,\n"
"QScrollBar::sub-page:horizontal {\n"
"    background: none;\n"
"}\n"
"QListWidget::item:selected:hover {\n"
"    background-color: #092F6B;\n"
"    color: white;\n"
"}")
        self.allFilesList.setFrameShape(QFrame.Shape.StyledPanel)
        self.allFilesList.setFrameShadow(QFrame.Shadow.Sunken)
        self.allFilesList.setLineWidth(1)
        self.allFilesList.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.rightMovebutton = QPushButton(self.centralwidget)
        self.rightMovebutton.setObjectName(u"rightMovebutton")
        self.rightMovebutton.setGeometry(QRect(280, 260, 21, 24))
        self.rightMovebutton.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.leftMoveBotton = QPushButton(self.centralwidget)
        self.leftMoveBotton.setObjectName(u"leftMoveBotton")
        self.leftMoveBotton.setGeometry(QRect(280, 310, 21, 24))
        self.leftMoveBotton.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.label_2 = QLabel(self.centralwidget)
        self.label_2.setObjectName(u"label_2")
        self.label_2.setGeometry(QRect(30, 95, 231, 16))
        self.label_2.setStyleSheet(u"QLabel{background-color: #E4E8EB; padding: 2px; color:black;}")
        self.label_2.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label_3 = QLabel(self.centralwidget)
        self.label_3.setObjectName(u"label_3")
        self.label_3.setGeometry(QRect(320, 95, 231, 16))
        self.label_3.setStyleSheet(u"QLabel{background-color: #E4E8EB; padding: 2px; color:black;}")
        self.label_3.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label_4 = QLabel(self.centralwidget)
        self.label_4.setObjectName(u"label_4")
        self.label_4.setGeometry(QRect(190, 560, 71, 16))
        self.label_5 = QLabel(self.centralwidget)
        self.label_5.setObjectName(u"label_5")
        self.label_5.setGeometry(QRect(480, 560, 71, 16))
        self.selectAllCheckBox = QCheckBox(self.centralwidget)
        self.selectAllCheckBox.setObjectName(u"selectAllCheckBox")
        self.selectAllCheckBox.setGeometry(QRect(30, 560, 79, 20))
        self.selectAllCheckBox.setFont(font)
        self.selectAllCheckBox.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.selectAllCheckBox.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.selectAllCheckBox.setStyleSheet(u"QCheckBox::indicator {\n"
"    width: 16px;\n"
"    height: 16px;\n"
"}\n"
"QCheckBox::indicator:checked{image:url(:/icons/icons/circle-check.png)}QCheckBox::indicator:unchecked{image:url(:/icons/icons/circle-dashed.png)}")
        self.selectAllCheckBox_2 = QCheckBox(self.centralwidget)
        self.selectAllCheckBox_2.setObjectName(u"selectAllCheckBox_2")
        self.selectAllCheckBox_2.setGeometry(QRect(320, 560, 79, 20))
        self.selectAllCheckBox_2.setFont(font)
        self.selectAllCheckBox_2.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.selectAllCheckBox_2.setStyleSheet(u"QCheckBox::indicator {\n"
"    width: 16px;\n"
"    height: 16px;\n"
"}\n"
"QCheckBox::indicator:checked{image:url(:/icons/icons/circle-check.png)}QCheckBox::indicator:unchecked{image:url(:/icons/icons/circle-dashed.png)}")
        self.probeTextEdit = QTextEdit(self.centralwidget)
        self.probeTextEdit.setObjectName(u"probeTextEdit")
        self.probeTextEdit.setGeometry(QRect(630, 70, 531, 241))
        self.probeTextEdit.setStyleSheet(u"QTextEdit{background-color:white;}\n"
"    QScrollBar:vertical {\n"
"        background: transparent;\n"
"        width: 6px;\n"
"        margin: 2px 0 2px 0;\n"
"    }\n"
"\n"
"    QScrollBar::handle:vertical {\n"
"        background: #1D83D5;\n"
"        border-radius: 3px;\n"
"        min-height: 20px;\n"
"    }\n"
"\n"
"    QScrollBar::add-line:vertical,\n"
"    QScrollBar::sub-line:vertical {\n"
"        height: 0px;\n"
"        background: none;\n"
"    }\n"
"\n"
"    QScrollBar::add-page:vertical,\n"
"    QScrollBar::sub-page:vertical {\n"
"        background: none;\n"
"    }")
        self.label_6 = QLabel(self.centralwidget)
        self.label_6.setObjectName(u"label_6")
        self.label_6.setGeometry(QRect(1080, 320, 81, 16))
        self.label_7 = QLabel(self.centralwidget)
        self.label_7.setObjectName(u"label_7")
        self.label_7.setGeometry(QRect(630, 50, 531, 20))
        self.label_7.setStyleSheet(u"QLabel{background-color: #E4E8EB; padding: 2px; color:black;}")
        self.label_7.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label_8 = QLabel(self.centralwidget)
        self.label_8.setObjectName(u"label_8")
        self.label_8.setGeometry(QRect(630, 350, 251, 16))
        self.label_9 = QLabel(self.centralwidget)
        self.label_9.setObjectName(u"label_9")
        self.label_9.setGeometry(QRect(630, 390, 281, 16))
        self.selectOutputButton = QPushButton(self.centralwidget)
        self.selectOutputButton.setObjectName(u"selectOutputButton")
        self.selectOutputButton.setGeometry(QRect(630, 440, 131, 24))
        self.selectOutputButton.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.outputPathLineEdit = QLineEdit(self.centralwidget)
        self.outputPathLineEdit.setObjectName(u"outputPathLineEdit")
        self.outputPathLineEdit.setGeometry(QRect(770, 440, 391, 22))
        self.outputPathLineEdit.setStyleSheet(u"QLineEdit{background-color: white;}")
        self.startAnalysisButton = QPushButton(self.centralwidget)
        self.startAnalysisButton.setObjectName(u"startAnalysisButton")
        self.startAnalysisButton.setGeometry(QRect(920, 500, 241, 41))
        font1 = QFont()
        font1.setPointSize(12)
        font1.setBold(True)
        self.startAnalysisButton.setFont(font1)
        self.startAnalysisButton.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.startAnalysisButton.setStyleSheet(u" QPushButton {\n"
"                background-color: #0D47A1;\n"
"                color: white;\n"
"                font-weight: bold;\n"
"                border: 1px solid #4183BC;\n"
"                border-radius: 6px;\n"
"            }")
        self.threasholdtextEdit = QDoubleSpinBox(self.centralwidget)
        self.threasholdtextEdit.setObjectName(u"threasholdtextEdit")
        self.threasholdtextEdit.setGeometry(QRect(890, 345, 70, 24))
        self.threasholdtextEdit.setMaximumSize(QSize(70, 16777215))
        self.threasholdtextEdit.setCursor(QCursor(Qt.CursorShape.ArrowCursor))
        self.threasholdtextEdit.setStyleSheet(u"QDoubleSpinBox{background-color: white; color:black;}")
        self.threasholdtextEdit.setDecimals(0)
        self.threasholdtextEdit.setMinimum(1.000000000000000)
        self.threasholdtextEdit.setValue(100.000000000000000)
        self.percentagetextEdit = QDoubleSpinBox(self.centralwidget)
        self.percentagetextEdit.setObjectName(u"percentagetextEdit")
        self.percentagetextEdit.setGeometry(QRect(930, 385, 70, 24))
        self.percentagetextEdit.setMaximumSize(QSize(70, 16777215))
        self.percentagetextEdit.setStyleSheet(u"QDoubleSpinBox{background-color: white; color:black;}")
        self.percentagetextEdit.setDecimals(0)
        self.percentagetextEdit.setMinimum(1.000000000000000)
        self.percentagetextEdit.setMaximum(100.000000000000000)
        self.percentagetextEdit.setValue(1.000000000000000)
        self.selectedFilesList = QListWidget(self.centralwidget)
        self.selectedFilesList.setObjectName(u"selectedFilesList")
        self.selectedFilesList.setGeometry(QRect(320, 110, 231, 441))
        self.selectedFilesList.setStyleSheet(u"QListWidget {\n"
"                background-color: white;\n"
"                color: #2c3e50;\n"
"                border: 1px solid #ccc;\n"
"            }\n"
"QListWidget::item:selected {\n"
"                background-color: #1663AE;\n"
"                color: white;\n"
"            }\n"
"QListWidget::item:hover {\n"
"                background-color: #EBF5FD;\n"
"                color : black;\n"
"            }\n"
"QScrollBar:vertical {\n"
"        background: transparent;\n"
"        width: 6px;\n"
"        margin: 2px 0 2px 0;\n"
"    }\n"
"\n"
"QScrollBar::handle:vertical {\n"
"        background: #1D83D5;\n"
"        border-radius: 3px;\n"
"        min-height: 20px;\n"
"    }\n"
"\n"
"QScrollBar::add-line:vertical,\n"
"QScrollBar::sub-line:vertical {\n"
"        height: 0px;\n"
"        background: none;\n"
"    }\n"
"\n"
"QScrollBar::add-page:vertical,\n"
"QScrollBar::sub-page:vertical {\n"
"        background: none;\n"
"    }\n"
"\n"
"QScrollBar:horizontal {\n"
"    background: transparent;\n"
"    h"
                        "eight: 6px;\n"
"    margin: 0 2px 0 2px;\n"
"}\n"
"\n"
"QScrollBar::handle:horizontal {\n"
"    background: #1D83D5;\n"
"    border-radius: 3px;\n"
"    min-width: 20px;\n"
"}\n"
"\n"
"QScrollBar::add-line:horizontal,\n"
"QScrollBar::sub-line:horizontal {\n"
"    width: 0px;\n"
"    background: none;\n"
"}\n"
"\n"
"QScrollBar::add-page:horizontal,\n"
"QScrollBar::sub-page:horizontal {\n"
"    background: none;\n"
"}\n"
"QListWidget::item:selected:hover {\n"
"    background-color: #092F6B;\n"
"    color: white;\n"
"}")
        self.selectedFilesList.setFrameShape(QFrame.Shape.StyledPanel)
        self.selectedFilesList.setFrameShadow(QFrame.Shadow.Sunken)
        self.selectedFilesList.setLineWidth(1)
        self.selectedFilesList.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.label_10 = QLabel(self.centralwidget)
        self.label_10.setObjectName(u"label_10")
        self.label_10.setGeometry(QRect(700, 500, 121, 16))
        font2 = QFont()
        font2.setPointSize(8)
        self.label_10.setFont(font2)
        self.parallelProcessing = QDoubleSpinBox(self.centralwidget)
        self.parallelProcessing.setObjectName(u"parallelProcessing")
        self.parallelProcessing.setGeometry(QRect(630, 495, 60, 24))
        self.parallelProcessing.setMaximumSize(QSize(60, 16777215))
        font3 = QFont()
        font3.setBold(False)
        self.parallelProcessing.setFont(font3)
        self.parallelProcessing.setCursor(QCursor(Qt.CursorShape.ArrowCursor))
        self.parallelProcessing.setStyleSheet(u"QDoubleSpinBox{background-color: white; color:black;}")
        self.parallelProcessing.setDecimals(0)
        self.parallelProcessing.setMinimum(1.000000000000000)
        self.parallelProcessing.setMaximum(61.000000000000000)
        self.parallelProcessing.setValue(4.000000000000000)
        self.parallelProcessingCheckbox = QCheckBox(self.centralwidget)
        self.parallelProcessingCheckbox.setObjectName(u"parallelProcessingCheckbox")
        self.parallelProcessingCheckbox.setGeometry(QRect(630, 530, 211, 20))
        palette = QPalette()
        brush = QBrush(QColor(0, 0, 0, 255))
        brush.setStyle(Qt.BrushStyle.SolidPattern)
        palette.setBrush(QPalette.ColorGroup.Active, QPalette.ColorRole.WindowText, brush)
        brush1 = QBrush(QColor(249, 250, 250, 255))
        brush1.setStyle(Qt.BrushStyle.SolidPattern)
        palette.setBrush(QPalette.ColorGroup.Active, QPalette.ColorRole.Button, brush1)
        palette.setBrush(QPalette.ColorGroup.Active, QPalette.ColorRole.Text, brush)
        palette.setBrush(QPalette.ColorGroup.Active, QPalette.ColorRole.ButtonText, brush)
        palette.setBrush(QPalette.ColorGroup.Active, QPalette.ColorRole.Base, brush1)
        palette.setBrush(QPalette.ColorGroup.Active, QPalette.ColorRole.Window, brush1)
        brush2 = QBrush(QColor(0, 0, 0, 128))
        brush2.setStyle(Qt.BrushStyle.SolidPattern)
#if QT_VERSION >= QT_VERSION_CHECK(5, 12, 0)
        palette.setBrush(QPalette.ColorGroup.Active, QPalette.ColorRole.PlaceholderText, brush2)
#endif
        palette.setBrush(QPalette.ColorGroup.Inactive, QPalette.ColorRole.WindowText, brush)
        palette.setBrush(QPalette.ColorGroup.Inactive, QPalette.ColorRole.Button, brush1)
        palette.setBrush(QPalette.ColorGroup.Inactive, QPalette.ColorRole.Text, brush)
        palette.setBrush(QPalette.ColorGroup.Inactive, QPalette.ColorRole.ButtonText, brush)
        palette.setBrush(QPalette.ColorGroup.Inactive, QPalette.ColorRole.Base, brush1)
        palette.setBrush(QPalette.ColorGroup.Inactive, QPalette.ColorRole.Window, brush1)
#if QT_VERSION >= QT_VERSION_CHECK(5, 12, 0)
        palette.setBrush(QPalette.ColorGroup.Inactive, QPalette.ColorRole.PlaceholderText, brush2)
#endif
        palette.setBrush(QPalette.ColorGroup.Disabled, QPalette.ColorRole.WindowText, brush)
        palette.setBrush(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Button, brush1)
        palette.setBrush(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text, brush)
        brush3 = QBrush(QColor(255, 255, 255, 135))
        brush3.setStyle(Qt.BrushStyle.SolidPattern)
        palette.setBrush(QPalette.ColorGroup.Disabled, QPalette.ColorRole.ButtonText, brush3)
        palette.setBrush(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Base, brush1)
        palette.setBrush(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Window, brush1)
#if QT_VERSION >= QT_VERSION_CHECK(5, 12, 0)
        palette.setBrush(QPalette.ColorGroup.Disabled, QPalette.ColorRole.PlaceholderText, brush2)
#endif
        self.parallelProcessingCheckbox.setPalette(palette)
        self.parallelProcessingCheckbox.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.parallelProcessingCheckbox.setStyleSheet(u"QCheckBox::indicator {\n"
"    width: 16px;\n"
"    height: 16px;\n"
"}\n"
"QCheckBox::indicator:checked{image:url(:/icons/icons/circle-check.png)}QCheckBox::indicator:unchecked{image:url(:/icons/icons/circle-dashed.png)}")
        self.label_11 = QLabel(self.centralwidget)
        self.label_11.setObjectName(u"label_11")
        self.label_11.setGeometry(QRect(1120, 0, 81, 71))
        self.label_11.setPixmap(QPixmap(u":/icons/icons/IDL_mark.PNG"))
        self.label_11.setScaledContents(True)
        In_silico_sequence_mining.setCentralWidget(self.centralwidget)
        self.label_11.raise_()
        self.label.raise_()
        self.selectFolderButton.raise_()
        self.inputPathLineEdit.raise_()
        self.rightMovebutton.raise_()
        self.leftMoveBotton.raise_()
        self.label_2.raise_()
        self.label_3.raise_()
        self.label_4.raise_()
        self.label_5.raise_()
        self.selectAllCheckBox.raise_()
        self.selectAllCheckBox_2.raise_()
        self.probeTextEdit.raise_()
        self.label_6.raise_()
        self.label_7.raise_()
        self.label_8.raise_()
        self.label_9.raise_()
        self.selectOutputButton.raise_()
        self.outputPathLineEdit.raise_()
        self.startAnalysisButton.raise_()
        self.threasholdtextEdit.raise_()
        self.percentagetextEdit.raise_()
        self.selectedFilesList.raise_()
        self.radioButton.raise_()
        self.radioButton_2.raise_()
        self.label_10.raise_()
        self.parallelProcessing.raise_()
        self.parallelProcessingCheckbox.raise_()
        self.allFilesList.raise_()

        self.retranslateUi(In_silico_sequence_mining)

        QMetaObject.connectSlotsByName(In_silico_sequence_mining)
    # setupUi

    def retranslateUi(self, In_silico_sequence_mining):
        In_silico_sequence_mining.setWindowTitle(QCoreApplication.translate("In_silico_sequence_mining", u"MainWindow", None))
        self.radioButton.setText(QCoreApplication.translate("In_silico_sequence_mining", u"FASTQ", None))
        self.radioButton_2.setText(QCoreApplication.translate("In_silico_sequence_mining", u"FASTA", None))
        self.label.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Select file format :", None))
        self.selectFolderButton.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Select Folder", None))
        self.rightMovebutton.setText(QCoreApplication.translate("In_silico_sequence_mining", u"\u25b6", None))
        self.leftMoveBotton.setText(QCoreApplication.translate("In_silico_sequence_mining", u"\u25c0", None))
        self.label_2.setText(QCoreApplication.translate("In_silico_sequence_mining", u"All FASTA/FASTQ files", None))
        self.label_3.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Selected FASTA/FASTQ files", None))
        self.label_4.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Total :", None))
        self.label_5.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Total :", None))
        self.selectAllCheckBox.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Select all", None))
        self.selectAllCheckBox_2.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Select all", None))
        self.label_6.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Probe : ", None))
        self.label_7.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Probe sequence for diverse targets (in FASTA format)", None))
        self.label_8.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Percentage of match to the probe sequence :", None))
        self.label_9.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Percentage of reads to extract from the total reads :", None))
        self.selectOutputButton.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Select Output Path", None))
        self.startAnalysisButton.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Start Analysis", None))
        self.label_10.setText(QCoreApplication.translate("In_silico_sequence_mining", u"(min : 1 ~ max : 61)", None))
        self.parallelProcessingCheckbox.setText(QCoreApplication.translate("In_silico_sequence_mining", u"Use custom parallel processing", None))
        self.label_11.setText("")
    # retranslateUi

