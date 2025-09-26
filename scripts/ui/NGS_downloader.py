# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'NGS_downloaderwbuMCO.ui'
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
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QComboBox, QFormLayout,
    QFrame, QHBoxLayout, QHeaderView, QLabel,
    QLineEdit, QMainWindow, QProgressBar, QPushButton,
    QSizePolicy, QSpacerItem, QTreeWidget, QTreeWidgetItem,
    QWidget)
from ui import resource_rc

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(1284, 749)
        MainWindow.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        MainWindow.setAnimated(True)
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.centralwidget.sizePolicy().hasHeightForWidth())
        self.centralwidget.setSizePolicy(sizePolicy)
        self.centralwidget.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.centralwidget.setStyleSheet(u"QWidget{background-color: #F3F4F5; color: black;}\n"
"     QPushButton {\n"
"                background-color: #547296;\n"
"                color: white;\n"
"                font-weight: bold;\n"
"                border: 1px solid #40546C;\n"
"                border-radius: 6px;\n"
"            }\n"
"            QPushButton:hover {\n"
"                background-color: #40546C;\n"
"            }\n"
"QLabel{background-color: #F3F4F5}")
        self.treeSearchResults = QTreeWidget(self.centralwidget)
        __qtreewidgetitem = QTreeWidgetItem()
        __qtreewidgetitem.setTextAlignment(8, Qt.AlignCenter);
        __qtreewidgetitem.setTextAlignment(7, Qt.AlignCenter);
        __qtreewidgetitem.setTextAlignment(6, Qt.AlignCenter);
        __qtreewidgetitem.setTextAlignment(5, Qt.AlignCenter);
        __qtreewidgetitem.setTextAlignment(4, Qt.AlignCenter);
        __qtreewidgetitem.setTextAlignment(3, Qt.AlignCenter);
        __qtreewidgetitem.setTextAlignment(2, Qt.AlignCenter);
        __qtreewidgetitem.setTextAlignment(1, Qt.AlignCenter);
        __qtreewidgetitem.setTextAlignment(0, Qt.AlignCenter);
        self.treeSearchResults.setHeaderItem(__qtreewidgetitem)
        self.treeSearchResults.setObjectName(u"treeSearchResults")
        self.treeSearchResults.setGeometry(QRect(10, 140, 1261, 192))
        self.treeSearchResults.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.treeSearchResults.setStyleSheet(u"QTreeWidget{background-color: white;\n"
"        border: 1px solid #ccc;\n"
"        border-radius: 3px\n"
"    }\n"
"QHeaderView::section {\n"
"    background-color: #2f3e51; \n"
"    color: #e0f0ff;                \n"
"    padding: 3px;\n"
"    border: 1px solid #1e2a38;   \n"
"    font-weight: bold;\n"
"}\n"
"QScrollBar:vertical {\n"
"        background: transparent;\n"
"        width: 6px;\n"
"        margin: 2px 0 2px 0;\n"
"    }\n"
"\n"
"QScrollBar::handle:vertical {\n"
"        background: #9E9E9E;\n"
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
"    height: 6px;\n"
"    margin: 0 2px 0 2px;\n"
"}\n"
"\n"
"QScrollBar::handle:horizontal {\n"
"    background: #"
                        "9E9E9E;\n"
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
"\n"
"QTreeWidget::item:selected {\n"
"    background-color: #1663AE; \n"
"    color: white;\n"
"}\n"
"\n"
"QTreeWidget::item:hover {\n"
"    background-color: #EBF5FD; \n"
"    color: black;\n"
"}\n"
"\n"
"QTreeWidget::item {\n"
"    border-bottom: 1px dotted #d0d0d0; \n"
"}")
        self.treeSearchResults.setLineWidth(1)
        self.treeSearchResults.setRootIsDecorated(False)
        self.treeSearchResults.header().setMinimumSectionSize(1)
        self.treeSearchResults.header().setDefaultSectionSize(139)
        self.treeDownloadList = QTreeWidget(self.centralwidget)
        __qtreewidgetitem1 = QTreeWidgetItem()
        __qtreewidgetitem1.setTextAlignment(8, Qt.AlignCenter);
        __qtreewidgetitem1.setTextAlignment(7, Qt.AlignCenter);
        __qtreewidgetitem1.setTextAlignment(6, Qt.AlignCenter);
        __qtreewidgetitem1.setTextAlignment(5, Qt.AlignCenter);
        __qtreewidgetitem1.setTextAlignment(4, Qt.AlignCenter);
        __qtreewidgetitem1.setTextAlignment(3, Qt.AlignCenter);
        __qtreewidgetitem1.setTextAlignment(2, Qt.AlignCenter);
        __qtreewidgetitem1.setTextAlignment(1, Qt.AlignCenter);
        __qtreewidgetitem1.setTextAlignment(0, Qt.AlignCenter);
        self.treeDownloadList.setHeaderItem(__qtreewidgetitem1)
        self.treeDownloadList.setObjectName(u"treeDownloadList")
        self.treeDownloadList.setGeometry(QRect(10, 420, 1261, 201))
        self.treeDownloadList.setStyleSheet(u"QTreeWidget{background-color: white;\n"
"        border: 1px solid #ccc;\n"
"        border-radius: 3px\n"
"    }\n"
"QHeaderView::section {\n"
"    background-color: #2f3e51; \n"
"    color: #e0f0ff;                \n"
"    padding: 3px;\n"
"    border: 1px solid #1e2a38;   \n"
"    font-weight: bold;\n"
"}\n"
"QScrollBar:vertical {\n"
"        background: transparent;\n"
"        width: 6px;\n"
"        margin: 2px 0 2px 0;\n"
"    }\n"
"\n"
"QScrollBar::handle:vertical {\n"
"        background: #9E9E9E;\n"
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
"    height: 6px;\n"
"    margin: 0 2px 0 2px;\n"
"}\n"
"\n"
"QScrollBar::handle:horizontal {\n"
"    background: #"
                        "9E9E9E;\n"
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
"\n"
"QTreeWidget::item:selected {\n"
"    background-color: #1663AE;\n"
"    color: white;\n"
"}\n"
"\n"
"QTreeWidget::item:hover {\n"
"    background-color: #EBF5FD; \n"
"    color: black;\n"
"}\n"
"\n"
"QTreeWidget::item {\n"
"    border-bottom: 1px dotted #d0d0d0; \n"
"}")
        self.treeDownloadList.setFrameShadow(QFrame.Shadow.Sunken)
        self.treeDownloadList.setAlternatingRowColors(True)
        self.treeDownloadList.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.treeDownloadList.setIndentation(20)
        self.treeDownloadList.setRootIsDecorated(False)
        self.treeDownloadList.header().setMinimumSectionSize(1)
        self.treeDownloadList.header().setDefaultSectionSize(139)
        self.treeDownloadList.header().setStretchLastSection(False)
        self.buttonDownloadStart = QPushButton(self.centralwidget)
        self.buttonDownloadStart.setObjectName(u"buttonDownloadStart")
        self.buttonDownloadStart.setGeometry(QRect(550, 700, 181, 41))
        self.buttonDownloadStart.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.buttonDownloadStart.setStyleSheet(u" QPushButton {\n"
"                background-color: #0D47A1;\n"
"                color: white;\n"
"                font-weight: bold;\n"
"                border: 1px solid #4183BC;\n"
"                border-radius: 6px;\n"
"            }")
        self.horizontalLayoutWidget_3 = QWidget(self.centralwidget)
        self.horizontalLayoutWidget_3.setObjectName(u"horizontalLayoutWidget_3")
        self.horizontalLayoutWidget_3.setGeometry(QRect(10, 670, 1261, 26))
        self.horizontalLayout_4 = QHBoxLayout(self.horizontalLayoutWidget_3)
        self.horizontalLayout_4.setSpacing(2)
        self.horizontalLayout_4.setObjectName(u"horizontalLayout_4")
        self.horizontalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.horizontalSpacer_10 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_4.addItem(self.horizontalSpacer_10)

        self.labelFormat = QLabel(self.horizontalLayoutWidget_3)
        self.labelFormat.setObjectName(u"labelFormat")
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.labelFormat.sizePolicy().hasHeightForWidth())
        self.labelFormat.setSizePolicy(sizePolicy1)
        self.labelFormat.setMinimumSize(QSize(90, 0))

        self.horizontalLayout_4.addWidget(self.labelFormat)

        self.comboFormat = QComboBox(self.horizontalLayoutWidget_3)
        self.comboFormat.addItem("")
        self.comboFormat.addItem("")
        self.comboFormat.setObjectName(u"comboFormat")
        self.comboFormat.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.comboFormat.setStyleSheet(u"QComboBox {\n"
"    border: 1px solid #5a5a5a; \n"
"    border-radius: 1px;             /* \ubaa8\uc11c\ub9ac \ub465\uae00\uac8c */\n"
"    padding: 2px;\n"
"    background-color: white;      /* \ubc30\uacbd\uc0c9 */\n"
"    color: black;                   /* \ud14d\uc2a4\ud2b8 \uc0c9\uc0c1 */\n"
"}\n"
"\n"
"QComboBox QAbstractItemView {\n"
"    background-color: white;           /* \ub4dc\ub86d\ub2e4\uc6b4 \ubc30\uacbd */\n"
"    color: black;                      /* \ud56d\ubaa9 \ud14d\uc2a4\ud2b8 \uc0c9 */\n"
"    selection-background-color: #448aff; /* \uc120\ud0dd\ub41c \ud56d\ubaa9 \ubc30\uacbd */\n"
"    selection-color: white;             /* \uc120\ud0dd\ub41c \ud56d\ubaa9 \ud14d\uc2a4\ud2b8 */\n"
"    border: 1px solid #555555;\n"
"}")

        self.horizontalLayout_4.addWidget(self.comboFormat)

        self.horizontalSpacer_11 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_4.addItem(self.horizontalSpacer_11)

        self.horizontalLayoutWidget_6 = QWidget(self.centralwidget)
        self.horizontalLayoutWidget_6.setObjectName(u"horizontalLayoutWidget_6")
        self.horizontalLayoutWidget_6.setGeometry(QRect(10, 630, 1261, 31))
        self.horizontalLayout_7 = QHBoxLayout(self.horizontalLayoutWidget_6)
        self.horizontalLayout_7.setSpacing(15)
        self.horizontalLayout_7.setObjectName(u"horizontalLayout_7")
        self.horizontalLayout_7.setContentsMargins(0, 0, 0, 0)
        self.horizontalSpacer_3 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_7.addItem(self.horizontalSpacer_3)

        self.buttonSelectDirectory = QPushButton(self.horizontalLayoutWidget_6)
        self.buttonSelectDirectory.setObjectName(u"buttonSelectDirectory")
        sizePolicy1.setHeightForWidth(self.buttonSelectDirectory.sizePolicy().hasHeightForWidth())
        self.buttonSelectDirectory.setSizePolicy(sizePolicy1)
        self.buttonSelectDirectory.setMinimumSize(QSize(110, 0))
        self.buttonSelectDirectory.setMaximumSize(QSize(150, 20))
        self.buttonSelectDirectory.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))

        self.horizontalLayout_7.addWidget(self.buttonSelectDirectory)

        self.lineEditDirectory = QLineEdit(self.horizontalLayoutWidget_6)
        self.lineEditDirectory.setObjectName(u"lineEditDirectory")
        sizePolicy1.setHeightForWidth(self.lineEditDirectory.sizePolicy().hasHeightForWidth())
        self.lineEditDirectory.setSizePolicy(sizePolicy1)
        self.lineEditDirectory.setMinimumSize(QSize(400, 0))
        self.lineEditDirectory.setMaximumSize(QSize(500, 20))
        self.lineEditDirectory.setStyleSheet(u"QLineEdit{background-color: white;}")
        self.lineEditDirectory.setMaxLength(1000)

        self.horizontalLayout_7.addWidget(self.lineEditDirectory)

        self.horizontalSpacer_7 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_7.addItem(self.horizontalSpacer_7)

        self.horizontalLayoutWidget_7 = QWidget(self.centralwidget)
        self.horizontalLayoutWidget_7.setObjectName(u"horizontalLayoutWidget_7")
        self.horizontalLayoutWidget_7.setGeometry(QRect(10, 340, 1261, 41))
        self.horizontalLayout_8 = QHBoxLayout(self.horizontalLayoutWidget_7)
        self.horizontalLayout_8.setObjectName(u"horizontalLayout_8")
        self.horizontalLayout_8.setContentsMargins(0, 0, 0, 0)
        self.horizontalSpacer_8 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_8.addItem(self.horizontalSpacer_8)

        self.buttonMoveToDownload = QPushButton(self.horizontalLayoutWidget_7)
        self.buttonMoveToDownload.setObjectName(u"buttonMoveToDownload")
        self.buttonMoveToDownload.setMinimumSize(QSize(180, 25))
        self.buttonMoveToDownload.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))

        self.horizontalLayout_8.addWidget(self.buttonMoveToDownload)

        self.buttonRemoveToDownloadList = QPushButton(self.horizontalLayoutWidget_7)
        self.buttonRemoveToDownloadList.setObjectName(u"buttonRemoveToDownloadList")
        self.buttonRemoveToDownloadList.setMinimumSize(QSize(180, 25))
        self.buttonRemoveToDownloadList.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))

        self.horizontalLayout_8.addWidget(self.buttonRemoveToDownloadList)

        self.horizontalSpacer_9 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_8.addItem(self.horizontalSpacer_9)

        self.horizontalLayoutWidget_8 = QWidget(self.centralwidget)
        self.horizontalLayoutWidget_8.setObjectName(u"horizontalLayoutWidget_8")
        self.horizontalLayoutWidget_8.setGeometry(QRect(10, 110, 281, 31))
        self.horizontalLayout_9 = QHBoxLayout(self.horizontalLayoutWidget_8)
        self.horizontalLayout_9.setObjectName(u"horizontalLayout_9")
        self.horizontalLayout_9.setContentsMargins(0, 0, 0, 0)
        self.buttonSelectAllinSearchList = QPushButton(self.horizontalLayoutWidget_8)
        self.buttonSelectAllinSearchList.setObjectName(u"buttonSelectAllinSearchList")
        sizePolicy.setHeightForWidth(self.buttonSelectAllinSearchList.sizePolicy().hasHeightForWidth())
        self.buttonSelectAllinSearchList.setSizePolicy(sizePolicy)
        self.buttonSelectAllinSearchList.setMinimumSize(QSize(80, 20))
        self.buttonSelectAllinSearchList.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))

        self.horizontalLayout_9.addWidget(self.buttonSelectAllinSearchList)

        self.labelSearchCount = QLabel(self.horizontalLayoutWidget_8)
        self.labelSearchCount.setObjectName(u"labelSearchCount")

        self.horizontalLayout_9.addWidget(self.labelSearchCount)

        self.layoutWidget = QWidget(self.centralwidget)
        self.layoutWidget.setObjectName(u"layoutWidget")
        self.layoutWidget.setGeometry(QRect(10, 390, 271, 31))
        self.horizontalLayout_11 = QHBoxLayout(self.layoutWidget)
        self.horizontalLayout_11.setObjectName(u"horizontalLayout_11")
        self.horizontalLayout_11.setContentsMargins(0, 0, 0, 0)
        self.buttonSelectAllinSearchList_2 = QPushButton(self.layoutWidget)
        self.buttonSelectAllinSearchList_2.setObjectName(u"buttonSelectAllinSearchList_2")
        sizePolicy.setHeightForWidth(self.buttonSelectAllinSearchList_2.sizePolicy().hasHeightForWidth())
        self.buttonSelectAllinSearchList_2.setSizePolicy(sizePolicy)
        self.buttonSelectAllinSearchList_2.setMinimumSize(QSize(80, 20))
        self.buttonSelectAllinSearchList_2.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))

        self.horizontalLayout_11.addWidget(self.buttonSelectAllinSearchList_2)

        self.labelDownloadCount = QLabel(self.layoutWidget)
        self.labelDownloadCount.setObjectName(u"labelDownloadCount")

        self.horizontalLayout_11.addWidget(self.labelDownloadCount)

        self.IDLlogo = QLabel(self.centralwidget)
        self.IDLlogo.setObjectName(u"IDLlogo")
        self.IDLlogo.setGeometry(QRect(1210, 10, 61, 51))
        self.IDLlogo.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.IDLlogo.setTextFormat(Qt.TextFormat.PlainText)
        self.IDLlogo.setPixmap(QPixmap(u":/icons/icons/IDL_mark.PNG"))
        self.IDLlogo.setScaledContents(True)
        self.filteringComboBox = QComboBox(self.centralwidget)
        self.filteringComboBox.setObjectName(u"filteringComboBox")
        self.filteringComboBox.setGeometry(QRect(930, 110, 341, 24))
        self.filteringComboBox.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.filteringComboBox.setStyleSheet(u"QComboBox {\n"
"    border: 1px solid #5a5a5a; \n"
"    border-radius: 1px;            \n"
"    padding: 2px;\n"
"    background-color: white;     \n"
"    color: black;                   \n"
"}\n"
"\n"
"QComboBox QAbstractItemView {\n"
"    background-color: white;          \n"
"    color: black;                     \n"
"    selection-background-color: #01579B;\n"
"    selection-color: white;             \n"
"    border: 1px solid #555555;\n"
"}\n"
"\n"
"QComboBox QAbstractItemView::item:hover {\n"
"    border-color: #1663AE;\n"
"    background-color: #EBF5FD; \n"
"    color: black;}\n"
"\n"
"\n"
"QScrollBar:vertical {\n"
"        background: transparent;\n"
"        width: 6px;\n"
"        margin: 2px 0 2px 0;\n"
"    }\n"
"\n"
"QScrollBar::handle:vertical {\n"
"        background: #9E9E9E;\n"
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
""
                        "QScrollBar::add-page:vertical,\n"
"QScrollBar::sub-page:vertical {\n"
"        background: none;\n"
"    }\n"
"\n"
"")
        self.Filteringlabel = QLabel(self.centralwidget)
        self.Filteringlabel.setObjectName(u"Filteringlabel")
        self.Filteringlabel.setGeometry(QRect(930, 90, 341, 20))
        self.Filteringlabel.setStyleSheet(u"QLabel{      background-color: #547296;\n"
"                color: white;\n"
"                font-weight: bold;\n"
"                border: 1px solid #40546C;\n"
"                border-radius: 1px;\n"
"           }")
        self.Filteringlabel.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.resetpushButton = QPushButton(self.centralwidget)
        self.resetpushButton.setObjectName(u"resetpushButton")
        self.resetpushButton.setGeometry(QRect(880, 110, 41, 21))
        self.resetpushButton.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.resetpushButton.setStyleSheet(u"     QPushButton {\n"
"                background-color: #9E9E9E;\n"
"                color: white;\n"
"                font-weight: bold;\n"
"                border: 1px solid black;\n"
"                border-radius: 6px;\n"
"            }\n"
"            QPushButton:hover {\n"
"                background-color: #616161;\n"
"            }")
        self.horizontalLayoutWidget_2 = QWidget(self.centralwidget)
        self.horizontalLayoutWidget_2.setObjectName(u"horizontalLayoutWidget_2")
        self.horizontalLayoutWidget_2.setGeometry(QRect(10, 10, 151, 26))
        self.horizontalLayout_3 = QHBoxLayout(self.horizontalLayoutWidget_2)
        self.horizontalLayout_3.setSpacing(0)
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.labelDatabase = QLabel(self.horizontalLayoutWidget_2)
        self.labelDatabase.setObjectName(u"labelDatabase")
        sizePolicy2 = QSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Preferred)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.labelDatabase.sizePolicy().hasHeightForWidth())
        self.labelDatabase.setSizePolicy(sizePolicy2)
        self.labelDatabase.setMinimumSize(QSize(90, 0))
        self.labelDatabase.setMaximumSize(QSize(90, 16777215))

        self.horizontalLayout_3.addWidget(self.labelDatabase)

        self.comboDatabase = QComboBox(self.horizontalLayoutWidget_2)
        self.comboDatabase.addItem("")
        self.comboDatabase.addItem("")
        self.comboDatabase.setObjectName(u"comboDatabase")
        sizePolicy.setHeightForWidth(self.comboDatabase.sizePolicy().hasHeightForWidth())
        self.comboDatabase.setSizePolicy(sizePolicy)
        self.comboDatabase.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.comboDatabase.setStyleSheet(u"QComboBox {\n"
"    border: 1px solid #5a5a5a; \n"
"    border-radius: 1px;             /* \ubaa8\uc11c\ub9ac \ub465\uae00\uac8c */\n"
"    padding: 2px;\n"
"    background-color: white;      /* \ubc30\uacbd\uc0c9 */\n"
"    color: black;                   /* \ud14d\uc2a4\ud2b8 \uc0c9\uc0c1 */\n"
"}\n"
"\n"
"QComboBox QAbstractItemView {\n"
"    background-color: white;           /* \ub4dc\ub86d\ub2e4\uc6b4 \ubc30\uacbd */\n"
"    color: black;                      /* \ud56d\ubaa9 \ud14d\uc2a4\ud2b8 \uc0c9 */\n"
"    selection-background-color: #448aff; /* \uc120\ud0dd\ub41c \ud56d\ubaa9 \ubc30\uacbd */\n"
"    selection-color: white;             /* \uc120\ud0dd\ub41c \ud56d\ubaa9 \ud14d\uc2a4\ud2b8 */\n"
"    border: 1px solid #555555;\n"
"}")

        self.horizontalLayout_3.addWidget(self.comboDatabase)

        self.formLayoutWidget_2 = QWidget(self.centralwidget)
        self.formLayoutWidget_2.setObjectName(u"formLayoutWidget_2")
        self.formLayoutWidget_2.setGeometry(QRect(440, 10, 391, 124))
        self.formLayout = QFormLayout(self.formLayoutWidget_2)
        self.formLayout.setObjectName(u"formLayout")
        self.formLayout.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        self.formLayout.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.formLayout.setContentsMargins(0, 0, 0, 0)
        self.lineEditEmail_2 = QLineEdit(self.formLayoutWidget_2)
        self.lineEditEmail_2.setObjectName(u"lineEditEmail_2")
        sizePolicy3 = QSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Fixed)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(self.lineEditEmail_2.sizePolicy().hasHeightForWidth())
        self.lineEditEmail_2.setSizePolicy(sizePolicy3)
        self.lineEditEmail_2.setMinimumSize(QSize(300, 0))
        self.lineEditEmail_2.setMaximumSize(QSize(255, 20))
        self.lineEditEmail_2.setStyleSheet(u"QLineEdit{background-color: white;}")

        self.formLayout.setWidget(0, QFormLayout.ItemRole.LabelRole, self.lineEditEmail_2)

        self.buttonSaveEmail_2 = QPushButton(self.formLayoutWidget_2)
        self.buttonSaveEmail_2.setObjectName(u"buttonSaveEmail_2")
        sizePolicy3.setHeightForWidth(self.buttonSaveEmail_2.sizePolicy().hasHeightForWidth())
        self.buttonSaveEmail_2.setSizePolicy(sizePolicy3)
        self.buttonSaveEmail_2.setMinimumSize(QSize(85, 20))
        self.buttonSaveEmail_2.setMaximumSize(QSize(80, 16777215))
        self.buttonSaveEmail_2.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))

        self.formLayout.setWidget(0, QFormLayout.ItemRole.FieldRole, self.buttonSaveEmail_2)

        self.lineEditQuery = QLineEdit(self.formLayoutWidget_2)
        self.lineEditQuery.setObjectName(u"lineEditQuery")
        self.lineEditQuery.setEnabled(True)
        sizePolicy4 = QSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Preferred)
        sizePolicy4.setHorizontalStretch(0)
        sizePolicy4.setVerticalStretch(0)
        sizePolicy4.setHeightForWidth(self.lineEditQuery.sizePolicy().hasHeightForWidth())
        self.lineEditQuery.setSizePolicy(sizePolicy4)
        self.lineEditQuery.setMinimumSize(QSize(300, 0))
        self.lineEditQuery.setMaximumSize(QSize(300, 20))
        self.lineEditQuery.setStyleSheet(u"QLineEdit{background-color: white;}")

        self.formLayout.setWidget(1, QFormLayout.ItemRole.LabelRole, self.lineEditQuery)

        self.buttonSearch = QPushButton(self.formLayoutWidget_2)
        self.buttonSearch.setObjectName(u"buttonSearch")
        sizePolicy.setHeightForWidth(self.buttonSearch.sizePolicy().hasHeightForWidth())
        self.buttonSearch.setSizePolicy(sizePolicy)
        self.buttonSearch.setMinimumSize(QSize(85, 20))
        self.buttonSearch.setMaximumSize(QSize(80, 16777215))
        self.buttonSearch.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))

        self.formLayout.setWidget(1, QFormLayout.ItemRole.FieldRole, self.buttonSearch)

        self.buttonLoadMore = QPushButton(self.formLayoutWidget_2)
        self.buttonLoadMore.setObjectName(u"buttonLoadMore")
        sizePolicy1.setHeightForWidth(self.buttonLoadMore.sizePolicy().hasHeightForWidth())
        self.buttonLoadMore.setSizePolicy(sizePolicy1)
        self.buttonLoadMore.setMinimumSize(QSize(85, 20))
        self.buttonLoadMore.setMaximumSize(QSize(85, 16777215))
        self.buttonLoadMore.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))

        self.formLayout.setWidget(3, QFormLayout.ItemRole.FieldRole, self.buttonLoadMore)

        self.SearchClearbutton = QPushButton(self.formLayoutWidget_2)
        self.SearchClearbutton.setObjectName(u"SearchClearbutton")
        sizePolicy3.setHeightForWidth(self.SearchClearbutton.sizePolicy().hasHeightForWidth())
        self.SearchClearbutton.setSizePolicy(sizePolicy3)
        self.SearchClearbutton.setMinimumSize(QSize(85, 20))
        self.SearchClearbutton.setMaximumSize(QSize(85, 16777215))
        self.SearchClearbutton.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.SearchClearbutton.setStyleSheet(u"     QPushButton {\n"
"                background-color: #9E9E9E;\n"
"                color: white;\n"
"                font-weight: bold;\n"
"                border: 1px solid black;\n"
"                border-radius: 6px;\n"
"            }\n"
"            QPushButton:hover {\n"
"                background-color: #616161;\n"
"            }")

        self.formLayout.setWidget(2, QFormLayout.ItemRole.FieldRole, self.SearchClearbutton)

        self.progressSearch = QProgressBar(self.formLayoutWidget_2)
        self.progressSearch.setObjectName(u"progressSearch")
        sizePolicy5 = QSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Minimum)
        sizePolicy5.setHorizontalStretch(1)
        sizePolicy5.setVerticalStretch(0)
        sizePolicy5.setHeightForWidth(self.progressSearch.sizePolicy().hasHeightForWidth())
        self.progressSearch.setSizePolicy(sizePolicy5)
        self.progressSearch.setMinimumSize(QSize(300, 20))
        self.progressSearch.setMaximumSize(QSize(1300, 15))
        self.progressSearch.setStyleSheet(u"QProgressBar {\n"
"    border: 1px solid #00CC66;\n"
"    border-radius: 0px;\n"
"    background-color: #E9F0F4;\n"
"    text-align: center;\n"
"    font: bold 12px;\n"
"    color: #006400; /* \ud14d\uc2a4\ud2b8: \uc9c4\ud55c \ucd08\ub85d (DarkGreen) */\n"
"}\n"
"\n"
"QProgressBar::chunk {\n"
"    background: qlineargradient(\n"
"        x1: 0, y1: 0, x2: 1, y2: 0,\n"
"        stop: 0 #00FF7F,\n"
"        stop: 1 #00CC66\n"
"    );\n"
"    border-radius: 0px;\n"
"    margin: 0px;\n"
"}")
        self.progressSearch.setMaximum(100)
        self.progressSearch.setValue(0)
        self.progressSearch.setTextVisible(False)

        self.formLayout.setWidget(3, QFormLayout.ItemRole.LabelRole, self.progressSearch)

        MainWindow.setCentralWidget(self.centralwidget)
        QWidget.setTabOrder(self.comboDatabase, self.treeSearchResults)
        QWidget.setTabOrder(self.treeSearchResults, self.treeDownloadList)
        QWidget.setTabOrder(self.treeDownloadList, self.buttonRemoveToDownloadList)
        QWidget.setTabOrder(self.buttonRemoveToDownloadList, self.buttonSelectDirectory)
        QWidget.setTabOrder(self.buttonSelectDirectory, self.lineEditDirectory)
        QWidget.setTabOrder(self.lineEditDirectory, self.comboFormat)
        QWidget.setTabOrder(self.comboFormat, self.buttonDownloadStart)

        self.retranslateUi(MainWindow)

        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"SRA to FASTQ Downloader", None))
        ___qtreewidgetitem = self.treeSearchResults.headerItem()
        ___qtreewidgetitem.setText(8, QCoreApplication.translate("MainWindow", u"Run Accession", None));
        ___qtreewidgetitem.setText(7, QCoreApplication.translate("MainWindow", u"Biosample", None));
        ___qtreewidgetitem.setText(6, QCoreApplication.translate("MainWindow", u"Bioproject", None));
        ___qtreewidgetitem.setText(5, QCoreApplication.translate("MainWindow", u"Library Strategy", None));
        ___qtreewidgetitem.setText(4, QCoreApplication.translate("MainWindow", u"Organism", None));
        ___qtreewidgetitem.setText(3, QCoreApplication.translate("MainWindow", u"Total Reads", None));
        ___qtreewidgetitem.setText(2, QCoreApplication.translate("MainWindow", u"Total Bases", None));
        ___qtreewidgetitem.setText(1, QCoreApplication.translate("MainWindow", u"Platform Model", None));
        ___qtreewidgetitem.setText(0, QCoreApplication.translate("MainWindow", u"Title", None));
        ___qtreewidgetitem1 = self.treeDownloadList.headerItem()
        ___qtreewidgetitem1.setText(8, QCoreApplication.translate("MainWindow", u"Run Accession", None));
        ___qtreewidgetitem1.setText(7, QCoreApplication.translate("MainWindow", u"Biosample", None));
        ___qtreewidgetitem1.setText(6, QCoreApplication.translate("MainWindow", u"Bioproject", None));
        ___qtreewidgetitem1.setText(5, QCoreApplication.translate("MainWindow", u"Library Strategy", None));
        ___qtreewidgetitem1.setText(4, QCoreApplication.translate("MainWindow", u"Organism", None));
        ___qtreewidgetitem1.setText(3, QCoreApplication.translate("MainWindow", u"Total Reads", None));
        ___qtreewidgetitem1.setText(2, QCoreApplication.translate("MainWindow", u"Total Bases", None));
        ___qtreewidgetitem1.setText(1, QCoreApplication.translate("MainWindow", u"Platform Model", None));
        ___qtreewidgetitem1.setText(0, QCoreApplication.translate("MainWindow", u"Title", None));
        self.buttonDownloadStart.setText(QCoreApplication.translate("MainWindow", u"Download Start", None))
        self.labelFormat.setText(QCoreApplication.translate("MainWindow", u"Select format:", None))
        self.comboFormat.setItemText(0, QCoreApplication.translate("MainWindow", u"FASTQ", None))
        self.comboFormat.setItemText(1, QCoreApplication.translate("MainWindow", u"FASTA", None))

        self.buttonSelectDirectory.setText(QCoreApplication.translate("MainWindow", u"Select Directory", None))
        self.buttonMoveToDownload.setText(QCoreApplication.translate("MainWindow", u"Move to Download List", None))
        self.buttonRemoveToDownloadList.setText(QCoreApplication.translate("MainWindow", u"Remove to Download List", None))
        self.buttonSelectAllinSearchList.setText(QCoreApplication.translate("MainWindow", u"Select all", None))
        self.labelSearchCount.setText(QCoreApplication.translate("MainWindow", u"Search List: 0 items     ", None))
        self.buttonSelectAllinSearchList_2.setText(QCoreApplication.translate("MainWindow", u"Select all", None))
        self.labelDownloadCount.setText(QCoreApplication.translate("MainWindow", u"Download List: 0 items     ", None))
        self.IDLlogo.setText("")
        self.Filteringlabel.setText(QCoreApplication.translate("MainWindow", u"Filtering on Organism", None))
        self.resetpushButton.setText(QCoreApplication.translate("MainWindow", u"Reset", None))
        self.labelDatabase.setText(QCoreApplication.translate("MainWindow", u"Select Database:", None))
        self.comboDatabase.setItemText(0, QCoreApplication.translate("MainWindow", u"SRA", None))
        self.comboDatabase.setItemText(1, QCoreApplication.translate("MainWindow", u"ENA", None))

        self.lineEditEmail_2.setText(QCoreApplication.translate("MainWindow", u"your Email for NCBI", None))
        self.buttonSaveEmail_2.setText(QCoreApplication.translate("MainWindow", u" Save Email ", None))
        self.lineEditQuery.setText(QCoreApplication.translate("MainWindow", u"Search Query", None))
        self.buttonSearch.setText(QCoreApplication.translate("MainWindow", u"Search", None))
        self.buttonLoadMore.setText(QCoreApplication.translate("MainWindow", u"Load More", None))
        self.SearchClearbutton.setText(QCoreApplication.translate("MainWindow", u"Clear search", None))
    # retranslateUi

