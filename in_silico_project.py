import sys
import os
import traceback
import psutil
import ctypes
import multiprocessing

project_root = os.path.dirname(os.path.abspath(__file__))
scripts_path = os.path.join(project_root, "scripts")
if scripts_path not in sys.path:
    sys.path.insert(0, scripts_path)

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QPushButton, QSizePolicy, QApplication
)
from PySide6.QtGui import QCursor
from PySide6.QtCore import Qt

from in_silico_sequence_mining.in_silico_sequence_mining import file_3_GUI
from NGS_data_downloader.launcher import start_gui_after_sratoolkit_check

def log_error(error_message):
    with open("error_log.txt", "a", encoding="utf-8") as log_file:
        log_file.write(error_message + "\n")

def apply_fixed_memory_limit_by_ratio(ratio=0.9):
    if os.name == "nt":
        total_mem = psutil.virtual_memory().total
        memory_limit_bytes = int(total_mem * ratio)

        min_memory = min(memory_limit_bytes // 2, 2 * 1024**3)
        max_memory = memory_limit_bytes

        kernel32 = ctypes.windll.kernel32
        handle = kernel32.OpenProcess(0x1F0FFF, False, os.getpid())
        success = kernel32.SetProcessWorkingSetSize(
            handle,
            ctypes.c_size_t(min_memory),
            ctypes.c_size_t(max_memory)
        )

        if not success:
            raise ctypes.WinError()

class LauncherWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("In silico Project")
        self.setFixedSize(350, 200)
        self.setWindowFlags(self.windowFlags() & ~Qt.WindowMaximizeButtonHint | Qt.WindowCloseButtonHint)
        self.setStyleSheet("""
            QWidget {
                background-color: #F3F4F5;
                color: #2c3e50;
                font-family: 'Segoe UI', sans-serif;
            }
            QPushButton {
                background-color: #1D83D5;
                color: white;
                font-weight: bold;
                border: 1px solid #40546C;
                border-radius: 8px;
                padding: 8px;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #306999;
            }
        """)

        central_widget = QWidget(self)
        layout = QVBoxLayout(central_widget)
        layout.setSpacing(25)
        layout.setContentsMargins(40, 40, 40, 40)
        layout.setAlignment(Qt.AlignCenter)

        ngs_button = QPushButton("NGS Data Downloader")
        ngs_button.setCursor(QCursor(Qt.PointingHandCursor))
        ngs_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        ngs_button.setMinimumHeight(40)
        ngs_button.clicked.connect(start_gui_after_sratoolkit_check)

        silico_button = QPushButton("In Silico Sequence Mining")
        silico_button.setCursor(QCursor(Qt.PointingHandCursor))
        silico_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        silico_button.setMinimumHeight(40)
        silico_button.clicked.connect(file_3_GUI)

        layout.addWidget(ngs_button)
        layout.addWidget(silico_button)

        self.setCentralWidget(central_widget)

if __name__ == "__main__":
    try:
        multiprocessing.freeze_support()
        apply_fixed_memory_limit_by_ratio(ratio=0.9)

        app = QApplication.instance()
        if not app:
            app = QApplication(sys.argv)

        window = LauncherWindow()
        window.show()
        sys.exit(app.exec())

    except FileNotFoundError as fnf_error:
        log_error(str(fnf_error))

    except ImportError as imp_error:
        log_error(str(imp_error))

    except Exception:
        log_error(traceback.format_exc())
        raise
