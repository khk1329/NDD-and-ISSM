import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from NGS_data_downloader.sratoolkit_installer import SraToolkitCheckDialog
from NGS_data_downloader.NGS_data_downloader import MainApp

from PySide6.QtWidgets import QApplication, QDialog

app = QApplication.instance()
if not app:
    app = QApplication(sys.argv)

main_window = None

def start_gui_after_sratoolkit_check():
    global main_window
    dialog = SraToolkitCheckDialog()

    if dialog.exec() == QDialog.Accepted:
        main_window = MainApp()
        main_window.setWindowTitle("NGS Data Downloader")
        main_window.show()

    app.exec()

if __name__ == "__main__":
    start_gui_after_sratoolkit_check()
