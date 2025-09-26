import subprocess
import os
import zipfile
import requests
import winreg
import ctypes
import sys
from PySide6.QtCore import QThread, Signal, QTimer, Qt
from PySide6.QtWidgets import QDialog, QVBoxLayout, QLabel, QProgressBar, QMessageBox, QPushButton, QApplication
from PySide6.QtGui import QCursor

def is_sratoolkit_installed(parent_dialog=None, on_success=None, toolkit_path=None):
    try:
        custom_env = os.environ.copy()
        if toolkit_path:
            custom_env["PATH"] += os.pathsep + toolkit_path

        subprocess.run(
            ["prefetch", "--help"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=custom_env
        )

        if parent_dialog:
            msg_box = QMessageBox(parent_dialog)
            msg_box.setIcon(QMessageBox.Information)
            msg_box.setWindowTitle("SRA Toolkit")
            msg_box.setText("✅ SRA Toolkit is already installed.")
            msg_box.setStandardButtons(QMessageBox.Ok)
            msg_box.setWindowFlag(Qt.WindowStaysOnTopHint)
            msg_box.setModal(True)

            ok_button = msg_box.button(QMessageBox.Ok)
            if ok_button:
                ok_button.setCursor(QCursor(Qt.PointingHandCursor))

            if msg_box.exec() == QMessageBox.Ok:
                if on_success:
                    QTimer.singleShot(100, on_success)
                QTimer.singleShot(200, parent_dialog.accept)

        return True

    except FileNotFoundError:
        if parent_dialog:
            parent_dialog.label.setText("❌ SRA Toolkit not found (file not found).")
            parent_dialog.progress.hide()
            parent_dialog.install_button.show()
        return False

    except subprocess.CalledProcessError:
        if parent_dialog:
            parent_dialog.label.setText("❌ SRA Toolkit not found (called process error).")
            parent_dialog.progress.hide()
            parent_dialog.install_button.show()
        return False


def is_admin():
    try:
        return ctypes.windll.shell32.IsUserAnAdmin() != 0
    except Exception:
        return False


def set_sratoolkit_env(toolkit_path):
    try:
        os.environ["PATH"] += os.pathsep + toolkit_path

        try:
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, r"Environment", 0, winreg.KEY_ALL_ACCESS) as reg_key:
                try:
                    current_path, _ = winreg.QueryValueEx(reg_key, "PATH")
                    current_paths = set(current_path.split(";"))
                except FileNotFoundError:
                    current_paths = set()

                if toolkit_path not in current_paths:
                    current_paths.add(toolkit_path)
                    new_path = ";".join(current_paths)
                    winreg.SetValueEx(reg_key, "PATH", 0, winreg.REG_EXPAND_SZ, new_path)
        except Exception:
            pass

        try:
            with winreg.OpenKey(
                winreg.HKEY_LOCAL_MACHINE,
                r"SYSTEM\CurrentControlSet\Control\Session Manager\Environment",
                0,
                winreg.KEY_ALL_ACCESS
            ) as reg_key:
                try:
                    current_path, _ = winreg.QueryValueEx(reg_key, "Path")
                    current_paths = set(current_path.split(";"))
                except FileNotFoundError:
                    current_paths = set()

                if toolkit_path not in current_paths:
                    current_paths.add(toolkit_path)
                    new_path = ";".join(current_paths)
                    winreg.SetValueEx(reg_key, "Path", 0, winreg.REG_EXPAND_SZ, new_path)
        except (PermissionError, Exception):
            pass

        send_environment_change_message()
        subprocess.run("taskkill /f /im explorer.exe", shell=True)
        subprocess.run("start explorer.exe", shell=True)

    except Exception:
        pass

def send_environment_change_message():
    HWND_BROADCAST = 0xFFFF
    WM_SETTINGCHANGE = 0x001A
    SMTO_ABORTIFHUNG = 0x0002
    ctypes.windll.user32.SendMessageTimeoutW(
        HWND_BROADCAST,
        WM_SETTINGCHANGE,
        0,
        "Environment",
        SMTO_ABORTIFHUNG,
        5000,
        None
    )

class SraToolkitInstaller(QThread):
    progress_update = Signal(str)
    installation_complete = Signal(bool, str)

    def run(self):
        toolkit_url = "https://ftp.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-win64.zip"
        zip_path = os.path.join(os.getcwd(), "sratoolkit.zip")
        install_root = os.path.join(os.getcwd(), "sratoolkit")

        try:
            self.progress_update.emit(" Downloading SRA Toolkit...")
            with requests.get(toolkit_url, stream=True) as response:
                response.raise_for_status()
                with open(zip_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)

            self.progress_update.emit(f"Extracting to {install_root}...")
            os.makedirs(install_root, exist_ok=True)
            with zipfile.ZipFile(zip_path, "r") as zip_ref:
                zip_ref.extractall(install_root)

            toolkit_dir = next((d for d in os.listdir(install_root) if d.startswith("sratoolkit")), None)
            if toolkit_dir:
                toolkit_path = os.path.join(install_root, toolkit_dir, "bin")
                set_sratoolkit_env(toolkit_path)
                if is_sratoolkit_installed():
                    self.installation_complete.emit(True, toolkit_path)
                else:
                    self.installation_complete.emit(False, "Toolkit installed, but prefetch is not recognized. A system reboot may be required.")
            else:
                self.installation_complete.emit(False, "Failed to locate extracted toolkit folder.")
        except Exception as e:
            self.installation_complete.emit(False, f"Error occurred during installation: {str(e)}")

class SraToolkitCheckDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SRA Toolkit Check")
        self.setFixedSize(420, 180)
        self.setStyleSheet("""
            QWidget{background-color: #F3F4F5; color: black;}
            QLabel {
                font-size: 13px;
                color: #2c3e50;
            }
            QProgressBar {
                border: 1px solid #00CC66;
                border-radius: 3px;
                background-color: #E9F0F4;
                height: 5px;
                text-align: center;
            }
            QProgressBar::chunk {
                background: qlineargradient(
                    x1: 0, y1: 0, x2: 1, y2: 0,
                    stop: 0 #00FF7F,
                    stop: 1 #00CC66
                );
                border-radius: 3px;
            }
            QPushButton {
                background-color: #1D83D5;
                color: white;
                font-weight: bold;
                border: 1px solid #40546C;
                border-radius: 8px;
                padding: 8px;
                font-size: 14px;
            }
            QPushButton:hover {
                background-color: #306999;
            }
        """)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(25, 25, 25, 25)
        layout.setSpacing(20)

        self.label = QLabel("Checking for SRA Toolkit installation...")
        self.label.setAlignment(Qt.AlignCenter)

        self.progress = QProgressBar()
        self.progress.setRange(0, 0)

        self.install_button = QPushButton("Install SRA Toolkit")
        self.install_button.setCursor(QCursor(Qt.PointingHandCursor))
        self.install_button.hide()
        self.install_button.clicked.connect(self.start_installation)

        layout.addWidget(self.label)
        layout.addWidget(self.progress)
        layout.addWidget(self.install_button)

        QTimer.singleShot(100, self.check_toolkit)

    def check_toolkit(self):
        if is_sratoolkit_installed(self):
            self.label.setText("✅ SRA Toolkit is already installed.")
            QTimer.singleShot(300, self.accept)
        else:
            self.label.setText("❌ SRA Toolkit is not installed.")
            self.progress.hide()
            self.install_button.show()

    def start_installation(self):
        self.install_button.setEnabled(False)
        self.progress.show()
        self.label.setText("⬇ Installing SRA Toolkit...")

        self.installer = SraToolkitInstaller()
        self.installer.progress_update.connect(self.label.setText)
        self.installer.installation_complete.connect(self.on_installation_complete)
        self.installer.start()

    def on_installation_complete(self, success, path_or_msg):
        if success and is_sratoolkit_installed(self):
            self.label.setText("✅ Installation successful. Restarting Explorer...")
            self.progress.hide()
            QTimer.singleShot(500, self.restart_explorer_and_notify)
        else:
            msg = QMessageBox(QMessageBox.Critical, "Installation Failed", path_or_msg, parent=self)
            msg.setWindowFlag(Qt.WindowStaysOnTopHint)
            msg.exec()
            self.reject()

    def restart_explorer_and_notify(self):
        try:
            subprocess.run("taskkill /f /im explorer.exe", shell=True)
            subprocess.run("start explorer.exe", shell=True)
            print("🔄 Explorer restarted.")
    
            msg = QMessageBox(parent=self)
            msg.setIcon(QMessageBox.Information)
            msg.setWindowTitle("SRA Toolkit Installation Complete")
            msg.setText("✅ SRA Toolkit has been successfully installed.")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.setWindowFlag(Qt.WindowStaysOnTopHint)
            msg.setModal(True)
    
            ok_button = msg.button(QMessageBox.Ok)
            ok_button.setCursor(QCursor(Qt.PointingHandCursor))
    
            if msg.exec() == QMessageBox.Ok:
                self.accept()
        except Exception as e:
            print(f"⚠️ Explorer restart failed: {e}")
            self.reject()
