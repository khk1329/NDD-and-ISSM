import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from NGS_data_downloader.sratoolkit_installer import SraToolkitCheckDialog
from NGS_data_downloader.NGS_data_downloader import MainApp

from PySide6.QtWidgets import QApplication, QDialog
import urllib.request
import tempfile

def _download_cacert(dest_path: str) -> str:
    url = "https://curl.se/ca/cacert.pem"
    with urllib.request.urlopen(url, timeout=30) as r, open(dest_path, "wb") as f:
        f.write(r.read())
    return dest_path


def ensure_ca_bundle() -> None:
    try:
        with urllib.request.urlopen("https://eutils.ncbi.nlm.nih.gov", timeout=5) as _:
            return
    except Exception:
        pass

    ca_path = None
    try:
        import certifi  
        ca_path = certifi.where()
    except Exception:
        cache_dir = os.path.join(os.path.expanduser("~"), ".ndd_certs")
        os.makedirs(cache_dir, exist_ok=True)
        ca_path = os.path.join(cache_dir, "cacert.pem")
        need_download = (
            not os.path.exists(ca_path)
            or os.path.getsize(ca_path) < 100_000  
        )
        if need_download:
            ca_path = _download_cacert(ca_path)

    os.environ["SSL_CERT_FILE"] = ca_path
    os.environ["REQUESTS_CA_BUNDLE"] = ca_path
    os.environ["CURL_CA_BUNDLE"] = ca_path

app = QApplication.instance()
if not app:
    app = QApplication(sys.argv)

main_window = None

def start_gui_after_sratoolkit_check():
    global main_window
    dialog = SraToolkitCheckDialog()

    if dialog.exec() == QDialog.Accepted:
        ensure_ca_bundle()

        main_window = MainApp()
        main_window.setWindowTitle("NGS Data Downloader")
        main_window.show()

    app.exec()

if __name__ == "__main__":
    start_gui_after_sratoolkit_check()
