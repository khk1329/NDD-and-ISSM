import sys
import os
import subprocess
import urllib.request

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from NGS_data_downloader.sratoolkit_installer import SraToolkitCheckDialog
from NGS_data_downloader.NGS_data_downloader import MainApp

from PySide6.QtWidgets import QApplication, QDialog


# ---------------------------
# CA bundle helpers
# ---------------------------
def _download_cacert(dest_path: str) -> str:
    """curl에서 배포하는 최신 Mozilla CA 번들을 내려받아 dest_path에 저장."""
    url = "https://curl.se/ca/cacert.pem"
    with urllib.request.urlopen(url, timeout=30) as r, open(dest_path, "wb") as f:
        f.write(r.read())
    return dest_path


def ensure_ca_bundle() -> str:
    """
    루트 CA 번들을 확보하고 관련 환경변수를 설정한 뒤, 최종 CA 파일 경로를 반환.
    certifi 우선 사용, 없으면 cacert.pem을 사용자 홈 하위에 저장.
    """
    # certifi 우선
    try:
        import certifi  # type: ignore
        ca_path = certifi.where()
    except Exception:
        cache_dir = os.path.join(os.path.expanduser("~"), ".ndd_certs")
        os.makedirs(cache_dir, exist_ok=True)
        ca_path = os.path.join(cache_dir, "cacert.pem")
        if not os.path.exists(ca_path) or os.path.getsize(ca_path) < 100_000:
            ca_path = _download_cacert(ca_path)

    # 환경변수 (Python/requests/biopython + SRA curl 공통)
    os.environ["SSL_CERT_FILE"] = ca_path
    os.environ["REQUESTS_CA_BUNDLE"] = ca_path
    os.environ["CURL_CA_BUNDLE"] = ca_path
    os.environ["SSL_CERT_DIR"] = os.path.dirname(ca_path)
    return ca_path


def configure_sra_toolkit_cafile(ca_path: str) -> None:
    """
    vdb-config에 CA 파일/디렉토리를 직접 기록하여 SRA Toolkit이 확실히 사용하도록 함.
    일부 환경에서 CURL_CA_BUNDLE만으론 무시될 수 있어 추가로 강제 주입.
    """
    if not ca_path or not os.path.exists(ca_path):
        return
    try:
        subprocess.run(
            ["vdb-config", "-s", f"/http/ca-file={ca_path}"],
            check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        subprocess.run(
            ["vdb-config", "-s", f"/http/ca-path={os.path.dirname(ca_path)}"],
            check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        # 프록시 환경변수 사용을 허용(있으면 자동 감지)
        subprocess.run(
            ["vdb-config", "-s", "/http/proxy/env=1"],
            check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
    except Exception:
        # 필요 시 print로 바꿔 로그 확인 가능
        pass


# ---------------------------
# QApplication helper
# ---------------------------
def get_app() -> QApplication:
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    return app


# ---------------------------
# Main flow
# ---------------------------
def start_gui_after_sratoolkit_check() -> None:
    app = get_app()  # ✅ 항상 여기서 보장
    dialog = SraToolkitCheckDialog()

    if dialog.exec() == QDialog.Accepted:
        ca_path = ensure_ca_bundle()          # ✅ CA 확보 + env 세팅
        configure_sra_toolkit_cafile(ca_path) # ✅ SRA Toolkit에 직접 주입

        main_window = MainApp()
        main_window.setWindowTitle("NGS Data Downloader")
        main_window.show()

    app.exec()


if __name__ == "__main__":
    start_gui_after_sratoolkit_check()
