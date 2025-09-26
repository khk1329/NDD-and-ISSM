import os
import requests
from Bio import Entrez
import subprocess
import threading
from datetime import datetime
import time
import csv
import json
import shutil
import sys
import stat
import tempfile 
import gzip
from concurrent.futures import ThreadPoolExecutor
import re
import ctypes
from urllib.parse import quote_plus
from collections import Counter
from PySide6.QtWidgets import (
    QWidget, QLabel, QPushButton, QVBoxLayout, QHBoxLayout, QTreeWidgetItem, QDockWidget, QListWidget, QListWidgetItem, QVBoxLayout, QWidget, QLineEdit, QSplitter,
    QMessageBox, QProgressBar, QFileDialog, QTreeWidget, QDialog, QApplication, QMainWindow, QAbstractItemView, QSizePolicy, QHeaderView
)
from PySide6.QtCore import Qt, QThread, Signal, QObject, QTimer
from PySide6.QtGui import QCursor
from functools import partial

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from ui.NGS_downloader import Ui_MainWindow
from queue import Queue, Empty

def get_resource_path(relative_path):
    if getattr(sys, 'frozen', False):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)

stop_flag = False

prefetch_retry_count = 3

log_file_name = None
ena_search_results = []
processes = []  
completed_items = []
failed_files = []
active_threads = []
log_lock = threading.Lock()  

def is_url_valid(url):
    try:
        response = requests.head(url)
        return response.status_code == 200
    except requests.RequestException:
        return False
    
def is_admin():
    try:
        return ctypes.windll.shell32.IsUserAnAdmin() != 0  
    except AttributeError:
        return False  
    
def format_number(value):
    try:
        return "{:,}".format(int(value))
    except ValueError:
        return "N/A"

class DownloadStatusWindow(QWidget):
    def __init__(self, download_items, controller, output_dir, file_format, start_button=None, database="SRA", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Download Status")
        self.setMinimumSize(800, 600)

        self.setStyleSheet("background-color: #F9FAFA;")

        self.download_items = download_items
        self.controller = controller
        self.output_dir = output_dir
        self.file_format = file_format
        self.start_button = start_button
        self.database = database
        self._canceled = False
        
        self.status_bars = {}
        self.status_labels = {}
        
        self.progress_queue = Queue()

        layout = QVBoxLayout(self)
        
        self.completed_count = 0
        self.failed_files = []      
        self.total_count = len(self.download_items)
        
        self._updated_files = set()
        self._completion_shown = False
              
        self.tree = QTreeWidget()
        self.tree.headerItem().setTextAlignment(0, Qt.AlignCenter)
        self.tree.headerItem().setTextAlignment(1, Qt.AlignCenter)
        self.tree.headerItem().setTextAlignment(2, Qt.AlignCenter)
        self.tree.setColumnCount(3)
        self.tree.setHeaderLabels(["Run Accession", "Progress", "Status"])
        self.tree.setRootIsDecorated(False)
        self.tree.setSelectionMode(QTreeWidget.NoSelection)
        self.tree.setFocusPolicy(Qt.NoFocus)
        self.tree.setStyleSheet("""
            QTreeWidget {
                font-size: 12px;
                color: #2c3e50;
                background-color: white;
                border: 1px solid #ccc;
            }
            QTreeWidget::item {
                padding: 8px;
                border-bottom: 1px dashed #ddd;
            }
            QHeaderView::section {
                background-color: #f0f0f0;
                color: #2c3e50;
                padding: 6px;
                font-weight: bold;
                border: none;
            }
            QScrollBar:vertical {
                background: transparent;
                width: 6px;
                margin: 2px 0 2px 0;
            }
        
        QScrollBar::handle:vertical {
                background: #9E9E9E;
                border-radius: 3px;
                min-height: 20px;
            }
        
        QScrollBar::add-line:vertical,
        QScrollBar::sub-line:vertical {
                height: 0px;
                background: none;
            }
        
        QScrollBar::add-page:vertical,
        QScrollBar::sub-page:vertical {
                background: none;
            }
        """)
        self.tree.header().setStretchLastSection(False)
        self.tree.header().setSectionResizeMode(0, QHeaderView.Fixed)
        self.tree.header().resizeSection(0, 150)
        self.tree.header().setSectionResizeMode(1, QHeaderView.Fixed)
        self.tree.header().resizeSection(1, 300)
        self.tree.header().setSectionResizeMode(2, QHeaderView.Stretch)
        self.tree.header().resizeSection(2, 300)

        layout.addWidget(self.tree)

        for acc, metadata in self.download_items:
            item = QTreeWidgetItem(self.tree)
            item.setText(0, acc)
            item.setText(2, "Waiting...")
            item.setTextAlignment(0, Qt.AlignCenter)

            bar = QProgressBar()
            bar.setValue(0)
            bar.setMaximum(100)
            bar.setFixedHeight(16)
            bar.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
            bar.setStyleSheet("""
                QProgressBar {
                    background-color: #E4E8EB;
                    border: none;
                    border-radius: 4px;
                    text-align: center;
                    font-size: 11px;
                    color: white;
                }
                QProgressBar::chunk {
                    background-color: #1D83D5;
                    border-radius: 4px;
                }
            """)
            self.tree.setItemWidget(item, 1, bar)
            self.status_bars[acc] = bar
            self.status_labels[acc] = item

        bottom_layout = QHBoxLayout()

        self.cancel_button = QPushButton("Cancel All")
        self.cancel_button.setFixedWidth(100)
        self.cancel_button.setCursor(Qt.PointingHandCursor)
        self.cancel_button.clicked.connect(self.handle_close)
        self.cancel_button.setStyleSheet("""
            QPushButton {
                background-color: #E53935;
                color: white;
                font-weight: bold;
                border: 1px solid #742220;
                border-radius: 6px;
                padding: 6px 12px;
            }
            QPushButton:hover {
                background-color: #742220;
            }
        """)
        bottom_layout.addWidget(self.cancel_button, alignment=Qt.AlignLeft)
        
        bottom_layout.addStretch(1)
        
        self.completed_label = QLabel(f"Completed: 0 / {self.total_count}  ")
        self.completed_label.setStyleSheet("font-size: 12px; color: #555;")
        bottom_layout.addWidget(self.completed_label, alignment=Qt.AlignRight)
        
        layout.addLayout(bottom_layout)
        
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.process_progress_queue)
        self.timer.start(200)

        self.start_downloads()

    def process_progress_queue(self):
        try:
            while True:
                data = self.progress_queue.get_nowait()
                acc = data.get("acc")
                status = data.get("status")
                progress = data.get("progress")

                if acc in self.status_labels:
                    self.status_labels[acc].setText(2, status)
                if acc in self.status_bars:
                    self.status_bars[acc].setValue(progress)
        except Empty:
            pass

    def start_downloads(self):
        self.max_concurrent = 4
        self.pending_queue = []
        self.active_workers = []
    
        for i, (acc, metadata) in enumerate(self.download_items):
    
            if self.database == "SRA":
                worker = SraDownloadWorker(
                    run_acc=acc,
                    metadata=metadata,
                    output_dir=self.output_dir,
                    selected_format=self.file_format,
                    progress_queue=self.progress_queue
                )
            else:
                worker = EnaDownloadWorker(
                    run_acc=acc,
                    metadata=metadata,
                    output_dir=self.output_dir,
                    progress_queue=self.progress_queue
                )
    
            worker.error.connect(lambda message, acc=acc: self.handle_worker_error(acc, message))
            worker.done.connect(lambda acc=acc: self.mark_worker_finished(acc))
            worker.done.connect(self.on_worker_done)
    
            self.pending_queue.append(worker)
    
        self.start_next_batch()

    def start_next_batch(self):
        while len(self.active_workers) < self.max_concurrent and self.pending_queue:
            worker = self.pending_queue.pop(0)
            self.active_workers.append(worker)
            worker.start()

    def on_worker_done(self):
        self.active_workers = [w for w in self.active_workers if w.isRunning()]
        self.start_next_batch()

    def handle_worker_error(self, acc, message, *args):
        self.failed_files.append(acc)
    
        if self.progress_queue:
            self.progress_queue.put({
                "acc": acc,
                "status": message,
                "progress": 0
            })
    
        self.mark_worker_finished(acc)
        self.on_worker_done()
        
    def mark_worker_finished(self, acc, *_):
        self.completed_count += 1
        self.completed_label.setText(f"Completed: {self.completed_count} / {self.total_count}  ")
    
        if self.completed_count == self.total_count and not self._completion_shown:
            self._completion_shown = True
            if self.start_button:
                self.start_button.setEnabled(True)
    
            if self._canceled:
                msg = "🛑 Download process was canceled."
            elif self.failed_files:
                failed_list = ", ".join(self.failed_files)
                msg = f"❌ Some files failed to download:\n{failed_list}"
            else:
                msg = "✅ All files downloaded successfully."
    
            msg_box = CustomMessageBox(title="Download Summary", message=msg, parent=self)
            msg_box.setWindowFlags(msg_box.windowFlags() | Qt.WindowStaysOnTopHint)
            msg_box.exec()
    
            self.force_close()
            
    def handle_close(self):
        confirm = CustomConfirmBox(
            title="Cancel Downloads",
            message="Do you really want to cancel?",
            on_close=self.mark_canceled,
            parent=self
        )
        confirm.setWindowFlags(confirm.windowFlags() | Qt.WindowStaysOnTopHint)
        confirm.exec()

    def mark_canceled(self):
        self._canceled = True
        
        for worker in self.active_workers:
            worker.stop()
            
            if self.progress_queue:
                self.progress_queue.put({
                    "acc": worker.run_acc,
                    "status": "Canceling...",
                    "progress": 0
                })
                
        for worker in self.pending_queue:
            if self.progress_queue:
                self.progress_queue.put({
                    "acc": worker.run_acc,
                    "status": "Canceled (Not started)",
                    "progress": 0
                })
            self.mark_worker_finished(worker.run_acc)
    
        self.pending_queue.clear()
        
    def force_close(self):
        self._force_closing = True
        self.close()
    
    def closeEvent(self, event):
        if getattr(self, '_force_closing', False):
            event.accept()
        else:
            self.handle_close()
            event.ignore()
                           
class SraDownloadWorker(QThread):
    progress = Signal(int)
    status = Signal(str)
    error = Signal(str)
    done = Signal(str)
    
    def __init__(self, run_acc, metadata, output_dir, selected_format="FASTQ", progress_queue=None):
        super().__init__()
        self.run_acc = run_acc
        self.metadata = metadata
        self.output_dir = output_dir
        self.selected_format = selected_format
        self.progress_queue = progress_queue 
        self.success = False
        self._stopped = False 
        self._output_created = False
        self._created_out_dir = None
        
    def stop(self):
        self._stopped = True
        
    def _put_progress(self, status, progress):
        if self.progress_queue:
            self.progress_queue.put({
                "acc": self.run_acc,
                "status": status,
                "progress": progress
            })

    def cleanup_temp_files(self, temp_path, temp_sra_root, run_tmp):
        if temp_path and os.path.exists(temp_path): shutil.rmtree(temp_path, ignore_errors=True)
        if temp_sra_root and os.path.exists(temp_sra_root): shutil.rmtree(temp_sra_root, ignore_errors=True)
        if run_tmp and os.path.exists(run_tmp): shutil.rmtree(run_tmp, ignore_errors=True)

    def _cleanup_out_dir_if_created(self):
        if not self._output_created or not self._created_out_dir:
            return
        try:
            if os.path.isdir(self._created_out_dir):
                shutil.rmtree(self._created_out_dir, ignore_errors=True)
        except Exception:
            pass
        finally:
            self._output_created = False
            self._created_out_dir = None

    def run(self):
        temp_path = None
        temp_sra_root = None
        run_tmp = None
        
        start_dt = datetime.now()
        start_str = start_dt.strftime("%Y-%m-%d %H:%M:%S")
        end_dt = None
        end_str = None
        elapsed_str = None
        
        try:
            if self._stopped:
                self._put_progress("❌ Canceled", 0)
                return
    
            self._put_progress("Prefetching SRA...", 10)
    
            run_tmp = os.path.join("C:/tmp_ndd", self.run_acc)
            try:
                shutil.rmtree(run_tmp, ignore_errors=True)
                os.makedirs(run_tmp, exist_ok=True)
            except Exception as e:
                self._put_progress(f"❌ Failed to prepare temp folder: {e}", 0)
                self.error.emit(f"Temp folder creation failed for {self.run_acc}: {e}")
                return
    
            try:
                temp_path = tempfile.mkdtemp(prefix="tmp_", dir=run_tmp)
            except Exception as e:
                self._put_progress(f"❌ Failed to create temp path: {e}", 0)
                self.error.emit(f"Temp path creation failed for {self.run_acc}: {e}")
                return
    
            temp_sra_root = os.path.join("C:/NDD_dummy_file", self.run_acc)
            os.makedirs(temp_sra_root, exist_ok=True)
    
            prefetch_worker = PrefetchWorker(
                run_acc_list=[self.run_acc],
                prefetch_output_dir=temp_sra_root,
                max_workers=1
            )
            prefetch_worker.run()
            
            if self._stopped:
                self._put_progress("❌ Canceled", 0)
                return
    
            candidate_paths = [
                os.path.join(temp_sra_root, self.run_acc, f"{self.run_acc}.sra"),
                os.path.join(temp_sra_root, f"{self.run_acc}.sra")
            ]
            sra_path = next((p for p in candidate_paths if os.path.exists(p)), None)

            if self._stopped:
                self._put_progress("❌ Canceled", 0)
                return

            if not sra_path:
                self._put_progress("❌ .sra file not found", 0)
                self.error.emit(f".sra file not found for {self.run_acc}")
                return
    
            self._put_progress("Downloading...", 50)
    
            out_dir = os.path.join(self.output_dir, self.run_acc)
            os.makedirs(out_dir, exist_ok=True)
            self._output_created = True
            self._created_out_dir = out_dir
    
            if self.selected_format.upper() == "FASTA":
                cmd = f'fastq-dump "{sra_path}" --outdir "{out_dir}" --split-files --fasta 0'
                fallback_cmd = None
            else:
                cmd = f'fasterq-dump "{sra_path}" -O "{out_dir}" --split-files --temp "{temp_path}"'
                fallback_cmd = f'fastq-dump "{sra_path}" -O "{out_dir}" --split-files --gzip'
                
            runner = ConversionRunner()
            runner.log.connect(lambda msg: self.status.emit(msg))
            
            retry_limit = 3
            for attempt in range(1, retry_limit + 1):
                self.success = runner.run(cmd, fallback_cmd)
                if self.success:
                    break
                else:
                    self.status.emit(f"🔁 Retry {attempt} failed for {self.run_acc}")
                    time.sleep(3)
            
            end_dt = datetime.now()
            end_str = end_dt.strftime("%Y-%m-%d %H:%M:%S")
            elapsed_str = str(end_dt - start_dt).split(".")[0]
            
            if self._stopped:
                self._put_progress("❌ Canceled", 0)
                self._cleanup_out_dir_if_created()
                return
            
            if not self.success:
                self._put_progress("❌ Conversion failed", 0)
                self.error.emit(f"Conversion failed for {self.run_acc}")
                return
    
            if self.selected_format.upper() == "FASTQ" and "fasterq-dump" in cmd:
                for fname in os.listdir(out_dir):
                    
                    if self._stopped:
                        return
                    
                    if fname.endswith(".fastq"):
                        try:
                            src = os.path.join(out_dir, fname)
                            dst = os.path.join(out_dir, fname + ".gz")
                            with open(src, "rb") as f_in, gzip.open(dst, "wb") as f_out:
                                shutil.copyfileobj(f_in, f_out)
                            os.remove(src)
                        except Exception as e:
                            self.status.emit(f"Compression error: {e}")
    
            try:
                log_metadata(
                    run_acc=self.run_acc,
                    output_dir=self.output_dir,
                    metadata=self.metadata or {},
                    elapsed_time=elapsed_str,
                    start_time=start_str,
                    end_time=end_str
                )
            except Exception as e:
                print(f"[WARN] log_metadata failed for {self.run_acc}: {e}")
    
            if self._stopped:
                self._put_progress("❌ Canceled", 0)
                self._cleanup_out_dir_if_created()
                return

            self._put_progress("✅ Download Completed", 100)
            self.progress.emit(1)
    
        except Exception as e:
            self._put_progress(f"❌ Error: {e}", 0)
            self._cleanup_out_dir_if_created()
            self.error.emit(str(e))
            
        finally:
            self.cleanup_temp_files(temp_path, temp_sra_root, run_tmp)
            self.done.emit(self.run_acc)
            
class EnaDownloadWorker(QThread):
    progress = Signal(int)
    status = Signal(str)
    error = Signal(str)
    done = Signal(str)
    
    def __init__(self, run_acc, metadata, output_dir, progress_queue=None):
        super().__init__()
        self.run_acc = run_acc
        self.metadata = metadata
        self.output_dir = output_dir
        self.progress_queue = progress_queue
        self._stopped = False
        self.success = False
        self._emitted_done = False
        self._output_created = False
        self._created_dir = None 

    def stop(self):
        self._stopped = True
        
    def _emit_done_once(self):
        if not self._emitted_done:
            self._emitted_done = True
            self.done.emit(self.run_acc)

    def _put_progress(self, status, progress):
        if self.progress_queue:
            self.progress_queue.put({
                "acc": self.run_acc,
                "status": status,
                "progress": progress
            })

    def _cleanup_output_dir(self):
        if not self._output_created or not self._created_dir:
            return
        try:
            if os.path.isdir(self._created_dir):
                def _onerror(func, path, exc_info):
                    try:
                        os.chmod(path, stat.S_IWRITE)
                        func(path)
                    except Exception:
                        pass
                for _ in range(2):
                    try:
                        shutil.rmtree(self._created_dir, onerror=_onerror)
                        break
                    except PermissionError:
                        time.sleep(0.1)
                    except FileNotFoundError:
                        break
        except Exception:
            pass
        finally:
            self._output_created = False
            self._created_dir = None

    def run(self):
        start_dt = datetime.now()
        start_str = start_dt.strftime("%Y-%m-%d %H:%M:%S")
        end_dt = None
        elapsed_str = None

        final_output_dir = None
        try:
            if self._stopped:
                self._put_progress("❌ Canceled", 0)
                self._cleanup_output_dir()
                return

            self._put_progress("Fetching ENA file info...", 5)

            api_url = (
                f"https://www.ebi.ac.uk/ena/portal/api/filereport"
                f"?accession={self.run_acc}&result=read_run&fields=fastq_ftp"
            )
            response = requests.get(api_url, timeout=60)
            response.raise_for_status()

            lines = response.text.strip().split("\n")
            if len(lines) < 2:
                raise FileNotFoundError(f"No FASTQ files found for {self.run_acc} in ENA.")

            fastq_ftp_urls = lines[1].strip().split("\t", 1)[-1]
            fastq_files = fastq_ftp_urls.split(";") if ";" in fastq_ftp_urls else [fastq_ftp_urls]
            fastq_files = [u.strip() for u in fastq_files if u.strip()]

            if not fastq_files:
                raise FileNotFoundError(f"No FASTQ files found for {self.run_acc} in ENA.")

            final_output_dir = os.path.abspath(os.path.join(self.output_dir, self.run_acc))
            os.makedirs(final_output_dir, exist_ok=True)
            self._output_created = True 
            self._created_dir = final_output_dir

            for file_url in fastq_files:
                if self._stopped:
                    self._put_progress("❌ Canceled", 0)
                    self._cleanup_output_dir()
                    return

                https_url = f"https://{file_url}"
                file_name = os.path.basename(file_url)
                output_file_path = os.path.join(final_output_dir, file_name)

                self._put_progress("Downloading...", 50)

                success = False
                for attempt in range(1, 4):
                    try:
                        with requests.get(https_url, stream=True, timeout=60) as r:
                            r.raise_for_status()
                            with open(output_file_path, "wb") as f:
                                for chunk in r.iter_content(chunk_size=8192):
                                    if self._stopped:
                                        self._put_progress("❌ Canceled", 0)
                                        self._cleanup_output_dir()
                                        return
                                    if chunk:
                                        f.write(chunk)
                        success = True
                        break
                    except Exception:
                        if self._stopped:
                            self._put_progress("❌ Canceled", 0)
                            self._cleanup_output_dir()
                            return
                        self._put_progress(f"Retry {attempt} failed for {file_name}", 40)
                        time.sleep(5)

                if not success:
                    self._put_progress("❌ Failed to download", 100)
                    self.error.emit(f"Failed to download {file_name} after 3 attempts.")
                    self._cleanup_output_dir()
                    return

            end_dt = datetime.now()
            elapsed_str = str(end_dt - start_dt).split(".")[0]
            total_size_mb = (
                sum(os.path.getsize(os.path.join(final_output_dir, f)) for f in os.listdir(final_output_dir))
                / (1024 * 1024)
                if os.path.exists(final_output_dir) else 0
            )

            try:
                log_metadata(
                    run_acc=self.run_acc,
                    output_dir=self.output_dir,
                    metadata=self.metadata or {},
                    elapsed_time=elapsed_str,
                    start_time=start_str,
                    end_time=end_dt.strftime("%Y-%m-%d %H:%M:%S"),
                    total_size=total_size_mb
                )
            except Exception as e:
                print(f"[WARN] log_metadata failed for {self.run_acc}: {e}")

            if self._stopped:
                self._put_progress("❌ Canceled", 0)
                self._cleanup_output_dir()
                return

            self.success = True
            self._put_progress("✅ Download Completed", 100)
            self.progress.emit(1)

        except Exception as e:
            self._put_progress(f"❌ Error: {e}", 0)
            self.error.emit(f"Exception for {self.run_acc}: {e}")
            self._cleanup_output_dir()

        finally:
            self._emit_done_once()
            
class CustomMessageBox(QDialog):
    def __init__(self, message, title="Notification", parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setWindowFlags(Qt.Dialog | Qt.WindowTitleHint | Qt.CustomizeWindowHint | Qt.WindowStaysOnTopHint)
        self.setAttribute(Qt.WA_StyledBackground, True)
        self.setFixedSize(320, 120)

        self.setStyleSheet("""
            QWidget {
                background-color: #F3F4F5; color: black;
            }
            QLabel {
                font-size: 14px;
                color: #2c3e50;
                padding: 4px;
            }
            QPushButton {
                background-color: #1D83D5;
                color: white;
                font-weight: bold;
                border: none;
                border-radius: 6px;
                padding: 6px 14px;
                font-size: 12px; 
            }
            QPushButton:hover {
                background-color: #306999;
            }
        """)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 10)
        layout.setSpacing(15)

        label = QLabel(message)
        label.setAlignment(Qt.AlignCenter)
        layout.addWidget(label)

        button_layout = QHBoxLayout()
        button_layout.addStretch()

        ok_button = QPushButton("Ok")
        ok_button.setCursor(QCursor(Qt.PointingHandCursor))
        ok_button.clicked.connect(self.accept)
        button_layout.addWidget(ok_button)

        button_layout.addStretch()
        layout.addLayout(button_layout)

    def showEvent(self, event):
        self.raise_()
        self.activateWindow()
        self.setWindowState(self.windowState() & ~Qt.WindowMinimized | Qt.WindowActive)
        QApplication.alert(self, 3000)
        super().showEvent(event)

    @staticmethod
    def warning(parent=None, title="Warning", message="Warning occurred."):
        dlg = CustomMessageBox(message, title, parent)
        dlg.exec()

    @staticmethod
    def info(parent=None, title="Info", message="Information"):
        dlg = CustomMessageBox(message, title, parent)
        dlg.exec()

    @staticmethod
    def error(parent=None, title="Error", message="An error occurred."):
        dlg = CustomMessageBox(message, title, parent)
        dlg.exec()

class CustomConfirmBox(QDialog):
    def __init__(self, title, message, on_close=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setFixedSize(300, 100)
        self.on_close = on_close
        self.setWindowFlags(Qt.Dialog | Qt.WindowTitleHint | Qt.CustomizeWindowHint)
        self.setAttribute(Qt.WA_StyledBackground, True)

        self.result = False

        self.setStyleSheet("""
            QWidget {
                background-color: #F3F4F5; color: black;
            }
            QLabel {
                background-color: #F3F4F5;
                font-size: 14px;
                color: black;
            }
            QPushButton {
                background-color: #1D83D5;
                color: white;
                border: none;
                border-radius: 6px;
                padding: 6px 14px;
                font-weight: bold;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #1565C0;
            }
        """)

        layout = QVBoxLayout(self)
        label = QLabel(message)
        label.setWordWrap(True)
        label.setAlignment(Qt.AlignCenter)
        layout.addWidget(label)

        button_layout = QHBoxLayout()
        yes_btn = QPushButton("Yes")
        no_btn = QPushButton("No")
        yes_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        no_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))

        yes_btn.clicked.connect(self.accept) 
        no_btn.clicked.connect(self.reject)

        button_layout.addStretch()
        button_layout.addWidget(yes_btn)
        button_layout.addWidget(no_btn)
        button_layout.addStretch()

        layout.addLayout(button_layout)

    def accept(self):
        if self.on_close:
            self.on_close()
        super().accept()

class SraSearchWorker(QThread):
    result_ready = Signal(list)
    error_occurred = Signal(str)

    def __init__(self, query, email, retmax=9999, retstart=0):
        super().__init__()
        self.query = query
        self.email = email
        self.retmax = retmax
        self.retstart = retstart

    def run(self):
        try:
            Entrez.email = self.email

            handle = Entrez.esearch(db="sra", term=self.query, retmax=self.retmax, retstart=self.retstart)
            record = Entrez.read(handle)
            handle.close()

            ids = ",".join(record["IdList"])
            summary_handle = Entrez.esummary(db="sra", id=ids, rettype="full")
            summary_records = Entrez.read(summary_handle)
            summary_handle.close()

            results = []

            for summary_record in summary_records:
                exp_xml = summary_record.get('ExpXml', '')
                fields = {
                    'Title': re.search(r'<Title>(.*?)</Title>', exp_xml),
                    'Platform': re.search(r'instrument_model="(.*?)"', exp_xml),
                    'Bases': re.search(r'total_bases="(\d+)"', exp_xml),
                    'Spots': re.search(r'total_spots="(\d+)"', exp_xml),
                    'Organism': re.search(r'ScientificName="(.*?)"', exp_xml),
                    'Strategy': re.search(r'<LIBRARY_STRATEGY>(.*?)</LIBRARY_STRATEGY>', exp_xml),
                    'Bioproject': re.search(r'<Bioproject>(.*?)</Bioproject>', exp_xml),
                    'Biosample': re.search(r'<Biosample>(.*?)</Biosample>', exp_xml),
                    'Run': re.search(r'acc="(.*?)"', summary_record.get('Runs', '')),
                }

                values = [
                    fields['Title'].group(1) if fields['Title'] else "",
                    fields['Platform'].group(1) if fields['Platform'] else "",
                    f"{int(fields['Bases'].group(1)):,}" if fields['Bases'] else "",
                    f"{int(fields['Spots'].group(1))*2:,}" if fields['Spots'] and 'library_layout="PAIRED"' in exp_xml else (f"{int(fields['Spots'].group(1)):,}" if fields['Spots'] else ""),
                    fields['Organism'].group(1) if fields['Organism'] else "",
                    fields['Strategy'].group(1) if fields['Strategy'] else "",
                    fields['Bioproject'].group(1) if fields['Bioproject'] else "",
                    fields['Biosample'].group(1) if fields['Biosample'] else "",
                    fields['Run'].group(1) if fields['Run'] else ""
                ]

                results.append(values)

            self.result_ready.emit(results)

        except Exception as e:
            self.error_occurred.emit(str(e))

class SearchHandler:
    def __init__(self, ui):
        self.ui = ui
        self.retstart = 0
        self.retmax = 9999
        self.total = None
        self.query = ""
        self.email = ""

    def start_search(self, reset=True):
        if reset:
            self.retstart = 0
            self.query = self.ui.lineEditQuery.text()
            self.email = self.ui.lineEditEmail_2.text()
            self.ui.treeSearchResults.clear()
        
            setattr(self, "_accum_results", [])
            if hasattr(self, "parent_app"):
                self.parent_app.reset_organism_combo() 
            
        if not self.query:
            CustomMessageBox.warning(None, "Warning", "Please enter a search query.")
            return
        if not self.email:
            CustomMessageBox.warning(None, "Warning","Please enter your email address.")
            return

        self.ui.buttonSearch.setEnabled(False)
        self.ui.progressSearch.setRange(0, 0) 

        self.worker = SraSearchWorker(self.query, self.email, retmax=self.retmax, retstart=self.retstart)
        self.worker.result_ready.connect(self.handle_results)
        self.worker.error_occurred.connect(self.handle_error)
        self.worker.start()

    def handle_results(self, results):
        self.ui.progressSearch.setRange(0, 1)
        self.ui.buttonSearch.setEnabled(True)
    
        for row in results:
            item_data = [
                row[8],  
                row[0],  
                row[1],  
                row[2],  
                row[3], 
                row[4],  
                row[5], 
                row[6],  
                row[7],  
            ]
            item = QTreeWidgetItem()
            for i, text in enumerate(item_data):
                item.setText(i, text)
                if i != 1: 
                    item.setTextAlignment(i, Qt.AlignCenter)
            self.ui.treeSearchResults.addTopLevelItem(item)
    
        count = self.ui.treeSearchResults.topLevelItemCount()
        self.ui.labelSearchCount.setText(f"Search List: {count:,} items")
        self.retstart += self.retmax
        
        results_for_org = []
        for row in results:
            run_acc = row[8] if len(row) > 8 else ""
            meta = {
                "organism": row[4] if len(row) > 4 else "",
                "title": row[0] if len(row) > 0 else "",
                "platform": row[1] if len(row) > 1 else "",
                "bases": row[2] if len(row) > 2 else "",
                "reads": row[3] if len(row) > 3 else "",
                "strategy": row[5] if len(row) > 5 else "",
                "bioproject": row[6] if len(row) > 6 else "",
                "biosample": row[7] if len(row) > 7 else "",
            }
            results_for_org.append((run_acc, meta))
        
        acc = getattr(self, "_accum_results", [])
        acc.extend(results_for_org)
        self._accum_results = acc
        
        if hasattr(self, "parent_app"):
            self.parent_app.update_organism_combo(self._accum_results)

    def handle_error(self, err):
        message = f"Search Failed:\n{err}"
        CustomMessageBox.error(None, "Failed", message)
        self.ui.progressSearch.setRange(0, 1)
        self.ui.buttonSearch.setEnabled(True)

    def load_more(self):
        if self.retstart == 0:
            CustomMessageBox.warning(None, "Warning", "No search started yet.")
            return
    
        if hasattr(self, "parent_app"):
            self.parent_app.clear_organism_filter_only()
    
        self.start_search(reset=False)


class EnaSearchWorker(QThread):
    result_ready = Signal(list)
    error_occurred = Signal(str)

    def __init__(self, query):
        super().__init__()
        self.query = query

    def run(self):
        try:
            raw = (self.query or "").strip()
            is_advanced = bool(re.search(r'\b(?:tax_eq|tax_tree|AND|OR|NOT|=)\b', raw, re.I))

            if is_advanced:
                query_expr = raw
            else:
                query_expr = f'study_title="{raw}"'

            url = "https://www.ebi.ac.uk/ena/portal/api/search"
            params = {
                "result": "read_run",
                "query": query_expr,
                "fields": "run_accession,study_accession,sample_accession,scientific_name,instrument_model,library_strategy",
                "format": "tsv",
                "limit": "0"
            }
            headers = {"User-Agent": "NGS-Downloader/1.0"}

            response = requests.get(url, params=params, headers=headers, timeout=30)
            response.raise_for_status()

            text = response.text.strip()
            if not text:
                self.result_ready.emit([])
                return

            lines = text.split("\n")
            if len(lines) <= 1:
                self.result_ready.emit([])
                return

            data_lines = lines[1:]
            results = [line.split("\t") for line in data_lines]
            self.result_ready.emit(results)

        except requests.HTTPError as e:
            msg = f"HTTP {e.response.status_code}: {e.response.text[:200]}"
            self.error_occurred.emit(msg)
        except Exception as e:
            self.error_occurred.emit(str(e))

class EnaSearchHandler:
    def __init__(self, ui):
        self.ui = ui

    def start_search(self):
        
        setattr(self, "_accum_results", [])
        if hasattr(self, "parent_app"):
            self.parent_app.clear_organism_filter_only()
            
        query = self.ui.lineEditQuery.text()
        tree_widget = self.ui.treeSearchResults
        progress_bar = self.ui.progressSearch
        button = self.ui.buttonSearch

        if not query:
            CustomMessageBox.warning(None, "Warning", "Please enter a search query.")
            return

        button.setEnabled(False)
        progress_bar.setRange(0, 0)

        self.worker = EnaSearchWorker(query)
        self.worker.result_ready.connect(lambda results: (
            self.update_tree_with_results(tree_widget, results),
            button.setEnabled(True),
            progress_bar.setRange(0, 1),
            self.ui.labelSearchCount.setText(f"Search List: {len(results):,} items")
        ))
        self.worker.error_occurred.connect(lambda err: (
            CustomMessageBox.error(None, "Search Failed", str(err)),
            button.setEnabled(True),
            progress_bar.setRange(0, 1)
        ))
        self.worker.start()

    def update_tree_with_results(self, tree_widget, results):
        tree_widget.clear()
        for row in results:
            item = QTreeWidgetItem()
            for i in range(len(row)):
                item.setText(i, row[i])
                if i != 3: 
                    item.setTextAlignment(i, Qt.AlignCenter)
            tree_widget.addTopLevelItem(item)
            
        results_for_org = []
        for row in results:
            run_acc = row[0] if len(row) > 0 else ""
            meta = {
                "organism": row[3] if len(row) > 3 else "",
                "instrument_model": row[4] if len(row) > 4 else "",
                "library_strategy": row[5] if len(row) > 5 else "",
                "study_accession": row[1] if len(row) > 1 else "",
                "sample_accession": row[2] if len(row) > 2 else "",
            }
            results_for_org.append((run_acc, meta))
    
        acc = getattr(self, "_accum_results", [])
        acc.extend(results_for_org)
        self._accum_results = acc
        
        if hasattr(self, "parent_app"):
            self.parent_app.update_organism_combo(self._accum_results)

class EmailConfigHandler:
    def __init__(self, ui):
        self.ui = ui
        self.config_path = os.path.join(os.getcwd(), "config.json")
        self.ui.buttonSaveEmail_2.clicked.connect(self.set_email_from_ui)
        self.load_email_to_ui()

    def is_valid_email(self, email):
        pattern = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
        return re.match(pattern, email)

    def save_email_to_config(self, email):
        with open(self.config_path, "w") as config_file:
            json.dump({"email": email}, config_file)

    def load_email_from_config(self):
        if os.path.exists(self.config_path):
            with open(self.config_path, "r") as config_file:
                data = json.load(config_file)
                return data.get("email", "")
        return ""

    def set_email_from_ui(self):
        email = self.ui.lineEditEmail_2.text().strip()
        if self.is_valid_email(email):
            Entrez.email = email
            self.save_email_to_config(email)
            CustomMessageBox.info(None, "Info", "Email saved successfully!")
        else:
            CustomMessageBox.error(None, "Error", "Please enter a valid email.")

    def load_email_to_ui(self):
        saved_email = self.load_email_from_config()
        if saved_email:
            self.ui.lineEditEmail_2.setText(saved_email)
            Entrez.email = saved_email

def apply_sra_column_styles(tree_widget):
    widths = [140, 300, 180, 100, 100, 150, 150, 150, 150]
    for i, width in enumerate(widths):
        tree_widget.setColumnWidth(i, width)
        tree_widget.headerItem().setTextAlignment(i, Qt.AlignCenter)
    tree_widget.updateGeometry()
    tree_widget.viewport().update()

def apply_ena_column_styles(tree_widget):
    widths = [180, 180, 180, 300, 180, 180]
    for i, width in enumerate(widths):
        tree_widget.setColumnWidth(i, width)
        tree_widget.headerItem().setTextAlignment(i, Qt.AlignCenter)                
    tree_widget.header().setStretchLastSection(True)
    tree_widget.updateGeometry()
    tree_widget.viewport().update()
            
class MainApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("NGS Data Downloader")
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
                
        self.move(QApplication.primaryScreen().availableGeometry().center() - self.rect().center())
        
        self.ui.treeSearchResults.setColumnWidth(0, 140) 
        self.ui.treeSearchResults.setColumnWidth(1, 300) 
        self.ui.treeSearchResults.setColumnWidth(2, 180)  
        self.ui.treeSearchResults.setColumnWidth(3, 100)  
        self.ui.treeSearchResults.setColumnWidth(4, 100)  
        self.ui.treeSearchResults.setColumnWidth(5, 150)  
        self.ui.treeSearchResults.setColumnWidth(6, 150)  
        self.ui.treeSearchResults.setColumnWidth(7, 150)  
        self.ui.treeSearchResults.setColumnWidth(8, 150)  
        
        self.ui.treeDownloadList.setColumnWidth(0, 140)
        self.ui.treeDownloadList.setColumnWidth(1, 300)
        self.ui.treeDownloadList.setColumnWidth(2, 180)
        self.ui.treeDownloadList.setColumnWidth(3, 100)
        self.ui.treeDownloadList.setColumnWidth(4, 100)
        self.ui.treeDownloadList.setColumnWidth(5, 150)
        self.ui.treeDownloadList.setColumnWidth(6, 150)
        self.ui.treeDownloadList.setColumnWidth(7, 150)
        self.ui.treeDownloadList.setColumnWidth(8, 150)

        self.email_handler = EmailConfigHandler(self.ui)
        self.current_search_slot = None
        self.current_search_handler = None
        self.ena_search_handler = None

        self.ui.treeSearchResults.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.ui.treeSearchResults.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.ui.treeDownloadList.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.ui.treeDownloadList.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._organism_filter_active = False

        self.connect_signals()

    def connect_signals(self):
        self.ui.buttonSelectDirectory.clicked.connect(self.select_directory)
        self.ui.buttonLoadMore.clicked.connect(lambda: self.current_search_handler.load_more() if self.current_search_handler else None)
        self.ui.SearchClearbutton.clicked.connect(lambda: SearchUIUtils.clear_search_results(self.ui))
        self.ui.filteringComboBox.currentIndexChanged.connect(self._apply_organism_combo_filter)
        self.ui.SearchClearbutton.clicked.connect(self._reset_top_org)
        self.ui.resetpushButton.clicked.connect(self.clear_organism_filter_only)
        self.ui.lineEditQuery.returnPressed.connect(lambda: self.ui.buttonSearch.click())
        self.ui.buttonDownloadStart.clicked.connect(self.start_download)
        self.ui.buttonMoveToDownload.clicked.connect(self.move_to_download)
        self.ui.buttonRemoveToDownloadList.clicked.connect(self.remove_from_download)
        self.ui.buttonSelectAllinSearchList.clicked.connect(lambda: select_all_items(self.ui.treeSearchResults))
        self.ui.buttonSelectAllinSearchList_2.clicked.connect(lambda: select_all_items(self.ui.treeDownloadList))
        self.ui.comboDatabase.currentTextChanged.connect(
            lambda _: on_database_change(
                self.ui, self.ui.comboDatabase,
                self.ui.buttonLoadMore, self.ui.treeSearchResults, self.ui.treeDownloadList, 
                self
            )
        )
        on_database_change(
            self.ui,
            self.ui.comboDatabase,
            self.ui.buttonLoadMore,
            self.ui.treeSearchResults,
            self.ui.treeDownloadList,
            self
        )
    def _reset_top_org(self):
        self.reset_organism_combo()

    def _apply_organism_combo_filter(self, _idx):
        data = self.ui.filteringComboBox.currentData()
        if not data:
            return
        self._organism_filter_active = True
        org = data.get("organism") if isinstance(data, dict) else str(data)
        handler = self.current_search_handler
        results = getattr(handler, "_accum_results", []) if handler else []
        if not results:
            return
        filtered = [(acc, meta) for acc, meta in results if (meta.get("organism") or "").strip() == org]
        self.populate_search_results(filtered)

    def populate_search_results(self, results):
        tree = self.ui.treeSearchResults
        tree.clear()

        db = self.ui.comboDatabase.currentText().strip().upper()

        if db == "SRA":
            headers = [
                "Run Accession", "Title", "Platform Model", "Total Bases",
                "Total Reads", "Organism", "Library Strategy", "Bioproject", "Biosample"
            ]
            tree.setColumnCount(len(headers))
            tree.setHeaderLabels(headers)
            apply_sra_column_styles(tree)

            for run_acc, meta in results:
                cols = [
                    run_acc,
                    meta.get("title", ""),
                    meta.get("platform", ""),
                    meta.get("bases", ""),
                    meta.get("reads", ""),
                    meta.get("organism", ""),
                    meta.get("strategy", ""),
                    meta.get("bioproject", ""),
                    meta.get("biosample", ""),
                ]
                item = QTreeWidgetItem()
                for i, txt in enumerate(cols):
                    item.setText(i, str(txt))
                    if i != 1:
                        item.setTextAlignment(i, Qt.AlignCenter)
                tree.addTopLevelItem(item)

        else: 
            headers = [
                "Run Accession", "Study Accession", "Sample Accession",
                "Organism", "Instrument Model", "Library Strategy"
            ]
            tree.setColumnCount(len(headers))
            tree.setHeaderLabels(headers)
            apply_ena_column_styles(tree)

            for run_acc, meta in results:
                cols = [
                    run_acc,
                    meta.get("study_accession", ""),
                    meta.get("sample_accession", ""),
                    meta.get("organism", ""),
                    meta.get("instrument_model", ""),
                    meta.get("library_strategy", ""),
                ]
                item = QTreeWidgetItem()
                for i, txt in enumerate(cols):
                    item.setText(i, str(txt))
                    if i != 3:  
                        item.setTextAlignment(i, Qt.AlignCenter)
                tree.addTopLevelItem(item)

        self.ui.labelSearchCount.setText(f"Search List: {tree.topLevelItemCount():,} items")
        tree.updateGeometry()
        tree.viewport().update()

    def update_organism_combo(self, results):
        c = Counter()
        for _, meta in (results or []):
            org = (meta.get("organism") or "").strip()
            if org:
                c[org] += 1
    
        total = sum(c.values()) or 1
        top = [(org, cnt, f"{(cnt/total)*100:.1f}%") for org, cnt in c.most_common(200)]
    
        cb = self.ui.filteringComboBox
        prev = cb.currentText() if cb.count() else ""
        cb.blockSignals(True)
        try:
            cb.clear()
            for org, cnt, pct in top:
                cb.addItem(f"{org} — {cnt:,} ({pct})", {"organism": org, "count": cnt, "percent": pct})
            if prev and any(prev == cb.itemText(i) for i in range(cb.count())):
                cb.setCurrentText(prev)
            else:
                cb.setCurrentIndex(-1)
        finally:
            cb.blockSignals(False)

    def clear_organism_filter_only(self):
        self._organism_filter_active = False
        cb = self.ui.filteringComboBox
        cb.blockSignals(True)
        try:
            cb.setCurrentIndex(-1)
        finally:
            cb.blockSignals(False)
    
        handler = self.current_search_handler
        results = getattr(handler, "_accum_results", []) if handler else []
        if results:
            self.populate_search_results(results)
    
    def reset_organism_combo(self):
        cb = self.ui.filteringComboBox
        cb.blockSignals(True)
        try:
            cb.clear()
        finally:
            cb.blockSignals(False)

    def select_directory(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_path:
            self.ui.lineEditDirectory.setText(os.path.normpath(dir_path))

    def prepare_download_data(self):
        selected_items = [self.ui.treeDownloadList.topLevelItem(i) for i in range(self.ui.treeDownloadList.topLevelItemCount())]
        if not selected_items:
            return None, "Please add items to download list."
        output_dir = self.ui.lineEditDirectory.text().strip()
        if not output_dir:
            return None, "Please select an output directory."
        selected_format = self.ui.comboFormat.currentText()
        download_data = []
        for item in selected_items:
            metadata = {
                "Run Accession": item.text(0),
                "Title": item.text(1),
                "Platform": item.text(2),
                "Bases": item.text(3),
                "Reads": item.text(4),
                "Organism": item.text(5),
                "Strategy": item.text(6),
                "Bioproject": item.text(7),
                "Biosample": item.text(8)
            }
            run_acc = metadata["Run Accession"]
            download_data.append((run_acc, metadata))
        return {
            "download_items": download_data,
            "output_dir": output_dir,
            "file_format": selected_format
        }, None

    def start_download(self):
        data, error = self.prepare_download_data()
        if error:
            CustomMessageBox.error(None, "error", error)
            return
        
        self.ui.buttonDownloadStart.setEnabled(False)
        status_window = DownloadStatusWindow(
            download_items=data["download_items"],
            controller=self,
            output_dir=data["output_dir"],
            file_format=data["file_format"],
            start_button=self.ui.buttonDownloadStart,
            database=self.ui.comboDatabase.currentText()
        )
        status_window.show()

    def move_to_download(self):
        move_to_download_list(
            self.ui.treeSearchResults,
            self.ui.treeDownloadList,
            self.ui.labelDownloadCount
        )

    def remove_from_download(self):
        remove_from_download_list(
            self.ui.treeDownloadList,
            self.ui.labelDownloadCount
        )

class SearchUIUtils:
    def clear_search_results(ui):
        ui.treeSearchResults.clear()
        ui.lineEditQuery.clear()
        ui.progressSearch.setRange(0, 1)
        ui.buttonSearch.setEnabled(True)
        ui.labelSearchCount.setText("Search List: 0 items")

class UIStateUtils:
    @staticmethod
    def disable_search_ui(ui):
        ui.buttonSearch.setEnabled(False)
        ui.progressSearch.setRange(0, 0)
        
    @staticmethod
    def enable_search_ui(ui):
        ui.buttonSearch.setEnabled(True)
        ui.progressSearch.setRange(0, 1)

    @staticmethod
    def disable_download_ui(ui):
        ui.buttonDownload.setEnabled(False)
        ui.progressDownload.setRange(0, 0)

    @staticmethod
    def enable_download_ui(ui):
        ui.buttonDownload.setEnabled(True)
        ui.progressDownload.setRange(0, 1)

    @staticmethod
    def select_directory(ui):
        directory = QFileDialog.getExistingDirectory()
        if directory:
            ui.labelSelectedDir.setText(directory)
        return directory

def run_command_with_retries(command, retries=3):
    for attempt in range(1, retries + 1):
        try:
            subprocess.run(command, shell=True, check=True)
            return True
        except subprocess.CalledProcessError:
            if attempt == retries:
                return False
            time.sleep(5)
    return False

def check_disk_space(directory, required_space_gb):
    try:
        total, used, free = shutil.disk_usage(directory)
        free_gb = free / (1024 ** 3)
        return free_gb >= required_space_gb
    except Exception:
        return False
    
class DiskMonitor(QThread):
    warning = Signal(str)

    def __init__(self, directory, required_gb, interval=60):
        super().__init__()
        self.directory = directory
        self.required_gb = required_gb
        self.interval = interval
        self._stop = False

    def run(self):
        while not self._stop:
            if not check_disk_space(self.directory, self.required_gb):
                self.warning.emit("🚨 Disk space is too low. Downloads will be stopped.")
                break

            for _ in range(self.interval):
                if self._stop:
                    break
                time.sleep(1)

    def stop(self):
        self._stop = True
        
def calculate_total_required_space(run_acc_list, prefetch_output_dir):
    total_required_space_gb = 0
    for run_acc in run_acc_list:
        sra_file_path = os.path.join(prefetch_output_dir, run_acc, f"{run_acc}.sra")
        if not os.path.exists(sra_file_path):
            continue
        sra_file_size_gb = os.path.getsize(sra_file_path) / (1024 ** 3)
        estimated_space = sra_file_size_gb * 5  
        total_required_space_gb += estimated_space
    return total_required_space_gb

class PrefetchWorker(QThread):
    finished = Signal()
    failed = Signal(str)

    def __init__(self, run_acc_list, prefetch_output_dir, max_workers=4):
        super().__init__()
        self.run_acc_list = run_acc_list
        self.output_dir = prefetch_output_dir
        self.max_workers = max_workers

    def run(self):
        os.makedirs(self.output_dir, exist_ok=True)
    
        def prefetch_one(run_acc):
        
            sra_file_path = os.path.join(self.output_dir, run_acc, f"{run_acc}.sra")
        
            if os.path.exists(sra_file_path):
                return True
        
            safe_output_dir = os.path.normpath(self.output_dir)
            cmd = f'prefetch "{run_acc}" --output-directory "{safe_output_dir}" --max-size 500G'
        
            success = run_command_with_retries(cmd)
        
            if not success or not os.path.exists(sra_file_path):
                return False
        
            return True

        failed_items = []

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {executor.submit(prefetch_one, run): run for run in self.run_acc_list}
            for future in futures:
                success = future.result()
                if not success:
                    failed_items.append(futures[future])

        if failed_items:
            self.failed.emit(", ".join(map(str, failed_items)))
        else:
            self.finished.emit()
            
def check_total_disk_space_and_prefetch(run_acc_list, prefetch_output_dir):
    os.makedirs(prefetch_output_dir, exist_ok=True)

    total_required_space_gb = calculate_total_required_space(run_acc_list, prefetch_output_dir)

    for run_acc in run_acc_list:
        sra_file_path = os.path.join(prefetch_output_dir, run_acc, f"{run_acc}.sra")
        if os.path.exists(sra_file_path):
            sra_file_size_gb = os.path.getsize(sra_file_path) / (1024 ** 3)
            estimated_space = sra_file_size_gb * 5
            total_required_space_gb += estimated_space

    if not check_disk_space(prefetch_output_dir, total_required_space_gb):
        raise RuntimeError(f"Not enough disk space. Required: {total_required_space_gb:.2f} GB")

    return total_required_space_gb

class ConversionRunner(QObject):
    finished = Signal(str)
    log = Signal(str)

    def __init__(self):
        super().__init__()
        self._stopped = False

    def stop(self):
        self._stopped = True

    def run(self, primary_cmd, fallback_cmd=None, retries=2):
        for attempt in range(1, retries + 1):
            if self._stopped:
                self.log.emit("⏹️ Conversion aborted by user.")
                self.finished.emit(False)
                return False  

            self.log.emit(f"▶️ Attempt {attempt}/{retries} for: {primary_cmd}")
            try:
                proc = subprocess.Popen(primary_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                if proc.returncode == 0:
                    self.log.emit("✅ Conversion succeeded.")
                    self.finished.emit(True)
                    return True  
                else:
                    self.log.emit(f"❌ Failed: {stderr.decode().strip()}")
            except Exception as e:
                self.log.emit(f"❌ Exception: {str(e)}")

            time.sleep(2)

        if fallback_cmd and not self._stopped:
            self.log.emit(f"🔁 Trying fallback: {fallback_cmd}")
            try:
                proc = subprocess.Popen(fallback_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()

                if proc.returncode == 0:
                    self.log.emit("✅ Fallback succeeded.")
                    self.finished.emit(True)
                    return True 
                else:
                    self.log.emit(f"❌ Fallback failed: {stderr.decode().strip()}")
            except Exception as e:
                self.log.emit(f"❌ Fallback exception: {str(e)}")

        self.finished.emit(False)
        return False 
    
def log_metadata(run_acc, output_dir, metadata, elapsed_time,
                 start_time=None, end_time=None, total_size=None):
    global log_file_name

    if log_file_name is None:
        log_file_name = f"NGS_data_downloader_{datetime.now().strftime('%Y%m%d_%H%M')}.csv"

    log_file_path = os.path.join(output_dir, log_file_name)

    with log_lock:
        file_exists = os.path.exists(log_file_path)

        with open(log_file_path, mode="a", newline="", encoding="utf-8") as log_file:
            writer = csv.writer(log_file)

            if not file_exists:
                writer.writerow([
                    "Run Accession", "Title", "Platform", "Bases", "Reads",
                    "Organism", "Strategy", "Bioproject", "Biosample",
                    "Start Time", "End Time", "Elapsed Time", "Size (MB)"
                ])

            if total_size is None:
                final_output_dir = os.path.abspath(os.path.join(output_dir, run_acc))
                if os.path.exists(final_output_dir):
                    total_size = sum(
                        os.path.getsize(os.path.join(final_output_dir, f))
                        for f in os.listdir(final_output_dir)
                        if os.path.isfile(os.path.join(final_output_dir, f))
                    ) / (1024 * 1024)
                else:
                    total_size = "N/A"

            formatted_size = f"{total_size:.2f} MB" if isinstance(total_size, (int, float)) else "N/A"

            writer.writerow([
                run_acc,
                metadata.get("Title", "N/A"),
                metadata.get("Platform", "N/A"),
                metadata.get("Bases", "N/A"),
                metadata.get("Reads", "N/A"),
                metadata.get("Organism", "N/A"),
                metadata.get("Strategy", "N/A"),
                metadata.get("Bioproject", "N/A"),
                metadata.get("Biosample", "N/A"),
                start_time or "N/A",
                end_time or "N/A",
                elapsed_time,
                formatted_size
            ])
           
def select_all_items(tree_widget: QTreeWidget):
    total_count = tree_widget.topLevelItemCount()
    selected_count = len(tree_widget.selectedItems())

    should_select = selected_count < total_count

    for i in range(total_count):
        item = tree_widget.topLevelItem(i)
        item.setSelected(should_select)

def update_download_count(label, tree_widget):
    count = tree_widget.topLevelItemCount()
    label.setText(f"Download List: {count:,} items")

def move_to_download_list(search_tree, download_tree, label_download_count):
    selected_items = search_tree.selectedItems()
    if not selected_items:
        return

    column_count = search_tree.columnCount()
    download_tree.setColumnCount(column_count)
    headers = [search_tree.headerItem().text(i) for i in range(column_count)]
    download_tree.setHeaderLabels(headers)

    existing_accessions = set()
    for i in range(download_tree.topLevelItemCount()):
        existing_accessions.add(download_tree.topLevelItem(i).text(0))

    for item in selected_items:
        run_accession = item.text(0)
        if run_accession in existing_accessions:
            continue
    
        new_item = QTreeWidgetItem()
        for i in range(item.columnCount()):
            text = item.text(i)
            new_item.setText(i, text)
            if column_count == 6 and i != 3:
                new_item.setTextAlignment(i, Qt.AlignCenter)
            elif column_count == 9 and i != 1:
                new_item.setTextAlignment(i, Qt.AlignCenter)
    
        download_tree.addTopLevelItem(new_item)
        existing_accessions.add(run_accession)

    label_download_count.setText(f"Download List: {download_tree.topLevelItemCount():,} items")
    
    if column_count == 6:
        apply_ena_column_styles(download_tree)
    elif column_count == 9:
        apply_sra_column_styles(download_tree)

def remove_from_download_list(download_tree, download_list_label):
    selected_items = download_tree.selectedItems()
    for item in selected_items:
        index = download_tree.indexOfTopLevelItem(item)
        if index != -1:
            download_tree.takeTopLevelItem(index)
    update_download_count(download_list_label, download_tree)

def on_database_change(ui, combo, load_more_button, tree_widget, download_widget, main_app=None):
    if main_app:
        main_app.reset_organism_combo()

    selected = combo.currentText().strip().upper()

    _toggle_email_area(ui, visible=(selected != "ENA"))

    tree_widget.clear()
    ui.treeDownloadList.clear()
    ui.labelDownloadCount.setText(f"Download List: {0:,} items")
    ui.labelSearchCount.setText(f"Search List: {0:,} items")
    load_more_button.setVisible(selected == "SRA")

    if main_app and main_app.current_search_slot:
        ui.buttonSearch.clicked.disconnect(main_app.current_search_slot)
        main_app.current_search_slot = None

    if selected == "ENA":
        ena_headers = [
            "Run Accession",
            "Study Accession",
            "Sample Accession",
            "Organism",
            "Instrument Model",
            "Library Strategy"
        ]
        
        tree_widget.setColumnCount(len(ena_headers))
        tree_widget.setHeaderLabels(ena_headers)
        download_widget.setColumnCount(len(ena_headers))
        download_widget.setHeaderLabels(ena_headers)
    
        apply_ena_column_styles(tree_widget)
        apply_ena_column_styles(download_widget)
        
        handler = EnaSearchHandler(ui)
        if main_app:
            main_app.current_search_handler = handler
            main_app.ena_search_handler = handler
            handler.parent_app = main_app
        main_app.current_search_slot = handler.start_search
        ui.buttonSearch.clicked.connect(main_app.current_search_slot)
    
    elif selected == "SRA":
        sra_headers = [
            "Run Accession",
            "Title",
            "Platform Model",
            "Total Bases",
            "Total Reads",
            "Organism",
            "Library Strategy",
            "Bioproject",
            "Biosample"
        ]
        
        tree_widget.setColumnCount(len(sra_headers))
        tree_widget.setHeaderLabels(sra_headers)
        download_widget.setColumnCount(len(sra_headers))
        download_widget.setHeaderLabels(sra_headers)
    
        apply_sra_column_styles(tree_widget)
        apply_sra_column_styles(download_widget)
    
        handler = SearchHandler(ui)
        if main_app:
            main_app.current_search_handler = handler
            
            handler.parent_app = main_app
        
        def sra_start():
            handler.start_search(reset=True)
    
        main_app.current_search_slot = sra_start
        ui.buttonSearch.clicked.connect(main_app.current_search_slot)
        
def _toggle_email_area(ui, visible: bool):
    widgets = [
        getattr(ui, "labelEmail_2", None),     
        ui.lineEditEmail_2,                     
        getattr(ui, "buttonSaveEmail_2", None)  
    ]
    for w in widgets:
        if w:
            w.setVisible(visible)

def toggle_load_more(button, selected_db):
    if selected_db == "SRA":
        button.setEnabled(True)
    else:
        button.setEnabled(False)

def download_selected(download_tree, output_dir, controller, format_selector, email_box, start_button=None):
    selected_items = [download_tree.topLevelItem(i) for i in range(download_tree.topLevelItemCount())]
    if not selected_items:
        CustomMessageBox("No items selected for download.").exec()
        return

    email = email_box.text().strip()
    if not email:
        CustomMessageBox("Please enter your email before downloading.").exec()
        return

    selected_format = format_selector.currentText()

    download_data = []
    for item in selected_items:
        metadata = {
            "Run Accession": item.text(0),
            "Title": item.text(1),
            "Platform": item.text(2),
            "Bases": item.text(3),
            "Reads": item.text(4),
            "Organism": item.text(5),
            "Strategy": item.text(6),
            "Bioproject": item.text(7),
            "Biosample": item.text(8)
        }
        run_acc = metadata["Run Accession"]
        download_data.append((run_acc, metadata))

    progress_window = DownloadStatusWindow(
        download_items=download_data,
        controller=controller,
        output_dir=output_dir,
        file_format=selected_format,
        start_button=start_button
    )
    progress_window.show()