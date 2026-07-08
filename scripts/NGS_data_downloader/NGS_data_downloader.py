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
import math
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

def _bytes_to_gb(n):
    return n / (1024**3)

def _fmt_gb(n):
    return f"{_bytes_to_gb(n):.1f} GB"

def _parse_mb_from_metadata(meta):

    if not isinstance(meta, dict):
        return None
    candidates = [
        "File capacity (MB)", "size_MB", "estimated_size_MB",
        "Size(MB)", "size", "size_mb", "Size_MB"
    ]
    for k in candidates:
        if k in meta and meta[k] is not None:
            try:
                return float(str(meta[k]).replace(",", "").strip())
            except:
                pass
    return None

def _ensure_free_space_or_fail(progress_put, error_emit, path_label, path_for_usage, required_bytes, run_acc):
    try:
        usage = shutil.disk_usage(path_for_usage)
        free_b = usage.free
    except Exception as e:
        progress_put(f"❌ {path_label} Free space check failed: {e}", 0)
        error_emit(f"Disk check failed on {path_label} for {run_acc}: {e}")
        return False

    if free_b < required_bytes:
        need = _fmt_gb(required_bytes)
        free = _fmt_gb(free_b)
        msg = f"❌ Insufficient disk space\nRequired space: {need} / {path_label} Free on: {free}"
        progress_put(msg, 0)

        error_emit(f"Not enough free space on {path_label} for {run_acc}")
        return False
    return True


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
        self.tree.setHeaderLabels(["Run Accession", "Activity", "Status"])
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
            
            bar.setRange(0, 1)
            bar.setValue(0)
            
            bar.setTextVisible(False)
            
            bar.setFixedHeight(16)
            bar.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
            bar.setStyleSheet("""
                QProgressBar {
                    background-color: #E4E8EB;
                    border: none;
                    border-radius: 4px;
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
                status = data.get("status", "")
    
                if acc in self.status_labels:
                    self.status_labels[acc].setText(2, status)
    
                bar = self.status_bars.get(acc)
                if bar is None:
                    continue
    
                normalized_status = status.lower()
    
                if "completed" in normalized_status:
                    bar.setRange(0, 1)
                    bar.setValue(1)
    
                elif (
                    "canceling" in normalized_status
                    or "canceled" in normalized_status
                    or "cancelled" in normalized_status
                ):
                    bar.setRange(0, 1)
                    bar.setValue(0)
    
                elif (
                    "error" in normalized_status
                    or "failed" in normalized_status
                    or "not found" in normalized_status
                    or "not available" in normalized_status
                    or "insufficient" in normalized_status
                ):
                    bar.setRange(0, 1)
                    bar.setValue(0)
    
                elif "waiting" in normalized_status:
                    bar.setRange(0, 1)
                    bar.setValue(0)

                else:
                    bar.setRange(0, 0)
    
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
    
            worker.error.connect(
                lambda message, acc=acc: self.handle_worker_error(acc, message)
            )
            
            worker.finished.connect(
                partial(
                    self.on_worker_thread_finished,
                    worker,
                    acc
                )
            )
    
            self.pending_queue.append(worker)
    
        self.start_next_batch()

    def start_next_batch(self):
        if self._canceled:
            return
    
        while (
            len(self.active_workers) < self.max_concurrent
            and self.pending_queue
        ):
            worker = self.pending_queue.pop(0)
            self.active_workers.append(worker)
            worker.start()

    def on_worker_thread_finished(self, worker, acc):
        if worker in self.active_workers:
            self.active_workers.remove(worker)
    
        self.mark_worker_finished(acc)
    
        if not self._canceled:
            self.start_next_batch()

    def handle_worker_error(self, acc, message, *args):
        if acc not in self.failed_files:
            self.failed_files.append(acc)
    
        if self.progress_queue:
            self.progress_queue.put({
                "acc": acc,
                "status": message,
                "progress": 0
            })
        
    def mark_worker_finished(self, acc, *_):
        if acc in self._updated_files:
            return
    
        self._updated_files.add(acc)
    
        self.completed_count = len(self._updated_files)
    
        self.completed_label.setText(
            f"Completed: {self.completed_count} / {self.total_count}  "
        )
    
        if (
            self.completed_count == self.total_count
            and not self._completion_shown
        ):
            self._completion_shown = True
    
            if self.start_button:
                self.start_button.setEnabled(True)
    
            if self._canceled:
                msg = "🛑 Download process was canceled."
    
            elif self.failed_files:
                msg = "❌ Some files failed to download."
    
            else:
                msg = "✅ All files downloaded successfully."
    
            msg_box = CustomMessageBox(
                title="Download Summary",
                message=msg,
                parent=self
            )
            msg_box.setWindowFlags(
                msg_box.windowFlags()
                | Qt.WindowStaysOnTopHint
            )
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
        
            acc = worker.run_acc
        
            if acc in self.status_labels:
                self.status_labels[acc].setText(
                    2,
                    "Canceling..."
                )
        
            if acc in self.status_bars:
                bar = self.status_bars[acc]
                bar.setRange(0, 1)
                bar.setValue(0)
        
            if self.progress_queue:
                self.progress_queue.put({
                    "acc": acc,
                    "status": "Canceling...",
                    "progress": 0
                })
                
        for worker in self.pending_queue:
            acc = worker.run_acc
            
            if acc in self.status_bars:
                bar = self.status_bars[acc]
                bar.setRange(0, 1)
                bar.setValue(0)            
            
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

class SraDownloadCancelled(Exception):
    pass
                           
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
        
        self._stop_event = threading.Event()
        
        self._process_lock = threading.Lock()
        self._active_process = None
        
        self._working_dir = None

    def _cleanup_working_dir(self):
        working_dir = self._working_dir
        self._working_dir = None
    
        if not working_dir or not os.path.exists(working_dir):
            return
    
        def handle_remove_error(func, path, exc_info):
            try:
                os.chmod(path, stat.S_IWRITE)
                func(path)
            except Exception:
                pass
    
        for _ in range(3):
            try:
                shutil.rmtree(
                    working_dir,
                    onerror=handle_remove_error
                )
                break
    
            except PermissionError:
                time.sleep(0.1)
    
            except FileNotFoundError:
                break
    
            except Exception:
                break
            
    def _commit_working_dir(self, final_output_dir):
        self._raise_if_cancelled()
    
        if not self._working_dir:
            raise RuntimeError(
                "SRA working directory is not available."
            )
    
        backup_dir = None
    
        if os.path.exists(final_output_dir):
            backup_dir = (
                f"{final_output_dir}.backup_"
                f"{int(time.time() * 1000)}"
            )
    
            os.replace(
                final_output_dir,
                backup_dir
            )
    
        try:
            os.replace(
                self._working_dir,
                final_output_dir
            )
    
            self._working_dir = None
    
        except Exception:
            if (
                backup_dir
                and os.path.exists(backup_dir)
                and not os.path.exists(final_output_dir)
            ):
                os.replace(
                    backup_dir,
                    final_output_dir
                )
    
            raise
    
        else:
            if backup_dir and os.path.exists(backup_dir):
                shutil.rmtree(
                    backup_dir,
                    ignore_errors=True
                )

    def stop(self):
        if self._stop_event.is_set():
            return
    
        self._stopped = True
        self._stop_event.set()
    
        with self._process_lock:
            process = self._active_process
    
        self._terminate_process_tree(process)
        
    def _cancel_requested(self):
        return self._stop_event.is_set()
    
    def _raise_if_cancelled(self):
        if self._cancel_requested():
            raise SraDownloadCancelled()

    def _set_active_process(self, process):
        with self._process_lock:
            self._active_process = process
    
    def _clear_active_process(self, process=None):
        with self._process_lock:
            if (
                process is None
                or self._active_process is process
            ):
                self._active_process = None

    def _terminate_process_tree(self, process):
        if process is None:
            return
    
        try:
            if process.poll() is not None:
                return
    
            if os.name == "nt":
                subprocess.run(
                    [
                        "taskkill",
                        "/PID",
                        str(process.pid),
                        "/T",
                        "/F"
                    ],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    check=False
                )
            else:
                process.terminate()
    
                try:
                    process.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    process.kill()
    
        except Exception:
            try:
                process.kill()
            except Exception:
                pass

    def _run_command_cancelable(self, command):
        self._raise_if_cancelled()
    
        process = None
    
        try:
            process = subprocess.Popen(
                command,
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
    
            self._set_active_process(process)
    
            while True:
                return_code = process.poll()
    
                if return_code is not None:
                    return return_code == 0
    
                if self._cancel_requested():
                    self._terminate_process_tree(process)
                    raise SraDownloadCancelled()
    
                time.sleep(0.2)
    
        finally:
            self._clear_active_process(process)
    
            if process is not None and process.poll() is None:
                self._terminate_process_tree(process)

    def _compress_fastq_cancelable(self, source_path, destination_path):
        try:
            with open(source_path, "rb") as source_file:
                with gzip.open(destination_path, "wb") as gzip_file:
    
                    while True:
                        self._raise_if_cancelled()
    
                        chunk = source_file.read(
                            1024 * 1024
                        )
    
                        if not chunk:
                            break
    
                        gzip_file.write(chunk)
    
            self._raise_if_cancelled()
    
        except SraDownloadCancelled:
            try:
                if os.path.exists(destination_path):
                    os.remove(destination_path)
            except Exception:
                pass
    
            raise
    
        except Exception:
            try:
                if os.path.exists(destination_path):
                    os.remove(destination_path)
            except Exception:
                pass
    
            raise

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
            self._raise_if_cancelled()
            est_mb = _parse_mb_from_metadata(self.metadata)
            est_sra_bytes = int((est_mb if est_mb else 2048) * 1024**2)

            FASTQ_FACTOR = 4.0
            PREFETCH_OVERHEAD = 1 * 1024**3

            prefetch_root = os.path.join("C:/NDD_dummy_file", self.run_acc)
            output_root = self.output_dir

            required_prefetch = est_sra_bytes + PREFETCH_OVERHEAD
            required_output   = int(est_sra_bytes * FASTQ_FACTOR)

            prefetch_drive_path = os.path.splitdrive(prefetch_root)[0] + os.sep
            output_drive_path   = os.path.splitdrive(output_root)[0] + os.sep

            if not _ensure_free_space_or_fail(self._put_progress, self.error.emit,
                                              f"{prefetch_drive_path} (Temporary storage)",
                                              prefetch_drive_path,
                                              required_prefetch, self.run_acc):
                return

            if not _ensure_free_space_or_fail(self._put_progress, self.error.emit,
                                              f"{output_drive_path} (Output folder)",
                                              output_drive_path,
                                              required_output, self.run_acc):
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
    
            self._raise_if_cancelled()
            
            safe_temp_sra_root = os.path.normpath(
                temp_sra_root
            )
            
            prefetch_cmd = (
                f'prefetch "{self.run_acc}" '
                f'--output-directory "{safe_temp_sra_root}" '
                f'--max-size 500G'
            )
            
            prefetch_success = self._run_command_cancelable(
                prefetch_cmd
            )
            
            self._raise_if_cancelled()
            
            if not prefetch_success:
                self._put_progress(
                    "❌ Prefetch failed",
                    0
                )
                self.error.emit(
                    f"Prefetch failed for {self.run_acc}"
                )
                return
            
            self._raise_if_cancelled()
    
            candidate_paths = [
                os.path.join(temp_sra_root, self.run_acc, f"{self.run_acc}.sra"),
                os.path.join(temp_sra_root, f"{self.run_acc}.sra")
            ]
            sra_path = next((p for p in candidate_paths if os.path.exists(p)), None)

            self._raise_if_cancelled()
            
            if not sra_path:
                self._put_progress("❌ .sra file not found", 0)
                self.error.emit(f".sra file not found for {self.run_acc}")
                return
            
            try:
                real_sra_bytes = os.path.getsize(sra_path)
            except Exception:
                real_sra_bytes = est_sra_bytes

            final_required_output = int(real_sra_bytes * FASTQ_FACTOR)

            if not _ensure_free_space_or_fail(self._put_progress, self.error.emit,
                                              f"{output_drive_path} (Output folder)",
                                              output_drive_path,
                                              final_required_output, self.run_acc):
                try:
                    if os.path.isfile(sra_path):
                        os.remove(sra_path)
                except Exception:
                    pass
                return
            
            self._put_progress("Downloading...", 50)
    
            final_output_dir = os.path.abspath(
                os.path.join(
                    self.output_dir,
                    self.run_acc
                )
            )
            
            os.makedirs(
                os.path.abspath(self.output_dir),
                exist_ok=True
            )
            
            self._working_dir = tempfile.mkdtemp(
                prefix=f".{self.run_acc}.partial_",
                dir=os.path.abspath(self.output_dir)
            )
            
            out_dir = self._working_dir
    
            if self.selected_format.upper() == "FASTA":
                cmd = f'fastq-dump "{sra_path}" --outdir "{out_dir}" --split-files --fasta 0'
                fallback_cmd = None
            else:
                cmd = f'fasterq-dump "{sra_path}" -O "{out_dir}" --split-files --temp "{temp_path}"'
                fallback_cmd = f'fastq-dump "{sra_path}" -O "{out_dir}" --split-files --gzip'
                
            retry_limit = 3
            self.success = False
            
            for attempt in range(1, retry_limit + 1):
                self._raise_if_cancelled()
            
                self._put_progress(
                    f"Converting {self.run_acc} "
                    f"(Attempt {attempt}/{retry_limit})...",
                    0
                )
            
                primary_success = self._run_command_cancelable(
                    cmd
                )
            
                self._raise_if_cancelled()
            
                if primary_success:
                    self.success = True
                    break
            
                if fallback_cmd:
                    self._put_progress(
                        f"Primary conversion failed. "
                        f"Trying fallback for {self.run_acc}...",
                        0
                    )
            
                    fallback_success = self._run_command_cancelable(
                        fallback_cmd
                    )
            
                    self._raise_if_cancelled()
            
                    if fallback_success:
                        self.success = True
                        break
            
                if attempt < retry_limit:
                    self._put_progress(
                        f"Conversion failed. Retrying "
                        f"{self.run_acc} "
                        f"({attempt + 1}/{retry_limit})...",
                        0
                    )
            
                    if self._stop_event.wait(3):
                        raise SraDownloadCancelled()

            self._raise_if_cancelled()
            
            if not self.success:
                self._put_progress("❌ Conversion failed", 0)
                self.error.emit(f"Conversion failed for {self.run_acc}")
                return
    
            if (
                self.selected_format.upper() == "FASTQ"
                and "fasterq-dump" in cmd
            ):
                fastq_files = [
                    file_name
                    for file_name in os.listdir(out_dir)
                    if file_name.endswith(".fastq")
                ]
            
                for file_name in fastq_files:
                    self._raise_if_cancelled()
            
                    source_path = os.path.join(
                        out_dir,
                        file_name
                    )
            
                    destination_path = (
                        source_path + ".gz"
                    )
            
                    self._put_progress(
                        f"Compressing {file_name}...",
                        0
                    )
            
                    try:
                        self._compress_fastq_cancelable(
                            source_path,
                            destination_path
                        )
            
                    except SraDownloadCancelled:
                        raise
            
                    except Exception as e:
                        self._put_progress(
                            f"❌ Compression failed: {file_name}",
                            0
                        )
            
                        self.error.emit(
                            f"Compression failed for "
                            f"{file_name}: {e}"
                        )
            
                        self.success = False
                        return
            
                    self._raise_if_cancelled()
            
                    try:
                        os.remove(source_path)
            
                    except Exception as e:
                        self._put_progress(
                            f"❌ Failed to remove uncompressed file: "
                            f"{file_name}",
                            0
                        )
            
                        self.error.emit(
                            f"Failed to remove uncompressed "
                            f"FASTQ file {file_name}: {e}"
                        )
            
                        self.success = False
                        return

            self._raise_if_cancelled()
            
            self._put_progress(
                "Finalizing download...",
                0
            )
            
            self._commit_working_dir(
                final_output_dir
            )
            
            out_dir = final_output_dir

            end_dt = datetime.now()
            end_str = end_dt.strftime("%Y-%m-%d %H:%M:%S")
            elapsed_str = str(end_dt - start_dt).split(".")[0]

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
                print(
                    f"[WARN] log_metadata failed "
                    f"for {self.run_acc}: {e}"
                )
            
            self.success = True
            self._put_progress(
                "✅ Download Completed",
                100
            )
            self.progress.emit(1)
    
        except SraDownloadCancelled:
            self.success = False
        
            self._put_progress(
                "🛑 Canceled",
                0
            )
        
            self._cleanup_working_dir()
        
        except Exception as e:
            self.success = False
            self._put_progress(
                f"❌ Error: {e}",
                0
            )
        
            self._cleanup_working_dir()
            self.error.emit(str(e))
            
        finally:
            self._cleanup_working_dir()
        
            self.cleanup_temp_files(
                temp_path,
                temp_sra_root,
                run_tmp
            )
        
            self.done.emit(self.run_acc)

class EnaDownloadCancelled(Exception):
    pass
        
class EnaDownloadWorker(QThread):
    progress = Signal(int)
    status = Signal(str)
    error = Signal(str)
    done = Signal(str)
    
    def __init__(
        self,
        run_acc,
        metadata,
        output_dir,
        progress_queue=None
    ):
        super().__init__()
    
        self.run_acc = run_acc
        self.metadata = metadata
        self.output_dir = output_dir
        self.progress_queue = progress_queue
    
        self._stopped = False
        self._stop_event = threading.Event()
        self._network_lock = threading.Lock()
        self._active_response = None
        self._session = None
    
        self.success = False
        self._emitted_done = False
        self._working_dir = None

    def stop(self):
        if self._stop_event.is_set():
            return
    
        self._stopped = True
        self._stop_event.set()
    
        with self._network_lock:
            response = self._active_response
            session = self._session
    
        if response is not None:
            try:
                response.close()
            except Exception:
                pass
    
        if session is not None:
            try:
                session.close()
            except Exception:
                pass
            
    def _cancel_requested(self):
        return self._stop_event.is_set()
    
    
    def _raise_if_cancelled(self):
        if self._cancel_requested():
            raise EnaDownloadCancelled()
    
    
    def _set_session(self, session):
        with self._network_lock:
            self._session = session
    
    
    def _set_active_response(self, response):
        with self._network_lock:
            self._active_response = response
    
    
    def _clear_active_response(self, response=None):
        with self._network_lock:
            if response is None or self._active_response is response:
                self._active_response = None
        
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

    def _cleanup_working_dir(self):
        working_dir = self._working_dir
        self._working_dir = None
    
        if not working_dir or not os.path.exists(working_dir):
            return
    
        def handle_remove_error(func, path, exc_info):
            try:
                os.chmod(path, stat.S_IWRITE)
                func(path)
            except Exception:
                pass
    
        for _ in range(3):
            try:
                shutil.rmtree(
                    working_dir,
                    onerror=handle_remove_error
                )
                break
    
            except PermissionError:
                time.sleep(0.1)
    
            except FileNotFoundError:
                break
    
            except Exception:
                break
            
    def _commit_working_dir(self, final_output_dir):
    
        self._raise_if_cancelled()
    
        if not self._working_dir:
            raise RuntimeError("ENA working directory is not available.")
    
        backup_dir = None
    
        if os.path.exists(final_output_dir):
            backup_dir = (
                f"{final_output_dir}.backup_"
                f"{int(time.time() * 1000)}"
            )
            os.replace(final_output_dir, backup_dir)
    
        try:
            os.replace(
                self._working_dir,
                final_output_dir
            )
            self._working_dir = None
    
        except Exception:
            if (
                backup_dir
                and os.path.exists(backup_dir)
                and not os.path.exists(final_output_dir)
            ):
                os.replace(
                    backup_dir,
                    final_output_dir
                )
    
            raise
    
        else:
            if backup_dir and os.path.exists(backup_dir):
                shutil.rmtree(
                    backup_dir,
                    ignore_errors=True
                )

    def run(self):
        start_dt = datetime.now()
        start_str = start_dt.strftime("%Y-%m-%d %H:%M:%S")
    
        final_output_dir = None
    
        session = requests.Session()
        session.headers.update({
            "User-Agent": "NGS-Data-Downloader/1.0"
        })
    
        self._set_session(session)
    
        try:
            self._raise_if_cancelled()
            self._put_progress(
                "Fetching ENA file info...",
                0
            )
    
            api_url = (
                "https://www.ebi.ac.uk/ena/portal/api/filereport"
                f"?accession={self.run_acc}"
                "&result=read_run"
                "&fields=fastq_ftp"
                "&format=tsv"
            )
    
            response = None
    
            try:
                self._raise_if_cancelled()
    
                response = session.get(
                    api_url,
                    timeout=(10, 10)
                )
    
                self._set_active_response(response)
    
                response.raise_for_status()
                self._raise_if_cancelled()
    
                response_text = response.text
    
            finally:
                self._clear_active_response(response)
    
                if response is not None:
                    try:
                        response.close()
                    except Exception:
                        pass
    
            self._raise_if_cancelled()
            lines = response_text.splitlines()
    
            if len(lines) < 2:
                raise FileNotFoundError(
                    f"No FASTQ file information was returned for "
                    f"{self.run_acc}."
                )
    
            headers = lines[0].split("\t")
            values = lines[1].split("\t")
    
            if len(values) < len(headers):
                values.extend(
                    [""] * (len(headers) - len(values))
                )
    
            row = dict(zip(headers, values))
    
            fastq_ftp_urls = row.get(
                "fastq_ftp",
                ""
            ).strip()
    
            if not fastq_ftp_urls:
                raise FileNotFoundError(
                    f"No downloadable FASTQ files are currently "
                    f"available for {self.run_acc} in ENA."
                )
    
            fastq_files = [
                url.strip()
                for url in fastq_ftp_urls.split(";")
                if url.strip()
            ]
    
            if not fastq_files:
                raise FileNotFoundError(
                    f"No downloadable FASTQ files are currently "
                    f"available for {self.run_acc} in ENA."
                )
    
            self._raise_if_cancelled()

            final_output_dir = os.path.abspath(
                os.path.join(
                    self.output_dir,
                    self.run_acc
                )
            )
    
            os.makedirs(
                os.path.abspath(self.output_dir),
                exist_ok=True
            )
    
            self._working_dir = tempfile.mkdtemp(
                prefix=f".{self.run_acc}.partial_",
                dir=os.path.abspath(self.output_dir)
            )

            for file_url in fastq_files:
                self._raise_if_cancelled()
    
                if file_url.startswith("ftp://"):
                    https_url = (
                        "https://"
                        + file_url[len("ftp://"):]
                    )
    
                elif file_url.startswith("https://"):
                    https_url = file_url
    
                elif file_url.startswith(
                    "ftp.sra.ebi.ac.uk/"
                ):
                    https_url = f"https://{file_url}"
    
                else:
                    raise ValueError(
                        f"Invalid ENA FASTQ URL returned for "
                        f"{self.run_acc}: {file_url}"
                    )
    
                file_name = os.path.basename(
                    https_url.split("?", 1)[0]
                )
    
                if not file_name:
                    raise ValueError(
                        f"Could not determine the FASTQ file name "
                        f"for {self.run_acc}: {https_url}"
                    )
    
                output_file_path = os.path.join(
                    self._working_dir,
                    file_name
                )
    
                success = False
                max_attempts = 3
                last_error = None

                for attempt in range(
                    1,
                    max_attempts + 1
                ):
                    response = None
    
                    try:
                        self._raise_if_cancelled()
    
                        self._put_progress(
                            f"Downloading {file_name} "
                            f"(Attempt {attempt}/{max_attempts})",
                            0
                        )
    
                        response = session.get(
                            https_url,
                            stream=True,
                            timeout=(15, 60)
                        )
    

                        self._set_active_response(
                            response
                        )
    
                        response.raise_for_status()
                        self._raise_if_cancelled()
    
                        with open(
                            output_file_path,
                            "wb"
                        ) as output_file:
    
                            for chunk in response.iter_content(
                                chunk_size=1024 * 1024
                            ):
                                self._raise_if_cancelled()
    
                                if chunk:
                                    output_file.write(
                                        chunk
                                    )
    
                        self._raise_if_cancelled()
    
                        success = True
                        break
    
                    except EnaDownloadCancelled:
                        raise
    
                    except Exception as e:
                        last_error = e
    

                        if self._cancel_requested():
                            raise EnaDownloadCancelled()
    
                        if attempt < max_attempts:
                            self._put_progress(
                                f"Download failed. Retrying "
                                f"{file_name} "
                                f"({attempt + 1}/{max_attempts})...",
                                0
                            )
    
                            if self._stop_event.wait(5):
                                raise EnaDownloadCancelled()
    
                    finally:
                        self._clear_active_response(
                            response
                        )
    
                        if response is not None:
                            try:
                                response.close()
                            except Exception:
                                pass
    
                if not success:
                    self.success = False
    
                    self._put_progress(
                        f"❌ Failed to download {file_name}",
                        0
                    )
    
                    error_detail = (
                        f": {last_error}"
                        if last_error is not None
                        else ""
                    )
    
                    self.error.emit(
                        f"Failed to download {file_name} "
                        f"after {max_attempts} attempts"
                        f"{error_detail}"
                    )
    
                    self._cleanup_working_dir()
                    return

            self._raise_if_cancelled()
    
            self._put_progress(
                "Finalizing download...",
                0
            )
    
            self._commit_working_dir(
                final_output_dir
            )
    
            end_dt = datetime.now()
            elapsed_str = str(
                end_dt - start_dt
            ).split(".")[0]
    
            total_size_mb = 0
    
            if os.path.isdir(final_output_dir):
                total_size_mb = (
                    sum(
                        os.path.getsize(
                            os.path.join(
                                final_output_dir,
                                file_name
                            )
                        )
                        for file_name in os.listdir(
                            final_output_dir
                        )
                        if os.path.isfile(
                            os.path.join(
                                final_output_dir,
                                file_name
                            )
                        )
                    )
                    / (1024 * 1024)
                )
    
            try:
                log_metadata(
                    run_acc=self.run_acc,
                    output_dir=self.output_dir,
                    metadata=self.metadata or {},
                    elapsed_time=elapsed_str,
                    start_time=start_str,
                    end_time=end_dt.strftime(
                        "%Y-%m-%d %H:%M:%S"
                    ),
                    total_size=total_size_mb
                )
    
            except Exception as e:
                print(
                    f"[WARN] log_metadata failed "
                    f"for {self.run_acc}: {e}"
                )
    
            self.success = True
    
            self._put_progress(
                "✅ Download Completed",
                100
            )
    
            self.progress.emit(1)

        except EnaDownloadCancelled:
            self.success = False
    
            self._put_progress(
                "🛑 Canceled",
                0
            )
    
            self._cleanup_working_dir()

        except Exception as e:
            self.success = False
    
            if self._cancel_requested():
                self._put_progress(
                    "🛑 Canceled",
                    0
                )
    
            else:
                self._put_progress(
                    f"❌ Error: {e}",
                    0
                )
    
                self.error.emit(
                    f"Exception for "
                    f"{self.run_acc}: {e}"
                )
    
            self._cleanup_working_dir()

        finally:
            self._clear_active_response()
    
            with self._network_lock:
                if self._session is session:
                    self._session = None
    
            try:
                session.close()
            except Exception:
                pass

            self._cleanup_working_dir()
    
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
        selected_items = [
            self.ui.treeDownloadList.topLevelItem(i)
            for i in range(
                self.ui.treeDownloadList.topLevelItemCount()
            )
        ]
    
        if not selected_items:
            return None, "Please add items to download list."
    
        output_dir = self.ui.lineEditDirectory.text().strip()
    
        if not output_dir:
            return None, "Please select an output directory."
    
        selected_format = self.ui.comboFormat.currentText()
        selected_database = (
            self.ui.comboDatabase.currentText()
            .strip()
            .upper()
        )
    
        download_data = []
    
        for item in selected_items:
            run_acc = item.text(0)
    
            if selected_database == "ENA":
                metadata = {
                    "Database": "ENA",
                    "Run Accession": run_acc,
                    "Study Accession": item.text(1),
                    "Sample Accession": item.text(2),
                    "Organism": item.text(3),
                    "Instrument Model": item.text(4),
                    "Library Strategy": item.text(5)
                }
    
            else:
                metadata = {
                    "Database": "SRA",
                    "Run Accession": run_acc,
                    "Title": item.text(1),
                    "Platform": item.text(2),
                    "Bases": item.text(3),
                    "Reads": item.text(4),
                    "Organism": item.text(5),
                    "Strategy": item.text(6),
                    "Bioproject": item.text(7),
                    "Biosample": item.text(8)
                }
    
            download_data.append(
                (run_acc, metadata)
            )
    
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
    env = os.environ.copy() 
    for attempt in range(1, retries + 1):
        try:
            subprocess.run(command, shell=True, check=True, env=env)  
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
    
def log_metadata(
    run_acc,
    output_dir,
    metadata,
    elapsed_time,
    start_time=None,
    end_time=None,
    total_size=None
):
    global log_file_name

    database = (
        metadata.get("Database", "SRA")
        .strip()
        .upper()
    )

    if log_file_name is None:
        log_file_name = (
            f"NGS_data_downloader_"
            f"{datetime.now().strftime('%Y%m%d_%H%M')}.csv"
        )

    base_name, extension = os.path.splitext(
        log_file_name
    )

    database_log_name = (
        f"{base_name}_{database}{extension}"
    )

    log_file_path = os.path.join(
        output_dir,
        database_log_name
    )

    with log_lock:
        file_exists = os.path.exists(
            log_file_path
        )

        if total_size is None:
            final_output_dir = os.path.abspath(
                os.path.join(
                    output_dir,
                    run_acc
                )
            )

            if os.path.isdir(final_output_dir):
                total_size = (
                    sum(
                        os.path.getsize(
                            os.path.join(
                                final_output_dir,
                                file_name
                            )
                        )
                        for file_name in os.listdir(
                            final_output_dir
                        )
                        if os.path.isfile(
                            os.path.join(
                                final_output_dir,
                                file_name
                            )
                        )
                    )
                    / (1024 * 1024)
                )
            else:
                total_size = None

        formatted_size = (
            f"{total_size:.2f} MB"
            if isinstance(
                total_size,
                (int, float)
            )
            else "N/A"
        )

        with open(
            log_file_path,
            mode="a",
            newline="",
            encoding="utf-8-sig"
        ) as log_file:
            writer = csv.writer(log_file)

            if database == "ENA":
                if not file_exists:
                    writer.writerow([
                        "Run Accession",
                        "Study Accession",
                        "Sample Accession",
                        "Organism",
                        "Instrument Model",
                        "Library Strategy",
                        "Start Time",
                        "End Time",
                        "Elapsed Time",
                        "Size (MB)"
                    ])

                writer.writerow([
                    run_acc,
                    metadata.get(
                        "Study Accession",
                        "N/A"
                    ),
                    metadata.get(
                        "Sample Accession",
                        "N/A"
                    ),
                    metadata.get(
                        "Organism",
                        "N/A"
                    ),
                    metadata.get(
                        "Instrument Model",
                        "N/A"
                    ),
                    metadata.get(
                        "Library Strategy",
                        "N/A"
                    ),
                    start_time or "N/A",
                    end_time or "N/A",
                    elapsed_time,
                    formatted_size
                ])

            else:
                if not file_exists:
                    writer.writerow([
                        "Run Accession",
                        "Title",
                        "Platform",
                        "Bases",
                        "Reads",
                        "Organism",
                        "Strategy",
                        "Bioproject",
                        "Biosample",
                        "Start Time",
                        "End Time",
                        "Elapsed Time",
                        "Size (MB)"
                    ])

                writer.writerow([
                    run_acc,
                    metadata.get(
                        "Title",
                        "N/A"
                    ),
                    metadata.get(
                        "Platform",
                        "N/A"
                    ),
                    metadata.get(
                        "Bases",
                        "N/A"
                    ),
                    metadata.get(
                        "Reads",
                        "N/A"
                    ),
                    metadata.get(
                        "Organism",
                        "N/A"
                    ),
                    metadata.get(
                        "Strategy",
                        "N/A"
                    ),
                    metadata.get(
                        "Bioproject",
                        "N/A"
                    ),
                    metadata.get(
                        "Biosample",
                        "N/A"
                    ),
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