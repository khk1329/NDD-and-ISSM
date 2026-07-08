import os
import gzip
import random
import threading
from Bio import SeqIO
from Bio.Seq import reverse_complement
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
from datetime import datetime
import itertools
import logging
import sys
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import seaborn as sns
import re
import psutil
import tempfile
from rapidfuzz import fuzz
import numpy as np
import aiofiles
import asyncio
import io
import subprocess
import shutil
from PySide6.QtCore import (
    QCoreApplication, QDate, QDateTime, QLocale, QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt, QTimer, QEventLoop, Slot, Q_ARG, QSettings, QThread, Signal
)
from PySide6.QtGui import (
    QBrush, QColor, QConicalGradient, QCursor, QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter, QPalette, QPixmap, QRadialGradient, QTransform
)
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QDoubleSpinBox, QFrame, QLabel, QLineEdit,
    QListWidget, QListWidgetItem, QMainWindow, QPushButton, QRadioButton, QSizePolicy,
    QTextEdit, QWidget, QFileDialog, QVBoxLayout, QHBoxLayout, QProgressBar, QComboBox,
    QSpinBox, QScrollArea, QMessageBox, QDialog, QGridLayout, QSpacerItem, QSizePolicy,
    QTreeWidget, QTreeWidgetItem, QHeaderView)
import queue
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED
from queue import Queue

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from ui.ISSM import Ui_In_silico_sequence_mining

logging.getLogger('matplotlib').setLevel(logging.WARNING)
analysis_cancelled = {"flag": False}

lock = threading.Lock()

_SEQKIT_PATH_CACHE = None
_SEQKIT_PATH_CHECKED = False
_SEQKIT_PATH_ERROR = None
_SEQKIT_PATH_LOCK = threading.Lock()

logging.basicConfig(
    level=logging.DEBUG,  
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]  
)

IUPAC_CODES = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
    "I": ["A", "C", "G", "T"]
}

def save_last_directory(key, path):
    settings = QSettings("MyCompany", "MyApp") 
    settings.setValue(key, path)

def load_last_directory(key):
    settings = QSettings("MyCompany", "MyApp")
    return settings.value(key, "")

def get_optimal_parallel_file_count(file_path=None):
    if file_path and os.path.exists(file_path):
        file_size = os.path.getsize(file_path)  
        file_size_gb = file_size / (1024**3)

        if file_size_gb > 10: 
            logging.warning(f"⚠️ Limiting the number of parallel processes due to the large file size ({file_size_gb:.2f}GB)")
            return 1

    cpu_count = os.cpu_count()
    mem = psutil.virtual_memory().total / (1024**3)

    if cpu_count >= 64 and mem >= 128:
        parallel_count = 4
    elif cpu_count >= 32 and mem >= 64:
        parallel_count = 3
    elif cpu_count >= 16 and mem >= 32:
        parallel_count = 2
    else:
        parallel_count = 1

    return min(parallel_count, 61)

def get_app_base_dir():
    if getattr(sys, "frozen", False):
        return getattr(sys, "_MEIPASS", os.path.dirname(sys.executable))

    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

def get_seqkit_path():
    global _SEQKIT_PATH_CACHE
    global _SEQKIT_PATH_CHECKED
    global _SEQKIT_PATH_ERROR

    with _SEQKIT_PATH_LOCK:
        if _SEQKIT_PATH_CHECKED:
            if _SEQKIT_PATH_CACHE:
                return _SEQKIT_PATH_CACHE

            raise FileNotFoundError(
                _SEQKIT_PATH_ERROR or
                "SeqKit executable was not found. Expected bundled tools/seqkit.exe or seqkit in PATH."
            )

        base_dir = get_app_base_dir()

        candidate_paths = [
            os.path.join(base_dir, "tools", "seqkit.exe"),
            os.path.join(os.path.dirname(base_dir), "tools", "seqkit.exe"),
            os.path.join(os.path.dirname(os.path.dirname(base_dir)), "tools", "seqkit.exe"),
        ]

        for seqkit_path in candidate_paths:
            if os.path.isfile(seqkit_path):
                _SEQKIT_PATH_CACHE = seqkit_path
                _SEQKIT_PATH_CHECKED = True
                logging.info(f"SeqKit executable found: {seqkit_path}")
                return _SEQKIT_PATH_CACHE

        path_seqkit = shutil.which("seqkit")
        if path_seqkit:
            _SEQKIT_PATH_CACHE = path_seqkit
            _SEQKIT_PATH_CHECKED = True
            logging.info(f"SeqKit executable found in PATH: {path_seqkit}")
            return _SEQKIT_PATH_CACHE

        _SEQKIT_PATH_ERROR = (
            "SeqKit executable was not found. "
            "Expected bundled tools/seqkit.exe or seqkit in PATH."
        )
        _SEQKIT_PATH_CHECKED = True

        raise FileNotFoundError(_SEQKIT_PATH_ERROR)

def run_seqkit_command(args, analysis_cancelled=None):
    creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    cmd = [get_seqkit_path()] + args

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        creationflags=creationflags
    )

    try:
        while proc.poll() is None:
            if analysis_cancelled and analysis_cancelled.get("flag"):
                proc.terminate()

                try:
                    proc.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()

                raise InterruptedError("SeqKit command cancelled by user.")

            time.sleep(0.1)

        stdout, stderr = proc.communicate()

        if proc.returncode != 0:
            raise subprocess.CalledProcessError(
                proc.returncode,
                cmd,
                output=stdout,
                stderr=stderr
            )

        return subprocess.CompletedProcess(
            cmd,
            proc.returncode,
            stdout,
            stderr
        )

    except Exception:
        if proc.poll() is None:
            try:
                proc.kill()
                proc.wait()
            except Exception:
                pass
        raise

def get_total_reads_with_seqkit(input_file, analysis_cancelled=None):
    result = run_seqkit_command(["stats", "-T", input_file], analysis_cancelled)

    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError("SeqKit stats returned no data.")

    header = lines[0].split("\t")
    values = lines[1].split("\t")

    if "num_seqs" not in header:
        raise ValueError("SeqKit stats output does not contain num_seqs column.")

    num_idx = header.index("num_seqs")
    return int(values[num_idx].replace(",", ""))

def run_seqkit_sample2(input_file, output_file, sample_size, analysis_cancelled=None):
    run_seqkit_command([
        "sample2",
        "-n", str(sample_size),
        "-2",
        "-r",
        "-o", output_file,
        input_file
    ], analysis_cancelled)
    
_ISSM_MATCHER_PATH_CACHE = None
_ISSM_MATCHER_PATH_CHECKED = False
_ISSM_MATCHER_PATH_ERROR = None
_ISSM_MATCHER_PATH_LOCK = threading.Lock()

def get_issm_matcher_path():
    global _ISSM_MATCHER_PATH_CACHE
    global _ISSM_MATCHER_PATH_CHECKED
    global _ISSM_MATCHER_PATH_ERROR

    with _ISSM_MATCHER_PATH_LOCK:
        if _ISSM_MATCHER_PATH_CHECKED:
            if _ISSM_MATCHER_PATH_CACHE:
                return _ISSM_MATCHER_PATH_CACHE
            raise FileNotFoundError(
                _ISSM_MATCHER_PATH_ERROR or
                "issm_matcher.exe was not found. Expected bundled tools/issm_matcher.exe or issm_matcher in PATH."
            )

        base_dir = get_app_base_dir()
        candidate_paths = [
            os.path.join(base_dir, "tools", "issm_matcher.exe"),
            os.path.join(os.path.dirname(base_dir), "tools", "issm_matcher.exe"),
            os.path.join(os.path.dirname(os.path.dirname(base_dir)), "tools", "issm_matcher.exe"),
        ]

        for matcher_path in candidate_paths:
            if os.path.isfile(matcher_path):
                _ISSM_MATCHER_PATH_CACHE = matcher_path
                _ISSM_MATCHER_PATH_CHECKED = True
                logging.info(f"ISSM matcher executable found: {matcher_path}")
                return _ISSM_MATCHER_PATH_CACHE

        path_matcher = shutil.which("issm_matcher")
        if path_matcher:
            _ISSM_MATCHER_PATH_CACHE = path_matcher
            _ISSM_MATCHER_PATH_CHECKED = True
            logging.info(f"ISSM matcher executable found in PATH: {path_matcher}")
            return _ISSM_MATCHER_PATH_CACHE

        _ISSM_MATCHER_PATH_ERROR = (
            "issm_matcher.exe was not found. "
            "Expected bundled tools/issm_matcher.exe or issm_matcher in PATH."
        )
        _ISSM_MATCHER_PATH_CHECKED = True
        raise FileNotFoundError(_ISSM_MATCHER_PATH_ERROR)

def normalize_matching_method(matching_method):
    method = str(matching_method or "fast_biological").strip().lower()
    if method in ("legacy", "legacy_fuzzy", "rapidfuzz", "fuzzy"):
        return "legacy_fuzzy"
    return "fast_biological"

def should_use_issm_matcher(threshold, sampling_percentage, matching_method="fast_biological"):
    try:
        thr = float(threshold)
    except Exception:
        return False

    if thr >= 100:
        return True

    return normalize_matching_method(matching_method) == "fast_biological"

def get_recommended_matcher_threads():
    cpu = os.cpu_count() or 1
    if cpu <= 2:
        return 1
    return max(1, min(cpu - 1, 8))

def write_temp_probe_fasta_for_matcher(probes, output_dir_path):
    os.makedirs(output_dir_path, exist_ok=True)
    fd, path = tempfile.mkstemp(prefix="issm_probes_", suffix=".fasta", dir=output_dir_path)
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        for idx, (name, seq) in enumerate(probes, start=1):
            safe_name = str(name).strip() or f"probe_{idx}"
            safe_seq = str(seq).replace(" ", "").replace("\t", "").upper()
            safe_seq = safe_seq.replace("I", "N")
            fh.write(f">{safe_name}\n{safe_seq}\n")
    return path

def parse_issm_matcher_probe_counts(result_file_path, probes):
    probe_counts = Counter({name: 0 for name, _ in probes})
    total_records = None
    total_matched = None
    in_probe_table = False

    with open(result_file_path, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if not line:
                if in_probe_table:
                    in_probe_table = False
                continue

            parts = line.split("\t")
            if len(parts) >= 2 and parts[0] in ("Total Selected Reads", "Total Selected Sequences"):
                try:
                    total_records = int(parts[1].replace(",", ""))
                except Exception:
                    pass
                continue
            if len(parts) >= 2 and parts[0] in ("Total Matched Reads", "Total Matched Sequences"):
                try:
                    total_matched = int(parts[1].replace(",", ""))
                except Exception:
                    pass
                continue

            if line.startswith("Probe Name\t"):
                in_probe_table = True
                continue
            if line.startswith("Matched Reads Information") or line.startswith("Matched Sequences Information"):
                in_probe_table = False
                continue

            if in_probe_table and len(parts) >= 3:
                probe_name = parts[0]
                try:
                    probe_counts[probe_name] = int(parts[2].replace(",", ""))
                except Exception:
                    probe_counts[probe_name] = 0

    return probe_counts, total_records, total_matched

def append_issm_matcher_summary(summary_file_path, original_filename, total_records, start_time, end_time):
    elapsed_time = str(end_time - start_time).split(".")[0]
    header_line = "Index\tFilename\tRecordCount\tStartTime\tEndTime\tElapsedTime\n"

    with lock:
        existing_lines = []
        if os.path.exists(summary_file_path):
            with open(summary_file_path, "r", encoding="utf-8") as fh:
                existing_lines = fh.read().strip().splitlines()

        if not existing_lines:
            existing_lines = [header_line.strip()]
        elif existing_lines[0] != header_line.strip():
            existing_lines = [header_line.strip()] + existing_lines

        index_counter = len(existing_lines) if len(existing_lines) > 1 else 1
        existing_lines.append(
            f"{index_counter}\t{os.path.basename(original_filename)}\t{total_records}\t{start_time}\t{end_time}\t{elapsed_time}"
        )

        with open(summary_file_path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(existing_lines) + "\n")

def set_busy_progress(update_queue, input_file, message):
    if update_queue is None:
        return

    update_queue.put(("progress_max", input_file, 0))
    update_queue.put(("progress_update", input_file, 0))
    update_queue.put(("status_update", input_file, message))

def emit_percent_progress(
    update_queue,
    input_file,
    phase,
    processed=0,
    total=None,
    percentage=None,
    unit="reads",
    reset=False
):

    if update_queue is None:
        return

    try:
        processed = max(0, int(processed))
    except (TypeError, ValueError):
        processed = 0

    try:
        valid_total = int(total) if total is not None else 0
    except (TypeError, ValueError):
        valid_total = 0

    if valid_total > 0:
        processed = min(processed, valid_total)
        pct = int((processed * 100) / valid_total)
        status = (
            f"{phase} | {pct}% | "
            f"{processed:,}/{valid_total:,} {unit}"
        )
    else:
        try:
            pct = int(percentage) if percentage is not None else 0
        except (TypeError, ValueError):
            pct = 0

        pct = max(0, min(100, pct))
        status = (
            f"{phase} | {pct}% | "
            f"{processed:,} {unit} processed"
        )

    pct = max(0, min(100, pct))

    if reset:
        update_queue.put(("progress_max", input_file, 100))

    update_queue.put(("progress_update", input_file, pct))
    update_queue.put(("status_update", input_file, status))

def set_terminal_progress(
    update_queue,
    input_file,
    message,
    progress_value=0
):
    if update_queue is None:
        return

    progress_value = max(0, min(100, int(progress_value)))
    update_queue.put(("progress_max", input_file, 100))
    update_queue.put(("progress_update", input_file, progress_value))
    update_queue.put(("status_update", input_file, message))


def run_issm_matcher_process(
    cmd,
    analysis_cancelled=None,
    update_queue=None,
    input_file=None,
    total_records=None,
    unit="reads"
):
    creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    stderr_lines = []
    line_queue = queue.Queue()
    last_emitted_percentage = -1
    last_emitted_processed = 0

    logging.info("ISSM matcher command started: " + " ".join(f'\"{c}\"' if " " in str(c) else str(c) for c in cmd))

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        creationflags=creationflags,
    )

    def _reader():
        try:
            for line in proc.stderr:
                line_queue.put(line.rstrip("\n"))
        except Exception as e:
            line_queue.put(f"__READER_ERROR__\t{e}")

    reader_thread = threading.Thread(target=_reader, daemon=True)
    reader_thread.start()

    try:
        while proc.poll() is None:
            if analysis_cancelled and analysis_cancelled.get("flag"):
                proc.terminate()
                try:
                    proc.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                raise InterruptedError("ISSM matcher cancelled by user.")

            while True:
                try:
                    line = line_queue.get_nowait()
                except queue.Empty:
                    break
                
                stderr_lines.append(line)
                logging.debug(f"ISSM matcher: {line}")

                if (
                    line.startswith("PROGRESS\t")
                    and update_queue is not None
                    and input_file is not None
                ):
                    try:
                        parts = line.split("\t")
                        processed = int(parts[1]) if len(parts) >= 2 else 0
                        bytes_read = int(parts[2]) if len(parts) >= 4 else 0
                        total_bytes = int(parts[3]) if len(parts) >= 4 else 0

                        if total_records and total_records > 0:
                            pct = max(
                                0,
                                min(
                                    100,
                                    int(
                                        min(processed, total_records)
                                        * 100
                                        / total_records
                                    ),
                                ),
                            )

                            if pct != last_emitted_percentage:
                                emit_percent_progress(
                                    update_queue=update_queue,
                                    input_file=input_file,
                                    phase="Matching",
                                    processed=processed,
                                    total=total_records,
                                    unit=unit,
                                )
                                last_emitted_percentage = pct
                                last_emitted_processed = processed

                        elif total_bytes > 0:
                            pct = max(
                                0,
                                min(
                                    100,
                                    int(
                                        min(bytes_read, total_bytes)
                                        * 100
                                        / total_bytes
                                    ),
                                ),
                            )

                            if pct != last_emitted_percentage:
                                emit_percent_progress(
                                    update_queue=update_queue,
                                    input_file=input_file,
                                    phase="Matching",
                                    processed=processed,
                                    total=None,
                                    percentage=pct,
                                    unit=unit,
                                )
                                last_emitted_percentage = pct
                                last_emitted_processed = processed

                        elif processed - last_emitted_processed >= 100000:
                            update_queue.put(
                                (
                                    "status_update",
                                    input_file,
                                    f"Matching | {processed:,} {unit} processed",
                                )
                            )
                            last_emitted_processed = processed

                    except Exception:
                        pass

            time.sleep(0.1)

        reader_thread.join(timeout=1)
        while True:
            try:
                line = line_queue.get_nowait()
            except queue.Empty:
                break
            stderr_lines.append(line)
            logging.debug(f"ISSM matcher: {line}")

        if proc.returncode != 0:
            tail = "\n".join(stderr_lines[-20:])
            raise subprocess.CalledProcessError(proc.returncode, cmd, stderr=tail)

        logging.info("ISSM matcher command completed.")
        return stderr_lines

    except Exception:
        if proc.poll() is None:
            try:
                proc.kill()
                proc.wait()
            except Exception:
                pass
        raise

def run_issm_matcher_exact_and_update(
    input_file,
    probes,
    total_records,
    output_file_path,
    output_dir_path,
    summary_file_path,
    results,
    selected_format,
    original_filename,
    start_time,
    analysis_cancelled,
    update_queue,
    matcher_threads=None,
    matcher_input_file=None,
    threshold=100,
    matching_method="fast_biological",
):
    probe_fasta_path = None
    try:
        actual_input_file = matcher_input_file or input_file
        if matcher_threads is None:
            matcher_threads = get_recommended_matcher_threads()
        try:
            matcher_threads = max(1, int(matcher_threads))
        except Exception:
            matcher_threads = get_recommended_matcher_threads()

        unit = "reads" if selected_format == "fastq" else "sequences"
        mode_label = (
            "Fast biological"
            if normalize_matching_method(matching_method) == "fast_biological"
            else "Matcher"
        )
        
        emit_percent_progress(
            update_queue=update_queue,
            input_file=input_file,
            phase="Matching",
            processed=0,
            total=total_records if total_records else None,
            percentage=0,
            unit=unit,
            reset=True
        )
        
        logging.info(
            f"{mode_label} matcher started with {matcher_threads} thread(s)."
        )

        probe_fasta_path = write_temp_probe_fasta_for_matcher(probes, output_dir_path)
        matched_fasta_path = os.path.join(
            output_dir_path,
            f"{os.path.splitext(os.path.basename(input_file))[0]}_matched_reads.fasta",
        )

        cmd = [
            get_issm_matcher_path(),
            "--input", actual_input_file,
            "--probes", probe_fasta_path,
            "--result", output_file_path,
            "--matched-fasta", matched_fasta_path,
            "--format", selected_format,
            "--threshold", str(float(threshold)),
            "--threads", str(matcher_threads),
            "--batch-size", "20000",
        ]

        run_issm_matcher_process(
            cmd,
            analysis_cancelled=analysis_cancelled,
            update_queue=update_queue,
            input_file=input_file,
            total_records=total_records,
            unit=unit
        )

        probe_counts, matcher_total_records, matcher_total_matched = (
            parse_issm_matcher_probe_counts(output_file_path, probes)
        )
        if matcher_total_records is not None:
            total_records = matcher_total_records

        end_time = datetime.now()
        append_issm_matcher_summary(
            summary_file_path,
            original_filename,
            total_records,
            start_time,
            end_time,
        )

        with lock:
            results["file_results"][original_filename] = {"sampled_records": total_records}
            for probe_name, _ in probes:
                count = probe_counts.get(probe_name, 0)
                results["file_results"][original_filename][probe_name] = count
                results["probe_results"].setdefault(probe_name, {})[original_filename] = count

        emit_percent_progress(
            update_queue=update_queue,
            input_file=input_file,
            phase="Matching",
            processed=total_records if total_records else 0,
            total=total_records if total_records else None,
            percentage=100,
            unit=unit
        )
        return probe_counts, total_records, matcher_total_matched

    finally:
        if probe_fasta_path and os.path.exists(probe_fasta_path):
            try:
                os.remove(probe_fasta_path)
            except Exception as e:
                logging.warning(f"Failed to remove temporary probe FASTA: {probe_fasta_path} — {e}")

def browse_folder(file_paths, folder_label_entry, all_files_list_box, selected_file_list_box, all_files_count_label, selected_files_count_label, selected_format):
    last_dir = load_last_directory("input_folder")
    folder_path = QFileDialog.getExistingDirectory(None, "Select Directory", last_dir)
    
    if not folder_path:
        return

    save_last_directory("input_folder", folder_path)

    list_files(folder_path, file_paths, all_files_list_box, selected_file_list_box, all_files_count_label, selected_files_count_label, selected_format)

    folder_label_entry.setReadOnly(False)  
    folder_label_entry.clear()
    folder_label_entry.setText(folder_path)
    folder_label_entry.setReadOnly(True)

def list_files(folder_path, file_paths, all_files_list_box, selected_file_list_box, all_files_count_label, selected_files_count_label, selected_format):

    fastq_extensions = (".fastq.gz", ".fastq", ".fq.gz", ".fq")
    fasta_extensions = (".fasta.gz", ".fasta", ".fa.gz", ".fa", ".fna.gz", ".fna", ".fas.gz", ".fas")

    if selected_format == "fastq":
        valid_extensions = fastq_extensions
    elif selected_format == "fasta":
        valid_extensions = fasta_extensions
    else:
        valid_extensions = fastq_extensions + fasta_extensions

    all_files_list_box.clear()
    file_paths.clear()

    for root, _, files in os.walk(folder_path):
        for file in files:
            file_lower = file.lower()
            if file_lower.endswith(valid_extensions):
                full_path = os.path.join(root, file)
                all_files_list_box.addItem(file)
                file_paths[file] = full_path

    update_file_count(all_files_list_box, selected_file_list_box, all_files_count_label, selected_files_count_label)

def update_file_count(all_files_list_box, selected_file_list_box, all_files_count_label, selected_files_count_label):
    all_files_count = all_files_list_box.count()
    selected_files_count = selected_file_list_box.count()
    all_files_count_label.setText(f"Total : {all_files_count}")
    selected_files_count_label.setText(f"Total : {selected_files_count}")

def toggle_all_selection(listbox, is_checked):
    if is_checked.isChecked():
        listbox.selectAll()
    else:
        listbox.clearSelection()

def update_select_all_checkbox(listbox, checkbox):
    total_items = listbox.count()
    selected_items = len(listbox.selectedItems())

    checkbox.blockSignals(True)

    if total_items == 0:
        checkbox.setChecked(False)
    elif total_items == selected_items:
        checkbox.setChecked(True)
    else:
        checkbox.setChecked(False)

    checkbox.blockSignals(False)

def move_selected_files(all_files_list_box, selected_file_list_box, file_paths, selected_file_paths, all_files_count_label, selected_files_count_label):
    selected_items = all_files_list_box.selectedItems()
    for item in selected_items:
        file_name = item.text()
        try:
            file_path = file_paths.pop(file_name)
            selected_file_paths[file_name] = file_path
            selected_file_list_box.addItem(file_name)
            all_files_list_box.takeItem(all_files_list_box.row(item))
        except KeyError:
            QMessageBox.warning(None, "Warning", f"The file '{file_name}' does not exist")

    update_file_count(all_files_list_box, selected_file_list_box, all_files_count_label, selected_files_count_label)
    
def move_back_files(all_files_list_box, selected_file_list_box, file_paths, selected_file_paths, all_files_count_label, selected_files_count_label):
    selected_items = selected_file_list_box.selectedItems()
    for item in selected_items:
        file_name = item.text()
        try:
            file_path = selected_file_paths.pop(file_name)
            file_paths[file_name] = file_path
            all_files_list_box.addItem(file_name)
            selected_file_list_box.takeItem(selected_file_list_box.row(item))
        except KeyError:
            QMessageBox.warning(None, "Warning", f"File '{file_name}' was not found in selected_file_paths")

    update_file_count(all_files_list_box, selected_file_list_box, all_files_count_label, selected_files_count_label)

def update_probe_count(probes_text_widget, probe_count_label):
    probes_text = probes_text_widget.toPlainText().strip()
    if not probes_text:
        probe_count_label.setText("Probe : 0")
        return

    lines = probes_text.splitlines()
    probes = [line for line in lines if line.startswith('>')]
    probe_count_label.setText(f"Probe : {len(probes)}")

def select_output_path(output_file_entry):
    settings = QSettings("MyCompany", "InSilicoTool")
    last_output = settings.value("lastOutputFolder", "")

    output_dir_path = QFileDialog.getExistingDirectory(None, "Select Output Directory", last_output)
    if not output_dir_path:
        return

    settings.setValue("lastOutputFolder", output_dir_path)
    output_file_entry.setText(output_dir_path)

def fuzzy_match_prepared(possible_sequences, read_sequence, threshold):
    if threshold >= 100:
        for seq in possible_sequences:
            if seq in read_sequence:
                return True, 100
        return False, 0

    for seq in possible_sequences:
        partial_ratio = fuzz.partial_ratio(seq, read_sequence)
        if partial_ratio >= threshold:
            return True, partial_ratio

    return False, 0

shared_probe_items = None

def init_worker(_probe_items):
    global shared_probe_items
    shared_probe_items = _probe_items

def init_worker_all(_probe_items, safe_cores=None):
    global shared_probe_items
    shared_probe_items = _probe_items

    if safe_cores is not None and os.name == 'nt':
        try:
            import psutil
            p = psutil.Process()
            p.cpu_affinity(safe_cores)
        except Exception as e:
            print(f"⚠️ Failed to set CPU affinity: {e}")

def match_batch_sequences(batch, threshold, batch_idx):
    matched = []
    probe_items = shared_probe_items or []

    for record in batch:
        if analysis_cancelled["flag"]:
            break

        read_id = record.id
        read_seq = str(record.seq)
        read_len = len(read_seq)

        for probe_name, probe_len, possible_sequences, possible_sequences_rc in probe_items:
            if read_len < probe_len:
                continue

            match_orig, score_orig = fuzzy_match_prepared(
                possible_sequences,
                read_seq,
                threshold
            )

            match_rc, score_rc = fuzzy_match_prepared(
                possible_sequences_rc,
                read_seq,
                threshold
            )

            if match_orig or match_rc:
                matched.append((read_id, probe_name, score_orig, score_rc, read_seq))

    batch_size = len(batch)
    return matched, batch_idx, batch_size

def match_batch_reads(records, threshold, batch_index):
    matched = []
    probe_items = shared_probe_items or []

    for read_id, read_seq in records:
        if analysis_cancelled["flag"]:
            break

        read_len = len(read_seq)

        for probe_name, probe_len, possible_sequences, possible_sequences_rc in probe_items:
            if read_len < probe_len:
                continue

            match_orig, score_orig = fuzzy_match_prepared(
                possible_sequences,
                read_seq,
                threshold
            )

            match_rc, score_rc = fuzzy_match_prepared(
                possible_sequences_rc,
                read_seq,
                threshold
            )

            if match_orig or match_rc:
                matched.append((read_id, probe_name, score_orig, score_rc, read_seq))

    batch_size = len(records)
    return matched, batch_index, batch_size

def run_analysis(
    probe_text: str,
    threshold: int,
    sampling_percentage: float,
    output_dir_path: str,
    selected_files: list,
    selected_format: str,
    analysis_cancelled,
    update_queue,
    matching_method: str = "fast_biological"
):
    logging.debug("Analysis started.")
    logging.debug(f"Selected format: {selected_format}")
    logging.debug(f"Output directory: {output_dir_path}")
    logging.debug(f"Files ({len(selected_files)}): {selected_files[:5]}{'...' if len(selected_files)>5 else ''}")
    matching_method = normalize_matching_method(matching_method)

    probes = []
    for chunk in [c for c in probe_text.strip().split('>') if c.strip()]:
        lines = chunk.splitlines()
        name = lines[0].strip() if lines else ""
        seq = ''.join(lines[1:]).strip().replace(' ', '').upper()
        if not name or not seq:
            raise ValueError("Invalid probe entry (name/sequence is empty).")
        probes.append((name, seq))
    logging.debug(f"Probes parsed: {len(probes)}")
    if probes:
        logging.debug(f"Example probe — {probes[0][0]}: {probes[0][1][:30]}{'...' if len(probes[0][1])>30 else ''}")

    if should_use_issm_matcher(threshold, sampling_percentage, matching_method):
        probe_items = []
        logging.debug("Fast matcher mode selected; skipped Python IUPAC variant expansion.")
    else:
        probe_items = prepare_probe_match_items(probes)
        logging.debug(f"Probe match items prepared: {len(probe_items)} items.")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_path = os.path.join(output_dir_path, f"summary_result_{ts}.txt")
    logging.debug(f"Summary will be written to: {summary_path}")

    results = {"file_results": {}, "probe_results": {}}
    completed_files = {'count': 0}
    total_files = len(selected_files)

    t = threading.Thread(
        target=process_files,
        args=(
            selected_files, probes, probe_items, threshold, sampling_percentage,
            output_dir_path, selected_format,
            completed_files, total_files, summary_path, results,
            analysis_cancelled, update_queue,
            matching_method
        ),
        daemon=True
    )
    t.start()
    t.join()

    update_queue.put(("finished", True, f"Completed: processed {total_files} files."))
    
def process_files(
    selected_files, probes, probe_items, threshold, sampling_percentage,
    output_dir_path, selected_format,
    completed_files, total_files, summary_file_path, results,
    analysis_cancelled, update_queue,
    matching_method="fast_biological"
):

    logging.debug("🚀 ThreadPoolExecutor Start")
    matching_method = normalize_matching_method(matching_method)
    use_fast_matcher = should_use_issm_matcher(threshold, sampling_percentage, matching_method)
    
    matcher_threads = None

    if use_fast_matcher:
        max_parallel = 1
        matcher_threads = get_recommended_matcher_threads()
        logging.info(
            f"ISSM fast matcher mode: file-level parallel count fixed to 1; "
            f"auto matcher threads = {matcher_threads}; "
            f"sampling={sampling_percentage}%; method={matching_method}; threshold={threshold}"
        )
    else:
        max_parallel = get_optimal_parallel_file_count()
        logging.info(f"Legacy fuzzy mode auto parallel count: {max_parallel}")

    def run_file(input_file, output_file_specific):
        file_key = os.path.basename(input_file)

        try:
            if analysis_cancelled["flag"]:
                update_queue.put(("status_update", file_key, "🚫 Cancelled"))
                return
        
            detect_func = (
                FASTQ_detect_and_quantify_probe
                if selected_format == "fastq"
                else FASTA_detect_and_quantify_probe_batched
            )

            detect_func(
                input_file,
                probes,
                probe_items,
                threshold,
                sampling_percentage,
                output_file_specific,
                completed_files,
                total_files,
                summary_file_path,
                results,
                output_dir_path,
                analysis_cancelled,
                update_queue,
                matcher_threads=matcher_threads,
                matching_method=matching_method
            )

        except InterruptedError:
            logging.info(f"🛑 Processing cancelled for file: {input_file}")
            update_queue.put(("status_update", file_key, "🚫 Cancelled"))

        except Exception as e:
            logging.error(f"❌ {input_file} Error during processing: {e}")
            update_queue.put(("status_update", file_key, f"❌ Error: {e}"))

    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        futures = {}

        for input_file in selected_files:
            file_key = os.path.basename(input_file)

            if not os.path.exists(input_file):
                completed_files["had_error"] = True
            
                update_queue.put((
                    "status_update",
                    file_key,
                    "❌ Error: Input file was not found."
                ))
            
                with lock:
                    completed_files["count"] += 1
            
                continue
    
            output_file_specific = os.path.join(
                output_dir_path,
                f"{os.path.splitext(file_key)[0]}_results.txt"
            )
    
            update_queue.put(("progress_max", file_key, 100))
            update_queue.put(("progress_update", file_key, 0))
            update_queue.put(("status_update", file_key, "Stand-by"))
    
            if analysis_cancelled["flag"]:
                update_queue.put(("status_update", file_key, "🚫 Cancelled"))
                continue

            future = executor.submit(run_file, input_file, output_file_specific)
            futures[future] = file_key
    
        try:
            pending = set(futures.keys())

            while pending:
                if analysis_cancelled["flag"]:
                    for future in pending:
                        file_key = futures[future]

                        if not future.running() and not future.done():
                            if future.cancel():
                                update_queue.put(("status_update", file_key, "🚫 Cancelled"))
                            else:
                                update_queue.put(("status_update", file_key, "🛑 Cancelling... Please wait."))
                        else:
                            update_queue.put(("status_update", file_key, "🛑 Cancelling... Please wait."))

                    break

                done, pending = wait(
                    pending,
                    timeout=0.2,
                    return_when=FIRST_COMPLETED
                )

                for future in done:
                    future.result()
    
        except Exception as e:
            logging.error(f"❌ Error occurred during analysis: {e}")
            
    update_queue.put(("analysis_complete", "", "All analyses have been completed!"))

def prepare_probe_match_items(probes):
    probe_items = []

    for probe_name, probe_seq in probes:
        logging.debug(f"Preparing probe: {probe_name}")

        probe_seq = probe_seq.upper()
        possible_sequences = generate_possible_sequences(probe_seq)

        possible_sequences_rc = [
            str(reverse_complement(seq))
            for seq in possible_sequences
        ]

        probe_items.append((
            probe_name,
            len(probe_seq),
            possible_sequences,
            possible_sequences_rc
        ))

    return probe_items

def generate_possible_sequences(iupac_sequence):
    logging.debug("IUPAC code conversion started")
    sequence_lists = []

    for base in iupac_sequence:
        if base in IUPAC_CODES:
            sequence_lists.append(IUPAC_CODES[base])
        else:
            raise ValueError(f"Unrecognized base '{base}' in probe sequence.")

    logging.debug(f"IUPAC code conversion completed: {iupac_sequence}")
    return [''.join(p) for p in itertools.product(*sequence_lists)]

def validate_probes_before_analysis(probe_text: str) -> tuple[bool, set]:
    IUPAC_CODES_SET = set(IUPAC_CODES.keys())
    invalid_bases_found = set()
    seen_probe_names = set()

    lines = probe_text.strip().splitlines()
    is_in_sequence = False
    expecting_sequence = False

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            is_in_sequence = True
            expecting_sequence = True

            probe_name = line[1:].strip()
            if not probe_name:
                invalid_bases_found.add("Empty probe name after '>'")
            elif probe_name in seen_probe_names:
                invalid_bases_found.add(f"Duplicate probe name: {probe_name}")
            else:
                seen_probe_names.add(probe_name)

        else:
            if not is_in_sequence:
                invalid_bases_found.add("Invalid FASTA format: missing header")
                continue

            if not line.isalpha():
                invalid_bases_found.add("Invalid characters in sequence line")
                continue

            if expecting_sequence:
                expecting_sequence = False

            for base in line.upper():
                if base not in IUPAC_CODES_SET:
                    invalid_bases_found.add(base)

    if expecting_sequence:
        invalid_bases_found.add("Missing sequence after >header")

    return (len(invalid_bases_found) == 0, invalid_bases_found)

def validate_gzip_file(file_path, analysis_cancelled=None):
    try:
        if analysis_cancelled and analysis_cancelled.get("flag"):
            return False
        with gzip.open(file_path, 'rt') as f:
            f.read(1)
        return True
    except Exception as e:
        logging.error(f"Invalid gzip file: {file_path}, Error: {e}")
        return False

def validate_fastq_format(file_path, analysis_cancelled=None):
    try:
        with gzip.open(file_path, "rt") if file_path.endswith(".gz") else open(file_path, "rt") as f:
            line_count = 0
            while True:
                if analysis_cancelled and analysis_cancelled.get("flag"):
                    raise ValueError("Analysis cancelled by user")

                header = f.readline().strip()
                if not header:
                    break
                sequence = f.readline().strip()
                plus_line = f.readline().strip()
                quality = f.readline().strip()
                line_count += 4

                if not header.startswith("@"):
                    raise ValueError(f"Invalid FASTQ header at line {line_count - 3}: {header}")
                if not plus_line.startswith("+"):
                    raise ValueError(f"Invalid '+' line at line {line_count - 1}: {plus_line}")
                if len(sequence) != len(quality):
                    raise ValueError(f"Sequence and quality lengths do not match at line {line_count}")
        logging.info(f"FASTQ format validation passed for file: {file_path}")
    except Exception as e:
        logging.error(f"FASTQ format validation failed: {e}")
        raise

def open_sequence_text_stream(file_path):
    if file_path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(file_path, "rb"), encoding="utf-8")
    return open(file_path, "rt", encoding="utf-8")

def shutdown_process_pool_fast(executor, pending_futures=None):
    if executor is None:
        return

    if pending_futures:
        for future in list(pending_futures):
            try:
                future.cancel()
            except Exception:
                pass

    processes = []
    try:
        processes = list(getattr(executor, "_processes", {}).values())
    except Exception:
        processes = []

    for process in processes:
        try:
            if process.is_alive():
                process.terminate()
        except Exception:
            pass

    try:
        executor.shutdown(wait=False, cancel_futures=True)
    except TypeError:
        executor.shutdown(wait=False)
    except Exception:
        pass

def FASTA_sampling(input_file, sampling_percentage, update_queue=None, min_reads=1, analysis_cancelled=None):
    temp_output_file = None

    try:
        base_filename = (
            os.path.basename(input_file)
            .replace(".gz", "")
            .replace(".fasta", "")
            .replace(".fa", "")
            .replace(".fna", "")
        )
        unique_id = uuid.uuid4().hex[:8]
        temp_output_file = os.path.join(
            tempfile.gettempdir(),
            f"tmp_{base_filename}_{unique_id}.fasta.gz"
        )

        if analysis_cancelled and analysis_cancelled.get("flag"):
            raise InterruptedError("FASTA sampling cancelled before starting.")

        if input_file.endswith(".gz") and not validate_gzip_file(input_file, analysis_cancelled):
            raise ValueError(f"Invalid gzip file: {input_file}")

        try:
            set_busy_progress(
                update_queue,
                input_file,
                "Counting sequences..."
            )

            total_sequences = get_total_reads_with_seqkit(input_file, analysis_cancelled)
            
            sample_size = max(int(total_sequences * (sampling_percentage / 100)), min_reads)
            sample_size = min(sample_size, total_sequences)

            if sample_size <= 0:
                raise ValueError(f"No sequences available for sampling: {input_file}")

            set_busy_progress(
                update_queue,
                input_file,
                f"Sampling {sample_size:,} sequences..."
            )

            if analysis_cancelled and analysis_cancelled.get("flag"):
                raise InterruptedError("FASTA sampling cancelled before SeqKit sampling.")

            run_seqkit_sample2(input_file, temp_output_file, sample_size, analysis_cancelled)

            if analysis_cancelled and analysis_cancelled.get("flag"):
                raise InterruptedError("FASTA sampling cancelled after SeqKit sampling.")

            if not os.path.exists(temp_output_file) or os.path.getsize(temp_output_file) == 0:
                raise ValueError("SeqKit FASTA sampling produced an empty output file.")

            if update_queue:
                update_queue.put((
                        "status_update",
                        input_file,
                        "Sampling completed. Preparing matching..."
                ))

            logging.info(f"✅ SeqKit FASTA sampling completed: {sample_size} sequences → {temp_output_file}")
            return temp_output_file, sample_size

        except FileNotFoundError:
            logging.warning("SeqKit not found. Falling back to Python FASTA sampling.")
            if update_queue:
                update_queue.put(("status_update", input_file, "SeqKit not found. Using standard Python FASTA sampling."))

        except InterruptedError:
            raise

        except subprocess.CalledProcessError as e:
            logging.warning(f"SeqKit FASTA sampling failed. Falling back to Python sampling: {e.stderr}")
            if update_queue:
                update_queue.put(("status_update", input_file, "SeqKit failed. Using standard Python FASTA sampling."))

        except Exception as e:
            logging.warning(f"Unexpected SeqKit FASTA sampling error. Falling back to Python sampling: {e}")
            if update_queue:
                update_queue.put(("status_update", input_file, "SeqKit error. Using standard Python FASTA sampling."))

        set_busy_progress(
            update_queue,
            input_file,
            "Counting sequences...",
        )

        total_sequences = 0
        with open_sequence_text_stream(input_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if analysis_cancelled and analysis_cancelled.get("flag"):
                    raise InterruptedError("FASTA sampling cancelled during counting.")
                total_sequences += 1

        sample_size = max(int(total_sequences * (sampling_percentage / 100)), min_reads)
        sample_size = min(sample_size, total_sequences)

        if sample_size <= 0:
            raise ValueError(f"No sequences available for sampling: {input_file}")

        set_busy_progress(
            update_queue,
            input_file,
            f"Sampling {sample_size:,} sequences...",
        )

        selected_indices = set(sorted(random.sample(range(total_sequences), sample_size)))
        sampled_count = 0

        with open_sequence_text_stream(input_file) as handle, gzip.open(temp_output_file, "wt") as output_handle:
            for i, record in enumerate(SeqIO.parse(handle, "fasta")):
                if analysis_cancelled and analysis_cancelled.get("flag"):
                    raise InterruptedError("FASTA sampling cancelled during writing.")

                if i in selected_indices:
                    SeqIO.write(record, output_handle, "fasta")
                    sampled_count += 1

        if update_queue:
            update_queue.put((
                "status_update",
                input_file,
                "Sampling completed. Preparing matching...",
            ))

        logging.info(f"✅ Python FASTA sampling completed: {sample_size} sequences → {temp_output_file}")
        return temp_output_file, sample_size

    except InterruptedError:
        if temp_output_file and os.path.exists(temp_output_file):
            try:
                os.remove(temp_output_file)
                logging.warning(f"🛑 FASTA sampling cancelled. Temp file removed: {temp_output_file}")
            except Exception as e:
                logging.warning(f"⚠️ Failed to remove FASTA temp file during cancel: {e}")

        if update_queue:
            update_queue.put(("status_update", input_file, "🚫 Cancelled"))

        return None, None

    except Exception as e:
        if temp_output_file and os.path.exists(temp_output_file):
            try:
                os.remove(temp_output_file)
                logging.info(f"🧹 Removed FASTA temp file after sampling error: {temp_output_file}")
            except Exception as err:
                logging.warning(f"⚠️ Failed to remove FASTA temp file after error: {err}")

        logging.error(f"❌ Error during FASTA sampling: {e}")

        if update_queue:
            if analysis_cancelled and analysis_cancelled.get("flag"):
                update_queue.put(("status_update", input_file, "🚫 Cancelled"))
            else:
                update_queue.put(("status_update", input_file, f"❌ Error: {e}"))

        return None, None

def FASTA_detect_and_quantify_probe_batched(
    input_file, probes, probe_items, threshold, sampling_percentage, output_file_path,
    completed_files, total_files, summary_file_path,
    results, output_dir_path,
    analysis_cancelled, update_queue,
    matcher_threads=None,
    matching_method="fast_biological"
):
    def fasta_batch_generator(handle, batch_size):
        batch = []
        for record in SeqIO.parse(handle, "fasta"):
            batch.append(record)
            if len(batch) == batch_size:
                yield batch
                batch = []
        if batch:
            yield batch
            
    processed_so_far = 0
    total_records = 0
    cancelled_early = False
    sampled_file = None
    sampled_count = None
    had_error = False
    
    try:
        start_time = datetime.now()
        original_filename = input_file

        file_to_use = input_file
        
        if sampling_percentage != 100:
            sampled_file, sampled_count = FASTA_sampling(
                input_file,
                sampling_percentage,
                update_queue=update_queue,
                min_reads=1,
                analysis_cancelled=analysis_cancelled
            )
        
            if not sampled_file:
                if analysis_cancelled.get("flag"):
                    cancelled_early = True
                    update_queue.put((
                        "status_update",
                        input_file,
                        "🚫 Cancelled"
                    ))
                else:
                    had_error = True
                    completed_files["had_error"] = True
            
                return
        
            file_to_use = sampled_file
            total_records = sampled_count

        if analysis_cancelled["flag"]:
            cancelled_early = True
            return

        if should_use_issm_matcher(threshold, sampling_percentage, matching_method):
            known_total_records = (
                total_records
                if total_records and total_records > 0
                else 0
            )
            try:
                _, matcher_total_records, _ = run_issm_matcher_exact_and_update(
                    input_file=original_filename,
                    matcher_input_file=file_to_use,
                    probes=probes,
                    total_records=known_total_records,
                    output_file_path=output_file_path,
                    output_dir_path=output_dir_path,
                    summary_file_path=summary_file_path,
                    results=results,
                    selected_format="fasta",
                    original_filename=original_filename,
                    start_time=start_time,
                    analysis_cancelled=analysis_cancelled,
                    update_queue=update_queue,
                    matcher_threads=matcher_threads,
                    threshold=threshold,
                    matching_method=matching_method,
                )
                processed_so_far = matcher_total_records or known_total_records or 0
                return
            except InterruptedError:
                cancelled_early = True
                return
            except Exception as e:
                logging.warning(f"ISSM fast matcher FASTA mode failed. Falling back to Python matcher: {e}")
                update_queue.put(("status_update", input_file, f"ISSM fast matcher failed. Falling back to Python matcher."))
        
        if not probe_items:
            logging.info("Preparing Python probe variants for fallback/legacy FASTA matching.")
            probe_items = prepare_probe_match_items(probes)

        if sampling_percentage == 100:
            set_busy_progress(
                update_queue,
                input_file,
                "Counting sequences...",
            )
        
            try:
                total_records = get_total_reads_with_seqkit(
                    file_to_use,
                    analysis_cancelled
                )

            except InterruptedError:
                cancelled_early = True
                return

            except Exception as e:
                logging.warning(
                    f"SeqKit FASTA count failed. Falling back to Python count: {e}"
                )
                update_queue.put((
                    "status_update",
                    input_file,
                    "SeqKit count failed. Counting sequences with Python."
                ))

                total_records = 0
                with open_sequence_text_stream(file_to_use) as handle:
                    for _ in SeqIO.parse(handle, "fasta"):
                        if analysis_cancelled.get("flag"):
                            cancelled_early = True
                            return
                        total_records += 1

        if total_records <= 0:
            raise ValueError(f"❌ No valid sequences found in {input_file}")
        
        probe_names = np.array([name for name, _ in probes])
        probe_dict = {name: seq for name, seq in probes}

        BATCH_SIZE = 300

        emit_percent_progress(
            update_queue=update_queue,
            input_file=input_file,
            phase="Matching",
            processed=0,
            total=total_records,
            unit="sequences",
            reset=True,
        )
        
        probe_counts = Counter({name: 0 for name in probe_names})
        matched_sequences = set()
        read_to_probes = {}
        total_matched_sequences = 0

        cpu_count = os.cpu_count() or 1
        num_workers = max(1, min(int(cpu_count * 0.75), 61))
        safe_cores = list(range(num_workers))
        
        with open(output_file_path, 'w', encoding="utf-8") as result_file:
            result_file.write("Metric\tValue\n")
            result_file.write(f"Total Selected Sequences\t{total_records}\n")
            result_file.write("Total Matched Sequences\t0\n")
            result_file.write("Probe Name\tProbe Sequence\tMatched Sequences\tPercentage of Total Sequences\n")
            for probe_name in probe_counts:
                probe_seq = probe_dict.get(probe_name, "N/A")
                result_file.write(f"{probe_name}\t{probe_seq}\t0\t0.00%\n")
            result_file.write("\n")
            result_file.write("\nMatched Sequences Information\n")
            result_file.write("Read ID\tProbe Name\tScore Original\tScore Reverse Complement\n")

        fasta_path = os.path.join(
            output_dir_path,
            f"{os.path.splitext(os.path.basename(input_file))[0]}_matched_reads.fasta"
        )
        fasta_fh = None
        written_reads = set()

        last_progress_update = 0
        progress_update_interval = max(
            BATCH_SIZE,
            total_records // 500,
        )
        completed_batches = set()
        pending_futures = set()
        max_in_flight = max(4, min(num_workers * 2, 64))
        executor = None

        try:
            executor = ProcessPoolExecutor(
                max_workers=num_workers,
                initializer=init_worker_all,
                initargs=(probe_items, safe_cores)
            )

            with open_sequence_text_stream(file_to_use) as fasta_handle, \
                 open(output_file_path, 'a', encoding="utf-8") as result_file:

                batch_iterator = enumerate(
                    fasta_batch_generator(fasta_handle, BATCH_SIZE),
                    start=1
                )
                input_exhausted = False

                while pending_futures or not input_exhausted:
                    if analysis_cancelled.get("flag"):
                        cancelled_early = True
                        update_queue.put((
                            "status_update",
                            input_file,
                            "🛑 Cancelling..."
                        ))
                        break

                    while (
                        not input_exhausted
                        and len(pending_futures) < max_in_flight
                    ):
                        try:
                            batch_index, record_batch = next(batch_iterator)
                        except StopIteration:
                            input_exhausted = True
                            break

                        if analysis_cancelled.get("flag"):
                            cancelled_early = True
                            break

                        pending_futures.add(
                            executor.submit(
                                match_batch_sequences,
                                record_batch,
                                threshold,
                                batch_index
                            )
                        )

                    if cancelled_early:
                        break

                    if not pending_futures:
                        continue

                    done, _ = wait(
                        pending_futures,
                        timeout=0.1,
                        return_when=FIRST_COMPLETED
                    )

                    if not done:
                        continue

                    pending_futures.difference_update(done)

                    for future in done:
                        if analysis_cancelled.get("flag"):
                            cancelled_early = True
                            break

                        try:
                            matches, batch_idx, batch_size = future.result()
                        except Exception as e:
                            had_error = True
                            logging.exception(
                                f"❌ Exception during FASTA batch processing: {e}"
                            )
                            update_queue.put((
                                "status_update",
                                input_file,
                                f"❌ Batch error: {e}"
                            ))
                            break

                        if batch_idx not in completed_batches:
                            processed_so_far += batch_size
                            completed_batches.add(batch_idx)

                        for match in matches:
                            read_id, probe_name, score_orig, score_rc, read_seq = match
                            probe_counts[probe_name] += 1
                            matched_sequences.add(read_id)
                            read_to_probes.setdefault(read_id, set()).add(probe_name)
                            result_file.write(
                                f"{read_id}\t{probe_name}\t{score_orig}\t{score_rc}\n"
                            )

                            if read_id not in written_reads:
                                if fasta_fh is None:
                                    try:
                                        fasta_fh = open(
                                            fasta_path,
                                            "w",
                                            encoding="utf-8"
                                        )
                                    except Exception as e:
                                        logging.exception(
                                            f"Failed to open FASTA for writing: "
                                            f"{fasta_path} — {e}"
                                        )
                                        fasta_fh = None

                                if fasta_fh is not None:
                                    fasta_fh.write(f">{read_id}\n{read_seq}\n")
                                    written_reads.add(read_id)

                        if (
                            processed_so_far - last_progress_update
                            >= progress_update_interval
                        ):
                            emit_percent_progress(
                                update_queue=update_queue,
                                input_file=input_file,
                                phase="Matching",
                                processed=processed_so_far,
                                total=total_records,
                                unit="sequences",
                            )
                            last_progress_update = processed_so_far

                    if cancelled_early or had_error:
                        break

        finally:
            if executor is not None:
                if (
                    cancelled_early
                    or had_error
                    or analysis_cancelled.get("flag")
                ):
                    shutdown_process_pool_fast(
                        executor,
                        pending_futures
                    )
                else:
                    executor.shutdown(wait=True)

            if fasta_fh is not None:
                try:
                    fasta_fh.close()
                except Exception:
                    pass

            if not written_reads and os.path.exists(fasta_path):
                try:
                    os.remove(fasta_path)
                    logging.debug(
                        f"No matched reads — removed empty FASTA: {fasta_path}"
                    )
                except Exception as e:
                    logging.warning(
                        f"Failed to remove empty FASTA ({fasta_path}): {e}"
                    )

        total_matched_sequences = len(matched_sequences)

        if cancelled_early or analysis_cancelled.get("flag") or had_error:
            return

        with open(output_file_path, 'r', encoding="utf-8") as f:
            lines = f.readlines()

        with open(output_file_path, 'w', encoding="utf-8") as f:
            i = 0
            while i < len(lines):
                line = lines[i]

                if line.startswith("Total Matched Sequences"):
                    f.write(f"Total Matched Sequences\t{total_matched_sequences}\n")

                elif line.startswith("Probe Name\t"):
                    f.write(line)
                    i += 1
                    while i < len(lines) and not lines[i].startswith("Matched Sequences Information"):
                        i += 1

                    for probe_name in probe_counts:
                        count = probe_counts[probe_name]
                        probe_seq = probe_dict.get(probe_name, "N/A")
                        percentage = (count / total_records) * 100 if total_records > 0 else 0
                        f.write(f"{probe_name}\t{probe_seq}\t{count}\t{percentage:.2f}%\n")
                    continue

                else:
                    f.write(line)

                i += 1

        end_time = datetime.now()
        elapsed_time = str(end_time - start_time).split(".")[0]

        with lock:
            results["file_results"][original_filename] = {"sampled_records": total_records}
            for probe_name in probe_counts:
                count = probe_counts[probe_name]
                results["file_results"][original_filename][probe_name] = count
                results["probe_results"].setdefault(probe_name, {})[original_filename] = count

            if not os.path.exists(summary_file_path) or os.path.getsize(summary_file_path) == 0:
                with open(summary_file_path, 'w', encoding="utf-8") as summary_file:
                    summary_file.write("Index\tFilename\tRecordCount\tStartTime\tEndTime\tElapsedTime\n")

            with open(summary_file_path, 'r', encoding="utf-8") as f:
                existing_lines = f.readlines()
                index_counter = len(existing_lines) if len(existing_lines) > 1 else 1

            with open(summary_file_path, 'a', encoding="utf-8") as summary_file:
                summary_file.write(
                    f"{index_counter}\t{os.path.basename(original_filename)}\t{total_records}\t{start_time}\t{end_time}\t{elapsed_time}\n"
                )

    except Exception as e:
        had_error = True
        completed_files["had_error"] = True
        logging.exception(f"❌ Error processing FASTA file: {e}")
        update_queue.put(("status_update", input_file, f"❌ Error: {e}"))

    finally:
        if sampled_file and sampled_file != input_file and os.path.exists(sampled_file):
            try:
                os.remove(sampled_file)
                logging.info(f"Successfully removed temporary FASTA file: {sampled_file}")
            except Exception as e:
                logging.warning(f"Failed to remove temporary FASTA file: {sampled_file}, reason: {e}")
    
        with lock:
            completed_files['count'] += 1
            is_last_file = completed_files.get('count', 0) >= total_files

        if (
            total_records > 0
            and not cancelled_early
            and not analysis_cancelled.get("flag")
            and not had_error
        ):
            processed_so_far = min(processed_so_far, total_records)
            emit_percent_progress(
                update_queue=update_queue,
                input_file=input_file,
                phase="Matching",
                processed=processed_so_far,
                total=total_records,
                unit="sequences",
            )

        if cancelled_early or analysis_cancelled.get("flag"):
            completed_files["cancelled"] = True
            update_queue.put(("status_update", input_file, "🚫 Cancelled"))

        elif had_error:
            completed_files["had_error"] = True

        else:
            if is_last_file:
                update_queue.put(("status_update", input_file, "📊 Generating plots and summary..."))
            try_generate_plots_once(
                completed_files,
                lock,
                total_files,
                analysis_cancelled,
                results,
                output_dir_path,
                probes
            )
            update_queue.put(("status_update", input_file, "✅ Completed"))

def FASTQ_sampling(input_file, sampling_percentage, file_format, update_queue=None, min_reads=5000, analysis_cancelled=None):
    temp_output_file = None

    try:
        base_filename = (
            os.path.basename(input_file)
            .replace(".gz", "")
            .replace(".fastq", "")
            .replace(".fq", "")
            .replace(".fasta", "")
            .replace(".fa", "")
        )
        unique_id = uuid.uuid4().hex[:8]
        temp_output_file = os.path.join(
            tempfile.gettempdir(),
            f"tmp_{base_filename}_{unique_id}.fastq.gz"
        )

        def open_text_stream(file_path):
            if file_path.endswith(".gz"):
                return io.TextIOWrapper(gzip.open(file_path, "rb"), encoding="utf-8")
            return open(file_path, "rt", encoding="utf-8")

        if analysis_cancelled and analysis_cancelled.get("flag"):
            raise InterruptedError("Sampling cancelled before starting.")

        if input_file.endswith(".gz") and not validate_gzip_file(input_file):
            raise ValueError(f"Invalid gzip file: {input_file}")

        try:
            set_busy_progress(
                update_queue,
                input_file,
                "Counting reads..."
            )
            
            total_reads = get_total_reads_with_seqkit(input_file, analysis_cancelled)

            sample_size = max(int(total_reads * (sampling_percentage / 100)), min_reads)
            sample_size = min(sample_size, total_reads)

            if sample_size <= 0:
                raise ValueError(f"No reads available for sampling: {input_file}")

            set_busy_progress(
                update_queue,
                input_file,
                f"Sampling {sample_size:,} reads..."
            )

            if analysis_cancelled and analysis_cancelled.get("flag"):
                raise InterruptedError("Sampling cancelled before SeqKit sampling.")

            run_seqkit_sample2(input_file, temp_output_file, sample_size, analysis_cancelled)

            if analysis_cancelled and analysis_cancelled.get("flag"):
                raise InterruptedError("Sampling cancelled after SeqKit sampling.")

            if not os.path.exists(temp_output_file) or os.path.getsize(temp_output_file) == 0:
                raise ValueError("SeqKit sampling produced an empty output file.")

            if update_queue:
                update_queue.put((
                    "status_update",
                    input_file,
                    "Sampling completed. Preparing matching...",
                ))


            logging.info(f"✅ SeqKit sampling completed: {sample_size} reads → {temp_output_file}")

            return temp_output_file, sample_size

        except FileNotFoundError:
            logging.warning("SeqKit not found. Falling back to Python sampling.")
            if update_queue:
                update_queue.put(("status_update", input_file, "SeqKit not found. Using standard Python sampling."))

        except InterruptedError:
            raise

        except subprocess.CalledProcessError as e:
            logging.warning(f"SeqKit sampling failed. Falling back to Python sampling: {e.stderr}")
            if update_queue:
                update_queue.put(("status_update", input_file, "SeqKit failed. Using standard Python sampling."))

        except Exception as e:
            logging.warning(f"Unexpected SeqKit sampling error. Falling back to Python sampling: {e}")
            if update_queue:
                update_queue.put(("status_update", input_file, "SeqKit error. Using standard Python sampling."))

        if file_format == "fastq":
            validate_fastq_format(input_file)

        set_busy_progress(
            update_queue,
            input_file,
            "Counting reads..."
        )

        total_reads = 0
        with open_text_stream(input_file) as handle:
            for record in SeqIO.parse(handle, file_format):
                if analysis_cancelled and analysis_cancelled.get("flag"):
                    raise InterruptedError("Sampling cancelled during counting.")
                total_reads += 1

        sample_size = max(int(total_reads * (sampling_percentage / 100)), min_reads)
        sample_size = min(sample_size, total_reads)

        if sample_size <= 0:
            raise ValueError(f"No reads available for sampling: {input_file}")

        set_busy_progress(
            update_queue,
            input_file,
            f"Sampling {sample_size:,} reads..."
        )

        selected_indices = set(sorted(random.sample(range(total_reads), sample_size)))
        sampled_count = 0

        with open_text_stream(input_file) as handle, gzip.open(temp_output_file, "wt") as output_handle:
            for i, record in enumerate(SeqIO.parse(handle, file_format)):
                if analysis_cancelled and analysis_cancelled.get("flag"):
                    raise InterruptedError("Sampling cancelled during writing.")

                if i in selected_indices:
                    SeqIO.write(record, output_handle, file_format)
                    sampled_count += 1

        if update_queue:
            update_queue.put((
                "status_update",
                input_file,
                "Sampling completed. Preparing matching...",
            ))
            
        logging.info(f"✅ Python sampling completed: {sample_size} reads → {temp_output_file}")

        return temp_output_file, sample_size

    except InterruptedError:
        if temp_output_file and os.path.exists(temp_output_file):
            try:
                os.remove(temp_output_file)
                logging.warning(f"🛑 FASTQ sampling cancelled. Temp file removed: {temp_output_file}")
            except Exception as e:
                logging.warning(f"⚠️ Failed to remove FASTQ temp file during cancel: {e}")

        if update_queue:
            update_queue.put(("status_update", input_file, "🚫 Cancelled"))

        return None, None

    except Exception as e:
        if temp_output_file and os.path.exists(temp_output_file):
            try:
                os.remove(temp_output_file)
                logging.info(f"🧹 Removed FASTQ temp file after sampling error: {temp_output_file}")
            except Exception as err:
                logging.warning(f"⚠️ Failed to remove FASTQ temp file after error: {err}")

        logging.error(f"❌ Error during FASTQ sampling: {e}")

        if update_queue:
            if analysis_cancelled and analysis_cancelled.get("flag"):
                update_queue.put(("status_update", input_file, "🚫 Cancelled"))
            else:
                update_queue.put(("status_update", input_file, f"❌ Error: {e}"))

        return None, None

def FASTQ_detect_and_quantify_probe(
    input_file, probes, probe_items, threshold, sampling_percentage, output_file_path,
    completed_files, total_files,
    summary_file_path, results, output_dir_path, analysis_cancelled, update_queue,
    matcher_threads=None,
    matching_method="fast_biological"
):
    def batch_generator(handle, batch_size):
        batch = []
        for record in SeqIO.parse(handle, "fastq"):
            batch.append((record.id, str(record.seq)))
            if len(batch) == batch_size:
                yield batch
                batch = []
        if batch:
            yield batch

    processed_so_far = 0
    cancelled_early = False
    sampled_file = None
    sampled_count = None
    total_records = 0
    had_error = False

    try:
        BATCH_SIZE = 500
        start_time = datetime.now()
        original_filename = input_file
        file_format = "fastq"

        file_to_use = input_file
        
        if sampling_percentage != 100:
            sampled_file, sampled_count = FASTQ_sampling(
                input_file, sampling_percentage, file_format,
                update_queue, 5000, analysis_cancelled
            )
            
            if not sampled_file:
                if analysis_cancelled.get("flag"):
                    cancelled_early = True
                    update_queue.put((
                        "status_update",
                        input_file,
                        "🚫 Cancelled"
                    ))
                else:
                    had_error = True
                    completed_files["had_error"] = True
            
                return
        
            file_to_use = sampled_file
            total_records = sampled_count
        

        if should_use_issm_matcher(threshold, sampling_percentage, matching_method):
            known_total_records = (
                total_records
                if total_records and total_records > 0
                else 0
            )
            try:
                _, matcher_total_records, _ = run_issm_matcher_exact_and_update(
                    input_file=original_filename,
                    matcher_input_file=file_to_use,
                    probes=probes,
                    total_records=known_total_records,
                    output_file_path=output_file_path,
                    output_dir_path=output_dir_path,
                    summary_file_path=summary_file_path,
                    results=results,
                    selected_format="fastq",
                    original_filename=original_filename,
                    start_time=start_time,
                    analysis_cancelled=analysis_cancelled,
                    update_queue=update_queue,
                    matcher_threads=matcher_threads,
                    threshold=threshold,
                    matching_method=matching_method,
                )
                processed_so_far = matcher_total_records or 0
                return
            except InterruptedError:
                cancelled_early = True
                return
            except Exception as e:
                logging.warning(f"ISSM fast matcher FASTQ mode failed. Falling back to Python matcher: {e}")
                update_queue.put(("status_update", input_file, f"ISSM fast matcher failed. Falling back to Python matcher."))

        if not probe_items:
            logging.info("Preparing Python probe variants for fallback/legacy FASTQ matching.")
            probe_items = prepare_probe_match_items(probes)

        if sampling_percentage == 100:
            set_busy_progress(
                update_queue,
                input_file,
                "Counting reads...",
            )
                
            try:
                total_records = get_total_reads_with_seqkit(
                    file_to_use,
                    analysis_cancelled
                )

            except InterruptedError:
                cancelled_early = True
                return

            except Exception as e:
                logging.warning(
                    f"SeqKit count failed. Falling back to Python count: {e}"
                )
                update_queue.put((
                    "status_update",
                    input_file,
                    "SeqKit count failed. Counting reads with Python."
                ))

                total_records = 0
                open_func = gzip.open if file_to_use.endswith(".gz") else open
                with open_func(file_to_use, "rt") as handle:
                    for _ in SeqIO.parse(handle, "fastq"):
                        if analysis_cancelled.get("flag"):
                            cancelled_early = True
                            return
                        total_records += 1

        if total_records <= 0:
            raise ValueError(f"❌ No valid reads found in {input_file}")

        emit_percent_progress(
            update_queue=update_queue,
            input_file=input_file,
            phase="Matching",
            processed=0,
            total=total_records,
            unit="reads",
            reset=True,
        )

        probe_counts = Counter({name: 0 for name, _ in probes})
        matched_reads = set()
        read_to_probes = {}
        probe_dict = {name: seq for name, seq in probes}
        total_matched_reads = 0
        cpu_count = os.cpu_count() or 1
        num_workers = max(1, min(int(cpu_count * 0.75), 61))
        fasta_path = os.path.join(
            output_dir_path,
            f"{os.path.splitext(os.path.basename(input_file))[0]}_matched_reads.fasta"
        )
        fasta_fh = None
        written_reads = set() 

        with open(output_file_path, 'w', encoding="utf-8") as result_file:
            result_file.write("Metric\tValue\n")
            result_file.write(f"Total Selected Reads\t{total_records}\n")
            result_file.write("Total Matched Reads\t0\n")
            result_file.write("Probe Name\tProbe Sequence\tMatched Reads\tPercentage of Total Reads\n")
            for probe_name in probe_counts:
                probe_seq = probe_dict.get(probe_name, "N/A")
                result_file.write(f"{probe_name}\t{probe_seq}\t0\t0.00%\n")
            result_file.write("\n")
            result_file.write("\nMatched Reads Information\n")
            result_file.write("Read ID\tProbe Name\tScore Original\tScore Reverse Complement\n")

        last_progress_update = 0
        progress_update_interval = max(
            10000,
            total_records // 500,
        )
        completed_batches = set()
        pending_futures = set()
        max_in_flight = max(4, min(num_workers * 2, 64))
        executor = None

        try:
            executor = ProcessPoolExecutor(
                max_workers=num_workers,
                initializer=init_worker,
                initargs=(probe_items,)
            )

            open_func = gzip.open if file_to_use.endswith(".gz") else open
            with open_func(file_to_use, "rt") as handle, \
                 open(output_file_path, 'a', encoding="utf-8") as result_file:

                batch_iterator = enumerate(
                    batch_generator(handle, BATCH_SIZE),
                    start=1
                )
                input_exhausted = False

                while pending_futures or not input_exhausted:
                    if analysis_cancelled.get("flag"):
                        cancelled_early = True
                        update_queue.put((
                            "status_update",
                            input_file,
                            "🛑 Cancelling..."
                        ))
                        break

                    while (
                        not input_exhausted
                        and len(pending_futures) < max_in_flight
                    ):
                        try:
                            batch_index, record_batch = next(batch_iterator)
                        except StopIteration:
                            input_exhausted = True
                            break

                        if analysis_cancelled.get("flag"):
                            cancelled_early = True
                            break

                        pending_futures.add(
                            executor.submit(
                                match_batch_reads,
                                record_batch,
                                threshold,
                                batch_index
                            )
                        )

                    if cancelled_early:
                        break

                    if not pending_futures:
                        continue

                    done, _ = wait(
                        pending_futures,
                        timeout=0.1,
                        return_when=FIRST_COMPLETED
                    )

                    if not done:
                        continue

                    pending_futures.difference_update(done)

                    for future in done:
                        if analysis_cancelled.get("flag"):
                            cancelled_early = True
                            break

                        try:
                            matches, batch_idx, batch_size = future.result()
                        except Exception as e:
                            had_error = True
                            logging.exception(
                                f"❌ Exception during FASTQ batch processing: {e}"
                            )
                            update_queue.put((
                                "status_update",
                                input_file,
                                f"❌ Batch error: {e}"
                            ))
                            break

                        if batch_idx not in completed_batches:
                            processed_so_far += batch_size
                            completed_batches.add(batch_idx)

                        for match in matches:
                            read_id, probe_name, score_orig, score_rc, read_seq = match
                            probe_counts[probe_name] += 1
                            matched_reads.add(read_id)
                            read_to_probes.setdefault(read_id, set()).add(probe_name)
                            result_file.write(
                                f"{read_id}\t{probe_name}\t{score_orig}\t{score_rc}\n"
                            )

                            if read_id not in written_reads:
                                if fasta_fh is None:
                                    try:
                                        fasta_fh = open(
                                            fasta_path,
                                            "w",
                                            encoding="utf-8"
                                        )
                                    except Exception as e:
                                        logging.exception(
                                            f"Failed to open FASTA for writing: "
                                            f"{fasta_path} — {e}"
                                        )
                                        fasta_fh = None

                                if fasta_fh is not None:
                                    fasta_fh.write(f">{read_id}\n{read_seq}\n")
                                    written_reads.add(read_id)

                        if (
                            processed_so_far - last_progress_update
                            >= progress_update_interval
                        ):
                            emit_percent_progress(
                                update_queue=update_queue,
                                input_file=input_file,
                                phase="Matching",
                                processed=processed_so_far,
                                total=total_records,
                                unit="reads",
                            )
                            last_progress_update = processed_so_far

                    if cancelled_early or had_error:
                        break

        finally:
            if executor is not None:
                if (
                    cancelled_early
                    or had_error
                    or analysis_cancelled.get("flag")
                ):
                    shutdown_process_pool_fast(
                        executor,
                        pending_futures
                    )
                else:
                    executor.shutdown(wait=True)

        total_matched_reads = len(matched_reads)

        if cancelled_early or analysis_cancelled.get("flag") or had_error:
            if fasta_fh is not None:
                try:
                    fasta_fh.close()
                except Exception:
                    pass
            return

        if fasta_fh is not None:
            try:
                fasta_fh.close()
            except Exception:
                pass
            if len(written_reads) == 0 and os.path.exists(fasta_path):
                try:
                    os.remove(fasta_path)
                    logging.debug(f"No matched reads — removed empty FASTA: {fasta_path}")
                except Exception as e:
                    logging.warning(f"Failed to remove empty FASTA ({fasta_path}): {e}")

        with open(output_file_path, 'r', encoding="utf-8") as f:
            lines = f.readlines()

        with open(output_file_path, 'w', encoding="utf-8") as f:
            i = 0
            while i < len(lines):
                line = lines[i]

                if line.startswith("Total Matched Reads"):
                    f.write(f"Total Matched Reads\t{total_matched_reads}\n")

                elif line.startswith("Probe Name\t"):
                    f.write(line)

                    i += 1
                    while i < len(lines) and not lines[i].startswith("Matched Reads Information"):
                        i += 1

                    for probe_name in probe_counts:
                        count = probe_counts[probe_name]
                        probe_seq = probe_dict.get(probe_name, "N/A")
                        percentage = (count / total_records) * 100 if total_records > 0 else 0
                        f.write(f"{probe_name}\t{probe_seq}\t{count}\t{percentage:.2f}%\n")
                    continue
                
                else:
                    f.write(line)

                i += 1

        end_time = datetime.now()
        elapsed_time = str(end_time - start_time).split(".")[0]

        async def save_summary():
            header_line = "Index\tFilename\tRecordCount\tStartTime\tEndTime\tElapsedTime\n"
            
            if not os.path.exists(summary_file_path):
                async with aiofiles.open(summary_file_path, 'w', encoding="utf-8") as summary_file:
                    await summary_file.write(header_line)
        
            async with aiofiles.open(summary_file_path, 'r', encoding="utf-8") as summary_file:
                content = await summary_file.read()
                existing_lines = content.strip().splitlines()
        
            if not existing_lines or existing_lines[0] != header_line.strip():
                content = header_line + '\n'.join(existing_lines)
                async with aiofiles.open(summary_file_path, 'w', encoding="utf-8") as summary_file:
                    await summary_file.write(content)
                existing_lines = content.strip().splitlines()
        
            index_counter = len(existing_lines) if len(existing_lines) > 1 else 1
        
            async with aiofiles.open(summary_file_path, 'a', encoding="utf-8") as summary_file:
                await summary_file.write(
                    f"{index_counter}\t{os.path.basename(original_filename)}\t{total_records}\t{start_time}\t{end_time}\t{elapsed_time}\n"
                )

        asyncio.run(save_summary())

        with lock:
            results["file_results"][original_filename] = {"sampled_records": total_records}
            for probe_name in probe_counts:
                count = probe_counts[probe_name]
                results["file_results"][original_filename][probe_name] = count
                results["probe_results"].setdefault(probe_name, {})[original_filename] = count

    except Exception as e:
        had_error = True
        completed_files["had_error"] = True
        logging.exception(f"❌ Error processing file {input_file}: {e}")
        update_queue.put(("status_update", input_file, f"❌ Error: {e}"))

    finally:
        if sampled_file and sampled_file != input_file and os.path.exists(sampled_file):
            try:
                os.remove(sampled_file)
                logging.warning(f"Successfully remove temp file: {sampled_file}")
            except Exception as e:
                logging.warning(f"🧹 Failed to remove temp file: {sampled_file}, reason: {e}")

        with lock:
            completed_files['count'] += 1
            is_last_file = completed_files.get('count', 0) >= total_files
    
        if (
            total_records > 0
            and not cancelled_early
            and not analysis_cancelled.get("flag")
            and not had_error
        ):
            processed_so_far = min(processed_so_far, total_records)
            emit_percent_progress(
                update_queue=update_queue,
                input_file=input_file,
                phase="Matching",
                processed=processed_so_far,
                total=total_records,
                unit="reads",
            )

        if cancelled_early or analysis_cancelled.get("flag"):
            completed_files["cancelled"] = True
            update_queue.put(("status_update", input_file, "🚫 Cancelled"))

        elif had_error:
            completed_files["had_error"] = True

        else:
            if is_last_file:
                update_queue.put(("status_update", input_file, "📊 Generating plots and summary..."))
            try_generate_plots_once(
                completed_files,
                lock,
                total_files,
                analysis_cancelled,
                results,
                output_dir_path,
                probes
            )
            update_queue.put(("status_update", input_file, "✅ Completed"))
       
def sanitize_filename(filename):
    return re.sub(r'[<>:"/\\|?*]', '_', filename)

def generate_all_plots(results, output_dir):
    logging.info("Generating all plots...")
    for file_name, probes in results["file_results"].items():
        plot_fastq_results(file_name, probes, output_dir)

    for probe_name, files in results["probe_results"].items():
        plot_probe_results(probe_name, files, output_dir)

def plot_probe_results(probe_name, file_counts, output_dir):
    safe_probe_name = sanitize_filename(probe_name)

    if len(file_counts) > 50:
        logging.warning(f"Skipping plot for probe '{safe_probe_name}' due to too many X-axis labels.")
        return

    probe_dir = os.path.join(output_dir, "probes")
    os.makedirs(probe_dir, exist_ok=True)

    files = [os.path.basename(f) for f in file_counts.keys()]
    counts = list(file_counts.values())

    plt.figure(figsize=(8, 5))
    plt.bar(files, counts, color='skyblue')
    plt.xlabel('FASTQ Files')
    plt.ylabel('Matched Reads')
    plt.title(f'Probe: {safe_probe_name}')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    plot_path = os.path.join(probe_dir, f"probe_{safe_probe_name}_plot.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    logging.info(f"Plot saved: {plot_path}")

def plot_fastq_results(file_name, probe_counts, output_dir):
    file_name_only = os.path.basename(file_name) 

    filtered_counts = {sanitize_filename(k): v for k, v in probe_counts.items() if k != "sampled_records"}

    if len(filtered_counts) > 50:
        logging.warning(f"Skipping plot for file '{file_name_only}' due to too many X-axis labels.")
        return

    fastq_dir = os.path.join(output_dir, "input_files")
    os.makedirs(fastq_dir, exist_ok=True)

    plt.figure(figsize=(8, 5))
    probes, counts = zip(*filtered_counts.items())

    plt.bar(probes, counts, color='orange')
    plt.xlabel('Probe')
    plt.ylabel('Matched Reads')
    plt.title(f'File: {sanitize_filename(file_name_only)}') 
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    plot_path = os.path.join(fastq_dir, f"{sanitize_filename(file_name_only)}_plot.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    logging.info(f"Plot saved: {plot_path}")

def generate_combined_results(output_dir_path, probes, results):
    total_results_path = os.path.join(output_dir_path, "total_results.txt")
    logging.debug(f"Results data: {results}")

    try:
        with open(total_results_path, 'w', encoding='utf-8') as total_file:
            total_file.write("Probe Name\tInput File Name\tMatched Reads\tPercentage of Total Analyzed Reads\n")

            for probe_name, files in results["probe_results"].items():
                for full_path, matched_reads in files.items():
                    file_name = os.path.basename(full_path) 
                    sampled_records = results["file_results"].get(full_path, {}).get("sampled_records", 0)
                    percentage = (matched_reads / sampled_records * 100) if sampled_records > 0 else 0
                    total_file.write(f"{probe_name}\t{file_name}\t{matched_reads}\t{percentage:.2f}%\n")

        logging.info(f"Total results file generated: {total_results_path}")
        generate_heatmap(total_results_path)

    except Exception as e:
        logging.error(f"❌ Error generating combined results: {e}")

def generate_heatmap(file_path):
    try:
        df = pd.read_csv(file_path, sep="\t")

        df["Percentage of Total Analyzed Reads"] = (
            df["Percentage of Total Analyzed Reads"].astype(str).str.replace("%", "").astype(float)
        )

        df["Probe Name"] = pd.Categorical(df["Probe Name"], categories=df["Probe Name"].unique(), ordered=True)

        heatmap_data = df.pivot(index="Probe Name", columns="Input File Name", values="Percentage of Total Analyzed Reads")

        num_files = len(heatmap_data.columns)
        num_probes = len(heatmap_data.index)

        fig_width = min(120, max(12, num_files * 0.5))
        fig_height = max(8, num_probes * 0.5)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        sns.heatmap(
            heatmap_data,
            cmap="Blues",
            annot=False,
            linewidths=0.5,
            cbar=True,
            ax=ax
        )

        plt.title("Heatmap of Percentage of Total Analyzed Reads", fontsize=14)
        plt.xlabel("Input File Name", fontsize=12)
        plt.ylabel("Probe Name", fontsize=12)
        plt.xticks(rotation=90)

        if num_files >= 100:
            ax2 = ax.twinx()
            ax2.set_ylim(ax.get_ylim())
            ax2.set_yticks(ax.get_yticks())
            ax2.set_yticklabels(heatmap_data.index, fontsize=10)
            ax2.set_ylabel("")

        output_dir = os.path.dirname(file_path)
        output_path = os.path.join(output_dir, "Total_heatmap.png")
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
        logging.info(f"Heatmap saved: {output_path}")

    except Exception as e:
        logging.error(f"❌ Error generating heatmap: {e}")

def try_generate_plots_once(completed_files, lock, total_files, analysis_cancelled, results, output_dir_path, probes):
    with lock:
        if completed_files.get("plotted", False):
            return

        if completed_files.get("count", 0) < total_files:
            return

        completed_files["plotted"] = True

    if analysis_cancelled.get("flag") or completed_files.get("cancelled", False):
        logging.warning("⚠️ Plot generation skipped because analysis was cancelled.")
        return

    if completed_files.get("had_error", False):
        logging.warning("⚠️ Plot generation skipped because one or more files failed.")
        return

    if not results.get("file_results") or not results.get("probe_results"):
        logging.warning("⚠️ Plot generation skipped because there are no complete results.")
        return

    logging.info("📊 Generating plots and heatmaps...")
    generate_all_plots(results, output_dir_path)
    generate_combined_results(output_dir_path, probes, results)
    logging.info("✅ Plots and heatmaps generated successfully!")

class ProgressUpdateWorker(QThread):
    progress_signal = Signal(str, str, str)

    def __init__(self, update_queue, parent=None):
        super().__init__(parent)
        self.update_queue = update_queue
        self._running = True

    def run(self):
        while self._running:
            try:
                key, file_name, value = self.update_queue.get(timeout=0.1)
                self.progress_signal.emit(file_name, key, str(value))
            except queue.Empty:
                continue
            except Exception as e:
                logging.error(f"❌ ProgressUpdateWorker error: {e}")

    def stop(self):
        self._running = False
        self.quit()
        self.wait()

class CustomProgressDialog(QDialog):
    def __init__(self, file_list, parent_gui, update_queue, analysis_cancelled):
        super().__init__()
        self.update_queue = update_queue
        self.analysis_cancelled = analysis_cancelled
        self.setWindowTitle("Processing files")
        self.setMinimumSize(950, 600)
        self.resize(950, 600)
        self.setWindowFlags(Qt.Window | Qt.WindowMinimizeButtonHint)

        self.file_bars = {}
        self.file_labels = {}
        self.parent_gui = parent_gui

        self.completed_count = 0
        self.total_count = len(file_list)
        self._updated_files = set()
        self._completion_shown = False
        self._error_files = set()

        self.completed_label = QLabel(f"Completed: 0 / {self.total_count}")
        self.completed_label.setStyleSheet("font-size: 12px; color: #555;")
        self.completed_label.setAlignment(Qt.AlignLeft)

        self.setStyleSheet("background-color: #F9FAFA;")

        main_layout = QVBoxLayout(self)

        self.tree = QTreeWidget()
        self.tree.setColumnCount(3)
        self.tree.setHeaderLabels(["Filename", "Progress", "Status"])
        self.tree.setRootIsDecorated(False)
        self.tree.setAlternatingRowColors(False)
        
        self.tree.setSelectionMode(QTreeWidget.NoSelection)
        self.tree.setFocusPolicy(Qt.NoFocus)
        
        self.tree.setStyleSheet("""
            QTreeWidget {
                font-size: 12px;
                color: #2c3e50;
                background-color: white;
                border: 1px solid #ccc;
                border-radius: 0px;
            }
            QTreeWidget::item {
                padding: 8px;
                border-bottom: 1px dashed #ddd;
            }
            QTreeWidget::item:hover {
                background: transparent;
            }
            QHeaderView::section {
                background-color: #f0f0f0;
                color: #2c3e50;
                padding: 6px;
                font-weight: bold;
                border: none;
            }
            QScrollBar:vertical {
                background: white;
                width: 6px;
                margin: 2px 0;
            }
            QScrollBar::handle:vertical {
                background: #1D83D5;
                border-radius: 3px;
                min-height: 20px;
            }
            QScrollBar::add-line:vertical,
            QScrollBar::sub-line:vertical,
            QScrollBar::add-page:vertical,
            QScrollBar::sub-page:vertical {
                background: none;
                height: 0px;
            }
            QScrollBar:horizontal {
                height: 0px;
            }
        """)
        self.tree.header().setStretchLastSection(False)
        self.tree.header().setSectionResizeMode(0, QHeaderView.Fixed)
        self.tree.header().resizeSection(0, 300)
        self.tree.header().setSectionResizeMode(1, QHeaderView.Fixed)
        self.tree.header().resizeSection(1, 260)
        self.tree.header().setSectionResizeMode(2, QHeaderView.Stretch)

        for filepath in file_list:
            filename = os.path.basename(filepath)
        
            item = QTreeWidgetItem(self.tree)
            item.setText(0, filename)
            item.setText(2, "Stand-by")
        
            item.setSizeHint(0, QSize(0, 44))

            progress_bar = QProgressBar()
            progress_bar.setValue(0)
            progress_bar.setMaximum(100)
            progress_bar.setFixedHeight(18)
            progress_bar.setStyleSheet("""
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
            progress_bar.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

            self.tree.setItemWidget(item, 1, progress_bar)
            self.file_bars[filename] = progress_bar
            self.file_labels[filename] = item

        main_layout.addWidget(self.tree)

        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setFixedWidth(100)
        self.cancel_button.setCursor(Qt.PointingHandCursor)
        self.cancel_button.clicked.connect(self.on_cancel)
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

        bottom_layout = QHBoxLayout()
        bottom_layout.setContentsMargins(20, 10, 20, 10)
        
        left_container = QWidget()
        left_layout = QHBoxLayout(left_container)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(self.completed_label)
        bottom_layout.addWidget(left_container, stretch=1)

        center_container = QWidget()
        center_layout = QHBoxLayout(center_container)
        center_layout.setContentsMargins(0, 0, 0, 0)
        center_layout.setAlignment(Qt.AlignCenter)
        center_layout.addWidget(self.cancel_button)
        bottom_layout.addWidget(center_container, stretch=1)

        right_spacer = QWidget()
        bottom_layout.addWidget(right_spacer, stretch=1)
        
        main_layout.addLayout(bottom_layout)

        self.worker_thread = ProgressUpdateWorker(update_queue)
        self.worker_thread.progress_signal.connect(self.update_progress)
        self.worker_thread.start()

    def update_status_all(self, message):
        for item in self.file_labels.values():
            item.setText(2, message)

    def closeEvent(self, event):
        if hasattr(self, 'worker_thread'):
            self.worker_thread.stop()
        super().closeEvent(event)

    def on_cancel(self):
        confirm_box = CustomConfirmBox(
            "⚠️ Cancel Analysis",
            "Are you sure you want to cancel the analysis?",
            on_close=self.handle_cancel_confirmed,
            parent=self
        )
        confirm_box.exec_()
        
    def handle_cancel_confirmed(self):
        self.analysis_cancelled["flag"] = True
        self.update_status_all("🛑 Cancelling... Please wait.")
        
    def show_completion_message(
        self,
        cancelled=False,
        had_error=False
    ):
        if self._completion_shown:
            return
    
        self._completion_shown = True
    
        if cancelled:
            message = "🚫 Analysis was cancelled by user."
        elif had_error:
            message = (
                "⚠️ Analysis finished, but one or more files "
                "could not be processed."
            )
        else:
            message = "✅ All tasks are completed!"
    
        CustomMessageBox(message, self).exec_()
        self.close()

    @Slot(str, str, str)
    def update_progress(self, file_name, key, value):
        file_key = os.path.basename(file_name)

        if key == 'progress_update':
            if file_key in self.file_bars:
                self.file_bars[file_key].setValue(int(float(value)))

        elif key == 'progress_max':
            if file_key in self.file_bars:
                self.file_bars[file_key].setMaximum(int(float(value)))

        elif key == "status_update":
            if file_key in self.file_labels:
                item = self.file_labels[file_key]
                item.setText(2, value)
                item.setToolTip(2, value)
        
                status = value.strip()
        
                if status.startswith("❌ Error") or status.startswith("❌ Batch error"):
                    self._error_files.add(file_key)
        
                terminal_statuses = {
                    "✅ Completed",
                    "🚫 Cancelled",
                }
        
                is_terminal_status = (
                    status in terminal_statuses
                    or status.startswith("❌ Error")
                    or status.startswith("❌ Batch error")
                )
        
                if is_terminal_status:
                    if file_key not in self._updated_files:
                        self._updated_files.add(file_key)
                        self.completed_count = len(self._updated_files)
        
                        self.completed_label.setText(
                            f"Completed: "
                            f"{self.completed_count} / {self.total_count}"
                        )
        
                        if self.completed_count == self.total_count:
                            self.show_completion_message(
                                cancelled=self.analysis_cancelled["flag"],
                                had_error=bool(self._error_files)
                            )
                    
class CustomMessageBox(QDialog):
    def __init__(self, message, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Analysis Completed")
        self.setWindowFlags(Qt.Dialog | Qt.WindowTitleHint | Qt.CustomizeWindowHint)
        self.setAttribute(Qt.WA_StyledBackground, True)
        self.setFixedSize(300, 100)
        self.setStyleSheet("""
            QLabel {
                font-size: 14px;
                color: black;
            }
            QPushButton {
                background-color: #1D83D5;
                color: white;
                border: none;
                border-radius: 6px;
                padding: 4px 12px;
                font-weight: bold;
                font-size: 12px; 
            }
            QPushButton:hover {
                background-color: #306999;
            }
        """)

        layout = QVBoxLayout(self)

        label = QLabel(message)
        label.setAlignment(Qt.AlignCenter)
        layout.addWidget(label)

        button = QPushButton("Ok")
        button.clicked.connect(self.accept)
        button.setFixedHeight(28)
        button.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        layout.addWidget(button, alignment=Qt.AlignCenter)
        
    def showEvent(self, event):
        self.raise_()
        self.activateWindow()
        self.setWindowState(self.windowState() & ~Qt.WindowMinimized | Qt.WindowActive)
        QApplication.alert(self, 3000)
        super().showEvent(event)
        
class CustomConfirmBox(QDialog):
    def __init__(self, title, message, on_close=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setFixedSize(350, 100)
        self.on_close = on_close
        self.setWindowFlags(Qt.Dialog | Qt.WindowTitleHint | Qt.CustomizeWindowHint)
        self.setAttribute(Qt.WA_StyledBackground, True)

        self.result = False

        self.setStyleSheet("""
            QLabel {
                font-size: 14px;
                color: #2c3e50;
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

class SamplingPolicyDialog(QDialog):
    def __init__(self, sampling, parent=None):
        super().__init__(parent)

        self.choice = None
        self.sampling = float(sampling)

        self.setWindowTitle("Intermediate Sampling Ratio")
        self.setWindowFlags(
            Qt.Dialog
            | Qt.WindowTitleHint
            | Qt.CustomizeWindowHint
        )
        self.setAttribute(Qt.WA_StyledBackground, True)
        self.setModal(True)
        self.setFixedWidth(570)

        self.setStyleSheet("""
            SamplingPolicyDialog {
                background-color: #F9FAFA;
            }

            QLabel#titleLabel {
                color: #2C3E50;
                font-size: 15px;
                font-weight: bold;
            }

            QLabel#messageLabel {
                color: #2C3E50;
                font-size: 13px;
            }

            QLabel#recommendationLabel {
                color: #555555;
                font-size: 12px;
                background-color: white;
                border: 1px solid #D8DEE3;
                border-radius: 6px;
                padding: 10px;
            }

            QPushButton {
                min-height: 30px;
                padding: 5px 12px;
                border-radius: 6px;
                font-size: 12px;
                font-weight: bold;
            }

            QPushButton#primaryButton {
                background-color: #1D83D5;
                color: white;
                border: none;
            }

            QPushButton#primaryButton:hover {
                background-color: #1565C0;
            }

            QPushButton#secondaryButton {
                background-color: #5F6B73;
                color: white;
                border: none;
            }

            QPushButton#secondaryButton:hover {
                background-color: #465159;
            }

            QPushButton#advancedButton {
                background-color: #F39C12;
                color: white;
                border: none;
            }

            QPushButton#advancedButton:hover {
                background-color: #D68910;
            }

            QPushButton#cancelButton {
                background-color: white;
                color: #555555;
                border: 1px solid #AEB6BF;
            }

            QPushButton#cancelButton:hover {
                background-color: #ECEFF1;
            }
        """)

        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(22, 18, 22, 18)
        main_layout.setSpacing(12)

        title_label = QLabel("⚠️ Intermediate sampling ratio")
        title_label.setObjectName("titleLabel")
        title_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title_label)

        message_label = QLabel(
            f"The selected sampling ratio is {self.sampling:g}%.\n"
            "Sampling ratios between 51% and 99% are not recommended"
        )
        message_label.setObjectName("messageLabel")
        message_label.setAlignment(Qt.AlignCenter)
        message_label.setWordWrap(True)
        main_layout.addWidget(message_label)

        recommendation_label = QLabel(
            "Strict sampling must first determine the exact read count \n"
            "and may be slower than 100% full extraction.\n\n"
            "Recommended settings:\n"
            "• 1–50% for subset screening\n"
            "• 100% for full extraction"
        )
        recommendation_label.setObjectName("recommendationLabel")
        recommendation_label.setWordWrap(True)
        main_layout.addWidget(recommendation_label)

        instruction_label = QLabel("Choose how to proceed.")
        instruction_label.setObjectName("messageLabel")
        instruction_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(instruction_label)

        recommended_layout = QHBoxLayout()
        recommended_layout.setSpacing(8)

        use_100_btn = QPushButton("Use 100%")
        use_100_btn.setObjectName("primaryButton")
        use_100_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        use_100_btn.setDefault(True)
        use_100_btn.clicked.connect(
            lambda: self._select_choice("use_100")
        )

        use_50_btn = QPushButton("Use 50%")
        use_50_btn.setObjectName("secondaryButton")
        use_50_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        use_50_btn.clicked.connect(
            lambda: self._select_choice("use_50")
        )

        recommended_layout.addStretch()
        recommended_layout.addWidget(use_100_btn)
        recommended_layout.addWidget(use_50_btn)
        recommended_layout.addStretch()

        main_layout.addLayout(recommended_layout)

        lower_layout = QHBoxLayout()
        lower_layout.setSpacing(8)

        advanced_btn = QPushButton(
            f"Continue with {self.sampling:g}%"
        )
        advanced_btn.setObjectName("advancedButton")
        advanced_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        advanced_btn.clicked.connect(
            lambda: self._select_choice("continue")
        )

        cancel_btn = QPushButton("Cancel")
        cancel_btn.setObjectName("cancelButton")
        cancel_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        cancel_btn.clicked.connect(self.reject)

        lower_layout.addStretch()
        lower_layout.addWidget(advanced_btn)
        lower_layout.addWidget(cancel_btn)
        lower_layout.addStretch()

        main_layout.addLayout(lower_layout)

        self.adjustSize()

    def _select_choice(self, choice):
        self.choice = choice
        self.accept()

    def showEvent(self, event):
        self.raise_()
        self.activateWindow()
        QApplication.alert(self, 3000)
        super().showEvent(event)    

class CustomAlertBox(QDialog):
    def __init__(
        self,
        title,
        message,
        parent=None,
        alert_type="warning"
    ):
        super().__init__(parent)

        self.setWindowTitle(title)
        self.setWindowFlags(
            Qt.Dialog
            | Qt.WindowTitleHint
            | Qt.CustomizeWindowHint
        )
        self.setAttribute(Qt.WA_StyledBackground, True)
        self.setModal(True)
        self.setMinimumWidth(380)
        self.setMaximumWidth(600)

        alert_settings = {
            "warning": {
                "icon": "⚠️",
                "title_color": "#D98200",
            },
            "error": {
                "icon": "❌",
                "title_color": "#C62828",
            },
            "information": {
                "icon": "ℹ️",
                "title_color": "#1D83D5",
            },
        }

        settings = alert_settings.get(
            alert_type,
            alert_settings["warning"]
        )

        self.setStyleSheet("""
            CustomAlertBox {
                background-color: #F9FAFA;
            }

            QLabel#alertTitle {
                font-size: 15px;
                font-weight: bold;
            }

            QLabel#alertMessage {
                color: #2C3E50;
                font-size: 13px;
                background-color: white;
                border: 1px solid #D8DEE3;
                border-radius: 6px;
                padding: 12px;
            }

            QPushButton {
                background-color: #1D83D5;
                color: white;
                border: none;
                border-radius: 6px;
                padding: 6px 20px;
                font-size: 12px;
                font-weight: bold;
                min-width: 70px;
                min-height: 28px;
            }

            QPushButton:hover {
                background-color: #306999;
            }
        """)

        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(22, 18, 22, 18)
        main_layout.setSpacing(12)

        title_label = QLabel(
            f"{settings['icon']} {title}"
        )
        title_label.setObjectName("alertTitle")
        title_label.setStyleSheet(
            f"color: {settings['title_color']};"
        )
        title_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title_label)

        message_label = QLabel(message)
        message_label.setObjectName("alertMessage")
        message_label.setWordWrap(True)
        message_label.setTextInteractionFlags(
            Qt.TextSelectableByMouse
        )
        main_layout.addWidget(message_label)

        button_layout = QHBoxLayout()
        button_layout.addStretch()

        ok_button = QPushButton("OK")
        ok_button.setCursor(
            QCursor(Qt.CursorShape.PointingHandCursor)
        )
        ok_button.setDefault(True)
        ok_button.clicked.connect(self.accept)

        button_layout.addWidget(ok_button)
        button_layout.addStretch()

        main_layout.addLayout(button_layout)

        self.adjustSize()

    def showEvent(self, event):
        self.raise_()
        self.activateWindow()
        QApplication.alert(self, 3000)
        super().showEvent(event)

class SequenceMiningGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowFlags(
            Qt.Window |
            Qt.WindowMinimizeButtonHint |
            Qt.WindowCloseButtonHint
        )

        self.ui = Ui_In_silico_sequence_mining()
        self.ui.setupUi(self)
        self.setup_matching_method_controls()

        self.file_paths = {}
        self.selected_file_paths = {}
        self.analysis_cancelled = {"flag": False}
        self.update_queue = Queue()

        self.connect_signals()

        self.set_file_format("fastq")
        self.ui.radioButton.setChecked(True)

        try:
            self.ui.threasholdtextEdit.textChanged.connect(lambda *_: self.update_matching_method_controls_state())
        except Exception:
            pass
        
    def _show_alert(
        self,
        title,
        message,
        alert_type="warning"
    ):
        CustomAlertBox(
            title=title,
            message=message,
            parent=self,
            alert_type=alert_type
        ).exec()
        
    def setup_matching_method_controls(self):
        self._compact_exact_window_height = 620
        self._expanded_method_window_height = 760
        try:
            self.resize(1200, self._compact_exact_window_height)
            if hasattr(self.ui, "centralwidget"):
                self.ui.centralwidget.setMinimumHeight(600)
        except Exception:
            pass

        self.matchingMethodFrame = QFrame(self.ui.centralwidget)
        self.matchingMethodFrame.setObjectName("matchingMethodFrame")
        self.matchingMethodFrame.setGeometry(QRect(630, 590, 531, 130))
        self.matchingMethodFrame.setFrameShape(QFrame.Shape.StyledPanel)
        self.matchingMethodFrame.setStyleSheet("QFrame#matchingMethodFrame{background-color:#F8FAFC;border:1px solid #D8DEE9;border-radius:6px;}")

        self.matchingMethodTitle = QLabel("Matching method (Threshold < 100 only)", self.matchingMethodFrame)
        self.matchingMethodTitle.setGeometry(QRect(12, 6, 420, 18))
        self.matchingMethodTitle.setStyleSheet("font-weight:bold;color:#2c3e50;")

        self.exactModeNotice = QLabel(
            "Threshold 100: exact matching is used automatically; method selection is ignored.",
            self.matchingMethodFrame
        )
        self.exactModeNotice.setGeometry(QRect(12, 25, 500, 17))
        self.exactModeNotice.setStyleSheet("font-size:10px;color:#1565C0;")
        self.exactModeNotice.setToolTip(
            "When threshold is 100%, ISSM always uses exact matching.\n"
            "Fast/Legacy method choice only affects threshold values below 100%."
        )

        self.fastBiologicalRadioButton = QRadioButton("Fast biological matching (Recommended, faster)", self.matchingMethodFrame)
        self.fastBiologicalRadioButton.setGeometry(QRect(12, 48, 300, 20))
        self.fastBiologicalRadioButton.setChecked(True)
        self.fastBiologicalRadioButton.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.fastBiologicalRadioButton.setToolTip(
            "Uses probe-length percent identity and allowed mismatches.\n"
            "Threshold means matched bp percentage across the probe length.\n"
            "Example: 20 bp at 95% allows up to 1 mismatch."
        )

        self.legacyFuzzyRadioButton = QRadioButton("Legacy fuzzy matching (Slower, compatibility mode)", self.matchingMethodFrame)
        self.legacyFuzzyRadioButton.setGeometry(QRect(12, 92, 330, 20))
        self.legacyFuzzyRadioButton.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.legacyFuzzyRadioButton.setToolTip(
            "Uses the original RapidFuzz partial-ratio score.\n"
            "Threshold means fuzzy string similarity score, not matched bp percentage.\n"
            "Use this only when reproducing previous ISSM fuzzy results."
        )

        self.fastBiologicalDescription = QLabel(
            "Threshold = probe-length percent identity / allowed mismatches.",
            self.matchingMethodFrame
        )
        self.fastBiologicalDescription.setGeometry(QRect(32, 68, 480, 17))
        self.fastBiologicalDescription.setStyleSheet("font-size:10px;color:#455A64;")

        self.legacyFuzzyDescription = QLabel(
            "Threshold = RapidFuzz partial-ratio score, not exact matched-bp percentage.",
            self.matchingMethodFrame
        )
        self.legacyFuzzyDescription.setGeometry(QRect(32, 112, 485, 17))
        self.legacyFuzzyDescription.setStyleSheet("font-size:10px;color:#455A64;")

        self.matchingMethodFrame.raise_()
        self.matchingMethodFrame.setVisible(False)
        self.update_matching_method_controls_state()

    def update_matching_method_controls_state(self):
        try:
            threshold = float(
                self.ui.threasholdtextEdit.value()
            )
        except Exception:
            threshold = None
    
        show_method_controls = (
            threshold is not None
            and threshold < 100
        )
    
        self.matchingMethodFrame.setVisible(
            show_method_controls
        )
    
        self.adjust_window_height_for_matching_method(
            show_method_controls
        )
    
        if not show_method_controls:
            self.fastBiologicalRadioButton.setChecked(True)
            self.fastBiologicalRadioButton.setEnabled(False)
            self.legacyFuzzyRadioButton.setEnabled(False)
            return
    
        self.fastBiologicalRadioButton.setEnabled(True)
        self.legacyFuzzyRadioButton.setEnabled(True)
    
        self.exactModeNotice.setText(
            "Threshold < 100: choose percent-identity "
            "matching or legacy fuzzy score."
        )

    def adjust_window_height_for_matching_method(self, show_method_controls):
        try:
            target_height = (
                self._expanded_method_window_height
                if show_method_controls
                else self._compact_exact_window_height
            )
            central_min_height = 740 if show_method_controls else 600

            if hasattr(self.ui, "centralwidget"):
                self.ui.centralwidget.setMinimumHeight(central_min_height)
                self.ui.centralwidget.setMaximumHeight(16777215)

            self.setMinimumHeight(target_height)
            self.setMaximumHeight(16777215)

            if not self.isMaximized() and not self.isFullScreen():
                self._force_window_height(target_height)
                QTimer.singleShot(0, lambda h=target_height: self._force_window_height(h))
                QTimer.singleShot(80, lambda h=target_height: self._force_window_height(h))
        except Exception:
            pass

    def _force_window_height(self, target_height):
        try:
            if self.isMaximized() or self.isFullScreen():
                return
            current_width = max(self.width(), 1200)

            self.setMinimumHeight(target_height)
            self.setMaximumHeight(target_height)
            self.resize(current_width, target_height)
            self.setMaximumHeight(16777215)
        except Exception:
            pass

    def get_matching_method(self):
        if getattr(self, "legacyFuzzyRadioButton", None) and self.legacyFuzzyRadioButton.isChecked():
            return "legacy_fuzzy"
        return "fast_biological"

    def connect_signals(self):
        self.ui.radioButton.toggled.connect(lambda: self.set_file_format("fastq"))
        self.ui.radioButton_2.toggled.connect(lambda: self.set_file_format("fasta"))

        self.ui.selectFolderButton.clicked.connect(lambda: browse_folder(
            self.file_paths,
            self.ui.inputPathLineEdit,
            self.ui.allFilesList,
            self.ui.selectedFilesList,
            self.ui.label_4,
            self.ui.label_5,
            self.file_format_var
        ))

        self.ui.rightMovebutton.clicked.connect(lambda: move_selected_files(
            self.ui.allFilesList,
            self.ui.selectedFilesList,
            self.file_paths,
            self.selected_file_paths,
            self.ui.label_4,
            self.ui.label_5
        ))
        
        self.ui.leftMoveBotton.clicked.connect(lambda: move_back_files(
            self.ui.allFilesList,
            self.ui.selectedFilesList,
            self.file_paths,
            self.selected_file_paths,
            self.ui.label_4,
            self.ui.label_5
        ))

        self.ui.selectOutputButton.clicked.connect(lambda: select_output_path(self.ui.outputPathLineEdit))
        
        self.ui.selectAllCheckBox.toggled.connect(
            lambda checked: toggle_all_selection(self.ui.allFilesList, self.ui.selectAllCheckBox)
        )
        
        self.ui.selectAllCheckBox_2.toggled.connect(
            lambda checked: toggle_all_selection(self.ui.selectedFilesList, self.ui.selectAllCheckBox_2)
        )       
        
        self.ui.allFilesList.itemSelectionChanged.connect(
            lambda: update_select_all_checkbox(self.ui.allFilesList, self.ui.selectAllCheckBox)
        )
        self.ui.selectedFilesList.itemSelectionChanged.connect(
            lambda: update_select_all_checkbox(self.ui.selectedFilesList, self.ui.selectAllCheckBox_2)
        )

        self.ui.probeTextEdit.textChanged.connect(
            lambda: update_probe_count(self.ui.probeTextEdit, self.ui.label_6)
        )

        self.ui.startAnalysisButton.clicked.connect(self.start_analysis)

    def set_file_format(self, format_str):
        self.file_format_var = format_str
        self.ui.allFilesList.clear()
        self.ui.selectedFilesList.clear()
        self.file_paths.clear()
        self.selected_file_paths.clear()
        self.ui.inputPathLineEdit.clear()
        
    def _restore_after_analysis(self):
        self.ui.startAnalysisButton.setEnabled(True)
        self.ui.startAnalysisButton.setText("Start Analysis")


    def _resolve_sampling_policy(self, sampling):
        if 0 < sampling <= 50 or sampling == 100:
            return sampling, True
    
        if 50 < sampling < 100:
            dialog = SamplingPolicyDialog(
                sampling=sampling,
                parent=self
            )
            dialog.exec()
    
            if dialog.choice == "use_100":
                self.ui.percentagetextEdit.setValue(100.0)
                return 100.0, True
    
            if dialog.choice == "use_50":
                self.ui.percentagetextEdit.setValue(50.0)
                return 50.0, True
    
            if dialog.choice == "continue":
                return sampling, True
    
            return sampling, False
    
        return sampling, True

    def start_analysis(self):
        probe_text = self.ui.probeTextEdit.toPlainText().strip()
        output_dir = self.ui.outputPathLineEdit.text().strip()
        file_format = "fastq" if self.ui.radioButton.isChecked() else "fasta"
        matching_method = self.get_matching_method()
    
        selected_files = list(self.selected_file_paths.values())
    
        if not selected_files:
            self._show_alert(
                "Input Required",
                "No files have been selected.\n\nPlease select at least one FASTQ or FASTA file."
            )
            return
        
        if not probe_text:
            self._show_alert(
                "Probe Required",
                "Please enter at least one probe sequence before starting the analysis."
            )
            return
    
        is_valid, invalid_bases = validate_probes_before_analysis(
            probe_text
        )
        
        if not is_valid:
            error_details = "\n".join(
                f"• {item}"
                for item in sorted(invalid_bases)
            )
        
            self._show_alert(
                title="Invalid Probe Format",
                message=(
                    "The following probe sequence errors were found:\n\n"
                    f"{error_details}\n\n"
                    "Please correct the probe input before continuing."
                ),
                alert_type="error"
            )
            return

        try:
            threshold = int(round(float(self.ui.threasholdtextEdit.value())))
        except Exception:
            QMessageBox.warning(self, "Warning", "Threshold must be an integer.")
            return

        try:
            sampling = float(
                self.ui.percentagetextEdit.value()
            )
        except Exception:
            self._show_alert(
                title="Invalid Sampling Percentage",
                message=(
                    "The sampling percentage must be a valid number."
                ),
                alert_type="error"
            )
            return
        
        if not (0 < sampling <= 100):
            self._show_alert(
                title="Invalid Sampling Percentage",
                message=(
                    "The sampling percentage must be greater than 0 "
                    "and less than or equal to 100."
                ),
                alert_type="error"
            )
            return

        sampling, sampling_ok = self._resolve_sampling_policy(sampling)
        if not sampling_ok:
            return
    
        if not output_dir:
            self._show_alert(
                title="Output Directory Required",
                message=(
                    "No output directory has been selected.\n\n"
                    "Please select a folder in which to save the results."
                )
            )
            return
        
        if not os.path.isdir(output_dir):
            self._show_alert(
                title="Invalid Output Directory",
                message=(
                    "The selected output directory does not exist:\n\n"
                    f"{output_dir}"
                ),
                alert_type="error"
            )
            return
    
        self.analysis_cancelled["flag"] = False

        self.ui.startAnalysisButton.setEnabled(False)
        self.ui.startAnalysisButton.setText("Processing...")

        self.completed_files = {"count": 0}
        self.total_files = len(selected_files)
    
        self.progress_window = CustomProgressDialog(
            file_list=self.selected_file_paths.keys(),
            parent_gui=self,
            update_queue=self.update_queue,
            analysis_cancelled=self.analysis_cancelled
        )
        for file_path in selected_files:
            self.update_queue.put(("status_update", file_path, "Stand-by"))
        self.progress_window.setModal(False)
        self.progress_window.show()
        self.progress_window.finished.connect(
            lambda _=None: self._restore_after_analysis()
        )
        self.progress_window.destroyed.connect(
            lambda _=None: self._restore_after_analysis()
        )

        threading.Thread(
            target=run_analysis,   
            args=(
                probe_text,           
                threshold,           
                sampling,            
                output_dir,           
                selected_files,       
                file_format,          
                self.analysis_cancelled,
                self.update_queue,
                matching_method
            ),
            daemon=True
        ).start()
            
def file_3_GUI():
    app = QApplication.instance()
    created = False

    if not app:
        app = QApplication(sys.argv)
        created = True

    gui = SequenceMiningGUI()
    gui.show()

    if created:
        def keep_alive():
            pass

        timer = QTimer()
        timer.timeout.connect(keep_alive)
        timer.start(1000)

        sys.exit(app.exec())