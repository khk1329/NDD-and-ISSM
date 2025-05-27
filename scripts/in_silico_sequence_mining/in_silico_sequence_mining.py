import os
import gzip
import random
import threading
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import reverse_complement
from collections import Counter
import matplotlib
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
from concurrent.futures import ProcessPoolExecutor, as_completed

from queue import Queue
from ui.ISSM import Ui_In_silico_sequence_mining
from ui import resource_rc

matplotlib.use('Agg')
logging.getLogger('matplotlib').setLevel(logging.WARNING)
analysis_cancelled = {"flag": False}

lock = threading.Lock()

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
    "N": ["A", "C", "G", "T"]
}

invalid_base_warning_shown = False

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

    if selected_format == "fastq":
        valid_extensions = (".fastq.gz", ".fastq", ".fq.gz", ".fq")
    elif selected_format == "fasta":
        valid_extensions = (".fasta.gz", ".fasta", ".fa.gz", ".fa", ".fna")
    else:
        valid_extensions = (".fastq.gz", ".fastq", ".fq.gz", ".fq", ".fasta.gz", ".fasta", ".fa.gz", ".fa", ".fna")

    all_files_list_box.clear()
    file_paths.clear()

    for root, _, files in os.walk(folder_path):
        for file in files:
            if file.endswith(valid_extensions):  
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

def fuzzy_match(cached_sequences, target_name, read_sequence, threshold):
    possible_sequences = cached_sequences.get(target_name, [])
    for seq in possible_sequences:
        partial_ratio = fuzz.partial_ratio(seq, read_sequence)
        if partial_ratio >= threshold:
            return True, partial_ratio
    return False, 0

def init_worker(_cached_sequences):
    global shared_cached_sequences
    shared_cached_sequences = _cached_sequences

def init_worker_all(shared_sequences, safe_cores=None):
    global shared_cached_sequences
    shared_cached_sequences = shared_sequences

    if safe_cores is not None and os.name == 'nt':
        try:
            import psutil
            p = psutil.Process()
            p.cpu_affinity(safe_cores)
        except Exception as e:
            print(f"⚠️ Failed to set CPU affinity: {e}")

def match_batch_sequences(batch, probes, threshold, batch_idx, cached_sequences):
    matched = []

    for record in batch: 
        if analysis_cancelled["flag"]:
            break
        read_id = record.id
        read_seq = str(record.seq)

        for probe_name, probe_seq in probes:
            probe_seq = probe_seq.upper()
            
            if len(read_seq) < len(probe_seq):
                continue

            match_orig, score_orig = fuzzy_match(cached_sequences, probe_name, read_seq, threshold)
            match_rc, score_rc = fuzzy_match(cached_sequences, probe_name, reverse_complement(read_seq), threshold)

            if match_orig or match_rc:
                matched.append((read_id, probe_name, score_orig, score_rc, read_seq))

    batch_size = len(batch)
    return matched, batch_idx, batch_size

def match_batch_reads(records, probes, threshold, batch_index):
    matched = []
    for read_id, read_seq in records:
        if analysis_cancelled["flag"]:
            break
        for probe_name, probe_seq in probes:
            probe_seq = probe_seq.upper()

            if len(read_seq) < len(probe_seq):
                continue

            match_orig, score_orig = fuzzy_match(shared_cached_sequences, probe_name, read_seq, threshold)
            match_rc, score_rc = fuzzy_match(shared_cached_sequences, probe_name, reverse_complement(read_seq), threshold)

            if match_orig or match_rc:
                matched.append((read_id, probe_name, score_orig, score_rc, read_seq))
                
    batch_size = len(records)
    
    return matched, batch_index, batch_size

def run_analysis(
    probes_text_widget, threshold_entry, percentage_entry, output_file_entry,
    selected_file_list_box, selected_file_paths, file_format_var,
    analysis_cancelled, update_queue,
    use_custom_parallel, custom_parallel_value 
):
    logging.debug("🔬 Analysis started")

    if hasattr(file_format_var, 'text'):
        selected_format = file_format_var.text().lower()
    elif hasattr(file_format_var, 'currentText'):
        selected_format = file_format_var.currentText().lower()
    else:
        selected_format = str(file_format_var).lower()

    probes_text = probes_text_widget.toPlainText().strip()
    if not probes_text:
        QMessageBox.warning(None, "Warning", "Please enter probe sequences")
        return

    try:
        probes = [
            (lines.split('\n')[0], ''.join(lines.split('\n')[1:]))
            for lines in probes_text.strip().split('>') if lines.strip()
        ]
        cached_sequences = prepare_probe_sequences(probes)
    except ValueError as e:
        QMessageBox.critical(None, "Error", str(e))
        return

    try:
        threshold = int(threshold_entry.text())
        sampling_percentage = float(percentage_entry.text())
    except ValueError:
        QMessageBox.critical(None, "Error", "Threshold / Sampling value is invalid")
        return

    output_dir_path = output_file_entry.text().strip()
    if not output_dir_path:
        QMessageBox.warning(None, "Warning", "Please select an output directory")
        return

    selected_files = [
        selected_file_paths.get(selected_file_list_box.item(i).text())
        for i in range(selected_file_list_box.count())
    ]
    selected_files = [f for f in selected_files if f]

    if not selected_files:
        QMessageBox.warning(None, "Warning", "Please select files to analyze")
        return
    
    completed_files = {'count': 0}
    total_files = len(selected_files)
    timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_file_path = os.path.join(output_dir_path, f"summary_result_{timestamp_str}.txt")
    results = {"file_results": {}, "probe_results": {}}
    
    threading.Thread(
        target=process_files,
        args=(
            selected_files, probes, threshold, sampling_percentage,
            output_dir_path, selected_format, cached_sequences,
            completed_files, total_files, summary_file_path, results,
            analysis_cancelled, update_queue,
            use_custom_parallel, custom_parallel_value
        ),
        daemon=True
    ).start()

    logging.debug("🔄 Background analysis thread started.")
    
def process_files(
    selected_files, probes, threshold, sampling_percentage,
    output_dir_path, selected_format, cached_sequences,
    completed_files, total_files, summary_file_path, results,
    analysis_cancelled, update_queue,
    use_custom_parallel, custom_parallel_value
):

    logging.debug("🚀 ThreadPoolExecutor Start")
    
    if use_custom_parallel:
        max_parallel = int(custom_parallel_value)
        logging.info(f"User-defined parallel count: {max_parallel}")
    else:
        max_parallel = get_optimal_parallel_file_count()
        logging.info(f"Auto-detected parallel count: {max_parallel}")

    def run_file(input_file, output_file_specific):
        try:
            file_key = os.path.basename(input_file)

            if analysis_cancelled["flag"]:
                update_queue.put(("status_update", file_key, "🚫 Cancelled"))
                return
        
            detect_func = (
                FASTQ_detect_and_quantify_probe
                if selected_format == "fastq"
                else FASTA_detect_and_quantify_probe_batched
            )

            detect_func(
                input_file, probes, threshold, sampling_percentage, output_file_specific,
                completed_files, total_files, summary_file_path, cached_sequences,
                results, output_dir_path, analysis_cancelled, update_queue
            )

        except Exception as e:
            logging.error(f"❌ {input_file} Error during processing: {e}")
            update_queue.put(("status_update", file_key, f"❌ Error: {e}"))

    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        futures = {}
        for input_file in selected_files:
            
            file_key = os.path.basename(input_file)

            if not os.path.exists(input_file):
                continue
    
            output_file_specific = os.path.join(
                output_dir_path,
                f"{os.path.splitext(file_key)[0]}_results.txt"
            )
    
            update_queue.put(("progress_max", file_key, 100))
            update_queue.put(("progress_update", file_key, 0))
            update_queue.put(("status_update", file_key, "Stand-by"))
    
            future = executor.submit(run_file, input_file, output_file_specific)
            futures[future] = file_key
    
        try:
            for future in as_completed(futures):
                if analysis_cancelled["flag"]:
                    for f in futures:
                        if not f.running() and not f.done():
                            f.cancel()
                            file_key = futures[f]
                            update_queue.put(("status_update", file_key, "🚫 Cancelled"))
                    break
    
                future.result()
    
        except Exception as e:
            logging.error(f"❌ Error occurred during analysis: {e}")
            
    update_queue.put(("analysis_complete", "", "All analyses have been completed!"))

def prepare_probe_sequences(probes):
    cached_sequences = {}

    for name, sequence in probes:
        logging.debug(f"Converting target: {name}")
        cached_sequences[name] = generate_possible_sequences(sequence.upper())

    return cached_sequences

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

def FASTA_sampling(input_file, sampling_percentage, update_queue=None, min_reads=1, analysis_cancelled=None):
    try:
        if analysis_cancelled and analysis_cancelled.get("flag"):
            raise InterruptedError("Sampling was cancelled before starting.")

        if update_queue:
            update_queue.put(("status_update", input_file, "Ready to sampling"))

        with open(input_file, "r") as handle:
            fasta_records = []
            for record in SeqIO.parse(handle, "fasta"):
                if analysis_cancelled and analysis_cancelled["flag"]:
                    return []
                fasta_records.append(record)

        total_reads = len(fasta_records)

        if total_reads == 0:
            logging.warning(f"⚠ No sequences found in {input_file}")
            if update_queue:
                update_queue.put(("status_update", input_file, "⚠ No sequences found"))
            return []

        sample_size = max(int(total_reads * (sampling_percentage / 100)), min_reads)
        sample_size = min(sample_size, total_reads)

        if analysis_cancelled and analysis_cancelled.get("flag"):
            raise InterruptedError("Sampling cancelled before sampling step.")

        sampled_records = random.sample(fasta_records, sample_size)
        logging.info(f"✅ Sampled {sample_size} sequences from {total_reads} total sequences in {input_file}")

        if update_queue:
            update_queue.put(("progress_max", input_file, sample_size))
            update_queue.put(("progress_update", input_file, sample_size))
            update_queue.put(("status_update", input_file, f"✅ Sampled {sample_size:,} sequences"))

        return sampled_records

    except InterruptedError as ie:
        logging.warning(f"🛑 Sampling cancelled for {input_file}: {ie}")
        if update_queue:
            update_queue.put(("status_update", input_file, "🚫 Cancelled"))
        return []

    except Exception as e:
        logging.error(f"❌ Error during FASTA sampling for {input_file}: {e}")
        if update_queue:
            if analysis_cancelled and analysis_cancelled.get("flag"):
                update_queue.put(("status_update", input_file, "🚫 Cancelled"))
            else:
                update_queue.put(("status_update", input_file, f"❌ Error: {e}"))
        return []

def FASTA_detect_and_quantify_probe_batched(
    input_file, probes, threshold, sampling_percentage, output_file_path,
    completed_files, total_files, summary_file_path,
    cached_sequences, results, output_dir_path,
    analysis_cancelled, update_queue
):
    def fasta_batch_generator(records, batch_size):
        batch = []
        for record in records:
            batch.append(record)
            if len(batch) == batch_size:
                yield batch
                batch = []
        if batch:
            yield batch
            
    processed_so_far = 0
    total_records = 0
    cancelled_early = False
    
    try:
        start_time = datetime.now()
        original_filename = input_file

        fasta_records = FASTA_sampling(
            input_file,
            sampling_percentage,
            update_queue=update_queue,
            min_reads=1,
            analysis_cancelled=analysis_cancelled
        )

        if not fasta_records:
            raise ValueError(f"❌ No valid sequences found in {input_file}")

        probe_names = np.array([name for name, _ in probes])
        probe_dict = {name: seq for name, seq in probes}

        if analysis_cancelled["flag"]:
            cancelled_early = True
            return

        total_records = len(fasta_records)
        BATCH_SIZE = 100

        update_queue.put(("progress_max", input_file, total_records))
        update_queue.put(("progress_update", input_file, 0))
        update_queue.put(("status_update", input_file, f"Matching sequences (0/{total_records:,})"))

        probe_counts = Counter({name: 0 for name in probe_names})
        matched_reads = set()
        read_to_probes = {}
        total_matched_reads = 0

        cpu_count = os.cpu_count() or 1
        num_workers = min(int(cpu_count * 0.75), 61)
        safe_cores = list(range(num_workers))
        
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
            result_file.write("Read ID\tProbe Name\tScore Original\tScore Reverse Complement\tRead Sequence\n")

        with open(output_file_path, 'a', encoding="utf-8") as result_file, \
             ProcessPoolExecutor(max_workers = num_workers, initializer=init_worker_all, initargs=(cached_sequences, safe_cores)) as executor:

            futures = []
            for batch_index, record_batch in enumerate(fasta_batch_generator(fasta_records, BATCH_SIZE), start=1):
                if analysis_cancelled["flag"]:
                    cancelled_early = True
                    break
                futures.append(executor.submit(match_batch_sequences, record_batch, probes, threshold, batch_index, cached_sequences))

            completed_batches = set()
            
            for future in as_completed(futures):
                if analysis_cancelled["flag"]:
                    cancelled_early = True
                    break
            
                try:
                    matches, batch_idx, batch_size = future.result()
            
                    if batch_idx not in completed_batches:
                        processed_so_far += batch_size
                        completed_batches.add(batch_idx)
            
                    for match in matches:
                        read_id, probe_name, score_orig, score_rc, read_seq = match
                        probe_counts[probe_name] += 1
                        matched_reads.add(read_id)
                        read_to_probes.setdefault(read_id, set()).add(probe_name)
                        result_file.write(f"{read_id}\t{probe_name}\t{score_orig}\t{score_rc}\t{read_seq}\n")
            
                    update_queue.put(("progress_update", input_file, processed_so_far))
                    update_queue.put(("status_update", input_file, f"Matching sequences ({processed_so_far:,}/{total_records:,})"))
            
                except Exception as e:
                    logging.exception(f"❌ Exception during FASTA batch processing: {e}")
                    update_queue.put(("status_update", input_file, f"❌ Batch error: {e}"))
                    
            total_matched_reads = len(matched_reads)

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
        logging.exception(f"❌ Error processing FASTA file: {e}")
        update_queue.put(("status_update", input_file, f"❌ Error: {e}"))

    finally:
        with lock:
            completed_files['count'] += 1

        update_queue.put(("progress_update", input_file, processed_so_far))
        update_queue.put(("status_update", input_file, f"Matching sequences ({processed_so_far:,}/{total_records:,})"))
        time.sleep(0.4)
        
        try_generate_plots_once(completed_files, lock, total_files, analysis_cancelled, results, output_dir_path, probes)

        if cancelled_early or (analysis_cancelled["flag"] and processed_so_far < total_records):
            update_queue.put(("status_update", input_file, "🚫 Cancelled"))
        else:
            update_queue.put(("status_update", input_file, "✅ Completed"))

def FASTQ_sampling(input_file, sampling_percentage, file_format, update_queue=None, min_reads=5000, analysis_cancelled=None):
    temp_output_file = None

    try:
        base_filename = os.path.basename(input_file).replace(".gz", "").replace(".fastq", "").replace(".fasta", "").replace(".fa", "")
        unique_id = uuid.uuid4().hex[:8]
        temp_output_file = os.path.join(tempfile.gettempdir(), f"tmp_{base_filename}_{unique_id}.fastq.gz")

        def open_text_stream(file_path):
            return io.TextIOWrapper(gzip.open(file_path, "rb"), encoding="utf-8") if file_path.endswith(".gz") else open(file_path, "rt", encoding="utf-8")

        if analysis_cancelled and analysis_cancelled.get("flag"):
            raise InterruptedError("Sampling cancelled before starting")

        if input_file.endswith(".gz") and not validate_gzip_file(input_file):
            raise ValueError(f"Invalid gzip file: {input_file}")

        if file_format == "fastq":
            validate_fastq_format(input_file)

        if update_queue:
            update_queue.put(("status_update", input_file, "Ready to sampling"))

        total_reads = 0
        with open_text_stream(input_file) as handle:
            for record in SeqIO.parse(handle, file_format):
                if analysis_cancelled and analysis_cancelled.get("flag"):
                    raise InterruptedError("Sampling cancelled during counting.")
                total_reads += 1

        sample_size = max(int(total_reads * (sampling_percentage / 100)), min_reads)
        sample_size = min(sample_size, total_reads)

        if update_queue:
            update_queue.put(("status_update", input_file, f"Sampling {sample_size:,} reads..."))
            update_queue.put(("progress_max", input_file, sample_size))
            update_queue.put(("progress_update", input_file, 0))

        selected_indices = set(sorted(random.sample(range(total_reads), sample_size)))
        sampled_count = 0

        with open_text_stream(input_file) as handle, gzip.open(temp_output_file, "wt") as output_handle:
            for i, record in enumerate(SeqIO.parse(handle, file_format)):
                if analysis_cancelled and analysis_cancelled.get("flag"):
                    raise InterruptedError("Sampling cancelled during writing.")
                if i in selected_indices:
                    SeqIO.write(record, output_handle, file_format)
                    sampled_count += 1
                    if update_queue and sampled_count % 10000 == 0:
                        update_queue.put(("progress_update", input_file, sampled_count))

        if update_queue:
            update_queue.put(("progress_update", input_file, sample_size))
            update_queue.put(("status_update", input_file, "✅ Sampling completed"))

        logging.info(f"✅ Sampling completed: {sample_size} reads → {temp_output_file}")
        return temp_output_file, total_reads

    except InterruptedError:
        if temp_output_file and os.path.exists(temp_output_file):
            try:
                os.remove(temp_output_file)
                logging.warning(f"🛑 Sampling cancelled. Temp file removed: {temp_output_file}")
            except Exception as e:
                logging.warning(f"⚠️ Failed to remove temp file during cancel: {e}")
        if update_queue:
            update_queue.put(("status_update", input_file, "🚫 Cancelled"))
        return None, None

    except Exception as e:
        if temp_output_file and os.path.exists(temp_output_file):
            try:
                os.remove(temp_output_file)
                logging.info(f"🧹 Removed temp file after sampling error: {temp_output_file}")
            except Exception as err:
                logging.warning(f"⚠️ Failed to remove temp file after error: {err}")

        logging.error(f"❌ Error during sampling: {e}")
        if update_queue:
            if analysis_cancelled and analysis_cancelled.get("flag"):
                update_queue.put(("status_update", input_file, "🚫 Cancelled"))
            else:
                update_queue.put(("status_update", input_file, f"❌ Error: {e}"))
        return None, None

def FASTQ_detect_and_quantify_probe(
    input_file, probes, threshold, sampling_percentage, output_file_path,
    completed_files, total_files,
    summary_file_path, cached_sequences, results, output_dir_path, analysis_cancelled, update_queue
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

    try:
        BATCH_SIZE = 500
        start_time = datetime.now()
        original_filename = input_file
        file_format = "fastq"

        file_to_use = input_file
        if sampling_percentage != 100:
            sampled_file, _ = FASTQ_sampling(
                input_file, sampling_percentage, file_format,
                update_queue, 5000, analysis_cancelled
            )
            
            if analysis_cancelled["flag"] or not sampled_file:
                update_queue.put(("status_update", input_file, "🚫 Cancelled"))
                cancelled_early = True
                return
            file_to_use = sampled_file


        update_queue.put(("status_update", input_file, "🔍 Counting reads..."))

        with gzip.open(file_to_use, "rt") if file_to_use.endswith(".gz") else open(file_to_use, "rt") as handle:
            total_records = sum(1 for _ in SeqIO.parse(handle, "fastq"))

        update_queue.put(("progress_max", input_file, total_records))
        update_queue.put(("progress_update", input_file, 0))
        update_queue.put(("status_update", input_file, f"Matching reads (0/{total_records:,})"))

        probe_counts = Counter({name: 0 for name, _ in probes})
        matched_reads = set()
        read_to_probes = {}
        probe_dict = {name: seq for name, seq in probes}
        total_matched_reads = 0
        cpu_count = os.cpu_count() or 1
        num_workers = min(int(cpu_count * 0.75), 61)
        futures_queue = queue.Queue()

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
            result_file.write("Read ID\tProbe Name\tScore Original\tScore Reverse Complement\tRead Sequence\n")

        with gzip.open(file_to_use, "rt") if file_to_use.endswith(".gz") else open(file_to_use, "rt") as handle, \
             open(output_file_path, 'a', encoding="utf-8") as result_file, \
             ProcessPoolExecutor(max_workers=num_workers, initializer=init_worker, initargs=(cached_sequences,)) as executor:

            def submit_batches():
                for batch_index, record_batch in enumerate(batch_generator(handle, BATCH_SIZE), 1):
                    if analysis_cancelled["flag"]:
                        logging.info(f"🛑 Skipping batch submit due to cancellation: {input_file}")
                        break
                    future = executor.submit(match_batch_reads, record_batch, probes, threshold, batch_index)
                    futures_queue.put(future)

            submit_thread = threading.Thread(target=submit_batches)
            submit_thread.start()

            last_progress_update = 0
            completed_batches = set()
            
            while submit_thread.is_alive() or not futures_queue.empty():
                if analysis_cancelled["flag"]:
                    cancelled_early = True
                    logging.info(f"🛑 Cancelled during processing loop: {input_file}")
                    break                
                try:
                    future = futures_queue.get(timeout=0.2)
                    matches, batch_idx, batch_size = future.result()
            
                    if batch_idx not in completed_batches:
                        processed_so_far += batch_size
                        completed_batches.add(batch_idx)
            
                    for match in matches:
                        read_id, probe_name, score_orig, score_rc, read_seq = match
                        probe_counts[probe_name] += 1
                        matched_reads.add(read_id)
                        read_to_probes.setdefault(read_id, set()).add(probe_name)
                        result_file.write(f"{read_id}\t{probe_name}\t{score_orig}\t{score_rc}\t{read_seq}\n")
            
                    if processed_so_far - last_progress_update >= 10000:
                        update_queue.put(("progress_update", input_file, processed_so_far))
                        update_queue.put(("status_update", input_file, f"Matching reads ({processed_so_far:,}/{total_records:,})"))

                        last_progress_update = processed_so_far
            
                except queue.Empty:
                    continue
                except Exception as e:
                    logging.exception(f"❌ Exception during FASTQ batch processing: {e}")
                    update_queue.put(("status_update", input_file, f"❌ Batch error: {e}"))

            submit_thread.join()
            
            total_matched_reads = len(matched_reads)
            
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
        logging.exception(f"❌ Error processing file {input_file}: {e}")
        update_queue.put(("status_update", input_file, f"❌ Error: {e}"))

    finally:
        if sampled_file and sampling_percentage < 100:
            try:
                os.remove(sampled_file)
                logging.warning(f"Successfully remove temp file: {sampled_file}")
            except Exception as e:
                logging.warning(f"🧹 Failed to remove temp file: {sampled_file}, reason: {e}")

        completed_files['count'] += 1
    
        processed_so_far = min(processed_so_far, total_records)
        update_queue.put(("progress_update", input_file, processed_so_far))
        update_queue.put(("status_update", input_file, f"Matching reads ({processed_so_far:,}/{total_records:,})"))
        time.sleep(0.4)
        
        try_generate_plots_once(completed_files, lock, total_files, analysis_cancelled, results, output_dir_path, probes)
    
        if cancelled_early or (analysis_cancelled["flag"] and processed_so_far < total_records):
            update_queue.put(("status_update", input_file, "🚫 Cancelled"))
        else:
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
    import matplotlib.pyplot as plt
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
            cmap="coolwarm",
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
        if (
            completed_files['count'] == total_files and
            not completed_files.get("plotted", False)
        ):
            completed_files["plotted"] = True

            if analysis_cancelled["flag"]:
                logging.warning("⚠️ Plot generation skipped due to user cancellation.")
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
        self.setMinimumSize(800, 600)
        self.resize(800, 600)
        self.setWindowFlags(Qt.Window | Qt.WindowMinimizeButtonHint)

        self.file_bars = {}
        self.file_labels = {}
        self.parent_gui = parent_gui

        self.completed_count = 0
        self.total_count = len(file_list)
        self._updated_files = set()
        self._completion_shown = False

        self.completed_label = QLabel("Completed: 0 / 0")
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
        self.tree.header().resizeSection(0, 250)
        self.tree.header().setSectionResizeMode(1, QHeaderView.Stretch)
        self.tree.header().setSectionResizeMode(2, QHeaderView.Fixed)
        self.tree.header().resizeSection(2, 250)

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
        
    def show_completion_message(self, cancelled=False):
        if self._completion_shown:
            return 
        self._completion_shown = True

        if cancelled:
            CustomMessageBox("🚫 Analysis was cancelled by user.", self).exec_()
        else:
            CustomMessageBox("✅ All tasks are completed!", self).exec_()
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

        elif key == 'status_update':
            if file_key in self.file_labels:
                item = self.file_labels[file_key]
                item.setText(2, value)
        
                if value.strip() in ("✅ Completed", "🚫 Cancelled"):
                    if file_key not in self._updated_files:
                        self._updated_files.add(file_key)
                        self.completed_count = len(self._updated_files)
                        self.completed_label.setText(f"Completed: {self.completed_count} / {self.total_count}")
        
                        if self.completed_count == self.total_count:
                            self.show_completion_message(cancelled=self.analysis_cancelled["flag"])
    
class CustomMessageBox(QDialog):
    def __init__(self, message, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Analysis Completed")
        self.setWindowFlags(Qt.Dialog | Qt.WindowTitleHint | Qt.CustomizeWindowHint)
        self.setAttribute(Qt.WA_StyledBackground, True)
        self.setFixedSize(300, 100)
        self.setStyleSheet("""
            QDialog {
                background-color: #ffffff;
            }
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
            QDialog {
                background-color: #ffffff;
            }
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

        self.file_paths = {}
        self.selected_file_paths = {}
        self.analysis_cancelled = {"flag": False}
        self.update_queue = Queue()

        self.connect_signals()

        self.set_file_format("fastq")
        self.ui.radioButton.setChecked(True)

        self.ui.parallelProcessing.setVisible(False)
        self.ui.label_10.setVisible(False)

        self.ui.parallelProcessingCheckbox.stateChanged.connect(self.toggle_parallel_controls)

    def toggle_parallel_controls(self, state):
        is_checked = bool(state)
        self.ui.parallelProcessing.setVisible(is_checked)
        self.ui.label_10.setVisible(is_checked)

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

    def start_analysis(self):
        if not self.selected_file_paths:
            QMessageBox.warning(self, "Warning", "No files selected!")
            return
    
        probe_text = self.ui.probeTextEdit.toPlainText()
        is_valid, invalid_bases = validate_probes_before_analysis(probe_text)
    
        if not is_valid:
            QMessageBox.warning(
                self,
                "Invalid Base or format Found",
                f"The following invalid bases or formatting errors were found in probe sequences:\n"
                f"{', '.join(invalid_bases)}\n\n"
                "Please correct them before continuing. Analysis will be aborted."
            )
            return

        file_format = "fastq" if self.ui.radioButton.isChecked() else "fasta"
        self.analysis_cancelled["flag"] = False
        
        self.completed_files = {"count": 0}
        self.total_files = len(self.selected_file_paths)

        self.progress_window = CustomProgressDialog(
            file_list=self.selected_file_paths.keys(),
            parent_gui=self,
            update_queue=self.update_queue,
            analysis_cancelled=self.analysis_cancelled
        )
        
        for file_path in self.selected_file_paths.values():
            self.update_queue.put(("status_update", file_path, "Stand-by"))

        self.progress_window.setModal(False)
        self.progress_window.show()

        threading.Thread(
            target=run_analysis,
            args=(
                self.ui.probeTextEdit,
                self.ui.threasholdtextEdit,
                self.ui.percentagetextEdit,
                self.ui.outputPathLineEdit,
                self.ui.selectedFilesList,
                self.selected_file_paths,
                file_format,
                self.analysis_cancelled,
                self.update_queue,
                self.ui.parallelProcessingCheckbox.isChecked(),
                self.ui.parallelProcessing.value()  
            ),
            daemon=True
        ).start()
            
def file_3_GUI():
    app = QApplication(sys.argv)
    gui = SequenceMiningGUI()
    gui.show()

    def keep_alive():
        pass
    timer = QTimer()
    timer.timeout.connect(keep_alive)
    timer.start(1000)

    sys.exit(app.exec())
