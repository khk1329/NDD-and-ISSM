import os
import requests
import tarfile
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, Label, Toplevel
from Bio import Entrez
import subprocess
import threading
from datetime import datetime
import time
import csv
from bs4 import BeautifulSoup
import zipfile
import json
import shutil
import sys
import tempfile 
import gzip
import traceback
from concurrent.futures import ThreadPoolExecutor
import re
from urllib.parse import quote
from ftplib import FTP
from PIL import Image, ImageTk
import winreg
import ctypes
import psutil
import logging

def get_resource_path(relative_path):
    if getattr(sys, 'frozen', False):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)

accession_tree = None
download_list_tree = None
stop_flag = False

search_state = {
    'query': '',
    'retmax': 50,
    'retstart': 0,
    'total': 0
}

prefetch_retry_count = 3

log_file_name = None
ena_search_results = []
processes = []  
completed_items = []
failed_files = []
active_threads = []
log_lock = threading.Lock()  

global start_button  

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

def set_sratoolkit_env(toolkit_path):
    try:
        os.environ["PATH"] += os.pathsep + toolkit_path

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
                print("✅ Updating User PATH completed")
            else:
                print("⚠️ The path is already registered in the user environment variables.")

        subprocess.run("taskkill /f /im explorer.exe", shell=True)
        subprocess.run("start explorer.exe", shell=True)

        print(f"✅ Environment variable successfully applied. (User + System): {toolkit_path}")

    except Exception as e:
        print(f"❌  Failed to set environment variable : {e}")

def check_and_install_sratoolkit(top):
    def check_sratoolkit_installed():
        try:
            subprocess.run(["prefetch", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            print("✅ SRA Toolkit successfully detected.")
            return True
        except FileNotFoundError:
            print("❌ 'prefetch' command not found! The environment variable may not be applied.")
            return False
        except subprocess.CalledProcessError as e:
            print(f"❌  An error occurred while running the SRA Toolkit: {e}")
            return False

    def install_sratoolkit():
        try:
            install_progress_label.config(text="Downloading SRA Toolkit...")
            top.update_idletasks()

            download_url = "https://ftp.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-win64.zip"
            download_path = os.path.join(os.getcwd(), "sratoolkit.zip")

            with requests.get(download_url, stream=True) as response:
                response.raise_for_status()
                with open(download_path, "wb") as file:
                    for chunk in response.iter_content(chunk_size=8192):
                        file.write(chunk)

            install_progress_label.config(text="Extracting SRA Toolkit...")
            top.update_idletasks()

            with zipfile.ZipFile(download_path, "r") as zip_ref:
                zip_ref.extractall(os.getcwd())

            toolkit_dir = next((d for d in os.listdir() if d.startswith("sratoolkit")), None)
            if toolkit_dir:
                toolkit_path = os.path.abspath(os.path.join(toolkit_dir, "bin"))

                set_sratoolkit_env(toolkit_path)

                messagebox.showinfo("Success", f"SRA Toolkit installation completed and environment variable applied.:\n{toolkit_path}")

                def close_install_window():
                    completion_window.destroy()
                    install_window.destroy()
                    top.destroy()
                    SRA_to_FASTQ_downloader_GUI()

                completion_window = tk.Toplevel(install_window)
                completion_window.title("Installation Complete")
                completion_window.geometry("350x100")

                tk.Label(completion_window, text="SRA Toolkit installation completed successfully.").pack(pady=15)
                tk.Button(completion_window, text="Ok", command=close_install_window, cursor="hand2").pack(pady=5)

                completion_window.update_idletasks()
                width = completion_window.winfo_width()
                height = completion_window.winfo_height()
                x = (completion_window.winfo_screenwidth() // 2) - (width // 2)
                y = (completion_window.winfo_screenheight() // 2) - (height // 2)
                completion_window.geometry(f'{width}x{height}+{x}+{y}')
            
            else:
                raise FileNotFoundError("Extracted toolkit directory not found.")

        except requests.RequestException as e:
            messagebox.showerror("Error", f"Download failed: {e}")
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error has occurred: {e}")

    if not check_sratoolkit_installed():
        install_window = tk.Toplevel(top)
        install_window.title("SRA Toolkit Installation")
        install_window.geometry("400x150")

        install_progress_label = tk.Label(install_window, text="Checking SRA Toolkit...")
        install_progress_label.pack(pady=20)

        install_button = tk.Button(
            install_window,
            text="Install SRA Toolkit", cursor="hand2",
            command=lambda: install_sratoolkit()
        )
        install_button.pack(pady=10)
    else:
        def close_info_window():
            info_window.destroy()
            top.destroy()
            SRA_to_FASTQ_downloader_GUI()
                    
        info_window = tk.Toplevel(top)
        info_window.title("Info")
        info_window.geometry("350x100")

        tk.Label(info_window, text="✅ SRA Toolkit is already installed and recognized.").pack(pady=15)
        tk.Button(info_window, text="Ok", command=close_info_window, cursor="hand2").pack(pady=5)

        info_window.update_idletasks()
        width = info_window.winfo_width()
        height = info_window.winfo_height()
        x = (info_window.winfo_screenwidth() // 2) - (width // 2)
        y = (info_window.winfo_screenheight() // 2) - (height // 2)
        info_window.geometry(f'{width}x{height}+{x}+{y}')

def format_number(value):
    try:
        return "{:,}".format(int(value))
    except ValueError:
        return "N/A"

def update_status(status_labels, index, message):
    status_labels[index].config(text=message)

def notify_completion(top, message, status_window):
    global start_button, log_file_name 

    completion_window = tk.Toplevel(top)
    completion_window.title("Download Completed")
    completion_window.geometry("350x100")

    tk.Label(completion_window, text=message).pack(pady=15)

    def close_windows():
        global log_file_name

        completion_window.destroy()
        status_window.destroy()

        log_file_name = None
        print("✅ log_file_name reset to None after completion.")

        if 'start_button' in globals() and start_button and start_button.winfo_exists():
            print("✅ Re-enabling Start button after completion.")
            start_button.config(state="normal")
        else:
            print("⚠️ start_button is not available! Trying alternative method.")

            for widget in top.winfo_children():
                if isinstance(widget, tk.Button) and widget.cget("text") == "Download Start":
                    print("✅ Found start_button! Re-enabling now.")
                    widget.config(state="normal")
                    break

    tk.Button(completion_window, text="OK", command=close_windows, cursor="hand2").pack(pady=5)

    completion_window.attributes("-topmost", True)

    completion_window.update_idletasks()
    width = completion_window.winfo_width()
    height = completion_window.winfo_height()
    x = (completion_window.winfo_screenwidth() // 2) - (width // 2)
    y = (completion_window.winfo_screenheight() // 2) - (height // 2)
    completion_window.geometry(f'{width}x{height}+{x}+{y}')

def open_download_status_window(top, download_items, output_dir, download_list_tree, selected_format, selected_db, start_button):
    global stop_flag, status_window
    
    if not download_items:
        messagebox.showerror("Error", "No items to download!")
        return
    
    if 'start_button' in globals() and start_button and start_button.winfo_exists():
        top.after(0, lambda: start_button.config(state="normal"))
        
    if 'status_window' in globals() and status_window.winfo_exists():
        print("⚠️ Download window is already open! Skipping new instance.")
        return

    stop_flag = False  

    status_window = tk.Toplevel(top)
    status_window.title("Download Status")
    status_window.geometry("700x400")

    if start_button and start_button.winfo_exists():
        print("✅ Disabling Start button...")
        start_button["state"] = "disabled"

    canvas = tk.Canvas(status_window)
    scrollbar = ttk.Scrollbar(status_window, orient="vertical", command=canvas.yview)
    scrollable_frame = ttk.Frame(canvas)

    scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    canvas.configure(yscrollcommand=scrollbar.set)

    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

    status_labels = []
    progress_bars = []

    for i, run_acc in enumerate(download_items, start=0):
        frame = tk.Frame(scrollable_frame)
        frame.pack(fill="x", pady=5)

        tk.Label(frame, text=run_acc, width=15, anchor="w").pack(side="left", padx=2)

        progress = ttk.Progressbar(frame, mode="indeterminate", length=250)
        progress.pack(side="left", padx=5)

        status_label = tk.Label(frame, text="Waiting...", width=40, anchor="w", justify="left")
        status_label.pack(side="left", padx=5, fill="x", expand=True)

        status_labels.append(status_label)
        progress_bars.append(progress)

    def check_completion():
        global completed_items, failed_files, stop_flag
    
        if stop_flag:
            print("🚨 Download stopped. Exiting check_completion().")
            return
    
        total_completed = len(completed_items) + len(failed_files)
        print(f"DEBUG: {total_completed}/{len(download_items)} downloads completed.")
    
        if total_completed >= len(download_items):
            print("DEBUG: All downloads completed!")
    
            if failed_files:
                message = f"Download completed with errors.\nFailed files: {', '.join(failed_files)}"
                top.after(0, lambda: messagebox.showwarning("Warning", message))  
            else:
                top.after(0, lambda: notify_completion(top, "✅ All files downloaded successfully.", status_window))
    
            completed_items.clear()
            failed_files.clear()
            active_threads.clear()
            stop_flag = True  
    
            if start_button and start_button.winfo_exists():
                top.after(0, lambda: start_button.config(state="normal"))
        else:
            print("DEBUG: Waiting for remaining downloads... Retrying in 15 seconds.")
            top.after(10000, check_completion)

    def force_delete_file(file_path):
        for proc in psutil.process_iter(['pid', 'name']):
            try:
                for file in proc.open_files():
                    if file.path == file_path:
                        print(f"⚠️ File {file_path} is in use by process {proc.info['name']} (PID: {proc.info['pid']}). Terminating process...")
                        proc.terminate()
                        proc.wait()
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
                print(f"🗑 Successfully deleted: {file_path}")
            except Exception as e:
                print(f"⚠️ Failed to delete {file_path}: {e}")
    
    def on_close(start_button):
        global completed_items, failed_files, active_threads, stop_flag, processes
    
        confirm_close = messagebox.askyesno(
            "Confirm Exit",
            "Downloads are still in progress.\nDo you really want to stop and close the window?"
        )
    
        if not confirm_close:
            print("⏳ Download window close cancelled by user.")
            return  
    
        print("🚨 Download window closed! Stopping active downloads and cleaning up...")
    
        stop_flag = True  
    
        for thread in active_threads:
            if thread.is_alive():
                print(f"⚠️ Stopping thread for {thread.name}")
                thread.join(timeout=1)
    
        active_threads.clear()
    
        terminate_processes()
    
        if 'start_button' in globals() and start_button:
            print("✅ Re-enabling Start button after download window closed.")
            start_button.config(state="normal")
        else:
            print("⚠️ start_button is not available! Trying alternative method.")
            
            for widget in top.winfo_children():
                if isinstance(widget, tk.Button) and widget.cget("text") == "Download Start":
                    print("✅ Found start_button! Re-enabling now.")
                    widget.config(state="normal")
                    break
    
        status_window.update_idletasks()
        status_window.destroy()
    
        print("✅ Download window closed successfully!")

    def terminate_processes():
        global processes
        for process in processes:
            if process.poll() is None:  
                process.terminate()
                print(f"🚨 Terminated process {process.pid}")
        processes.clear()

    status_window.protocol("WM_DELETE_WINDOW", lambda: on_close(start_button))

    if selected_db == "SRA":
        threading.Thread(
            target=download_items_worker,
            args=(download_items, output_dir, status_labels, progress_bars, top, status_window, selected_format)
        ).start()
    elif selected_db == "ENA":
        for index, run_acc in enumerate(download_items):
            threading.Thread(
                target=download_and_convert_ena_fastq,
                args=(run_acc, output_dir, status_labels, progress_bars, index, top, len(download_items), status_window)
            ).start()

    top.after(1000, check_completion)

def download_items_worker(download_items, output_dir, status_labels, progress_bars, top, status_window, selected_format):
    global completed_items, failed_files
    prefetch_output_dir = r"C:\SRA_Virtual_Folder\prefetch_output"
    os.makedirs(prefetch_output_dir, exist_ok=True)

    total_items = len(download_items)

    for index, run_acc in enumerate(download_items):
        threading.Thread(
            target=download_and_convert_sra_to_fastq,
            args=(
                run_acc, output_dir, status_labels, progress_bars, index,
                top, total_items, status_window,
                None, None, None, selected_format
            )
        ).start()
       
def update_accession_list_count(accession_tree, search_count_label):
    count = len(accession_tree.get_children())
    search_count_label.config(text=f"Search List: {count} items")

def search_sra(top, search_box, email_box, accession_tree, search_button, search_progress, load_more_button):
    query = search_box.get().strip()
    if not query:
        messagebox.showerror("Error", "Please enter a search query.")
        return

    disable_search_ui(search_button, search_progress)

    global search_state
    search_state['query'] = query
    search_state['retstart'] = 0
    search_state['total'] = 0

    search_thread = threading.Thread(
        target=search_and_display,
        args=(
            query,                
            email_box,            
            accession_tree,        
            0,                    
            50,                   
            True,
            search_button,
            search_progress,
            load_more_button       
        )
    )
    search_thread.daemon = True
    search_thread.start()

def search_ena(query, accession_tree, search_button, search_progress, load_more_button):
    
    def run_search():
        global ena_search_results
        ena_search_results = [] 
        
        try:
            disable_search_ui(search_button, search_progress)
            search_progress.start()

            encoded_query = quote(query)
            search_url = (
                f"https://www.ebi.ac.uk/ena/portal/api/search?result=read_run"
                f"&query=study_title%3D{encoded_query}"
                "&fields=run_accession,study_accession,sample_accession,scientific_name,instrument_model,library_strategy"
            )

            print(f"DEBUG: ENA Search URL: {search_url}")

            response = requests.get(search_url)
            response.raise_for_status()

            results = response.text.strip().split("\n")
            print(f"DEBUG: ENA API Response (First 500 chars): {response.text[:500]}")

            if len(results) > 1:
                accession_tree.delete(*accession_tree.get_children())

                for line in results[1:]:
                    fields = line.split("\t")
                    if len(fields) >= 6:
                        metadata = {
                            "Run Accession": fields[0],
                            "Title": fields[3],  
                            "Platform Model": fields[4],
                            "Total Bases": "N/A",
                            "Total Reads": "N/A",
                            "Organism": fields[3],  
                            "Library Strategy": fields[5],
                            "Bioproject": fields[1],
                            "Biosample": fields[2]
                        }

                        ena_search_results.append(metadata)

                        accession_tree.insert("", "end", values=(
                            fields[3], fields[4], "", "", fields[3], fields[5], fields[1], fields[2], fields[0]
                        ))
            else:
                messagebox.showwarning("Warning", "No results found in ENA database.")

        except requests.exceptions.RequestException as e:
            messagebox.showerror("Error", f"ENA search failed: {e}")

        finally:
            search_progress.stop()
            update_accession_list_count(accession_tree, search_count_label)
            enable_search_ui(search_button, search_progress)

    threading.Thread(target=run_search).start()
    
def save_email_to_config(email):
    config_path = os.path.join(os.getcwd(), "config.json")
    try:
        with open(config_path, "w") as config_file:
            json.dump({"email": email}, config_file)
    except Exception as e:
        print(f"Failed to save email to config: {e}")

def load_email_from_config():
    config_path = os.path.join(os.getcwd(), "config.json")
    if os.path.exists(config_path):
        try:
            with open(config_path, "r") as config_file:
                data = json.load(config_file)
                return data.get("email", "")
        except Exception as e:
            print(f"Failed to load email from config: {e}")
    return ""    

search_results = []  

def is_valid_email(email):
    pattern = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
    return re.match(pattern, email)

def set_email(email_box):
    email = email_box.get()
    if email:
        Entrez.email = email
        save_email_to_config(email)
        messagebox.showinfo("Info", "Email saved successfully!")
    else:
        messagebox.showerror("Error", "Please enter a valid email.")

def search_ngs_data(top, search_box, email_box, accession_tree, search_button, search_progress, load_more_button, database_combo, tree_frame):
    query = search_box.get().strip()
    database = database_combo.get() 

    if not query:
        messagebox.showerror("Error", "Please enter a search query.")
        return

    print(f"DEBUG: Searching in {database} with query: {query}")

    retry_count = 5
    while retry_count > 0:
        if accession_tree and accession_tree.winfo_exists():
            break
        print(f"WARNING: accession_tree is not ready yet. Retrying... ({retry_count})")
        time.sleep(0.1)
        retry_count -= 1

    if not accession_tree or not accession_tree.winfo_exists():
        print("ERROR: accession_tree does not exist! Trying to find it in GUI...")

        for widget in tree_frame.winfo_children():
            if isinstance(widget, ttk.Treeview):
                print("DEBUG: Found missing accession_tree. Reassigning...")
                accession_tree = widget
                break

    if not accession_tree or not accession_tree.winfo_exists():
        messagebox.showerror("Error", "Search list is not available. Please retry.")
        print("ERROR: accession_tree does not exist! Aborting search.")
        return

    try:
        accession_tree.delete(*accession_tree.get_children())
        update_accession_list_count(accession_tree, search_count_label)
    except Exception as e:
        print(f"WARNING: Unable to clear accession_tree. Reason: {e}")

    if database == "SRA":
        search_sra(top, search_box, email_box, accession_tree, search_button, search_progress, load_more_button)
    elif database == "ENA":
        search_ena(query, accession_tree, search_button, search_progress, load_more_button)
    else:
        messagebox.showerror("Error", "Selected database is not supported.")    
    
is_searching = False  

def clear_search_results(search_button, search_progress, accession_tree, load_more_button, search_box):
    global is_searching, search_state, search_results

    is_searching = False  
    search_results.clear()  

    enable_search_ui(search_button, search_progress)
    search_progress.stop()  

    try:
        accession_tree.delete(*accession_tree.get_children()) 
        update_accession_list_count(accession_tree, search_count_label)
        print("DEBUG: Search results cleared successfully!")
    except Exception as e:
        print(f"WARNING: Failed to clear accession_tree. Error: {e}")

    search_box.delete(0, tk.END)
    print("DEBUG: Search entry cleared.")

def search_and_display(query, email_box, accession_tree, retstart, retmax, is_new_search, search_button, search_progress, load_more_button):
    global search_state, search_results, is_searching
    is_searching = True

    user_email = email_box.get().strip()
    if not user_email:
        messagebox.showerror("Error", "Please enter a valid email before searching.")
        enable_search_ui(search_button, search_progress)
        return

    Entrez.email = user_email

    try:
        handle = Entrez.esearch(db="sra", term=query, retmax=retmax, retstart=retstart)
        record = Entrez.read(handle)
        handle.close()

        search_state['total'] = int(record['Count'])
        search_state['retstart'] += len(record['IdList'])

        if is_new_search:
            accession_tree.delete(*accession_tree.get_children())
            search_results.clear()

        ids = ",".join(record["IdList"])
        summary_handle = Entrez.esummary(db="sra", id=ids, rettype="full")
        summary_records = Entrez.read(summary_handle)
        summary_handle.close()

        results_batch = []

        for summary_record in summary_records:
            exp_xml = summary_record.get('ExpXml', '')
            platform_model = re.search(r'instrument_model="(.*?)"', exp_xml)
            total_bases = re.search(r'total_bases="(\d+)"', exp_xml)
            total_spots = re.search(r'total_spots="(\d+)"', exp_xml)
            experiment_name = re.search(r'<Title>(.*?)</Title>', exp_xml)
            organism_name = re.search(r'ScientificName="(.*?)"', exp_xml)
            library_strategy = re.search(r'<LIBRARY_STRATEGY>(.*?)</LIBRARY_STRATEGY>', exp_xml)
            bioproject = re.search(r'<Bioproject>(.*?)</Bioproject>', exp_xml)
            biosample = re.search(r'<Biosample>(.*?)</Biosample>', exp_xml)
            runs_acc = re.search(r'acc="(.*?)"', summary_record.get('Runs', ''))

            platform_model = platform_model.group(1) if platform_model else ''
            total_bases = format_number(total_bases.group(1)) if total_bases else ''
            total_reads = format_number(int(total_spots.group(1)) * 2) if total_spots and 'library_layout="PAIRED"' in exp_xml else format_number(int(total_spots.group(1))) if total_spots else ''
            experiment_name = experiment_name.group(1) if experiment_name else ''
            organism_name = organism_name.group(1) if organism_name else ''
            library_strategy = library_strategy.group(1) if library_strategy else ''
            bioproject = bioproject.group(1) if bioproject else ''
            biosample = biosample.group(1) if biosample else ''
            runs_acc = runs_acc.group(1) if runs_acc else ''

            results_batch.append((
                experiment_name, platform_model, total_bases, total_reads, organism_name,
                library_strategy, bioproject, biosample, runs_acc
            ))

            search_results.append({
                "Run Accession": runs_acc, "Title": experiment_name, "Platform Model": platform_model,
                "Total Bases": total_bases, "Total Reads": total_reads, "Organism": organism_name,
                "Library Strategy": library_strategy, "Bioproject": bioproject, "Biosample": biosample
            })

        for result in results_batch:
            accession_tree.insert("", "end", values=result)

    except Exception as e:
        messagebox.showerror("Error", f"Failed to search SRA database: {e}")
    
    update_accession_list_count(accession_tree, search_count_label)
    enable_search_ui(search_button, search_progress)

def load_more_results(accession_tree, search_button, search_progress, load_more_button, email_box, selected_db):
    global search_state

    if search_state['retstart'] < search_state['total']:
        disable_search_ui(search_button, search_progress)

        if selected_db == "SRA":
            threading.Thread(
                target=search_and_display,
                args=(
                    search_state['query'],     
                    email_box,                 
                    accession_tree,             
                    search_state['retstart'],   
                    search_state['retmax'],     
                    False,                     
                    search_button,             
                    search_progress,            
                    load_more_button           
                )
            ).start()
    else:
        messagebox.showinfo("Info", "No more results.")

def on_load_more_clicked(tree_frame, search_button, search_progress, load_more_button, email_box, database_combo):
    global accession_tree

    for widget in tree_frame.winfo_children():
        if isinstance(widget, ttk.Treeview):
            accession_tree = widget
            print("DEBUG: Found latest accession_tree for Load More.")
            break

    if not accession_tree or not accession_tree.winfo_exists():
        messagebox.showerror("Error", "Search results are not available. Please retry.")
        print("ERROR: accession_tree does not exist! Aborting load more.")
        return

    selected_db = database_combo.get()
    if selected_db == "SRA":
        load_more_results(accession_tree, search_button, search_progress, load_more_button, email_box, selected_db)
    else:
        messagebox.showinfo("Info", "ENA does not support Load More functionality.")

    update_accession_list_count(accession_tree, search_count_label)

def disable_search_ui(search_button, search_progress):
    search_button["state"] = "disabled"
    search_progress.start()

def enable_search_ui(search_button, search_progress):
    search_button["state"] = "normal"
    search_progress.stop()

def select_directory(selected_directory_label):
    output_dir = filedialog.askdirectory()
    selected_directory_label.config(text=f"{output_dir}", bg="#F6F7F8")

def disable_download_ui(download_button, download_progress):
    download_button["state"] = "disabled"
    download_progress.start()

def enable_download_ui(download_button, download_progress):
    download_button["state"] = "normal"
    download_progress.stop()

def run_command_with_retries(command, retries=3):
    for attempt in range(1, retries + 1):
        try:
            print(f"Attempt {attempt}/{retries}: Running command: {command}")
            subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            return True
        except subprocess.CalledProcessError as e:
            print(f"Attempt {attempt} failed: {e}")
            if attempt == retries:
                print(f"Command failed after {retries} attempts: {command}")
                return False
            time.sleep(5)  
    return False

def check_disk_space(directory, required_space_gb):
    try:
        total, used, free = shutil.disk_usage(directory)
        free_gb = free / (1024 ** 3)
        print(f"Available: {free_gb:.2f} GB, Required: {required_space_gb:.2f} GB")
        return free_gb >= required_space_gb
    except Exception as e:
        print(f"Error checking disk space: {e}")
        return False
    
def check_disk_space_periodically(directory, required_space_gb, interval=60):
    global stop_flag

    while not stop_flag:
        if not check_disk_space(directory, required_space_gb):
            print("🚨 Disk space is running low! Stopping downloads.")
            stop_flag = True  
            return  

        for _ in range(interval):
            if stop_flag:
                print("⏹ Disk space check stopped.")
                return
            time.sleep(5) 
        
def calculate_total_required_space(run_acc_list, prefetch_output_dir):
    total_required_space_gb = 0
    for run_acc in run_acc_list:
        sra_file_path = os.path.join(prefetch_output_dir, run_acc, f"{run_acc}.sra")
        if not os.path.exists(sra_file_path):
            print(f"File {sra_file_path} not found after prefetch.")
            continue
        try:
            sra_file_size_gb = os.path.getsize(sra_file_path) / (1024 ** 3)
            estimated_space = sra_file_size_gb * 5  
            total_required_space_gb += estimated_space
            print(f"File size of {sra_file_path}: {sra_file_size_gb:.2f} GB (Estimated space: {estimated_space:.2f} GB)")
        except Exception as e:
            print(f"Error calculating size for {sra_file_path}: {e}")
    return total_required_space_gb

def prefetch_sra_files(run_acc_list, prefetch_output_dir):
    os.makedirs(prefetch_output_dir, exist_ok=True)
    
    def prefetch_single(run_acc):
        sra_file_path = os.path.join(prefetch_output_dir, run_acc, f"{run_acc}.sra")
        if os.path.exists(sra_file_path):
            print(f"{run_acc}.sra already exists. Skipping prefetch.")
            return True
        
        prefetch_cmd = f'prefetch "{run_acc}" --output-directory "{prefetch_output_dir}" --max-size 500G'
        return run_command_with_retries(prefetch_cmd) 
        print(f"Running prefetch: {prefetch_cmd}")
        try:
            subprocess.run(prefetch_cmd, shell=True, check=True)
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error prefetching {run_acc}: {e}")
            return False
    
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(prefetch_single, run_acc) for run_acc in run_acc_list]
        for future in futures:
            if not future.result():
                print("Prefetch failed for a run.")
            
def check_total_disk_space_and_prefetch(top, run_acc_list, prefetch_output_dir):
    global stop_flag
    prefetch_sra_files(run_acc_list, prefetch_output_dir)
    
    if stop_flag:  
        print("🚨 Download stopped by user. Skipping prefetch retries.")
        return False

    total_required_space_gb = calculate_total_required_space(run_acc_list, prefetch_output_dir)
    if not check_disk_space(prefetch_output_dir, total_required_space_gb):
        raise RuntimeError(f"Not enough disk space. Required: {total_required_space_gb:.2f} GB")

    for run_acc in run_acc_list:
        sra_file_path = os.path.join(prefetch_output_dir, run_acc, f"{run_acc}.sra")
        if os.path.exists(sra_file_path):
            sra_file_size_gb = os.path.getsize(sra_file_path) / (1024 ** 3)
            estimated_space = sra_file_size_gb * 5
            total_required_space_gb += estimated_space
    return total_required_space_gb

def run_conversion_command(primary_cmd, fallback_cmd=None, retries=2):
    global stop_flag, processes

    for attempt in range(1, retries + 1):
        if stop_flag:  
            print("🚨 Download stopped by user. Aborting conversion.")
            return False

        print(f"Attempt {attempt}/{retries} for command: {primary_cmd}")

        try:
            process = subprocess.Popen(primary_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            processes.append(process)  
            stdout, stderr = process.communicate()

            if process.returncode == 0:
                print(f"✅ Primary command succeeded: {primary_cmd}")
                return True
            else:
                print(f"❌ Primary command failed on attempt {attempt}: {stderr.decode()}")

        except Exception as e:
            print(f"❌ Error running primary command on attempt {attempt}: {e}")

        if attempt < retries:
            print("⏳ Retrying primary command in 10 seconds...")
            time.sleep(3)  

    if fallback_cmd:
        for fallback_attempt in range(1, retries + 1):
            if stop_flag:  
                print("🚨 Download stopped by user. Aborting conversion.")
                return False

            print(f"Attempt {fallback_attempt}/{retries} for fallback command: {fallback_cmd}")

            try:
                process = subprocess.Popen(fallback_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                processes.append(process) 
                stdout, stderr = process.communicate()

                if process.returncode == 0:
                    print(f"✅ Fallback command succeeded: {fallback_cmd}")
                    return True
                else:
                    print(f"❌ Fallback command failed on attempt {fallback_attempt}: {stderr.decode()}")

            except Exception as e:
                print(f"❌ Error running fallback command on attempt {fallback_attempt}: {e}")

            if fallback_attempt < retries:
                print("⏳ Retrying fallback command in 10 seconds...")
                time.sleep(3)

    print(f"❌ Both primary and fallback commands failed after {retries} attempts.")
    return False

def log_metadata(run_acc, output_dir, metadata, elapsed_time, total_size=None):
    global log_file_name  

    if log_file_name is None:
        log_file_name = f"NGS_data_downloader_{datetime.now().strftime('%Y%m%d_%H%M')}.csv"

    log_file_path = os.path.join(output_dir, log_file_name)

    with log_lock:  
        file_exists = os.path.exists(log_file_path)

        with open(log_file_path, mode="a", newline="") as log_file:
            writer = csv.writer(log_file)

            if not file_exists:
                writer.writerow([
                    "Run Accession", "Title", "Platform Model", "Total Bases", "Total Reads",
                    "Organism", "Library Strategy", "Bioproject", "Biosample", "Elapsed Time", "Size (MB)"
                ])

            formatted_metadata = {
                "Title": metadata.get("Title", "N/A"),
                "Platform Model": metadata.get("Platform Model", "N/A"),
                "Total Bases": metadata.get("Total Bases", "N/A"),
                "Total Reads": metadata.get("Total Reads", "N/A"),
                "Organism": metadata.get("Organism", "N/A"),
                "Library Strategy": metadata.get("Library Strategy", "N/A"),
                "Bioproject": metadata.get("Bioproject", "N/A"),
                "Biosample": metadata.get("Biosample", "N/A"),
            }

            if total_size is None:
                final_output_dir = os.path.abspath(os.path.join(output_dir, run_acc))
                if os.path.exists(final_output_dir):
                    total_size = sum(
                        os.path.getsize(os.path.join(final_output_dir, f)) for f in os.listdir(final_output_dir)
                    ) / (1024 * 1024)  
                else:
                    total_size = "N/A"  

            formatted_size = f"{total_size:.2f} MB" if isinstance(total_size, (int, float)) else "N/A"

            writer.writerow([
                run_acc,
                formatted_metadata["Title"], formatted_metadata["Platform Model"],
                formatted_metadata["Total Bases"], formatted_metadata["Total Reads"],
                formatted_metadata["Organism"], formatted_metadata["Library Strategy"],
                formatted_metadata["Bioproject"], formatted_metadata["Biosample"],
                elapsed_time, formatted_size
            ])

    print(f"✅ Metadata logged successfully for {run_acc}")
            
def download_and_convert_sra_to_fastq(run_acc, output_dir, status_labels, progress_bars, index, 
                                      top, total_items, status_window, download_log, start_time_str, 
                                      download_list_tree, selected_format):
    global completed_items, failed_files, active_threads, stop_flag, start_button  

    if run_acc in completed_items:
        print(f"⚠️ {run_acc} is already in completed items. Skipping re-download.")
        return

    if stop_flag:
        print("🚨 Download stopped by user. Aborting SRA conversion.")
        return

    if not start_time_str:
        start_time_str = datetime.now().strftime("%Y%m%d_%H%M%S")

    task_start_time = time.time()
    print(f"Task started for {run_acc} at {datetime.now().strftime('%m-%d %H:%M:%S')}")

    prefetch_output_dir = r"C:\SRA_Virtual_Folder\prefetch_output"
    sra_folder_path = os.path.join(prefetch_output_dir, run_acc)
    sra_file_path = os.path.join(sra_folder_path, f"{run_acc}.sra")

    thread = threading.current_thread()
    active_threads.append(thread)

    def update_progress(progress, message):
        if status_window.winfo_exists():
            top.after(0, lambda: progress_bars[index].config(mode="determinate", value=progress))
            top.after(0, lambda: status_labels[index].config(text=message))
            
    try:
        update_progress(5, "Checking disk space and prefetching SRA file...")

        if stop_flag:
            print("🚨 Download stopped by user. Aborting.")
            return False

        if not check_total_disk_space_and_prefetch(top, [run_acc], prefetch_output_dir):
            raise RuntimeError(f"Prefetch failed or insufficient disk space for {run_acc}")

        if stop_flag:
            print("🚨 Download stopped by user. Aborting SRA conversion.")
            return False

        if not os.path.exists(sra_file_path):
            raise FileNotFoundError(f"Prefetch completed, but SRA file does not exist at: {sra_file_path}")

        print(f"Using SRA file: {sra_file_path}")

        final_output_dir = os.path.abspath(os.path.join(output_dir, run_acc))
        os.makedirs(final_output_dir, exist_ok=True)
        if not os.access(final_output_dir, os.W_OK):
            raise PermissionError(f"Cannot write to output directory: {final_output_dir}")

        update_progress(10, "Preparing directories...")

        conversion_success = False
        
        if selected_format == "FASTQ":
            update_progress(30, "Running fastq-dump...")

            fastq_cmd = f'fastq-dump --split-files --gzip "{sra_file_path}" --outdir "{final_output_dir}"'
            conversion_success = run_conversion_command(fastq_cmd)

        elif selected_format == "FASTA":
            update_progress(30, "Running fastq-dump for FASTA format...")
            fasta_cmd = f'fastq-dump "{sra_file_path}" --outdir "{final_output_dir}" --split-files --fasta 0'
            conversion_success = run_conversion_command(fasta_cmd)

        if not conversion_success:
            raise RuntimeError(f"fastq-dump failed for {run_acc}")

        if stop_flag:
            print("🚨 Download stopped by user. Aborting SRA conversion.")
            return False

        update_progress(80, "Writing metadata to log...")

        metadata = next((item for item in search_results if item["Run Accession"] == run_acc), {})

        log_metadata(run_acc, output_dir, metadata, f"{(time.time() - task_start_time) // 60}m")

        update_progress(100, "Task completed successfully!")

    except Exception as e:
        error_message = str(e)
        print(f"Error for {run_acc}: {error_message}")
        print(traceback.format_exc())
        failed_files.append(run_acc)
        update_progress(0, "Error")

        if status_window.winfo_exists():
            top.after(0, lambda: messagebox.showerror("Error", f"An error occurred for {run_acc}:\n{error_message}"))

    finally:
        completed_items.append(run_acc)
        print(f"✅ DEBUG: {run_acc} added to completed_items. Now: {len(completed_items)} completed.")

        if thread in active_threads:
            active_threads.remove(thread)

        try:
            shutil.rmtree(sra_folder_path)
            print(f"🗑 Deleted SRA folder: {sra_folder_path}")
        except Exception as e:
            print(f"⚠️ Failed to delete SRA folder {sra_folder_path}: {e}")

        if 'start_button' in globals() and start_button and start_button.winfo_exists():
            top.after(0, lambda: start_button.config(state="normal"))
        else:
            print("⚠️ start_button is not available! Skipping state change.")

def get_ena_fastq_url(run_acc):
    api_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={run_acc}&result=read_run&fields=fastq_ftp"
    response = requests.get(api_url)
    if response.status_code == 200 and response.text.strip():
        fastq_ftp_urls = response.text.strip().split("\n")[1]  
        return fastq_ftp_urls.split(";")  
    return []

def download_and_convert_ena_fastq(run_acc, output_dir, status_labels, progress_bars, index, top, total_items, status_window):
    global completed_items, failed_files, ena_search_results

    task_start_time = time.time()
    print(f"Task started for {run_acc} at {datetime.now().strftime('%m-%d %H:%M:%S')}")

    if index >= len(progress_bars) or index >= len(status_labels):
        print(f"⚠️ IndexError 방지: index={index}, progress_bars={len(progress_bars)}, status_labels={len(status_labels)}")
        return

    def update_progress(progress, message):
        if status_labels[index].winfo_exists():
            top.after(0, lambda: progress_bars[index].config(mode="determinate", value=progress))
            top.after(0, lambda: status_labels[index].config(text=message))

    elapsed_time = 0  
    error_message = ""  
    final_output_dir = os.path.abspath(os.path.join(output_dir, run_acc))
    os.makedirs(final_output_dir, exist_ok=True)  

    try:
        update_progress(10, "Fetching ENA file paths...")
        
        api_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={run_acc}&result=read_run&fields=fastq_ftp"
        response = requests.get(api_url)
        response.raise_for_status()

        lines = response.text.strip().split("\n")
        if len(lines) < 2:
            raise FileNotFoundError(f"⚠️ No FASTQ files found for {run_acc} in ENA.")

        fastq_ftp_urls = lines[1].strip().split("\t", 1)[-1]
        fastq_files = fastq_ftp_urls.split(";") if ";" in fastq_ftp_urls else [fastq_ftp_urls]

        if not fastq_files or fastq_files[0].strip() == "":
            raise FileNotFoundError(f"⚠️ No FASTQ files found for {run_acc} in ENA.")

        update_progress(30, "Preparing to download...")

        for file_url in fastq_files:
            https_url = f"https://{file_url.strip()}"
            file_name = os.path.basename(file_url.strip())
            output_file_path = os.path.join(final_output_dir, file_name)

            update_progress(50, "Downloading...")
            
            retry_count = 3
            for attempt in range(retry_count):
                try:
                    with requests.get(https_url, stream=True) as http_response:
                        http_response.raise_for_status()
                        with open(output_file_path, "wb") as file:
                            for chunk in http_response.iter_content(chunk_size=8192):
                                file.write(chunk)
                    break
                except requests.RequestException as err:
                    print(f"❌ Download attempt {attempt + 1} failed: {err}")
                    if attempt == retry_count - 1:
                        error_message = f"❌ {run_acc}: Failed to download {file_name} after {retry_count} attempts."
                        failed_files.append(run_acc)
                        update_progress(0, "Download Failed")
                        top.after(0, lambda: messagebox.showerror("Error", error_message))
                        return
                    time.sleep(5)  

        update_progress(100, "Download completed!")
        completed_items.append(run_acc)

    except Exception as e:
        error_message = f"An error occurred for {run_acc}:\n{str(e)}"
        print(f"❌ {error_message}")
        traceback.print_exc()
        failed_files.append(run_acc)
        update_progress(0, "Error")

    finally:
        elapsed_time = f"{(time.time() - task_start_time) // 60}m"

        if os.path.exists(final_output_dir):
            total_size = sum(
                os.path.getsize(os.path.join(final_output_dir, f)) for f in os.listdir(final_output_dir)
            ) / (1024 * 1024)  
        else:
            total_size = "N/A"

        metadata = next((item for item in ena_search_results if item["Run Accession"] == run_acc), {})

        print(f"✅ Download completed for {run_acc} in {elapsed_time}")

        log_metadata(run_acc, output_dir, metadata, elapsed_time, total_size)

        if error_message:
            top.after(0, lambda: messagebox.showerror("Error", error_message))

def select_all_items(tree):
    for item in tree.get_children():
        tree.selection_add(item)

def SRA_to_FASTQ_downloader_GUI():
    global accession_tree, download_list_tree, search_count_label, start_button  
   
    top = tk.Toplevel()
    top.title("SRA to FASTQ Downloader")
    top.geometry("1500x750")
    top.resizable(False, False)
    top.configure(bg="#F6F7F8")
    
    image_path = get_resource_path("IDL_mark.png")

    if os.path.exists(image_path):  
        image = Image.open(image_path)
        image = image.resize((80, 80)) 
        photo = ImageTk.PhotoImage(image)

        image_label = tk.Label(top, image=photo, bg="#F6F7F8")
        image_label.image = photo  
        image_label.place(relx=0.94, rely=0.02)
    else:
        image_label = tk.Label(top, text="", font=("", 8, ""), bg="#F6F7F8")
        image_label.place(relx=0.9, rely=-0)
    
    email_label = tk.Label(top, text="Enter your email:", bg="#F6F7F8")
    email_label.grid(row=0, column=0, padx=10, pady=10, sticky="w")
    email_label.place(relx=0.34, rely=0.01)

    email_box = tk.Entry(top)
    email_box.grid(row=0, column=1, sticky="w")
    email_box.place(relx=0.42, rely=0.01, relwidth=0.13)
    email_box.insert(0, load_email_from_config())  

    def set_email():
        email = email_box.get()
        if email:
            Entrez.email = email
            save_email_to_config(email)
            messagebox.showinfo("Info", "Email saved successfully!")
        else:
            messagebox.showerror("Error", "Please enter a valid email.")

    save_email_button = tk.Button(top, cursor="hand2", bg="#E0EDF9", font=("Helvetica", 9), borderwidth=0.5, text="Save Email", command=set_email, width=10)
    save_email_button.grid(row=0, column=2)
    save_email_button.place(relx=0.56, rely=0.01)
    
    format_label = tk.Label(top, text="Select format:", bg="#F6F7F8")
    format_label.grid(row=1, column=0, padx=10, pady=10, sticky="w")
    format_label.place(relx=0.42, rely=0.88)

    format_combo = ttk.Combobox(top, values=["FASTQ", "FASTA"], state="readonly", width=10, cursor="hand2")
    format_combo.grid(row=1, column=1, padx=10, pady=10, sticky="w")
    format_combo.set("FASTQ")  # Set default value
    format_combo.place(relx=0.48, rely=0.88)
    
    def download_selected(download_list_frame, selected_directory_label, download_button, download_label, 
                          top, format_combo, database_combo, start_button):
        global download_list_tree
    
        latest_download_list_tree = None
        for widget in download_list_frame.winfo_children():
            if isinstance(widget, ttk.Treeview):
                latest_download_list_tree = widget
                print("DEBUG: Found latest download_list_tree for download.")
                break
    
        if not latest_download_list_tree or not latest_download_list_tree.winfo_exists():
            messagebox.showerror("Error", "Download list is not available. Please retry.")
            print("ERROR: download_list_tree does not exist! Aborting download.")
            return
    
        output_dir = selected_directory_label.cget("text").split(": ")[0].strip()
        if output_dir == 'None' or not os.path.isdir(output_dir):
            messagebox.showerror("Error", "Please select a valid output directory.")
            return
    
        selected_items = latest_download_list_tree.get_children()
        if not selected_items:
            messagebox.showerror("Error", "No accession selected!")
            return
    
        download_items = [latest_download_list_tree.item(item)['values'][-1] for item in selected_items]
    
        if any(acc == "N/A" or not acc for acc in download_items):
            messagebox.showerror("Error", "One or more items have an invalid Run Accession!")
            return
    
        selected_format = format_combo.get()
        selected_db = database_combo.get()  
    
        print(f"DEBUG: Calling open_download_status_window with {len(download_items)} items, format: {selected_format}, DB: {selected_db}")
    
        if 'start_button' in globals() and start_button and start_button.winfo_exists():
            print("✅ Re-enabling Start button before download.")
            start_button.config(state="normal")
    
        open_download_status_window(top, download_items, output_dir, latest_download_list_tree, selected_format, selected_db, start_button)
    
    database_label = tk.Label(top, text="Select Database:", bg="#F6F7F8")
    database_label.place(relx=0.18, rely=0.01)

    database_combo = ttk.Combobox(top, values=["SRA", "ENA"], state="readonly", width=5, cursor="hand2")
    database_combo.place(relx=0.25, rely=0.01)
    database_combo.set("SRA")  
    
    search_label = tk.Label(top, text="Enter search query:", bg="#F6F7F8")
    search_label.place(relx=0.34, rely=0.06)
    
    search_box = tk.Entry(top)
    search_box.place(relx=0.42, rely=0.06, relwidth=0.13)
    
    search_button = tk.Button(
        top, text="Search", command=lambda: search_ngs_data(
            top, search_box, email_box, accession_tree, search_button, search_progress, load_more_button, database_combo, tree_frame  # 🔹 tree_frame 추가
        ), cursor="hand2", bg="#E0EDF9", font=("Helvetica", 9), borderwidth=0.5, width=10, height=1
    )

    search_button.place(relx=0.56, rely=0.06)
    
    search_box.bind('<Return>', lambda event: search_button.invoke())
    
    clear_button = tk.Button(
        top, text="Clear Search", 
        command=lambda: clear_search_results(search_button, search_progress, accession_tree, load_more_button, search_box),
        cursor="hand2", bg="#E0EDF9", font=("Helvetica", 9), borderwidth=0.5, width=12, height=1
        )
    clear_button.place(relx=0.65, rely=0.06)

    search_progress = ttk.Progressbar(top, mode="indeterminate", maximum=50) 
    search_progress.place(relx=0.42, rely=0.1, relwidth=0.13)
    
    load_more_button = tk.Button(
        top, cursor="hand2", bg="#E0EDF9", font=("Helvetica", 9), borderwidth=0.5,
        text="Load More",
        command=lambda: on_load_more_clicked(tree_frame, search_button, search_progress, load_more_button, email_box, database_combo),
        width=10
    )
    load_more_button.place(relx=0.56, rely=0.1)
    
    def toggle_load_more(database):
        if database == "SRA":
            load_more_button["state"] = "normal"
            load_more_button.place(relx=0.56, rely=0.1)  
        elif database == "ENA":
            load_more_button["state"] = "disabled"
            load_more_button.place_forget() 

        print(f"DEBUG: Load More button state set for {database}.")
        if len(accession_tree.get_children()) == 0:
            print("DEBUG: Search results still empty after toggle_load_more.")
        else:
            print("DEBUG: Search results NOT empty after toggle_load_more.")
    
    def on_database_change(event):
        global accession_tree, download_list_tree, is_searching, search_state, search_results

        clear_search_results(search_button, search_progress, accession_tree, load_more_button, search_box)  

        selected_db = database_combo.get()
        print(f"DEBUG: on_database_change triggered, Database changed to {selected_db}")

        toggle_load_more(selected_db)

        try:
            if download_list_tree and download_list_tree.winfo_exists():
                download_list_tree.delete(*download_list_tree.get_children())
                download_count_label.config(text="Download List: 0 items", bg="#F6F7F8")
            else:
                raise NameError
        except NameError:
            print("DEBUG: download_list_tree is not available! Recreating...")
            for widget in download_list_frame.winfo_children():
                widget.destroy()

            columns = ("Title", "Platform Model", "Total Bases", "Total Reads", "Organism", 
                       "Library Strategy", "Bioproject", "Biosample", "Runs Accession")

            download_list_tree = ttk.Treeview(download_list_frame, columns=columns, show="headings", height=10)

            for col in columns:
                download_list_tree.heading(col, text=col)
                if col in ["Total Bases", "Total Reads"]:
                    download_list_tree.column(col, anchor="e")
                else:
                    download_list_tree.column(col, anchor="center")

            scrollbar_dl = ttk.Scrollbar(download_list_frame, orient="vertical", command=download_list_tree.yview)
            download_list_tree.configure(yscrollcommand=scrollbar_dl.set)
            scrollbar_dl.pack(side='right', fill='y')
            download_list_tree.pack(side='left', fill='both', expand=True)

        try:
            if accession_tree and accession_tree.winfo_exists():
                accession_tree.delete(*accession_tree.get_children())
                update_accession_list_count(accession_tree, search_count_label)
            else:
                raise NameError
        except NameError:
            print("DEBUG: accession_tree is not available! Recreating...")
            for widget in tree_frame.winfo_children():
                widget.destroy()

            accession_tree = ttk.Treeview(tree_frame, columns=columns, show="headings")

            for col in columns:
                accession_tree.heading(col, text=col)
                if col in ["Total Bases", "Total Reads"]:
                    accession_tree.column(col, anchor="e")
                else:
                    accession_tree.column(col, anchor="center")

            scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=accession_tree.yview)
            accession_tree.configure(yscrollcommand=scrollbar.set)
            scrollbar.pack(side='right', fill='y')
            accession_tree.pack(side='left', fill='both', expand=True)

        completed_items.clear()  
        failed_files.clear()  
        search_results.clear()
        update_accession_list_count(accession_tree, search_count_label)
        print("DEBUG: Reset completed!")
        
        top.update_idletasks()

    database_combo.bind("<<ComboboxSelected>>", on_database_change)
    
    tree_frame = tk.Frame(top)
    tree_frame.place(relx=0.05, rely=0.15, relwidth=0.9, relheight=0.35)
    
    search_count_label = tk.Label(top, text="Search List: 0 items", bg="#F6F7F8")
    search_count_label.place(relx=0.05, rely=0.12)

    select_all_accession_button = tk.Button(
        top, text="Select All", bg="#E0EDF9", cursor="hand2",
        command=lambda: select_all_items(accession_tree)
        )
    select_all_accession_button.place(relx=0.05, rely=0.08)

    columns = ("Title", "Platform Model", "Total Bases", "Total Reads", "Organism", "Library Strategy", "Bioproject", "Biosample", "Runs Accession")
    accession_tree = ttk.Treeview(tree_frame, columns=columns, show="headings")
    
    accession_tree.column("Title", width=250)
    accession_tree.column("Platform Model", width=150, anchor="center")
    accession_tree.column("Total Bases", width=120, anchor="e")
    accession_tree.column("Total Reads", width=120, anchor="e")
    accession_tree.column("Organism", width=150, anchor="center")
    accession_tree.column("Library Strategy", width=120, anchor="center")
    accession_tree.column("Bioproject", width=115, anchor="center")
    accession_tree.column("Biosample", width=115, anchor="center")
    accession_tree.column("Runs Accession", width=120, anchor="center")
    
    for col in columns:
        accession_tree.heading(col, text=col)
        if col in ["Total Bases", "Total Reads"]:
            accession_tree.column(col, anchor="e")
        else:
            accession_tree.column(col, anchor="center")
    
    scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=accession_tree.yview)
    accession_tree.configure(yscrollcommand=scrollbar.set)
    scrollbar.pack(side='right', fill='y')
    accession_tree.pack(side='left', fill='both', expand=True)
    
    update_accession_list_count(accession_tree, search_count_label)
    
    move_button = tk.Button(
        top, cursor="hand2", bg="#E0EDF9", font=("Helvetica", 9), borderwidth=0.5,
        text="Move to Download List",
        command=lambda: move_to_download_list(tree_frame, download_list_frame),  # ✅ download_list_frame 추가
        width=25
    )
    move_button.place(relx=0.35, rely=0.53)

    download_list_frame = tk.Frame(top)
    download_list_frame.place(relx=0.05, rely=0.6, relwidth=0.9, relheight=0.2)
        
    download_list_tree = ttk.Treeview(download_list_frame, columns=columns, show="headings", height=10)
    download_list_tree.column("Title", width=250)
    download_list_tree.column("Platform Model", width=150, anchor="center")
    download_list_tree.column("Total Bases", width=120, anchor="e")
    download_list_tree.column("Total Reads", width=120, anchor="e")
    download_list_tree.column("Organism", width=150, anchor="center")
    download_list_tree.column("Library Strategy", width=120, anchor="center")
    download_list_tree.column("Bioproject", width=115, anchor="center")
    download_list_tree.column("Biosample", width=115, anchor="center")
    download_list_tree.column("Runs Accession", width=120, anchor="center")
    
    for col in columns:
        download_list_tree.heading(col, text=col)
        if col in ["Total Bases", "Total Reads"]:
            download_list_tree.column(col, anchor="e")
        else:
            download_list_tree.column(col, anchor="center")
    
    scrollbar_dl = ttk.Scrollbar(download_list_frame, orient="vertical", command=download_list_tree.yview)
    download_list_tree.configure(yscrollcommand=scrollbar_dl.set)
    scrollbar_dl.pack(side='right', fill='y')
    download_list_tree.pack(side='left', fill='both', expand=True)

    download_count_label = tk.Label(top, text="Download List: 0 items", bg="#F6F7F8")
    download_count_label.place(relx=0.05, rely=0.57)
    
    select_all_download_button = tk.Button(
        top, text="Select All", bg="#E0EDF9", cursor="hand2",
        command=lambda: select_all_items(download_list_tree)
        )
    select_all_download_button.place(relx=0.05, rely=0.53)    

    def update_download_count(download_list_frame):
        latest_download_list_tree = None
        for widget in download_list_frame.winfo_children():
            if isinstance(widget, ttk.Treeview):
                latest_download_list_tree = widget
                print("DEBUG: Found latest download_list_tree for counting.")
                break
    
        if not latest_download_list_tree or not latest_download_list_tree.winfo_exists():
            print("ERROR: Cannot update count. download_list_tree does not exist!")
            return
    
        count = len(latest_download_list_tree.get_children())
        download_count_label.config(text=f"Download List: {count} items")
    
    def move_to_download_list(tree_frame, download_list_frame):
        global accession_tree, download_list_tree
    
        for widget in tree_frame.winfo_children():
            if isinstance(widget, ttk.Treeview):
                accession_tree = widget
                print("DEBUG: Found latest accession_tree.")
                break
    
        if not accession_tree or not accession_tree.winfo_exists():
            messagebox.showerror("Error", "Search list is not available. Please retry.")
            print("ERROR: accession_tree does not exist! Aborting move.")
            return
    
        latest_download_list_tree = None
        for widget in download_list_frame.winfo_children():
            if isinstance(widget, ttk.Treeview):
                latest_download_list_tree = widget
                print("DEBUG: Found latest download_list_tree.")
                break
    
        if not latest_download_list_tree or not latest_download_list_tree.winfo_exists():
            messagebox.showerror("Error", "Download list is not available. Please retry.")
            print("ERROR: download_list_tree does not exist! Aborting move.")
            return
    
        selected_items = accession_tree.selection()
        if not selected_items:
            messagebox.showerror("Error", "No accession selected!")
            return
    
        for item in selected_items:
            values = accession_tree.item(item)['values']
            accession_tree.delete(item)
            latest_download_list_tree.insert("", "end", values=values)
    
        update_download_count(download_list_frame) 
        update_accession_list_count(accession_tree, search_count_label)  

    def remove_from_download_list(download_list_frame):
        latest_download_list_tree = None
        for widget in download_list_frame.winfo_children():
            if isinstance(widget, ttk.Treeview):
                latest_download_list_tree = widget
                print("DEBUG: Found latest download_list_tree for removal.")
                break
    
        if not latest_download_list_tree or not latest_download_list_tree.winfo_exists():
            messagebox.showerror("Error", "Download list is not available. Please retry.")
            print("ERROR: download_list_tree does not exist! Aborting removal.")
            return
    
        selected_items = latest_download_list_tree.selection()
        if not selected_items:
            messagebox.showerror("Error", "No item selected!")
            return
    
        for item in selected_items:
            latest_download_list_tree.delete(item)
    
        update_download_count(download_list_frame)
        update_accession_list_count(accession_tree, search_count_label) 

    remove_button = tk.Button(
        top, cursor="hand2", bg="#E0EDF9", font=("Helvetica", 9), borderwidth=0.5,
        text="Remove from Download List",
        command=lambda: remove_from_download_list(download_list_frame),
        width=25
    )
    remove_button.place(relx=0.5, rely=0.53)
    
    file_select_button = tk.Button(top, cursor="hand2", bg="#E0EDF9", font=("Helvetica", 9), borderwidth=0.5, text="Select Directory", command=lambda: select_directory(selected_directory_label), width=25)
    file_select_button.place(relx=0.42, rely=0.81)
            
    selected_directory_label = tk.Label(top, text=": None", bg="#F6F7F8", font=("",9,"italic"))
    selected_directory_label.place(relx=0.55, rely=0.81)
    
    download_label = tk.Label(top, text="")
    
    download_button = tk.Button(
        top, font=("Helvetica", 12, "bold"), relief="ridge", cursor="hand2", bg="#4183BC", fg="white",
        activebackground="#4183BC", activeforeground="white", borderwidth=3, text="Download Start", 
        command=lambda: download_selected(download_list_frame, selected_directory_label, download_button, 
                                          download_label, top, format_combo, database_combo, download_button),  
        width=20
    )


    download_button.place(relx=0.41, rely=0.92)
    
    top.mainloop() 
    
def start_gui_after_sratoolkit_check():
    top = tk.Tk()
    top.withdraw() 

    check_and_install_sratoolkit(top)
