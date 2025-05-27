import sys
import os
import traceback
import tkinter as tk
from tkinter import Button
import threading
import psutil
import time
import ctypes
import gc
import multiprocessing

from ui.ISSM import Ui_In_silico_sequence_mining
from ui import resource_rc
from in_silico_sequence_mining.in_silico_sequence_mining import file_3_GUI
from NGS_data_downloader.NGS_data_downloader import start_gui_after_sratoolkit_check

stop_event = threading.Event()

def get_resource_path(relative_path):
    if getattr(sys, 'frozen', False):  
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)

def log_error(error_message):
    with open("error_log.txt", "a", encoding="utf-8") as log_file:
        log_file.write(error_message + "\n")

def apply_fixed_memory_limit_by_ratio(ratio=0.9):
    if os.name == "nt":
        try:
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

            print(f"✅ Fixed memory limit applied: {min_memory / 1e9:.2f}GB ~ {max_memory / 1e9:.2f}GB")
        except Exception as e:
            print(f"❌ Memory limit setting failed: {e}")

def Program_start_GUI():
    root = tk.Tk()
    root.title("In silico project")
    root.geometry("300x150")
    root.configure(bg="#FEFEFE")
    
    def on_closing():
        print("On closing... Set up the flag to stop threads")
        stop_event.set()
        root.destroy()

    browse_button1 = tk.Button(
        root, text="NGS data downloader", command=start_gui_after_sratoolkit_check,
        font=("Helvetica", 10), relief="ridge", cursor="hand2",
        bg="#F0F5F9", fg="black", activebackground="#F5F5F5", borderwidth=3
    )
    browse_button1.place(relx=0.5, rely=0.3, anchor="center", relwidth=0.6)

    browse_button2 = tk.Button(
        root, text="In silico sequence mining", command=file_3_GUI,
        font=("Helvetica", 10), relief="ridge", cursor="hand2",
        bg="#F0F5F9", fg="black", activebackground="#F5F5F5", borderwidth=3
    )
    browse_button2.place(relx=0.5, rely=0.7, anchor="center", relwidth=0.6)

    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.mainloop()

if __name__ == "__main__":
    try:
        multiprocessing.freeze_support()
        
        apply_fixed_memory_limit_by_ratio(ratio=0.9)

        Program_start_GUI()

    except FileNotFoundError as fnf_error:
        log_error(str(fnf_error))
        print(fnf_error)

    except ImportError as imp_error:
        log_error(str(imp_error))
        print(imp_error)

    except Exception:
        log_error(traceback.format_exc())
        raise
