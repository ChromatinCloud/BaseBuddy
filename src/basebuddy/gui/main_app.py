import customtkinter
import tkinter.filedialog
from pathlib import Path
import threading
import queue # For thread communication
from typing import Dict, Any, Optional, Callable, List # Added List
import copy # For deepcopying params for runners
import socket # For IGV communication
import re # For basic range validation

# Assuming runner and utils are accessible via PYTHONPATH
from basebuddy import runner as src_bb_runner # Updated import
# Ensure signature_utils is imported if apply_signature_to_fasta needs it for path resolution
from basebuddy import signature_utils as sig_utils # Updated import
from basebuddy.utils import BaseBuddyInputError, BaseBuddyToolError, BaseBuddyFileError, BaseBuddyConfigError # Updated import

# Hardcoded for now, ideally from a shared config or bb_runners
KNOWN_ART_PROFILES = {
    "illumina": ["HS25", "HSXt", "HSXn", "MSv1", "MSv3", "MiS", "NS50"],
}
DEFAULT_ILLUMINA_PROFILE = "HS25"
KNOWN_NANOSIM_MODELS = [
    "dna_r9.4.1_e8.1", "dna_r10.3_e8.2", "rna_r9.4.1_e8.1_cdna", "rna_r9.4.1_e8.1_native", "GRCh38_wg_R9_guppy4"
]
DEFAULT_NANOSIM_MODEL = "dna_r9.4.1_e8.1"
BUNDLED_SIG_IDS = ["SBS1", "SBS5", "Custom Path..."] # Added "Custom Path..."


class BaseBuddyGUI(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        self.title("BaseBuddy GUI")
        self.geometry("800x850") # Increased height slightly for new tab

        customtkinter.set_appearance_mode("System")
        customtkinter.set_default_color_theme("blue")

        self.tab_view = customtkinter.CTkTabview(self)
        self.tab_view.pack(padx=10, pady=10, fill="both", expand=True)

        self.short_read_tab = self.tab_view.add("Short Read Sim")
        self.spike_tab = self.tab_view.add("Variant Spiking")
        self.long_read_tab = self.tab_view.add("Long Read Sim")
        self.apply_sig_tab = self.tab_view.add("Apply Signature") # New Tab

        self.thread_queue = queue.Queue()
        self.spike_variants_results: Optional[Dict[str, Any]] = None

        self.status_label = customtkinter.CTkLabel(self, text="Status:")
        self.status_label.pack(padx=20, pady=(5,0), anchor="w")
        self.status_textbox = customtkinter.CTkTextbox(self, height=120, state="disabled", wrap="word")
        self.status_textbox.pack(padx=20, pady=(0,10), fill="x", expand=False)

        self.create_short_read_tab()
        self.create_spike_variants_tab()
        self.create_long_read_tab()
        self.create_apply_signature_tab() # Create the new tab

        self.poll_thread()

    def update_status(self, message: str, is_error: bool = False, clear_first: bool = False):
        self.status_textbox.configure(state="normal")
        if clear_first:
            self.status_textbox.delete("1.0", "end")
        self.status_textbox.insert("end", message + "\n")
        self.status_textbox.configure(state="disabled")
        self.status_textbox.see("end")


    def poll_thread(self):
        try:
            message_type, data, button_to_enable = self.thread_queue.get_nowait()
            if message_type == "completion":
                result_dict, success = data

                if success:
                    if button_to_enable == self.run_spike_button:
                        self.spike_variants_results = result_dict
                        if self.spike_variants_results.get("igv_session_file"):
                            self.igv_button_spike.configure(state="normal")

                run_name = result_dict.get("run_name", "Unknown run")
                output_dir = result_dict.get("output_directory", "N/A")
                msg = f"Run '{run_name}' completed {'successfully' if success else 'with issues'}.\n"
                msg += f"Output Directory: {output_dir}\n"

                if "output_files" in result_dict and result_dict["output_files"]: # Check if list is not empty
                    msg += "Output Files:\n"
                    for f_info in result_dict["output_files"]:
                        msg += f"  - {f_info['name']} ({f_info['type']}): {f_info['path']}\n"

                # Specific outputs for apply_signature_to_fasta
                if "output_modified_fasta_path" in result_dict:
                    msg += f"Mutated FASTA: {result_dict['output_modified_fasta_path']}\n"
                if result_dict.get("num_mutations_applied") is not None:
                     msg += f"Mutations Applied: {result_dict['num_mutations_applied']}\n"

                if "manifest_path" in result_dict: msg += f"Manifest: {result_dict['manifest_path']}\n"
                if "output_bam" in result_dict: msg += f"Output BAM: {result_dict['output_bam']}\n"
                if result_dict.get("output_bam_index"): msg += f"Output BAI: {result_dict['output_bam_index']}\n"
                if "igv_session_file" in result_dict: msg += f"IGV Session: {result_dict['igv_session_file']}\n"

                self.update_status(msg, is_error=not success, clear_first=True)

            elif message_type == "error":
                error_message = data
                self.update_status(f"ERROR: {error_message}", is_error=True, clear_first=True)

            elif message_type == "info":
                self.update_status(data, is_error=False, clear_first=False)

            if button_to_enable:
                button_to_enable.configure(state="normal")
                self._re_enable_all_run_buttons(except_button=button_to_enable)

        except queue.Empty:
            pass
        finally:
            self.after(100, self.poll_thread)

    def _create_path_entry(self, parent, label_text, row, dialog_func, is_file=True, entry_var=None): # Added entry_var
        customtkinter.CTkLabel(master=parent, text=label_text).grid(row=row, column=0, padx=10, pady=5, sticky="w")
        # Use entry_var if provided, otherwise create a new one
        entry_widget = customtkinter.CTkEntry(master=parent, width=350, textvariable=entry_var)
        entry_widget.grid(row=row, column=1, padx=10, pady=5, sticky="ew")

        def browse_callback():
            path = tkinter.filedialog.askopenfilename() if is_file else tkinter.filedialog.askdirectory()
            if path:
                if entry_var:
                    entry_var.set(path)
                else:
                    entry_widget.delete(0, "end"); entry_widget.insert(0, path)

        customtkinter.CTkButton(master=parent, text="Browse", width=80, command=browse_callback).grid(row=row, column=2, padx=5, pady=5)
        return entry_widget # Return the widget itself if no var passed, or if needed

    def create_short_read_tab(self):
        tab = self.short_read_tab
        tab.grid_columnconfigure(1, weight=1)
        row_idx = 0

        self.sr_ref_fasta_entry = self._create_path_entry(tab, "Reference FASTA:", row_idx, tkinter.filedialog.askopenfilename)
        row_idx += 1
        self.sr_output_root_entry = self._create_path_entry(tab, "Output Root Dir:", row_idx, tkinter.filedialog.askdirectory, is_file=False)
        row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Run Name (optional):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sr_run_name_entry = customtkinter.CTkEntry(master=tab); self.sr_run_name_entry.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2); row_idx += 1

        customtkinter.CTkLabel(master=tab, text="ART Platform:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sr_art_platform_var = customtkinter.StringVar(value="illumina")
        customtkinter.CTkOptionMenu(master=tab, variable=self.sr_art_platform_var, values=["illumina"]).grid(row=row_idx, column=1, padx=10, pady=5, sticky="w"); row_idx += 1

        customtkinter.CTkLabel(master=tab, text="ART Profile:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sr_art_profile_var = customtkinter.StringVar(value=DEFAULT_ILLUMINA_PROFILE)
        customtkinter.CTkOptionMenu(master=tab, variable=self.sr_art_profile_var, values=KNOWN_ART_PROFILES["illumina"]).grid(row=row_idx, column=1, padx=10, pady=5, sticky="w"); row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Read Length:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); self.sr_read_length_entry = customtkinter.CTkEntry(master=tab); self.sr_read_length_entry.insert(0, "150"); self.sr_read_length_entry.grid(row=row_idx, column=1, sticky="ew", columnspan=2, padx=10, pady=5); row_idx+=1
        customtkinter.CTkLabel(master=tab, text="Depth:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); self.sr_depth_entry = customtkinter.CTkEntry(master=tab); self.sr_depth_entry.insert(0, "50"); self.sr_depth_entry.grid(row=row_idx, column=1, sticky="ew", columnspan=2, padx=10, pady=5); row_idx+=1
        customtkinter.CTkLabel(master=tab, text="Mean Fragment:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); self.sr_mean_frag_entry = customtkinter.CTkEntry(master=tab); self.sr_mean_frag_entry.insert(0, "400"); self.sr_mean_frag_entry.grid(row=row_idx, column=1, sticky="ew", columnspan=2, padx=10, pady=5); row_idx+=1
        customtkinter.CTkLabel(master=tab, text="Std Dev Fragment:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); self.sr_std_dev_frag_entry = customtkinter.CTkEntry(master=tab); self.sr_std_dev_frag_entry.insert(0, "50"); self.sr_std_dev_frag_entry.grid(row=row_idx, column=1, sticky="ew", columnspan=2, padx=10, pady=5); row_idx+=1

        customtkinter.CTkLabel(master=tab, text="Genomic Ranges (optional, e.g., chr1:1k-2k):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.short_sim_ranges_textbox = customtkinter.CTkTextbox(master=tab, height=80, wrap="word")
        self.short_sim_ranges_textbox.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2)
        self.short_sim_ranges_textbox.insert("1.0", "chr1:1000-2000\nchrM:1-16569") # Example placeholder
        row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Timeout (s):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); self.sr_timeout_entry = customtkinter.CTkEntry(master=tab); self.sr_timeout_entry.insert(0, "3600.0"); self.sr_timeout_entry.grid(row=row_idx, column=1, sticky="ew", columnspan=2, padx=10, pady=5); row_idx+=1

        self.sr_single_end_var = customtkinter.BooleanVar(value=False)
        customtkinter.CTkCheckBox(master=tab, text="Single-End Simulation", variable=self.sr_single_end_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); row_idx += 1
        self.sr_overwrite_var = customtkinter.BooleanVar(value=False)
        customtkinter.CTkCheckBox(master=tab, text="Overwrite Output", variable=self.sr_overwrite_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); row_idx += 1
        self.sr_auto_index_var = customtkinter.BooleanVar(value=True)
        customtkinter.CTkCheckBox(master=tab, text="Auto-Index FASTA", variable=self.sr_auto_index_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); row_idx += 1

        self.run_short_sim_button = customtkinter.CTkButton(master=tab, text="Run Simulation", command=self.run_short_simulation_thread)
        self.run_short_sim_button.grid(row=row_idx, column=0, columnspan=3, padx=10, pady=(10,0))

    def run_short_simulation_thread(self):
        if hasattr(self, 'igv_button_spike'): self.igv_button_spike.configure(state="disabled")
        self._start_runner_thread(
            button_to_disable=self.run_short_sim_button,
            runner_func=src_bb_runner.simulate_short,
            param_extractor_func=self._get_short_sim_params,
            result_handler_func=self._handle_short_sim_results
        )

    def _get_short_sim_params(self) -> Dict[str, Any]:
        ref_fasta = self.sr_ref_fasta_entry.get(); output_root = self.sr_output_root_entry.get()
        if not ref_fasta:
            self.update_status("Reference FASTA path is required.", is_error=True, clear_first=True)
            raise ValueError("Reference FASTA path is required.")
        if not output_root:
            self.update_status("Output root directory is required.", is_error=True, clear_first=True)
            raise ValueError("Output root directory is required.")
        run_name_str = self.sr_run_name_entry.get() or None

        # Parse genomic ranges from textbox
        ranges_text = self.short_sim_ranges_textbox.get("1.0", "end-1c").strip()
        parsed_genomic_ranges = None
        if ranges_text:
            parsed_genomic_ranges = [line.strip() for line in ranges_text.splitlines() if line.strip()]
            if not parsed_genomic_ranges: # All lines were whitespace
                parsed_genomic_ranges = None
            else:
                # Basic validation for format (example: chr:start-end)
                range_pattern = re.compile(r"^[\w.-]+:\d+-\d+$")
                for r_str in parsed_genomic_ranges:
                    if not range_pattern.fullmatch(r_str):
                        msg = f"Malformed genomic range: '{r_str}'. Expected format like 'chr1:1000-2000'."
                        self.update_status(msg, is_error=True, clear_first=True)
                        raise ValueError(msg)
                self.update_status(f"Parsed {len(parsed_genomic_ranges)} genomic ranges.", clear_first=True)

        command_params_for_manifest = {
            "id_prefix": "gui_sim_reads", "no_aln_output": False,
            "reference_fasta": ref_fasta, "depth": int(self.sr_depth_entry.get()),
            "read_length": int(self.sr_read_length_entry.get()), "art_profile": self.sr_art_profile_var.get(),
            "mean_fragment_length": int(self.sr_mean_frag_entry.get()),
            "std_dev_fragment_length": int(self.sr_std_dev_frag_entry.get()),
            "is_paired_end": not self.sr_single_end_var.get(),
            "art_platform": self.sr_art_platform_var.get(),
            "overwrite_output": self.sr_overwrite_var.get(),
            "output_root_dir": output_root, "run_name": run_name_str,
            "auto_index_fasta": self.sr_auto_index_var.get(),
            "genomic_ranges_from_gui": parsed_genomic_ranges # For manifest
        }
        return {
            "output_root_dir": Path(output_root), "run_name": run_name_str,
            "command_params": command_params_for_manifest, "reference_fasta": ref_fasta,
            "depth": int(self.sr_depth_entry.get()), "read_length": int(self.sr_read_length_entry.get()),
            "art_profile": self.sr_art_profile_var.get(),
            "mean_fragment_length": int(self.sr_mean_frag_entry.get()),
            "std_dev_fragment_length": int(self.sr_std_dev_frag_entry.get()),
            "is_paired_end": not self.sr_single_end_var.get(),
            "overwrite_output": self.sr_overwrite_var.get(),
            "art_platform": self.sr_art_platform_var.get(),
            "timeout": float(self.sr_timeout_entry.get()),
            "auto_index_fasta": self.sr_auto_index_var.get(),
            "variants_list": None, # Not implemented in this tab's GUI yet
            "genomic_ranges": parsed_genomic_ranges
        }

    def _handle_short_sim_results(self, result_dict: Dict[str, Any]):
        pass

    def create_spike_variants_tab(self):
        tab = self.spike_tab
        tab.grid_columnconfigure(1, weight=1)
        row_idx = 0

        self.sv_ref_fasta_entry = self._create_path_entry(tab, "Reference FASTA:", row_idx, tkinter.filedialog.askopenfilename)
        row_idx += 1
        self.sv_input_bam_entry = self._create_path_entry(tab, "Input BAM:", row_idx, tkinter.filedialog.askopenfilename)
        row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Variants to Spike (chr pos ref alt):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.spike_variants_textbox = customtkinter.CTkTextbox(master=tab, height=100, wrap="word")
        self.spike_variants_textbox.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2)
        self.spike_variants_textbox.insert("1.0", "chr1 12345 A T\nchr2 67890 C G") # Example placeholder
        row_idx += 1

        self.sv_output_root_entry = self._create_path_entry(tab, "Output Root Dir:", row_idx, tkinter.filedialog.askdirectory, is_file=False)
        row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Run Name (optional):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sv_run_name_entry = customtkinter.CTkEntry(master=tab); self.sv_run_name_entry.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2); row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Output BAM Prefix:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sv_output_prefix_entry = customtkinter.CTkEntry(master=tab); self.sv_output_prefix_entry.insert(0, "spiked_output"); self.sv_output_prefix_entry.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2); row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Timeout (s):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); self.sv_timeout_entry = customtkinter.CTkEntry(master=tab); self.sv_timeout_entry.insert(0, "7200.0"); self.sv_timeout_entry.grid(row=row_idx, column=1, sticky="ew", columnspan=2, padx=10, pady=5); row_idx+=1

        self.sv_overwrite_var = customtkinter.BooleanVar(value=False)
        customtkinter.CTkCheckBox(master=tab, text="Overwrite Output", variable=self.sv_overwrite_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); row_idx += 1
        self.sv_auto_index_input_bam_var = customtkinter.BooleanVar(value=True)
        customtkinter.CTkCheckBox(master=tab, text="Auto-Index Input BAM", variable=self.sv_auto_index_input_bam_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); row_idx += 1
        self.sv_auto_index_fasta_var = customtkinter.BooleanVar(value=True)
        customtkinter.CTkCheckBox(master=tab, text="Auto-Index FASTA (Ref)", variable=self.sv_auto_index_fasta_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); row_idx += 1

        button_frame = customtkinter.CTkFrame(master=tab, fg_color="transparent")
        button_frame.grid(row=row_idx, column=0, columnspan=3, pady=(10,0))

        self.run_spike_button = customtkinter.CTkButton(master=button_frame, text="Run Spiking", command=self.run_spike_variants_thread)
        self.run_spike_button.pack(side="left", padx=5)

        self.igv_button_spike = customtkinter.CTkButton(master=button_frame, text="Open in IGV", command=self.on_igv_button_spike_click, state="disabled")
        self.igv_button_spike.pack(side="left", padx=5)

    def run_spike_variants_thread(self):
        if hasattr(self, 'igv_button_spike'): self.igv_button_spike.configure(state="disabled")
        self._start_runner_thread(
            button_to_disable=self.run_spike_button,
            runner_func=src_bb_runner.spike_variants,
            param_extractor_func=self._get_spike_variants_params,
            result_handler_func=self._handle_spike_results
        )

    def _get_spike_variants_params(self) -> Optional[Dict[str, Any]]:
        ref_fasta = self.sv_ref_fasta_entry.get()
        input_bam_str = self.sv_input_bam_entry.get() # Assuming single BAM input for now from GUI
        output_root = self.sv_output_root_entry.get()

        if not all([ref_fasta, input_bam_str, output_root]):
            self.update_status("Reference FASTA, Input BAM, and Output Root Dir are required.", is_error=True, clear_first=True)
            return None

        output_bam_prefix = self.sv_output_prefix_entry.get()
        if not output_bam_prefix:
            self.update_status("Output BAM prefix is required.", is_error=True, clear_first=True)
            return None

        run_name_str = self.sv_run_name_entry.get() or None

        variants_text = self.spike_variants_textbox.get("1.0", "end-1c").strip()
        parsed_variants_list = []
        if not variants_text:
            self.update_status("No variants provided in the textbox.", is_error=True, clear_first=True)
            return None

        lines = variants_text.splitlines()
        for i, line in enumerate(lines):
            line = line.strip()
            if not line: continue
            parts = line.split()
            if len(parts) != 4:
                self.update_status(f"Error parsing variant line {i+1}: '{line}'. Expected 4 components (chr pos ref alt).", is_error=True, clear_first=True)
                return None
            chrom, pos_str, ref, alt = parts
            if not pos_str.isdigit():
                self.update_status(f"Error parsing variant line {i+1}: Position '{pos_str}' is not an integer.", is_error=True, clear_first=True)
                return None
            parsed_variants_list.append({
                "chromosome": chrom, "position": int(pos_str),
                "ref_allele": ref, "alt_allele": alt
            })

        if not parsed_variants_list:
            self.update_status("No valid variants parsed from the textbox.", is_error=True, clear_first=True)
            return None

        self.update_status(f"Parsed {len(parsed_variants_list)} variants from input.", clear_first=True)

        vaf_val = 0.5
        seed_val = 0

        command_params_for_manifest = {
            "reference_fasta": ref_fasta,
            "input_bams_gui": [input_bam_str],
            "num_variants_from_gui": len(parsed_variants_list),
            "output_prefix_for_bam": output_bam_prefix,
            "overwrite_output": self.sv_overwrite_var.get(),
            "auto_index_input_bam": self.sv_auto_index_input_bam_var.get(),
            "auto_index_fasta": self.sv_auto_index_fasta_var.get(),
            "timeout": float(self.sv_timeout_entry.get()),
            "output_root_dir": output_root, "run_name": run_name_str,
            "vaf": vaf_val, "seed": seed_val
        }

        return {
            "output_root_dir": Path(output_root), "run_name": run_name_str,
            "command_params": command_params_for_manifest,
            "reference_fasta": ref_fasta,
            "input_bams": [input_bam_str],
            "variants_list": parsed_variants_list,
            "output_prefix_for_bam": output_bam_prefix,
            "overwrite_output": self.sv_overwrite_var.get(),
            "auto_index_input_bam": self.sv_auto_index_input_bam_var.get(),
            "auto_index_fasta": self.sv_auto_index_fasta_var.get(),
            "timeout": float(self.sv_timeout_entry.get())
        }

    def _handle_spike_results(self, result_dict: Dict[str, Any]):
        self.spike_variants_results = result_dict
        if self.spike_variants_results.get("igv_session_file"):
            self.igv_button_spike.configure(state="normal")

    def on_igv_button_spike_click(self):
        if not self.spike_variants_results:
            self.update_status("No variant spiking results available to send.", is_error=True, clear_first=True)
            return
        session_file = self.spike_variants_results.get("igv_session_file")
        if session_file and Path(session_file).exists():
            self.update_status(f"Attempting to load IGV session: {session_file}", clear_first=True)
            self.send_to_igv_desktop(session_file_path=session_file)
        else:
            self.update_status("Could not find IGV session file in results.", is_error=True, clear_first=True)

    def send_to_igv_desktop(self, session_file_path: Optional[str] = None,
                            reference_fasta_path: Optional[str] = None,
                            bam_file_path: Optional[str] = None):
        try:
            self.update_status("Connecting to IGV on localhost:60151...", clear_first=True)
            s = socket.create_connection(("localhost", 60151), timeout=3)
            if session_file_path and Path(session_file_path).exists():
                abs_session_path = str(Path(session_file_path).resolve())
                cmd_str = f"load \"{abs_session_path}\"\n"
                s.sendall(cmd_str.encode())
                self.update_status(f"Sent session file to IGV: {abs_session_path}")
            else:
                self.update_status("No valid session file provided for IGV.", is_error=True); s.close(); return
            s.close()
            self.update_status("Commands sent to IGV. Check IGV application.", is_error=False)
        except socket.timeout: self.update_status("Error: IGV connection timed out.", is_error=True)
        except socket.error as e: self.update_status(f"Error connecting to IGV: {e}.", is_error=True)
        except Exception as e: self.update_status(f"Error sending to IGV: {e}", is_error=True)

    def create_long_read_tab(self):
        tab = self.long_read_tab
        tab.grid_columnconfigure(1, weight=1)
        row_idx = 0

        self.lr_ref_fasta_entry = self._create_path_entry(tab, "Reference FASTA:", row_idx, tkinter.filedialog.askopenfilename)
        row_idx += 1
        self.lr_output_root_entry = self._create_path_entry(tab, "Output Root Dir:", row_idx, tkinter.filedialog.askdirectory, is_file=False)
        row_idx += 1
        customtkinter.CTkLabel(master=tab, text="Run Name (optional):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.lr_run_name_entry = customtkinter.CTkEntry(master=tab); self.lr_run_name_entry.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2); row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Depth:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); self.lr_depth_entry = customtkinter.CTkEntry(master=tab); self.lr_depth_entry.insert(0, "30"); self.lr_depth_entry.grid(row=row_idx, column=1, sticky="ew", columnspan=2, padx=10, pady=5); row_idx+=1
        customtkinter.CTkLabel(master=tab, text="NanoSim Model:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w");
        self.lr_nanosim_model_var = customtkinter.StringVar(value=DEFAULT_NANOSIM_MODEL)
        customtkinter.CTkOptionMenu(master=tab, variable=self.lr_nanosim_model_var, values=KNOWN_NANOSIM_MODELS).grid(row=row_idx, column=1, padx=10, pady=5, sticky="w"); row_idx += 1
        customtkinter.CTkLabel(master=tab, text="Timeout (s):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); self.lr_timeout_entry = customtkinter.CTkEntry(master=tab); self.lr_timeout_entry.insert(0, "3600.0"); self.lr_timeout_entry.grid(row=row_idx, column=1, sticky="ew", columnspan=2, padx=10, pady=5); row_idx+=1

        self.lr_overwrite_var = customtkinter.BooleanVar(value=False)
        customtkinter.CTkCheckBox(master=tab, text="Overwrite Output", variable=self.lr_overwrite_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); row_idx += 1
        self.lr_auto_index_fasta_var = customtkinter.BooleanVar(value=True)
        customtkinter.CTkCheckBox(master=tab, text="Auto-Index FASTA", variable=self.lr_auto_index_fasta_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); row_idx += 1

        self.run_long_sim_button = customtkinter.CTkButton(master=tab, text="Run Simulation", command=self.run_long_simulation_thread)
        self.run_long_sim_button.grid(row=row_idx, column=0, columnspan=3, padx=10, pady=20)

    def run_long_simulation_thread(self):
        if hasattr(self, 'igv_button_spike'): self.igv_button_spike.configure(state="disabled")
        self._start_runner_thread(
            button_to_disable=self.run_long_sim_button,
            runner_func=src_bb_runner.simulate_long,
            param_extractor_func=self._get_long_sim_params,
            result_handler_func=None
        )

    def _get_long_sim_params(self) -> Dict[str, Any]:
        ref_fasta = self.lr_ref_fasta_entry.get(); output_root = self.lr_output_root_entry.get()
        if not ref_fasta: raise ValueError("Reference FASTA path is required.")
        if not output_root: raise ValueError("Output root directory is required.")
        run_name_str = self.lr_run_name_entry.get() or None
        depth_val = int(self.lr_depth_entry.get())
        model_val = self.lr_nanosim_model_var.get()
        timeout_val = float(self.lr_timeout_entry.get())
        overwrite_val = self.lr_overwrite_var.get()
        auto_index_val = self.lr_auto_index_fasta_var.get()

        command_params_for_runner = {
            "reference_fasta": ref_fasta, "depth": depth_val, "model": model_val,
            "overwrite_output": overwrite_val, "auto_index_fasta": auto_index_val,
            "timeout": timeout_val, "output_root_dir": output_root, "run_name": run_name_str,
        }
        return {
            "output_root_dir": Path(output_root), "run_name": run_name_str,
            "command_params": command_params_for_runner, "reference_fasta": ref_fasta,
            "depth": depth_val, "model": model_val,
            "overwrite_output": overwrite_val, "auto_index_fasta": auto_index_val,
            "timeout": timeout_val
            # variants_list and genomic_ranges will be passed as None by default by _execute_runner if not in dict
        }

    def _start_runner_thread(self, button_to_disable: customtkinter.CTkButton, runner_func: Callable,
                             param_extractor_func: Callable, result_handler_func: Optional[Callable] = None):
        # Disable all run buttons
        if hasattr(self, 'run_short_sim_button'): self.run_short_sim_button.configure(state="disabled")
        if hasattr(self, 'run_spike_button'): self.run_spike_button.configure(state="disabled")
        if hasattr(self, 'run_long_sim_button'): self.run_long_sim_button.configure(state="disabled")
        if hasattr(self, 'run_apply_sig_button'): self.run_apply_sig_button.configure(state="disabled")

        button_to_disable.configure(state="disabled")


        self.update_status("Processing...", clear_first=True)

        try:
            runner_args = param_extractor_func()
            if runner_args is None: # Param extractor indicated an error
                 self._re_enable_all_run_buttons() # Re-enable all since it was a GUI-side validation
                 # The specific button_to_disable should already be re-enabled by the error path in param_extractor or here.
                 button_to_disable.configure(state="normal")
                 return

            thread = threading.Thread(target=self._execute_runner, args=(runner_func, runner_args, button_to_disable, result_handler_func))
            thread.start()
        except ValueError as ve:
            self.thread_queue.put(("error", str(ve), button_to_disable))
        except Exception as e:
             self.thread_queue.put(("error", f"Unexpected setup error: {str(e)}", button_to_disable))


    def _re_enable_all_run_buttons(self, except_button: Optional[customtkinter.CTkButton] = None):
        buttons = []
        if hasattr(self, 'run_short_sim_button'): buttons.append(self.run_short_sim_button)
        if hasattr(self, 'run_spike_button'): buttons.append(self.run_spike_button)
        if hasattr(self, 'run_long_sim_button'): buttons.append(self.run_long_sim_button)
        if hasattr(self, 'run_apply_sig_button'): buttons.append(self.run_apply_sig_button)

        for btn in buttons:
            if btn is not except_button: # Don't re-enable the one that just finished if it's handled by poll_thread
                btn.configure(state="normal")


    def _execute_runner(self, runner_func: Callable, runner_args: Dict[str, Any],
                        button_to_disable: customtkinter.CTkButton, result_handler_func: Optional[Callable]):
        try:
            # Ensure all expected optional args are present if runner_func expects them
            # For example, simulate_short and simulate_long now expect variants_list and genomic_ranges
            if runner_func == src_bb_runner.simulate_short or runner_func == src_bb_runner.simulate_long:
                runner_args.setdefault("variants_list", None)
                runner_args.setdefault("genomic_ranges", None)
            elif runner_func == src_bb_runner.apply_signature_to_fasta:
                 runner_args.setdefault("target_regions", None)


            result = runner_func(**runner_args)
            if result_handler_func:
                result_handler_func(result)
            self.thread_queue.put(("completion", (result, True), button_to_disable))
        except (BaseBuddyInputError, BaseBuddyToolError, BaseBuddyFileError, BaseBuddyConfigError, NameError) as e:
            self.thread_queue.put(("error", str(e), button_to_disable))
        except Exception as e:
            tb_str = traceback.format_exc()
            self.thread_queue.put(("error", f"An unexpected error occurred: {str(e)}\nTraceback:\n{tb_str}", button_to_disable))


def run_gui():
    app = BaseBuddyGUI()
    app.mainloop()

if __name__ == "__main__":
    run_gui()

[end of src/basebuddy/gui/main_app.py]
