import customtkinter
import tkinter.filedialog
from pathlib import Path
import threading
import queue # For thread communication
from typing import Dict, Any, Optional, Callable, List # Added List
import copy # For deepcopying params for runners
import socket # For IGV communication
import re # For basic range validation
import sys # For platform check in on_open_output_dir_click
import os # For startfile in on_open_output_dir_click
import subprocess # For Popen in on_open_output_dir_click


# Assuming runner and utils are accessible via PYTHONPATH
from src.basebuddy import runner as src_bb_runner
# Ensure signature_utils is imported if apply_signature_to_fasta needs it for path resolution
from src.basebuddy import signature_utils as sig_utils
from bb_utils import BaseBuddyInputError, BaseBuddyToolError, BaseBuddyFileError, BaseBuddyConfigError

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
        self.germline_sim_tab = self.tab_view.add("Germline Simulation") # New Germline Tab

        self.thread_queue = queue.Queue()
        self.spike_variants_results: Optional[Dict[str, Any]] = None

        self.status_label = customtkinter.CTkLabel(self, text="Status:")
        self.status_label.pack(padx=20, pady=(5,0), anchor="w")
        self.status_textbox = customtkinter.CTkTextbox(self, height=120, state="disabled", wrap="word")
        self.status_textbox.pack(padx=20, pady=(0,10), fill="x", expand=False)

        self.open_output_dir_button = customtkinter.CTkButton(self, text="Open Output Directory", command=self.on_open_output_dir_click, state="disabled")
        self.open_output_dir_button.pack(padx=20, pady=(0,10), anchor="w")
        self.last_successful_output_dir: Optional[str] = None

        self.create_short_read_tab()
        self.create_spike_variants_tab()
        self.create_long_read_tab()
        self.create_apply_signature_tab() # Create the new tab
        self.create_germline_sim_tab() # Create Germline Sim Tab

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
                        self.spike_variants_results = result_dict # Store results
                        if self.spike_variants_results.get("igv_session_file"):
                            self.igv_button_spike.configure(state="normal")
                    elif button_to_disable == self.run_germline_sim_button:
                        self.germline_sim_results = result_dict # Store results
                        modified_fasta = result_dict.get("final_modified_fasta_path")
                        sim_results = result_dict.get("simulated_reads_results", {})
                        # Check if sim_results has 'output_files' and it's not empty, or has a BAM path directly
                        has_data_files = "output_files" in sim_results and sim_results["output_files"]

                        if modified_fasta and Path(modified_fasta).exists() and has_data_files:
                            self.igv_button_germline.configure(state="normal")
                        else:
                            self.igv_button_germline.configure(state="disabled")
                            logger.warning("Germline sim IGV button not enabled. FASTA or data files missing/invalid in results.")
                            logger.debug(f"Germline results for IGV check: FASTA='{modified_fasta}', HasData='{has_data_files}', SimResultsKeys='{sim_results.keys()}'")


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

                if success and result_dict.get("output_directory"):
                    self.last_successful_output_dir = result_dict.get("output_directory")
                    if hasattr(self, 'open_output_dir_button'): self.open_output_dir_button.configure(state="normal")


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

        customtkinter.CTkLabel(master=tab, text="Input BAM(s) (one per line):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="nw")
        self.sv_input_bams_textbox = customtkinter.CTkTextbox(master=tab, height=80, wrap="word")
        self.sv_input_bams_textbox.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2)
        row_idx += 1

        self.sv_snp_vcf_entry = self._create_path_entry(tab, "SNP VCF File (optional):", row_idx, tkinter.filedialog.askopenfilename)
        row_idx += 1
        self.sv_indel_vcf_entry = self._create_path_entry(tab, "Indel VCF File (optional):", row_idx, tkinter.filedialog.askopenfilename)
        row_idx += 1
        self.sv_picard_jar_entry = self._create_path_entry(tab, "Picard JAR Path (optional):", row_idx, tkinter.filedialog.askopenfilename)
        row_idx += 1

        customtkinter.CTkLabel(master=tab, text="VAF (Allele Freq):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sv_vaf_entry = customtkinter.CTkEntry(master=tab)
        self.sv_vaf_entry.insert(0, "0.05") # Default VAF
        self.sv_vaf_entry.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2)
        row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Seed (for VAF):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sv_seed_entry = customtkinter.CTkEntry(master=tab)
        self.sv_seed_entry.insert(0, "0") # Default seed
        self.sv_seed_entry.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2)
        row_idx += 1

        self.sv_output_root_entry = self._create_path_entry(tab, "Output Root Dir:", row_idx, tkinter.filedialog.askdirectory, is_file=False)
        row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Run Name (optional):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sv_run_name_entry = customtkinter.CTkEntry(master=tab); self.sv_run_name_entry.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2); row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Output File Prefix:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w") # Changed label from "Output BAM Prefix"
        self.sv_output_prefix_entry = customtkinter.CTkEntry(master=tab); self.sv_output_prefix_entry.insert(0, "spiked_output"); self.sv_output_prefix_entry.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2); row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Timeout (s, for sub-tools):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); self.sv_timeout_entry = customtkinter.CTkEntry(master=tab); self.sv_timeout_entry.insert(0, "7200.0"); self.sv_timeout_entry.grid(row=row_idx, column=1, sticky="ew", columnspan=2, padx=10, pady=5); row_idx+=1

        self.sv_overwrite_var = customtkinter.BooleanVar(value=False)
        customtkinter.CTkCheckBox(master=tab, text="Overwrite Output", variable=self.sv_overwrite_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); row_idx += 1
        self.sv_auto_index_input_bam_var = customtkinter.BooleanVar(value=True)
        customtkinter.CTkCheckBox(master=tab, text="Auto-Index Input BAM(s)", variable=self.sv_auto_index_input_bam_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); row_idx += 1 # Changed label
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
        output_root = self.sv_output_root_entry.get()

        if not ref_fasta:
            self.update_status("Reference FASTA path is required.", is_error=True, clear_first=True); return None
        if not output_root:
            self.update_status("Output Root Directory is required.", is_error=True, clear_first=True); return None

        input_bams_text = self.sv_input_bams_textbox.get("1.0", "end-1c").strip()
        input_bam_paths_str = [line.strip() for line in input_bams_text.splitlines() if line.strip()]
        if not input_bam_paths_str:
            self.update_status("At least one Input BAM path is required.", is_error=True, clear_first=True)
            return None

        snp_vcf_str = self.sv_snp_vcf_entry.get().strip() or None
        indel_vcf_str = self.sv_indel_vcf_entry.get().strip() or None

        if not snp_vcf_str and not indel_vcf_str:
            self.update_status("At least one of SNP VCF or Indel VCF must be provided.", is_error=True, clear_first=True)
            return None

        picard_jar_str = self.sv_picard_jar_entry.get().strip() or None

        output_bam_prefix = self.sv_output_prefix_entry.get()
        if not output_bam_prefix:
            self.update_status("Output File Prefix is required.", is_error=True, clear_first=True)
            return None

        run_name_str = self.sv_run_name_entry.get() or None

        vaf_val_str = self.sv_vaf_entry.get()
        seed_val_str = self.sv_seed_entry.get()
        try:
            vaf_val = float(vaf_val_str)
            if not (0 < vaf_val <= 1.0):
                 self.update_status("VAF must be a number between 0 (exclusive) and 1.0.", is_error=True, clear_first=True); return None
        except ValueError:
            self.update_status(f"Invalid VAF value: '{vaf_val_str}'. Must be a number.", is_error=True, clear_first=True); return None
        try:
            seed_val = int(seed_val_str)
        except ValueError:
            self.update_status(f"Invalid Seed value: '{seed_val_str}'. Must be an integer.", is_error=True, clear_first=True); return None

        command_params_for_manifest = {
            "reference_fasta": ref_fasta,
            "input_bams_gui": input_bam_paths_str,
            "output_prefix_for_bam": output_bam_prefix,
            "overwrite_output": self.sv_overwrite_var.get(),
            "auto_index_input_bam": self.sv_auto_index_input_bam_var.get(),
            "auto_index_fasta": self.sv_auto_index_fasta_var.get(),
            "timeout_per_bam": float(self.sv_timeout_entry.get()),
            "output_root_dir": output_root,
            "run_name": run_name_str,
            "vaf": vaf_val,
            "seed": seed_val,
            "picard_jar": picard_jar_str,
            "snp_vcf_path": snp_vcf_str,
            "indel_vcf_path": indel_vcf_str
        }

        return {
            "output_root_dir": Path(output_root),
            "run_name": run_name_str,
            "command_params": command_params_for_manifest,
            "reference_fasta": ref_fasta,
            "input_bams": input_bam_paths_str,
            "snp_vcf_path": snp_vcf_str,
            "indel_vcf_path": indel_vcf_str,
            "output_prefix_for_bam": output_bam_prefix,
            "overwrite_output": self.sv_overwrite_var.get(),
            "auto_index_input_bam": self.sv_auto_index_input_bam_var.get(),
            "auto_index_fasta": self.sv_auto_index_fasta_var.get()
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
        # Reference path might be useful if session doesn't load genome, or for context.
        ref_fasta = self.spike_variants_results.get("reference_fasta_used")

        if session_file and Path(session_file).exists():
            self.update_status(f"Attempting to load IGV session: {session_file}", clear_first=True)
            # Pass reference as well, in case session needs it or for future use
            self.send_to_igv_desktop(session_file_path=session_file, reference_fasta_path=ref_fasta)
        else:
            self.update_status("Could not find IGV session file in results.", is_error=True, clear_first=True)

    def send_to_igv_desktop(self, session_file_path: Optional[str] = None,
                            reference_fasta_path: Optional[str] = None,
                            data_file_paths: Optional[List[str]] = None,
                            locus: Optional[str] = None): # Added locus
        commands_to_send = []
        # Resolve and check reference FASTA path first if provided, as it might be needed by both scenarios
        abs_ref_path = None
        if reference_fasta_path:
            abs_ref_path = str(Path(reference_fasta_path).resolve())
            if not Path(abs_ref_path).exists():
                self.update_status(f"Error: Reference FASTA for IGV not found: {abs_ref_path}", is_error=True); return

        if session_file_path:
            abs_session_path = str(Path(session_file_path).resolve())
            if not Path(abs_session_path).exists():
                self.update_status(f"Error: Session file not found: {abs_session_path}", is_error=True); return
            # If session specifies a genome, this might be redundant or override. IGV handles this.
            if abs_ref_path:
                 commands_to_send.append(f"genome \"{abs_ref_path}\"")
            commands_to_send.append(f"load \"{abs_session_path}\"")
        elif abs_ref_path and data_file_paths:
            commands_to_send.append(f"genome \"{abs_ref_path}\"")
            for data_file in data_file_paths:
                abs_data_file_path = str(Path(data_file).resolve())
                if not Path(abs_data_file_path).exists():
                    self.update_status(f"Error: Data file for IGV not found: {abs_data_file_path}", is_error=True); return
                commands_to_send.append(f"load \"{abs_data_file_path}\"")
        else:
            self.update_status("Error: Insufficient information for IGV. Provide session or reference + data files.", is_error=True); return

        if locus: # Add goto command if locus is specified
            commands_to_send.append(f"goto {locus}")

        if not commands_to_send: # Should not happen if logic above is correct
            self.update_status("No commands to send to IGV.", is_error=True); return

        try:
            self.update_status(f"Connecting to IGV on localhost:60151...", clear_first=True)
            s = socket.create_connection(("localhost", 60151), timeout=3)
            full_command_str = "\n".join(commands_to_send) + "\n" # Ensure trailing newline
            s.sendall(full_command_str.encode())
            s.close()
            self.update_status(f"Sent commands to IGV: {'; '.join(commands_to_send)}. Check IGV application.", is_error=False)
        except socket.timeout: self.update_status("Error: IGV connection timed out (is IGV running and listening?).", is_error=True)
        except socket.error as e: self.update_status(f"Error connecting to IGV: {e}. (Is IGV running and listening on port 60151?).", is_error=True)
        except Exception as e: self.update_status(f"An unexpected error occurred sending commands to IGV: {e}", is_error=True)

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
        if hasattr(self, 'run_germline_sim_button'): self.run_germline_sim_button.configure(state="disabled") # Ensure this is also disabled

        # Disable open output dir button at start of new run
        if hasattr(self, 'open_output_dir_button'): self.open_output_dir_button.configure(state="disabled")
        if hasattr(self, 'last_successful_output_dir'): self.last_successful_output_dir = None


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
        if hasattr(self, 'run_germline_sim_button'): buttons.append(self.run_germline_sim_button)

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

    def create_germline_sim_tab(self):
        tab = self.germline_sim_tab
        tab.grid_columnconfigure(1, weight=1) # Allow input column to expand
        current_row_idx = 0

        # Core inputs
        self.gs_ref_fasta_entry = self._create_path_entry(tab, "Reference FASTA:", current_row_idx, tkinter.filedialog.askopenfilename)
        current_row_idx += 1
        self.gs_germline_vcf_entry = self._create_path_entry(tab, "Germline VCF:", current_row_idx, tkinter.filedialog.askopenfilename)
        current_row_idx += 1
        self.gs_output_root_entry = self._create_path_entry(tab, "Output Root Dir:", current_row_idx, tkinter.filedialog.askdirectory, is_file=False)
        current_row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Run Name (optional):").grid(row=current_row_idx, column=0, padx=10, pady=5, sticky="w")
        self.gs_run_name_entry = customtkinter.CTkEntry(master=tab)
        self.gs_run_name_entry.grid(row=current_row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2); current_row_idx += 1

        # Read simulation type choice
        customtkinter.CTkLabel(master=tab, text="Read Simulation Type:").grid(row=current_row_idx, column=0, padx=10, pady=5, sticky="w")
        self.gs_read_sim_type_var = customtkinter.StringVar(value="Short Reads")
        sim_type_menu = customtkinter.CTkOptionMenu(master=tab, variable=self.gs_read_sim_type_var,
                                                     values=["Short Reads", "Long Reads"],
                                                     command=self._on_gs_sim_type_change)
        sim_type_menu.grid(row=current_row_idx, column=1, padx=10, pady=5, sticky="w");
        self.gs_dynamic_params_start_row = current_row_idx + 1 # Save for placing conditional frames
        current_row_idx += 1


        # Frames for conditional parameters (will be gridded by _on_gs_sim_type_change)
        self.gs_short_read_params_frame = customtkinter.CTkFrame(master=tab, fg_color="transparent")
        self.gs_long_read_params_frame = customtkinter.CTkFrame(master=tab, fg_color="transparent")

        self._populate_gs_short_read_params(self.gs_short_read_params_frame)
        self._populate_gs_long_read_params(self.gs_long_read_params_frame)

        # Common options (place them after where dynamic frames will go)
        current_row_idx = self.gs_dynamic_params_start_row + 1

        self.gs_overwrite_var = customtkinter.BooleanVar(value=False)
        customtkinter.CTkCheckBox(master=tab, text="Overwrite Output", variable=self.gs_overwrite_var).grid(row=current_row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); current_row_idx += 1

        self.gs_auto_index_fasta_var = customtkinter.BooleanVar(value=True)
        customtkinter.CTkCheckBox(master=tab, text="Auto-Index FASTAs", variable=self.gs_auto_index_fasta_var).grid(row=current_row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); current_row_idx += 1

        self.gs_keep_intermediates_var = customtkinter.BooleanVar(value=False)
        customtkinter.CTkCheckBox(master=tab, text="Keep Intermediate Files (SimuG outputs)", variable=self.gs_keep_intermediates_var).grid(row=current_row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); current_row_idx += 1

        self.run_germline_sim_button = customtkinter.CTkButton(master=tab, text="Run Germline Workflow", command=self.run_germline_simulation_thread)
        self.run_germline_sim_button.grid(row=current_row_idx, column=0, columnspan=3, padx=10, pady=(20,10))

        self._on_gs_sim_type_change(self.gs_read_sim_type_var.get()) # Initial setup

    def _on_gs_sim_type_change(self, selection):
        dynamic_frame_row = self.gs_dynamic_params_start_row

        if selection == "Short Reads":
            self.gs_long_read_params_frame.grid_forget()
            self.gs_short_read_params_frame.grid(row=dynamic_frame_row, column=0, columnspan=3, sticky="ew", padx=5, pady=5)
        else: # Long Reads
            self.gs_short_read_params_frame.grid_forget()
            self.gs_long_read_params_frame.grid(row=dynamic_frame_row, column=0, columnspan=3, sticky="ew", padx=5, pady=5)

    def _populate_gs_short_read_params(self, parent_frame):
        sr_row_idx = 0
        parent_frame.grid_columnconfigure(1, weight=1)

        customtkinter.CTkLabel(master=parent_frame, text="ART Platform:").grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w")
        self.gs_sr_art_platform_var = customtkinter.StringVar(value="illumina")
        customtkinter.CTkOptionMenu(master=parent_frame, variable=self.gs_sr_art_platform_var, values=["illumina"]).grid(row=sr_row_idx, column=1, padx=5, pady=2, sticky="ew"); sr_row_idx +=1

        customtkinter.CTkLabel(master=parent_frame, text="ART Profile:").grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w")
        self.gs_sr_art_profile_var = customtkinter.StringVar(value=DEFAULT_ILLUMINA_PROFILE)
        customtkinter.CTkOptionMenu(master=parent_frame, variable=self.gs_sr_art_profile_var, values=KNOWN_ART_PROFILES["illumina"]).grid(row=sr_row_idx, column=1, padx=5, pady=2, sticky="ew"); sr_row_idx += 1

        customtkinter.CTkLabel(master=parent_frame, text="Read Length:").grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_sr_read_length_entry = customtkinter.CTkEntry(master=parent_frame); self.gs_sr_read_length_entry.insert(0, "150");
        self.gs_sr_read_length_entry.grid(row=sr_row_idx, column=1, sticky="ew", columnspan=2, padx=5, pady=2); sr_row_idx+=1

        customtkinter.CTkLabel(master=parent_frame, text="Depth:").grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_sr_depth_entry = customtkinter.CTkEntry(master=parent_frame); self.gs_sr_depth_entry.insert(0, "50");
        self.gs_sr_depth_entry.grid(row=sr_row_idx, column=1, sticky="ew", columnspan=2, padx=5, pady=2); sr_row_idx+=1

        customtkinter.CTkLabel(master=parent_frame, text="Mean Fragment:").grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_sr_mean_frag_entry = customtkinter.CTkEntry(master=parent_frame); self.gs_sr_mean_frag_entry.insert(0, "400");
        self.gs_sr_mean_frag_entry.grid(row=sr_row_idx, column=1, sticky="ew", columnspan=2, padx=5, pady=2); sr_row_idx+=1

        customtkinter.CTkLabel(master=parent_frame, text="Std Dev Fragment:").grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_sr_std_dev_frag_entry = customtkinter.CTkEntry(master=parent_frame); self.gs_sr_std_dev_frag_entry.insert(0, "50");
        self.gs_sr_std_dev_frag_entry.grid(row=sr_row_idx, column=1, sticky="ew", columnspan=2, padx=5, pady=2); sr_row_idx+=1

        self.gs_sr_is_paired_end_var = customtkinter.BooleanVar(value=True) # Default to paired-end
        customtkinter.CTkCheckBox(master=parent_frame, text="Paired-End Reads", variable=self.gs_sr_is_paired_end_var).grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w", columnspan=3); sr_row_idx += 1

    def _populate_gs_long_read_params(self, parent_frame):
        lr_row_idx = 0
        parent_frame.grid_columnconfigure(1, weight=1)

        customtkinter.CTkLabel(master=parent_frame, text="NanoSim Model:").grid(row=lr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_lr_nanosim_model_var = customtkinter.StringVar(value=DEFAULT_NANOSIM_MODEL)
        customtkinter.CTkOptionMenu(master=parent_frame, variable=self.gs_lr_nanosim_model_var, values=KNOWN_NANOSIM_MODELS).grid(row=lr_row_idx, column=1, padx=5, pady=2, sticky="ew"); lr_row_idx += 1

        customtkinter.CTkLabel(master=parent_frame, text="Depth:").grid(row=lr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_lr_depth_entry = customtkinter.CTkEntry(master=parent_frame); self.gs_lr_depth_entry.insert(0, "30");
        self.gs_lr_depth_entry.grid(row=lr_row_idx, column=1, sticky="ew", columnspan=2, padx=5, pady=2); lr_row_idx+=1

        customtkinter.CTkLabel(master=parent_frame, text="Num Reads (optional, overrides depth):").grid(row=lr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_lr_num_reads_entry = customtkinter.CTkEntry(master=parent_frame);
        self.gs_lr_num_reads_entry.grid(row=lr_row_idx, column=1, sticky="ew", columnspan=2, padx=5, pady=2); lr_row_idx+=1

    def _get_germline_sim_params(self) -> Optional[Dict[str, Any]]:
        ref_fasta = self.gs_ref_fasta_entry.get()
        germline_vcf = self.gs_germline_vcf_entry.get()
        output_root = self.gs_output_root_entry.get()

        if not all([ref_fasta, germline_vcf, output_root]):
            self.update_status("Reference FASTA, Germline VCF, and Output Root Dir are required.", is_error=True, clear_first=True)
            return None

        run_name_str = self.gs_run_name_entry.get() or None
        sim_type_selected = self.gs_read_sim_type_var.get()
        read_sim_params_dict = {}

        if sim_type_selected == "Short Reads":
            read_sim_type_for_runner = "short"
            try:
                read_sim_params_dict = {
                    "art_platform": self.gs_sr_art_platform_var.get(),
                    "art_profile": self.gs_sr_art_profile_var.get(),
                    "read_length": int(self.gs_sr_read_length_entry.get()),
                    "depth": int(self.gs_sr_depth_entry.get()),
                    "mean_fragment_length": int(self.gs_sr_mean_frag_entry.get()),
                    "std_dev_fragment_length": int(self.gs_sr_std_dev_frag_entry.get()),
                    "is_paired_end": self.gs_sr_is_paired_end_var.get()
                }
            except ValueError as e:
                self.update_status(f"Invalid short read parameter: {e}", is_error=True, clear_first=True); return None
        else: # Long Reads
            read_sim_type_for_runner = "long"
            num_reads_str = self.gs_lr_num_reads_entry.get()
            try:
                read_sim_params_dict = {
                    "model": self.gs_lr_nanosim_model_var.get(),
                    "depth": int(self.gs_lr_depth_entry.get()),
                    "num_reads": int(num_reads_str) if num_reads_str else None
                }
            except ValueError as e:
                self.update_status(f"Invalid long read parameter: {e}", is_error=True, clear_first=True); return None

        command_params_for_manifest = {
            "reference_fasta_original": ref_fasta, "germline_vcf_path": germline_vcf,
            "read_sim_type_selected": sim_type_selected,
            "read_sim_params_from_gui": copy.deepcopy(read_sim_params_dict),
            "overwrite_output": self.gs_overwrite_var.get(),
            "auto_index_fasta": self.gs_auto_index_fasta_var.get(),
            "keep_intermediate_files": self.gs_keep_intermediates_var.get(),
            "output_root_dir": output_root, "run_name": run_name_str
        }

        return {
            "output_root_dir": Path(output_root),
            "reference_fasta_path": ref_fasta,
            "germline_vcf_path": germline_vcf,
            "read_sim_type": read_sim_type_for_runner,
            "read_sim_params": read_sim_params_dict,
            "run_name": run_name_str,
            "command_params": command_params_for_manifest,
            "overwrite_output": self.gs_overwrite_var.get(),
            "auto_index_fasta": self.gs_auto_index_fasta_var.get(),
            "keep_intermediate_files": self.gs_keep_intermediates_var.get()
        }

    def run_germline_simulation_thread(self):
        if hasattr(self, 'igv_button_spike'): self.igv_button_spike.configure(state="disabled")
        self._start_runner_thread(
            button_to_disable=self.run_germline_sim_button,
            runner_func=src_bb_runner.run_germline_simulation_workflow,
            param_extractor_func=self._get_germline_sim_params,
            result_handler_func=None
        )

    def create_germline_sim_tab(self):
        tab = self.germline_sim_tab
        tab.grid_columnconfigure(1, weight=1) # Allow input column to expand
        current_row_idx = 0

        # Core inputs
        self.gs_ref_fasta_entry = self._create_path_entry(tab, "Reference FASTA:", current_row_idx, tkinter.filedialog.askopenfilename)
        current_row_idx += 1
        self.gs_germline_vcf_entry = self._create_path_entry(tab, "Germline VCF:", current_row_idx, tkinter.filedialog.askopenfilename)
        current_row_idx += 1
        self.gs_output_root_entry = self._create_path_entry(tab, "Output Root Dir:", current_row_idx, tkinter.filedialog.askdirectory, is_file=False)
        current_row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Run Name (optional):").grid(row=current_row_idx, column=0, padx=10, pady=5, sticky="w")
        self.gs_run_name_entry = customtkinter.CTkEntry(master=tab)
        self.gs_run_name_entry.grid(row=current_row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2); current_row_idx += 1

        # Read simulation type choice
        customtkinter.CTkLabel(master=tab, text="Read Simulation Type:").grid(row=current_row_idx, column=0, padx=10, pady=5, sticky="w")
        self.gs_read_sim_type_var = customtkinter.StringVar(value="Short Reads")
        sim_type_menu = customtkinter.CTkOptionMenu(master=tab, variable=self.gs_read_sim_type_var,
                                                     values=["Short Reads", "Long Reads"],
                                                     command=self._on_gs_sim_type_change)
        sim_type_menu.grid(row=current_row_idx, column=1, padx=10, pady=5, sticky="w");
        self.gs_dynamic_params_start_row = current_row_idx + 1
        current_row_idx += 1 # Increment for the dynamic frame row itself


        # Frames for conditional parameters
        self.gs_short_read_params_frame = customtkinter.CTkFrame(master=tab, fg_color="transparent")
        # self.gs_short_read_params_frame.grid_columnconfigure(1, weight=1) # Configure inside populate
        self.gs_long_read_params_frame = customtkinter.CTkFrame(master=tab, fg_color="transparent")
        # self.gs_long_read_params_frame.grid_columnconfigure(1, weight=1) # Configure inside populate

        self._populate_gs_short_read_params(self.gs_short_read_params_frame)
        self._populate_gs_long_read_params(self.gs_long_read_params_frame)

        # Common options will be placed starting from current_row_idx,
        # which is now effectively gs_dynamic_params_start_row + 1 (or more if dynamic content is tall)
        # The _on_gs_sim_type_change will grid one of the frames at gs_dynamic_params_start_row
        # So, we set current_row_idx to be after the space potentially occupied by the dynamic frame.
        # This assumes the dynamic frames don't exceed 1 effective "row" in height for layout planning.
        # A more robust way would be to have a dedicated frame for options after dynamic part.
        # For now, assuming dynamic part takes up the equivalent of one row slot for options below it.
        current_row_idx = self.gs_dynamic_params_start_row + 1


        self.gs_overwrite_var = customtkinter.BooleanVar(value=False)
        customtkinter.CTkCheckBox(master=tab, text="Overwrite Output", variable=self.gs_overwrite_var).grid(row=current_row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); current_row_idx += 1

        self.gs_auto_index_fasta_var = customtkinter.BooleanVar(value=True)
        customtkinter.CTkCheckBox(master=tab, text="Auto-Index FASTAs", variable=self.gs_auto_index_fasta_var).grid(row=current_row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); current_row_idx += 1

        self.gs_keep_intermediates_var = customtkinter.BooleanVar(value=False)
        customtkinter.CTkCheckBox(master=tab, text="Keep Intermediate Files (SimuG outputs)", variable=self.gs_keep_intermediates_var).grid(row=current_row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); current_row_idx += 1

        self.run_germline_sim_button = customtkinter.CTkButton(master=tab, text="Run Germline Workflow", command=self.run_germline_simulation_thread)
        self.run_germline_sim_button.grid(row=current_row_idx, column=0, columnspan=3, padx=10, pady=(20,5)) # Reduced bottom pady
        current_row_idx +=1

        self.igv_button_germline = customtkinter.CTkButton(master=tab, text="Open Outputs in IGV", command=self.on_igv_button_germline_click, state="disabled")
        self.igv_button_germline.grid(row=current_row_idx, column=0, columnspan=3, padx=10, pady=(5,10)) # Added pady

        self._on_gs_sim_type_change(self.gs_read_sim_type_var.get()) # Initial setup

    def _on_gs_sim_type_change(self, selection):
        dynamic_frame_row = self.gs_dynamic_params_start_row

        if selection == "Short Reads":
            self.gs_long_read_params_frame.grid_forget()
            self.gs_short_read_params_frame.grid(row=dynamic_frame_row, column=0, columnspan=3, sticky="nsew", padx=5, pady=5)
        else: # Long Reads
            self.gs_short_read_params_frame.grid_forget()
            self.gs_long_read_params_frame.grid(row=dynamic_frame_row, column=0, columnspan=3, sticky="nsew", padx=5, pady=5)

    def _populate_gs_short_read_params(self, parent_frame):
        sr_row_idx = 0
        parent_frame.grid_columnconfigure(1, weight=1)

        customtkinter.CTkLabel(master=parent_frame, text="ART Platform:").grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w")
        self.gs_sr_art_platform_var = customtkinter.StringVar(value="illumina")
        customtkinter.CTkOptionMenu(master=parent_frame, variable=self.gs_sr_art_platform_var, values=["illumina"]).grid(row=sr_row_idx, column=1, padx=5, pady=2, sticky="ew"); sr_row_idx +=1

        customtkinter.CTkLabel(master=parent_frame, text="ART Profile:").grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w")
        self.gs_sr_art_profile_var = customtkinter.StringVar(value=DEFAULT_ILLUMINA_PROFILE)
        customtkinter.CTkOptionMenu(master=parent_frame, variable=self.gs_sr_art_profile_var, values=KNOWN_ART_PROFILES["illumina"]).grid(row=sr_row_idx, column=1, padx=5, pady=2, sticky="ew"); sr_row_idx += 1

        customtkinter.CTkLabel(master=parent_frame, text="Read Length:").grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_sr_read_length_entry = customtkinter.CTkEntry(master=parent_frame); self.gs_sr_read_length_entry.insert(0, "150");
        self.gs_sr_read_length_entry.grid(row=sr_row_idx, column=1, sticky="ew", padx=5, pady=2); sr_row_idx+=1 # Removed columnspan

        customtkinter.CTkLabel(master=parent_frame, text="Depth:").grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_sr_depth_entry = customtkinter.CTkEntry(master=parent_frame); self.gs_sr_depth_entry.insert(0, "50");
        self.gs_sr_depth_entry.grid(row=sr_row_idx, column=1, sticky="ew", padx=5, pady=2); sr_row_idx+=1

        customtkinter.CTkLabel(master=parent_frame, text="Mean Fragment:").grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_sr_mean_frag_entry = customtkinter.CTkEntry(master=parent_frame); self.gs_sr_mean_frag_entry.insert(0, "400");
        self.gs_sr_mean_frag_entry.grid(row=sr_row_idx, column=1, sticky="ew", padx=5, pady=2); sr_row_idx+=1

        customtkinter.CTkLabel(master=parent_frame, text="Std Dev Fragment:").grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_sr_std_dev_frag_entry = customtkinter.CTkEntry(master=parent_frame); self.gs_sr_std_dev_frag_entry.insert(0, "50");
        self.gs_sr_std_dev_frag_entry.grid(row=sr_row_idx, column=1, sticky="ew", padx=5, pady=2); sr_row_idx+=1

        self.gs_sr_is_paired_end_var = customtkinter.BooleanVar(value=True)
        customtkinter.CTkCheckBox(master=parent_frame, text="Paired-End Reads", variable=self.gs_sr_is_paired_end_var).grid(row=sr_row_idx, column=0, padx=5, pady=2, sticky="w", columnspan=2); sr_row_idx += 1

    def _populate_gs_long_read_params(self, parent_frame):
        lr_row_idx = 0
        parent_frame.grid_columnconfigure(1, weight=1)

        customtkinter.CTkLabel(master=parent_frame, text="NanoSim Model:").grid(row=lr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_lr_nanosim_model_var = customtkinter.StringVar(value=DEFAULT_NANOSIM_MODEL)
        customtkinter.CTkOptionMenu(master=parent_frame, variable=self.gs_lr_nanosim_model_var, values=KNOWN_NANOSIM_MODELS).grid(row=lr_row_idx, column=1, padx=5, pady=2, sticky="ew"); lr_row_idx += 1

        customtkinter.CTkLabel(master=parent_frame, text="Depth:").grid(row=lr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_lr_depth_entry = customtkinter.CTkEntry(master=parent_frame); self.gs_lr_depth_entry.insert(0, "30");
        self.gs_lr_depth_entry.grid(row=lr_row_idx, column=1, sticky="ew", padx=5, pady=2); lr_row_idx+=1

        customtkinter.CTkLabel(master=parent_frame, text="Num Reads (optional, overrides depth):").grid(row=lr_row_idx, column=0, padx=5, pady=2, sticky="w");
        self.gs_lr_num_reads_entry = customtkinter.CTkEntry(master=parent_frame);
        self.gs_lr_num_reads_entry.grid(row=lr_row_idx, column=1, sticky="ew", padx=5, pady=2); lr_row_idx+=1

    def _get_germline_sim_params(self) -> Optional[Dict[str, Any]]:
        ref_fasta = self.gs_ref_fasta_entry.get()
        germline_vcf = self.gs_germline_vcf_entry.get()
        output_root = self.gs_output_root_entry.get()

        if not all([ref_fasta, germline_vcf, output_root]):
            self.update_status("Reference FASTA, Germline VCF, and Output Root Dir are required.", is_error=True, clear_first=True)
            return None

        run_name_str = self.gs_run_name_entry.get() or None
        sim_type_selected = self.gs_read_sim_type_var.get()
        read_sim_params_dict = {}

        if sim_type_selected == "Short Reads":
            read_sim_type_for_runner = "short"
            try:
                read_sim_params_dict = {
                    "art_platform": self.gs_sr_art_platform_var.get(),
                    "art_profile": self.gs_sr_art_profile_var.get(),
                    "read_length": int(self.gs_sr_read_length_entry.get()),
                    "depth": int(self.gs_sr_depth_entry.get()),
                    "mean_fragment_length": int(self.gs_sr_mean_frag_entry.get()),
                    "std_dev_fragment_length": int(self.gs_sr_std_dev_frag_entry.get()),
                    "is_paired_end": self.gs_sr_is_paired_end_var.get()
                }
            except ValueError as e:
                self.update_status(f"Invalid short read parameter: {e}", is_error=True, clear_first=True); return None
        else: # Long Reads
            read_sim_type_for_runner = "long"
            num_reads_str = self.gs_lr_num_reads_entry.get()
            try:
                read_sim_params_dict = {
                    "model": self.gs_lr_nanosim_model_var.get(),
                    "depth": int(self.gs_lr_depth_entry.get()),
                    "num_reads": int(num_reads_str) if num_reads_str else None
                }
            except ValueError as e:
                self.update_status(f"Invalid long read parameter: {e}", is_error=True, clear_first=True); return None

        command_params_for_manifest = {
            "reference_fasta_original": ref_fasta, "germline_vcf_path": germline_vcf,
            "read_sim_type_selected": sim_type_selected,
            "read_sim_params_from_gui": copy.deepcopy(read_sim_params_dict),
            "overwrite_output": self.gs_overwrite_var.get(),
            "auto_index_fasta": self.gs_auto_index_fasta_var.get(),
            "keep_intermediate_files": self.gs_keep_intermediates_var.get(),
            "output_root_dir": output_root, "run_name": run_name_str
        }

        return {
            "output_root_dir": Path(output_root),
            "reference_fasta_path": ref_fasta,
            "germline_vcf_path": germline_vcf,
            "read_sim_type": read_sim_type_for_runner,
            "read_sim_params": read_sim_params_dict,
            "run_name": run_name_str,
            "command_params": command_params_for_manifest,
            "overwrite_output": self.gs_overwrite_var.get(),
            "auto_index_fasta": self.gs_auto_index_fasta_var.get(),
            "keep_intermediate_files": self.gs_keep_intermediates_var.get()
        }

    def run_germline_simulation_thread(self):
        if hasattr(self, 'igv_button_spike'): self.igv_button_spike.configure(state="disabled")
        # Disable other IGV buttons if they exist
        if hasattr(self, 'igv_button_germline'): self.igv_button_germline.configure(state="disabled")

        self._start_runner_thread(
            button_to_disable=self.run_germline_sim_button,
            runner_func=src_bb_runner.run_germline_simulation_workflow,
            param_extractor_func=self._get_germline_sim_params,
            result_handler_func=None # No specific handler for now, poll_thread handles generic completion
        )

    def on_igv_button_germline_click(self):
        if not hasattr(self, 'germline_sim_results') or not self.germline_sim_results: # Check instance variable
            self.update_status("No germline simulation results available.", is_error=True, clear_first=True)
            return

        ref_fasta = self.germline_sim_results.get("final_modified_fasta_path")
        sim_outputs = self.germline_sim_results.get("simulated_reads_results", {})

        data_files_to_load = []
        # output_directory in sim_outputs is the sub-run directory (e.g., .../germline_sim_run_XXXX/germline_sim_run_XXXX_short_reads)
        sub_run_output_dir = sim_outputs.get("output_directory")

        if sub_run_output_dir and "output_files" in sim_outputs:
            for f_info in sim_outputs["output_files"]:
                # f_info["path"] should be relative to sub_run_output_dir
                abs_path = Path(sub_run_output_dir) / f_info["path"]
                if abs_path.exists():
                    data_files_to_load.append(str(abs_path))
                else:
                    self.update_status(f"Warning: Data file for IGV not found: {abs_path}", is_error=True)

        if not ref_fasta or not Path(ref_fasta).exists():
            self.update_status("Error: Modified reference FASTA from germline simulation not found.", is_error=True); return

        # Check if there are any BAM files specifically, as IGV needs alignments
        bam_files_found = any(f.lower().endswith(".bam") for f in data_files_to_load)
        if not bam_files_found:
             self.update_status("No BAM files found in germline simulation results to load into IGV.", is_error=True); return

        self.update_status(f"Preparing to load into IGV: Ref: {Path(ref_fasta).name}, Data: {[Path(f).name for f in data_files_to_load]}", clear_first=True)
        self.send_to_igv_desktop(reference_fasta_path=ref_fasta, data_file_paths=data_files_to_load)

    def on_open_output_dir_click(self):
        if hasattr(self, 'last_successful_output_dir') and self.last_successful_output_dir:
            output_path = Path(self.last_successful_output_dir)
            if output_path.exists() and output_path.is_dir():
                self.update_status(f"Attempting to open output directory: {output_path}")
                try:
                    if sys.platform == "win32":
                        os.startfile(str(output_path))
                    elif sys.platform == "darwin": # macOS
                        subprocess.Popen(["open", str(output_path)])
                    else: # linux variants
                        subprocess.Popen(["xdg-open", str(output_path)])
                except FileNotFoundError:
                    self.update_status(f"Error: Could not find a command to open the directory. Please navigate manually: {output_path}", is_error=True)
                except Exception as e:
                    self.update_status(f"Error opening directory {output_path}: {e}", is_error=True)
            else:
                self.update_status(f"Error: Output directory not found or is not a directory: {output_path}", is_error=True)
        else:
            self.update_status("No output directory from a successful run is available to open.", is_error=True)

if __name__ == "__main__":
    run_gui()

[end of src/basebuddy/gui/main_app.py]
