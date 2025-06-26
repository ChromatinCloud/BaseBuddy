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
import traceback # For error handling
import logging # For logger
import json # For saving/loading settings

logger = logging.getLogger(__name__)

# Assuming runner and utils are accessible via PYTHONPATH
from basebuddy import runner as src_bb_runner
# Ensure signature_utils is imported if apply_signature_to_fasta needs it for path resolution
from basebuddy import signature_utils as sig_utils
from basebuddy import utils as bb_utils
from basebuddy.utils import BaseBuddyInputError, BaseBuddyToolError, BaseBuddyFileError, BaseBuddyConfigError

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

# Settings directory location
SETTINGS_DIR = Path("refs/gui_settings")

# Default reference paths for common genome builds
DEFAULT_REFERENCES = {
    "hg19/GRCh37": [
        "refs/GRCh37/hs37d5.fa",
        "refs/GRCh37/human_g1k_v37.fasta",
        "refs/GRCh37/GRCh37.fa"
    ],
    "hg38/GRCh38": [
        "refs/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "refs/GRCh38/GRCh38.fa",
        "refs/GRCh38/hg38.fa"
    ],
    "T2T-CHM13": [
        "refs/T2T/chm13v2.0.fa",
        "refs/T2T/T2T-CHM13.fa"
    ],
    "mm10": [
        "refs/mm10/mm10.fa",
        "refs/mm10/GRCm38.fa"
    ],
    "mm39": [
        "refs/mm39/mm39.fa",
        "refs/mm39/GRCm39.fa"
    ]
}

# Example genomic ranges for each build
EXAMPLE_RANGES = {
    "hg19/GRCh37": "7:140453136-140453236\nMT:1-16569",  # BRAF V600E in GRCh37 (no 'chr' prefix)
    "hg38/GRCh38": "chr7:140753336-140753436\nchrM:1-16569",  # BRAF V600E in GRCh38
    "T2T-CHM13": "chr7:140753336-140753436\nchrM:1-16569",
    "mm10": "chr7:45000000-45001000\nchrM:1-16299",
    "mm39": "chr7:45000000-45001000\nchrM:1-16299",
    "": "chr1:1000-2000\nchrM:1-16569"  # Default
}


class BaseBuddyGUI(customtkinter.CTk):
    def __init__(self):
        super().__init__()
        
        # Load saved settings
        self.settings = self.load_settings()

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

        # Create a frame to hold the status label and copy button
        status_frame = customtkinter.CTkFrame(self, fg_color="transparent")
        status_frame.pack(padx=20, pady=(5,0), fill="x")
        
        self.status_label = customtkinter.CTkLabel(status_frame, text="Status:")
        self.status_label.pack(side="left", anchor="w")
        
        self.copy_status_button = customtkinter.CTkButton(status_frame, text="Copy Status", command=self.copy_status_to_clipboard, width=100)
        self.copy_status_button.pack(side="left", padx=(10, 0))
        
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
        
        # Apply loaded settings to all tabs
        self.apply_settings()

        self.poll_thread()
        
        # Save settings when window is closed
        self.protocol("WM_DELETE_WINDOW", self.on_closing)
    
    def load_settings(self) -> Dict[str, Any]:
        """Load saved GUI settings from individual tab JSON files"""
        settings = {}
        
        # Ensure settings directory exists
        SETTINGS_DIR.mkdir(parents=True, exist_ok=True)
        
        # Map of tab names to their settings files
        tab_files = {
            "short_read": "short_read_sim.settings.json",
            "spike_variants": "variant_spiking.settings.json", 
            "long_read": "long_read_sim.settings.json",
            "apply_signature": "apply_signature.settings.json",
            "germline_sim": "germline_simulation.settings.json"
        }
        
        # Load each tab's settings
        for tab_name, filename in tab_files.items():
            settings_file = SETTINGS_DIR / filename
            try:
                if settings_file.exists():
                    with open(settings_file, 'r') as f:
                        settings[tab_name] = json.load(f)
                        logger.info(f"Loaded settings for {tab_name} from {filename}")
            except Exception as e:
                logger.warning(f"Failed to load settings for {tab_name}: {e}")
                
        return settings
    
    def save_settings(self):
        """Save current GUI settings to individual tab JSON files"""
        try:
            # Ensure directory exists
            SETTINGS_DIR.mkdir(parents=True, exist_ok=True)
            
            # Map of tab names to their settings files and getter methods
            tab_configs = {
                "short_read": ("short_read_sim.settings.json", self.get_short_read_settings),
                "spike_variants": ("variant_spiking.settings.json", self.get_spike_variants_settings),
                "long_read": ("long_read_sim.settings.json", self.get_long_read_settings),
                "apply_signature": ("apply_signature.settings.json", self.get_apply_signature_settings),
                "germline_sim": ("germline_simulation.settings.json", self.get_germline_sim_settings)
            }
            
            # Save each tab's settings to its own file
            for tab_name, (filename, getter_method) in tab_configs.items():
                settings_file = SETTINGS_DIR / filename
                try:
                    tab_settings = getter_method()
                    with open(settings_file, 'w') as f:
                        json.dump(tab_settings, f, indent=2)
                    logger.info(f"Saved settings for {tab_name} to {filename}")
                except Exception as e:
                    logger.error(f"Failed to save settings for {tab_name}: {e}")
                
        except Exception as e:
            logger.error(f"Failed to save settings: {e}")
    
    def on_closing(self):
        """Handle window close event"""
        self.save_settings()
        self.destroy()
    
    # Settings getter methods for each tab
    def get_short_read_settings(self) -> Dict[str, Any]:
        """Get current settings from short read tab"""
        return {
            "reference_fasta": self.sr_ref_fasta_entry.get(),
            "output_root": self.sr_output_root_entry.get(),
            "run_name": self.sr_run_name_entry.get(),
            "art_platform": self.sr_art_platform_var.get(),
            "art_profile": self.sr_art_profile_var.get(),
            "genome_build": self.sr_genome_build_var.get(),
            "custom_build": self.sr_custom_build_entry.get() if hasattr(self, "sr_custom_build_entry") else "",
            "read_length": self.sr_read_length_entry.get(),
            "depth": self.sr_depth_entry.get(),
            "mean_fragment": self.sr_mean_frag_entry.get(),
            "std_dev_fragment": self.sr_std_dev_frag_entry.get(),
            "genomic_ranges": self.short_sim_ranges_textbox.get("1.0", "end-1c"),
            "timeout": self.sr_timeout_entry.get(),
            "single_end": self.sr_single_end_var.get(),
            "overwrite": self.sr_overwrite_var.get(),
            "auto_align": self.sr_auto_align_var.get() if hasattr(self, "sr_auto_align_var") else False,
            "aligner": self.sr_aligner_var.get() if hasattr(self, "sr_aligner_var") else "bwa"
        }
    
    def get_spike_variants_settings(self) -> Dict[str, Any]:
        """Get current settings from spike variants tab"""
        return {
            "reference_fasta": self.sv_ref_fasta_entry.get() if hasattr(self, "sv_ref_fasta_entry") else "",
            "input_bam": self.sv_input_bam_entry.get() if hasattr(self, "sv_input_bam_entry") else "",
            "manual_variants": self.sv_manual_variants_textbox.get("1.0", "end-1c") if hasattr(self, "sv_manual_variants_textbox") else "",
            "vcf_file": self.sv_vcf_entry.get() if hasattr(self, "sv_vcf_entry") else "",
            "picard_jar": self.sv_picard_jar_entry.get() if hasattr(self, "sv_picard_jar_entry") else "",
            "vaf": self.sv_vaf_entry.get() if hasattr(self, "sv_vaf_entry") else "0.05",
            "seed": self.sv_seed_entry.get() if hasattr(self, "sv_seed_entry") else "0",
            "strand_bias": self.sv_strand_bias_entry.get() if hasattr(self, "sv_strand_bias_entry") else "0.5",
            "output_root": self.sv_output_root_entry.get() if hasattr(self, "sv_output_root_entry") else "",
            "run_name": self.sv_run_name_entry.get() if hasattr(self, "sv_run_name_entry") else "",
            "output_prefix": self.sv_output_prefix_entry.get() if hasattr(self, "sv_output_prefix_entry") else "spiked_output",
            "overwrite": self.sv_overwrite_var.get() if hasattr(self, "sv_overwrite_var") else False,
            "auto_index_bam": self.sv_auto_index_input_bam_var.get() if hasattr(self, "sv_auto_index_input_bam_var") else True,
            "auto_index_fasta": self.sv_auto_index_fasta_var.get() if hasattr(self, "sv_auto_index_fasta_var") else True
        }
    
    def get_long_read_settings(self) -> Dict[str, Any]:
        """Get current settings from long read tab"""
        return {
            "reference_fasta": self.lr_ref_fasta_entry.get(),
            "output_root": self.lr_output_root_entry.get(),
            "run_name": self.lr_run_name_entry.get(),
            "nanosim_model": self.lr_nanosim_model_var.get(),
            "depth": self.lr_depth_entry.get(),
            "num_reads": self.lr_num_reads_entry.get(),
            "perfect_reads": self.lr_perfect_reads_var.get(),
            "kmer_bias": self.lr_kmer_bias_var.get(),
            "overwrite": self.lr_overwrite_var.get()
        }
    
    def get_apply_signature_settings(self) -> Dict[str, Any]:
        """Get current settings from apply signature tab"""
        return {
            "reference_fasta": self.apply_sig_ref_fasta_entry.get() if hasattr(self, "apply_sig_ref_fasta_entry") else "",
            "output_fasta": self.apply_sig_output_fasta_entry.get() if hasattr(self, "apply_sig_output_fasta_entry") else "",
            "signature_id": self.apply_sig_id_var.get() if hasattr(self, "apply_sig_id_var") else "SBS1",
            "custom_path": self.apply_sig_custom_path_entry.get() if hasattr(self, "apply_sig_custom_path_entry") else "",
            "num_mutations": self.apply_sig_num_mutations_entry.get() if hasattr(self, "apply_sig_num_mutations_entry") else "100",
            "random_seed": self.apply_sig_random_seed_entry.get() if hasattr(self, "apply_sig_random_seed_entry") else "",
            "genomic_ranges": self.apply_sig_ranges_textbox.get("1.0", "end-1c") if hasattr(self, "apply_sig_ranges_textbox") else ""
        }
    
    def get_germline_sim_settings(self) -> Dict[str, Any]:
        """Get current settings from germline simulation tab"""
        settings = {
            "reference_fasta": self.gs_ref_fasta_entry.get() if hasattr(self, "gs_ref_fasta_entry") else "",
            "germline_vcf": self.gs_germline_vcf_entry.get() if hasattr(self, "gs_germline_vcf_entry") else "",
            "output_root": self.gs_output_root_entry.get() if hasattr(self, "gs_output_root_entry") else "",
            "run_name": self.gs_run_name_entry.get() if hasattr(self, "gs_run_name_entry") else "",
            "read_sim_type": self.gs_read_sim_type_var.get() if hasattr(self, "gs_read_sim_type_var") else "Short Reads",
            "overwrite": self.gs_overwrite_var.get() if hasattr(self, "gs_overwrite_var") else False,
            "auto_index_fasta": self.gs_auto_index_fasta_var.get() if hasattr(self, "gs_auto_index_fasta_var") else True,
            "keep_intermediates": self.gs_keep_intermediates_var.get() if hasattr(self, "gs_keep_intermediates_var") else False
        }
        
        # Add short read specific settings
        if hasattr(self, "gs_sr_art_platform_var"):
            settings["sr_art_platform"] = self.gs_sr_art_platform_var.get()
        if hasattr(self, "gs_sr_art_profile_var"):
            settings["sr_art_profile"] = self.gs_sr_art_profile_var.get()
        if hasattr(self, "gs_sr_read_length_entry"):
            settings["sr_read_length"] = self.gs_sr_read_length_entry.get()
        if hasattr(self, "gs_sr_depth_entry"):
            settings["sr_depth"] = self.gs_sr_depth_entry.get()
        if hasattr(self, "gs_sr_mean_frag_entry"):
            settings["sr_mean_fragment"] = self.gs_sr_mean_frag_entry.get()
        if hasattr(self, "gs_sr_std_dev_frag_entry"):
            settings["sr_std_dev_fragment"] = self.gs_sr_std_dev_frag_entry.get()
        if hasattr(self, "gs_sr_is_paired_end_var"):
            settings["sr_is_paired_end"] = self.gs_sr_is_paired_end_var.get()
            
        # Add long read specific settings
        if hasattr(self, "gs_lr_nanosim_model_var"):
            settings["lr_nanosim_model"] = self.gs_lr_nanosim_model_var.get()
        if hasattr(self, "gs_lr_depth_entry"):
            settings["lr_depth"] = self.gs_lr_depth_entry.get()
        if hasattr(self, "gs_lr_num_reads_entry"):
            settings["lr_num_reads"] = self.gs_lr_num_reads_entry.get()
            
        return settings
    
    # Settings application methods
    def apply_settings(self):
        """Apply loaded settings to all tabs"""
        if hasattr(self, 'settings') and self.settings:
            if 'short_read' in self.settings:
                self.apply_short_read_settings(self.settings['short_read'])
            if 'spike_variants' in self.settings:
                self.apply_spike_variants_settings(self.settings['spike_variants'])
            if 'long_read' in self.settings:
                self.apply_long_read_settings(self.settings['long_read'])
            if 'apply_signature' in self.settings:
                self.apply_signature_settings(self.settings['apply_signature'])
            if 'germline_sim' in self.settings:
                self.apply_germline_sim_settings(self.settings['germline_sim'])
    
    def apply_short_read_settings(self, settings: Dict[str, Any]):
        """Apply saved settings to short read tab"""
        if 'reference_fasta' in settings and hasattr(self, 'sr_ref_fasta_entry'):
            self.sr_ref_fasta_entry.delete(0, 'end')
            self.sr_ref_fasta_entry.insert(0, settings['reference_fasta'])
        if 'output_root' in settings and hasattr(self, 'sr_output_root_entry'):
            self.sr_output_root_entry.delete(0, 'end')
            self.sr_output_root_entry.insert(0, settings['output_root'])
        if 'run_name' in settings and hasattr(self, 'sr_run_name_entry'):
            self.sr_run_name_entry.delete(0, 'end')
            self.sr_run_name_entry.insert(0, settings['run_name'])
        if 'art_platform' in settings and hasattr(self, 'sr_art_platform_var'):
            self.sr_art_platform_var.set(settings['art_platform'])
        if 'art_profile' in settings and hasattr(self, 'sr_art_profile_var'):
            self.sr_art_profile_var.set(settings['art_profile'])
        if 'genome_build' in settings and hasattr(self, 'sr_genome_build_var'):
            self.sr_genome_build_var.set(settings['genome_build'])
        if 'read_length' in settings and hasattr(self, 'sr_read_length_entry'):
            self.sr_read_length_entry.delete(0, 'end')
            self.sr_read_length_entry.insert(0, settings['read_length'])
        if 'depth' in settings and hasattr(self, 'sr_depth_entry'):
            self.sr_depth_entry.delete(0, 'end')
            self.sr_depth_entry.insert(0, settings['depth'])
        if 'mean_fragment' in settings and hasattr(self, 'sr_mean_frag_entry'):
            self.sr_mean_frag_entry.delete(0, 'end')
            self.sr_mean_frag_entry.insert(0, settings['mean_fragment'])
        if 'std_dev_fragment' in settings and hasattr(self, 'sr_std_dev_frag_entry'):
            self.sr_std_dev_frag_entry.delete(0, 'end')
            self.sr_std_dev_frag_entry.insert(0, settings['std_dev_fragment'])
        if 'is_paired_end' in settings and hasattr(self, 'sr_is_paired_end_var'):
            self.sr_is_paired_end_var.set(settings['is_paired_end'])
        if 'overwrite' in settings and hasattr(self, 'sr_overwrite_var'):
            self.sr_overwrite_var.set(settings['overwrite'])
        if 'auto_index' in settings and hasattr(self, 'sr_auto_index_var'):
            self.sr_auto_index_var.set(settings['auto_index'])
        if 'auto_align' in settings and hasattr(self, 'sr_auto_align_var'):
            self.sr_auto_align_var.set(settings['auto_align'])
        if 'genomic_ranges' in settings and hasattr(self, 'sr_genomic_ranges_text'):
            self.sr_genomic_ranges_text.delete("1.0", "end")
            self.sr_genomic_ranges_text.insert("1.0", settings['genomic_ranges'])
        if 'aligner' in settings and hasattr(self, 'sr_aligner_var'):
            self.sr_aligner_var.set(settings['aligner'])
        if 'output_bam' in settings and hasattr(self, 'sr_output_bam_entry'):
            self.sr_output_bam_entry.delete(0, 'end')
            self.sr_output_bam_entry.insert(0, settings['output_bam'])
    
    def apply_spike_variants_settings(self, settings: Dict[str, Any]):
        """Apply saved settings to spike variants tab"""
        if 'reference_fasta' in settings and hasattr(self, 'sv_ref_fasta_entry'):
            self.sv_ref_fasta_entry.delete(0, 'end')
            self.sv_ref_fasta_entry.insert(0, settings['reference_fasta'])
        if 'input_bam' in settings and hasattr(self, 'sv_input_bam_entry'):
            self.sv_input_bam_entry.delete(0, 'end')
            self.sv_input_bam_entry.insert(0, settings['input_bam'])
        if 'manual_variants' in settings and hasattr(self, 'sv_manual_variants_textbox'):
            self.sv_manual_variants_textbox.delete("1.0", "end")
            self.sv_manual_variants_textbox.insert("1.0", settings['manual_variants'])
        if 'vcf_file' in settings and hasattr(self, 'sv_vcf_entry'):
            self.sv_vcf_entry.delete(0, 'end')
            self.sv_vcf_entry.insert(0, settings['vcf_file'])
        if 'picard_jar' in settings and hasattr(self, 'sv_picard_jar_entry'):
            self.sv_picard_jar_entry.delete(0, 'end')
            self.sv_picard_jar_entry.insert(0, settings['picard_jar'])
        if 'vaf' in settings and hasattr(self, 'sv_vaf_entry'):
            self.sv_vaf_entry.delete(0, 'end')
            self.sv_vaf_entry.insert(0, settings['vaf'])
        if 'seed' in settings and hasattr(self, 'sv_seed_entry'):
            self.sv_seed_entry.delete(0, 'end')
            self.sv_seed_entry.insert(0, settings['seed'])
        if 'strand_bias' in settings and hasattr(self, 'sv_strand_bias_entry'):
            self.sv_strand_bias_entry.delete(0, 'end')
            self.sv_strand_bias_entry.insert(0, settings['strand_bias'])
        if 'output_root' in settings and hasattr(self, 'sv_output_root_entry'):
            self.sv_output_root_entry.delete(0, 'end')
            self.sv_output_root_entry.insert(0, settings['output_root'])
        if 'run_name' in settings and hasattr(self, 'sv_run_name_entry'):
            self.sv_run_name_entry.delete(0, 'end')
            self.sv_run_name_entry.insert(0, settings['run_name'])
        if 'output_prefix' in settings and hasattr(self, 'sv_output_prefix_entry'):
            self.sv_output_prefix_entry.delete(0, 'end')
            self.sv_output_prefix_entry.insert(0, settings['output_prefix'])
        if 'overwrite' in settings and hasattr(self, 'sv_overwrite_var'):
            self.sv_overwrite_var.set(settings['overwrite'])
        if 'auto_index_bam' in settings and hasattr(self, 'sv_auto_index_input_bam_var'):
            self.sv_auto_index_input_bam_var.set(settings['auto_index_bam'])
        if 'auto_index_fasta' in settings and hasattr(self, 'sv_auto_index_fasta_var'):
            self.sv_auto_index_fasta_var.set(settings['auto_index_fasta'])
    
    def apply_long_read_settings(self, settings: Dict[str, Any]):
        """Apply saved settings to long read tab"""
        if 'reference_fasta' in settings and hasattr(self, 'lr_ref_fasta_entry'):
            self.lr_ref_fasta_entry.delete(0, 'end')
            self.lr_ref_fasta_entry.insert(0, settings['reference_fasta'])
        if 'output_root' in settings and hasattr(self, 'lr_output_root_entry'):
            self.lr_output_root_entry.delete(0, 'end')
            self.lr_output_root_entry.insert(0, settings['output_root'])
        if 'run_name' in settings and hasattr(self, 'lr_run_name_entry'):
            self.lr_run_name_entry.delete(0, 'end')
            self.lr_run_name_entry.insert(0, settings['run_name'])
        if 'nanosim_model' in settings and hasattr(self, 'lr_nanosim_model_var'):
            self.lr_nanosim_model_var.set(settings['nanosim_model'])
        if 'depth' in settings and hasattr(self, 'lr_depth_entry'):
            self.lr_depth_entry.delete(0, 'end')
            self.lr_depth_entry.insert(0, settings['depth'])
        if 'num_reads' in settings and hasattr(self, 'lr_num_reads_entry'):
            self.lr_num_reads_entry.delete(0, 'end')
            self.lr_num_reads_entry.insert(0, settings['num_reads'])
        if 'overwrite' in settings and hasattr(self, 'lr_overwrite_var'):
            self.lr_overwrite_var.set(settings['overwrite'])
        if 'auto_index' in settings and hasattr(self, 'lr_auto_index_var'):
            self.lr_auto_index_var.set(settings['auto_index'])
    
    def apply_signature_settings(self, settings: Dict[str, Any]):
        """Apply saved settings to apply signature tab"""
        if 'reference_fasta' in settings and hasattr(self, 'apply_sig_ref_fasta_entry'):
            self.apply_sig_ref_fasta_entry.delete(0, 'end')
            self.apply_sig_ref_fasta_entry.insert(0, settings['reference_fasta'])
        if 'output_fasta' in settings and hasattr(self, 'apply_sig_output_fasta_entry'):
            self.apply_sig_output_fasta_entry.delete(0, 'end')
            self.apply_sig_output_fasta_entry.insert(0, settings['output_fasta'])
        if 'num_mutations' in settings and hasattr(self, 'apply_sig_num_mutations_entry'):
            self.apply_sig_num_mutations_entry.delete(0, 'end')
            self.apply_sig_num_mutations_entry.insert(0, settings['num_mutations'])
        if 'signature_id' in settings and hasattr(self, 'apply_sig_id_var'):
            self.apply_sig_id_var.set(settings['signature_id'])
        if 'custom_path' in settings and hasattr(self, 'apply_sig_custom_path_entry'):
            self.apply_sig_custom_path_entry.delete(0, 'end')
            self.apply_sig_custom_path_entry.insert(0, settings['custom_path'])
        if 'random_seed' in settings and hasattr(self, 'apply_sig_random_seed_entry'):
            self.apply_sig_random_seed_entry.delete(0, 'end')
            self.apply_sig_random_seed_entry.insert(0, settings['random_seed'])
        if 'genomic_ranges' in settings and hasattr(self, 'apply_sig_ranges_textbox'):
            self.apply_sig_ranges_textbox.delete("1.0", "end")
            self.apply_sig_ranges_textbox.insert("1.0", settings['genomic_ranges'])
    
    def apply_germline_sim_settings(self, settings: Dict[str, Any]):
        """Apply saved settings to germline simulation tab"""
        if 'reference_fasta' in settings and hasattr(self, 'gs_ref_fasta_entry'):
            self.gs_ref_fasta_entry.delete(0, 'end')
            self.gs_ref_fasta_entry.insert(0, settings['reference_fasta'])
        if 'germline_vcf' in settings and hasattr(self, 'gs_germline_vcf_entry'):
            self.gs_germline_vcf_entry.delete(0, 'end')
            self.gs_germline_vcf_entry.insert(0, settings['germline_vcf'])
        if 'output_root' in settings and hasattr(self, 'gs_output_root_entry'):
            self.gs_output_root_entry.delete(0, 'end')
            self.gs_output_root_entry.insert(0, settings['output_root'])
        if 'run_name' in settings and hasattr(self, 'gs_run_name_entry'):
            self.gs_run_name_entry.delete(0, 'end')
            self.gs_run_name_entry.insert(0, settings['run_name'])
        if 'read_sim_type' in settings and hasattr(self, 'gs_read_sim_type_var'):
            self.gs_read_sim_type_var.set(settings['read_sim_type'])
            self._on_gs_sim_type_change(settings['read_sim_type'])
        # Apply short read specific settings if short reads selected
        if settings.get('read_sim_type') == 'Short Reads':
            if 'sr_art_platform' in settings and hasattr(self, 'gs_sr_art_platform_var'):
                self.gs_sr_art_platform_var.set(settings['sr_art_platform'])
            if 'sr_art_profile' in settings and hasattr(self, 'gs_sr_art_profile_var'):
                self.gs_sr_art_profile_var.set(settings['sr_art_profile'])
            if 'sr_read_length' in settings and hasattr(self, 'gs_sr_read_length_entry'):
                self.gs_sr_read_length_entry.delete(0, 'end')
                self.gs_sr_read_length_entry.insert(0, settings['sr_read_length'])
            if 'sr_depth' in settings and hasattr(self, 'gs_sr_depth_entry'):
                self.gs_sr_depth_entry.delete(0, 'end')
                self.gs_sr_depth_entry.insert(0, settings['sr_depth'])
            if 'sr_mean_fragment' in settings and hasattr(self, 'gs_sr_mean_frag_entry'):
                self.gs_sr_mean_frag_entry.delete(0, 'end')
                self.gs_sr_mean_frag_entry.insert(0, settings['sr_mean_fragment'])
            if 'sr_std_dev_fragment' in settings and hasattr(self, 'gs_sr_std_dev_frag_entry'):
                self.gs_sr_std_dev_frag_entry.delete(0, 'end')
                self.gs_sr_std_dev_frag_entry.insert(0, settings['sr_std_dev_fragment'])
            if 'sr_is_paired_end' in settings and hasattr(self, 'gs_sr_is_paired_end_var'):
                self.gs_sr_is_paired_end_var.set(settings['sr_is_paired_end'])
        # Apply long read specific settings if long reads selected
        elif settings.get('read_sim_type') == 'Long Reads':
            if 'lr_nanosim_model' in settings and hasattr(self, 'gs_lr_nanosim_model_var'):
                self.gs_lr_nanosim_model_var.set(settings['lr_nanosim_model'])
            if 'lr_depth' in settings and hasattr(self, 'gs_lr_depth_entry'):
                self.gs_lr_depth_entry.delete(0, 'end')
                self.gs_lr_depth_entry.insert(0, settings['lr_depth'])
            if 'lr_num_reads' in settings and hasattr(self, 'gs_lr_num_reads_entry'):
                self.gs_lr_num_reads_entry.delete(0, 'end')
                self.gs_lr_num_reads_entry.insert(0, settings['lr_num_reads'])
        if 'overwrite' in settings and hasattr(self, 'gs_overwrite_var'):
            self.gs_overwrite_var.set(settings['overwrite'])
        if 'auto_index_fasta' in settings and hasattr(self, 'gs_auto_index_fasta_var'):
            self.gs_auto_index_fasta_var.set(settings['auto_index_fasta'])
        if 'keep_intermediates' in settings and hasattr(self, 'gs_keep_intermediates_var'):
            self.gs_keep_intermediates_var.set(settings['keep_intermediates'])

    def update_status(self, message: str, is_error: bool = False, clear_first: bool = False):
        self.status_textbox.configure(state="normal")
        if clear_first:
            self.status_textbox.delete("1.0", "end")
        self.status_textbox.insert("end", message + "\n")
        self.status_textbox.configure(state="disabled")
        self.status_textbox.see("end")
    
    def copy_status_to_clipboard(self):
        """Copy the contents of the status textbox to clipboard"""
        status_content = self.status_textbox.get("1.0", "end-1c")
        if status_content:
            self.clipboard_clear()
            self.clipboard_append(status_content)
            # Briefly change button text to indicate success
            self.copy_status_button.configure(text="Copied!")
            self.after(1500, lambda: self.copy_status_button.configure(text="Copy Status"))


    def poll_thread(self):
        try:
            message_type, data, button_to_enable = self.thread_queue.get_nowait()
            if message_type == "completion":
                logger.info(f"Received completion message in poll_thread")
                result_dict, success = data

                if success:
                    if button_to_enable == self.run_spike_button:
                        self.spike_variants_results = result_dict # Store results
                        if self.spike_variants_results.get("igv_session_file"):
                            self.igv_button_spike.configure(state="normal")
                    elif button_to_enable == self.run_germline_sim_button:
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
                
                # Enhanced success message
                if success:
                    msg = f"SIMULATION SUCCESSFUL\n"
                    msg += f"{'='*50}\n"
                    msg += f"Run Name: {run_name}\n"
                    msg += f"Output Directory: {output_dir}\n"
                    msg += f"{'='*50}\n"
                else:
                    msg = f"Run '{run_name}' completed with issues.\n"
                    msg += f"Output Directory: {output_dir}\n"

                if "output_files" in result_dict and result_dict["output_files"]: # Check if list is not empty
                    msg += "\nGenerated Files:\n"
                    # Group files by type
                    fastq_files = [f for f in result_dict["output_files"] if f["type"] == "FASTQ"]
                    bam_files = [f for f in result_dict["output_files"] if f["type"] in ["BAM", "BAI"]]
                    other_files = [f for f in result_dict["output_files"] if f["type"] not in ["FASTQ", "BAM", "BAI"]]
                    
                    if fastq_files:
                        msg += "  FASTQ Files:\n"
                        for f_info in fastq_files:
                            msg += f"    - {f_info['path']}\n"
                    
                    if bam_files:
                        msg += "  Alignment Files:\n"
                        for f_info in bam_files:
                            msg += f"    - {f_info['path']} ({f_info['type']})\n"
                    
                    if other_files:
                        msg += "  Other Files:\n"
                        for f_info in other_files:
                            msg += f"    - {f_info['path']} ({f_info['type']})\n"

                # Specific outputs for apply_signature_to_fasta
                if "output_modified_fasta_path" in result_dict:
                    msg += f"Mutated FASTA: {result_dict['output_modified_fasta_path']}\n"
                if result_dict.get("num_mutations_applied") is not None:
                     msg += f"Mutations Applied: {result_dict['num_mutations_applied']}\n"

                # Add summary statistics if available
                if success:
                    msg += f"\n{'='*50}\n"
                    msg += "TIP: Click 'Open Output Directory' to view all files\n"
                    
                if "manifest_path" in result_dict: 
                    msg += f"\nManifest: {result_dict['manifest_path']}\n"
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
            # Use the provided dialog_func if it's callable
            if callable(dialog_func):
                path = dialog_func()
            else:
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
        # Configure columns for side-by-side layout
        tab.grid_columnconfigure(1, weight=1)
        tab.grid_columnconfigure(3, weight=1)
        row_idx = 0

        # Create custom callback for reference FASTA that auto-detects genome build
        def browse_reference():
            path = tkinter.filedialog.askopenfilename(
                title="Select Reference FASTA",
                filetypes=[("FASTA files", "*.fa *.fasta *.fna"), ("All files", "*.*")]
            )
            if path:
                self.sr_ref_fasta_entry.delete(0, "end")
                self.sr_ref_fasta_entry.insert(0, path)
                # Auto-detect genome build from filename
                self._auto_detect_genome_build(path)
            return path
        
        self.sr_ref_fasta_entry = self._create_path_entry(tab, "Reference FASTA:", row_idx, browse_reference)
        if saved_ref := self.settings.get("short_read", {}).get("reference_fasta"):
            self.sr_ref_fasta_entry.insert(0, saved_ref)
        row_idx += 1
        self.sr_output_root_entry = self._create_path_entry(tab, "Output Root Dir:", row_idx, tkinter.filedialog.askdirectory, is_file=False)
        if saved_output := self.settings.get("short_read", {}).get("output_root"):
            self.sr_output_root_entry.insert(0, saved_output)
        row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Run Name (optional):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sr_run_name_entry = customtkinter.CTkEntry(master=tab)
        if saved_run_name := self.settings.get("short_read", {}).get("run_name"):
            self.sr_run_name_entry.insert(0, saved_run_name)
        self.sr_run_name_entry.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2); row_idx += 1

        # ART Platform and Profile on same row
        customtkinter.CTkLabel(master=tab, text="ART Platform:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sr_art_platform_var = customtkinter.StringVar(value=self.settings.get("short_read", {}).get("art_platform", "illumina"))
        customtkinter.CTkOptionMenu(master=tab, variable=self.sr_art_platform_var, values=["illumina"], width=120).grid(row=row_idx, column=1, padx=10, pady=5, sticky="w")
        
        customtkinter.CTkLabel(master=tab, text="ART Profile:").grid(row=row_idx, column=2, padx=10, pady=5, sticky="w")
        self.sr_art_profile_var = customtkinter.StringVar(value=self.settings.get("short_read", {}).get("art_profile", DEFAULT_ILLUMINA_PROFILE))
        customtkinter.CTkOptionMenu(master=tab, variable=self.sr_art_profile_var, values=KNOWN_ART_PROFILES["illumina"], width=120).grid(row=row_idx, column=3, padx=10, pady=5, sticky="w"); row_idx += 1
        
        # Genome build selection
        customtkinter.CTkLabel(master=tab, text="Genome Build:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sr_genome_build_var = customtkinter.StringVar(value=self.settings.get("short_read", {}).get("genome_build", ""))
        self.sr_genome_build_menu = customtkinter.CTkOptionMenu(master=tab, variable=self.sr_genome_build_var, 
                                                               values=["", "hg19/GRCh37", "hg38/GRCh38", "T2T-CHM13", "mm10", "mm39", "Other"])
        self.sr_genome_build_menu.grid(row=row_idx, column=1, padx=10, pady=5, sticky="w")
        # Custom build entry (shown when "Other" is selected)
        self.sr_custom_build_entry = customtkinter.CTkEntry(master=tab, placeholder_text="Enter custom build")
        if saved_custom := self.settings.get("short_read", {}).get("custom_build"):
            self.sr_custom_build_entry.insert(0, saved_custom)
        self.sr_custom_build_entry.grid(row=row_idx, column=2, padx=10, pady=5, sticky="w")
        self.sr_custom_build_entry.grid_remove()  # Initially hidden
        
        self.sr_genome_build_var.trace_add("write", self._on_genome_build_change)
        row_idx += 1

        # Read Length and Depth on same row
        customtkinter.CTkLabel(master=tab, text="Read Length:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sr_read_length_entry = customtkinter.CTkEntry(master=tab, width=120)
        self.sr_read_length_entry.insert(0, self.settings.get("short_read", {}).get("read_length", "150"))
        self.sr_read_length_entry.grid(row=row_idx, column=1, sticky="w", padx=10, pady=5)
        
        customtkinter.CTkLabel(master=tab, text="Depth:").grid(row=row_idx, column=2, padx=10, pady=5, sticky="w")
        self.sr_depth_entry = customtkinter.CTkEntry(master=tab, width=120)
        self.sr_depth_entry.insert(0, self.settings.get("short_read", {}).get("depth", "50"))
        self.sr_depth_entry.grid(row=row_idx, column=3, sticky="w", padx=10, pady=5)
        row_idx += 1
        # Mean Fragment and Std Dev Fragment on same row
        customtkinter.CTkLabel(master=tab, text="Mean Fragment:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sr_mean_frag_entry = customtkinter.CTkEntry(master=tab, width=120)
        self.sr_mean_frag_entry.insert(0, self.settings.get("short_read", {}).get("mean_fragment", "400"))
        self.sr_mean_frag_entry.grid(row=row_idx, column=1, sticky="w", padx=10, pady=5)
        
        customtkinter.CTkLabel(master=tab, text="Std Dev Fragment:").grid(row=row_idx, column=2, padx=10, pady=5, sticky="w")
        self.sr_std_dev_frag_entry = customtkinter.CTkEntry(master=tab, width=120)
        self.sr_std_dev_frag_entry.insert(0, self.settings.get("short_read", {}).get("std_dev_fragment", "50"))
        self.sr_std_dev_frag_entry.grid(row=row_idx, column=3, sticky="w", padx=10, pady=5)
        row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Genomic Ranges:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.short_sim_ranges_textbox = customtkinter.CTkTextbox(master=tab, height=80, wrap="word")
        self.short_sim_ranges_textbox.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=3)
        # Use saved ranges, or default based on genome build
        saved_ranges = self.settings.get("short_read", {}).get("genomic_ranges", "")
        if saved_ranges:
            self.short_sim_ranges_textbox.insert("1.0", saved_ranges)
        else:
            # Use build-specific example if available
            build = self.sr_genome_build_var.get()
            default_example = EXAMPLE_RANGES.get(build, EXAMPLE_RANGES[""])
            self.short_sim_ranges_textbox.insert("1.0", default_example)
        row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Timeout (s):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w"); self.sr_timeout_entry = customtkinter.CTkEntry(master=tab); self.sr_timeout_entry.insert(0, self.settings.get("short_read", {}).get("timeout", "3600.0")); self.sr_timeout_entry.grid(row=row_idx, column=1, sticky="ew", columnspan=2, padx=10, pady=5); row_idx+=1

        self.sr_single_end_var = customtkinter.BooleanVar(value=self.settings.get("short_read", {}).get("single_end", False))
        customtkinter.CTkCheckBox(master=tab, text="Single-End Simulation", variable=self.sr_single_end_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); row_idx += 1
        self.sr_overwrite_var = customtkinter.BooleanVar(value=self.settings.get("short_read", {}).get("overwrite", False))
        customtkinter.CTkCheckBox(master=tab, text="Overwrite Output", variable=self.sr_overwrite_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); row_idx += 1
        self.sr_auto_index_var = customtkinter.BooleanVar(value=self.settings.get("short_read", {}).get("auto_index", True))
        customtkinter.CTkCheckBox(master=tab, text="Auto-Index FASTA", variable=self.sr_auto_index_var).grid(row=row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); row_idx += 1
        
        # Auto-align checkbox
        self.sr_auto_align_var = customtkinter.BooleanVar(value=self.settings.get("short_read", {}).get("auto_align", False))
        self.sr_auto_align_checkbox = customtkinter.CTkCheckBox(master=tab, text="Auto-align to BAM", variable=self.sr_auto_align_var, command=self._on_auto_align_toggle)
        self.sr_auto_align_checkbox.grid(row=row_idx, column=0, padx=10, pady=5, sticky="w", columnspan=3); row_idx += 1
        
        # Aligner selection (shown only when auto-align is checked)
        self.sr_aligner_label = customtkinter.CTkLabel(master=tab, text="Aligner:")
        self.sr_aligner_label.grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sr_aligner_label.grid_remove()  # Initially hidden
        
        self.sr_aligner_var = customtkinter.StringVar(value=self.settings.get("short_read", {}).get("aligner", "bwa"))
        self.sr_aligner_menu = customtkinter.CTkOptionMenu(master=tab, variable=self.sr_aligner_var, values=["bwa"])
        self.sr_aligner_menu.grid(row=row_idx, column=1, padx=10, pady=5, sticky="w")
        self.sr_aligner_menu.grid_remove()  # Initially hidden
        row_idx += 1

        # Create button frame for Run and IGV buttons
        button_frame = customtkinter.CTkFrame(master=tab, fg_color="transparent")
        button_frame.grid(row=row_idx, column=0, columnspan=4, padx=10, pady=(10,0))
        
        self.run_short_sim_button = customtkinter.CTkButton(master=button_frame, text="Run Simulation", command=self.run_short_simulation_thread)
        self.run_short_sim_button.pack(side="left", padx=5)
        
        self.igv_button_short = customtkinter.CTkButton(master=button_frame, text="Import into IGV", command=self.on_igv_button_short_click, state="disabled")
        self.igv_button_short.pack(side="left", padx=5)
        
        # Check if auto-align was saved as enabled
        if self.sr_auto_align_var.get():
            self._on_auto_align_toggle()
    
    def _auto_detect_genome_build(self, ref_path: str):
        """Auto-detect genome build from reference filename and set it in the dropdown"""
        ref_path_lower = ref_path.lower()
        
        # Check for common genome build patterns in filename
        if "grch38" in ref_path_lower or "hg38" in ref_path_lower:
            self.sr_genome_build_var.set("hg38/GRCh38")
            self.update_status("Auto-detected genome build: hg38/GRCh38", clear_first=False)
        elif "grch37" in ref_path_lower or "hg19" in ref_path_lower or "hs37d5" in ref_path_lower:
            self.sr_genome_build_var.set("hg19/GRCh37")
            self.update_status("Auto-detected genome build: hg19/GRCh37", clear_first=False)
        elif "t2t" in ref_path_lower or "chm13" in ref_path_lower:
            self.sr_genome_build_var.set("T2T-CHM13")
            self.update_status("Auto-detected genome build: T2T-CHM13", clear_first=False)
        elif "mm10" in ref_path_lower:
            self.sr_genome_build_var.set("mm10")
            self.update_status("Auto-detected genome build: mm10", clear_first=False)
        elif "mm39" in ref_path_lower:
            self.sr_genome_build_var.set("mm39")
            self.update_status("Auto-detected genome build: mm39", clear_first=False)
        else:
            # Could not detect - set to empty
            self.sr_genome_build_var.set("")
            self.update_status("Could not auto-detect genome build from filename", is_error=False, clear_first=False)
    
    def _suggest_reference_for_build(self, build: str):
        """Suggest and optionally set reference path based on genome build"""
        # Get the base directory where BaseBuddy is installed
        base_dir = Path(__file__).parent.parent.parent
        
        current_ref = self.sr_ref_fasta_entry.get()
        
        # Only suggest if field is empty
        if not current_ref:
            # Check for existing reference files
            for ref_path in DEFAULT_REFERENCES.get(build, []):
                full_path = base_dir / ref_path
                if full_path.exists():
                    # Found a matching reference file
                    self.sr_ref_fasta_entry.delete(0, "end")
                    self.sr_ref_fasta_entry.insert(0, str(full_path))
                    self.update_status(f"Reference auto-filled for {build}: {full_path.name}", clear_first=False)
                    break
        elif self._is_mismatched_reference(current_ref, build):
            # Warn about potential mismatch but don't change the user's selection
            self.update_status(f"Warning: Selected reference may not match {build} coordinates", is_error=True, clear_first=False)
    
    def _is_mismatched_reference(self, ref_path: str, build: str) -> bool:
        """Check if the current reference path doesn't match the selected build"""
        ref_path_lower = ref_path.lower()
        
        # Check for build mismatches
        if build == "hg19/GRCh37":
            return "grch38" in ref_path_lower or "hg38" in ref_path_lower
        elif build == "hg38/GRCh38":
            return "grch37" in ref_path_lower or "hg19" in ref_path_lower or "hs37d5" in ref_path_lower
        elif build == "mm10":
            return "mm39" in ref_path_lower
        elif build == "mm39":
            return "mm10" in ref_path_lower
        
        return False

    def _on_auto_align_toggle(self):
        """Show/hide aligner options based on auto-align checkbox"""
        if self.sr_auto_align_var.get():
            self.sr_aligner_label.grid()
            self.sr_aligner_menu.grid()
        else:
            self.sr_aligner_label.grid_remove()
            self.sr_aligner_menu.grid_remove()
    
    def _on_genome_build_change(self, *args):
        """Show/hide custom build entry and suggest reference paths based on genome build selection"""
        build = self.sr_genome_build_var.get()
        
        if build == "Other":
            self.sr_custom_build_entry.grid()
        else:
            self.sr_custom_build_entry.grid_remove()
            
        # Suggest reference path based on genome build
        if build and build != "Other" and build in DEFAULT_REFERENCES:
            self._suggest_reference_for_build(build)
            
        # Update genomic ranges example if the textbox is empty or has default values
        if hasattr(self, 'short_sim_ranges_textbox'):
            current_ranges = self.short_sim_ranges_textbox.get("1.0", "end-1c").strip()
            # Check if it's empty or contains a default example
            if not current_ranges or any(default in current_ranges for default in EXAMPLE_RANGES.values()):
                example = EXAMPLE_RANGES.get(build, EXAMPLE_RANGES[""])
                self.short_sim_ranges_textbox.delete("1.0", "end")
                self.short_sim_ranges_textbox.insert("1.0", example)

    def run_short_simulation_thread(self):
        if hasattr(self, 'igv_button_spike'): self.igv_button_spike.configure(state="disabled")
        if hasattr(self, 'igv_button_short'): self.igv_button_short.configure(state="disabled")
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
        
        # Get genome build
        genome_build = self.sr_genome_build_var.get()
        if genome_build == "Other":
            genome_build = self.sr_custom_build_entry.get().strip()
        elif not genome_build:
            genome_build = None

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
            "auto_align": self.sr_auto_align_var.get(),
            "aligner": self.sr_aligner_var.get() if self.sr_auto_align_var.get() else "bwa",
            "genome_build": genome_build,
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
            "auto_align": self.sr_auto_align_var.get(),
            "aligner": self.sr_aligner_var.get() if self.sr_auto_align_var.get() else "bwa",
            "genome_build": genome_build,
            "variants_list": None, # Not implemented in this tab's GUI yet
            "genomic_ranges": parsed_genomic_ranges
        }

    def _handle_short_sim_results(self, result_dict: Dict[str, Any]):
        # Store results for IGV button
        self.short_sim_results = result_dict
        
        # Enable IGV button if BAM file was created
        if result_dict.get("output_files"):
            bam_files = [f for f in result_dict["output_files"] if f["type"] == "BAM"]
            if bam_files:
                self.igv_button_short.configure(state="normal")
    
    def on_igv_button_short_click(self):
        if not hasattr(self, 'short_sim_results') or not self.short_sim_results:
            self.update_status("No simulation results available to send.", is_error=True, clear_first=True)
            return
        
        # Get BAM files from results and construct full paths
        output_files = self.short_sim_results.get("output_files", [])
        output_dir = self.short_sim_results.get("output_directory", "")
        
        # Construct full paths for BAM files
        bam_files = []
        for f in output_files:
            if f["type"] == "BAM":
                # If path is relative (just filename), join with output directory
                file_path = f["path"]
                if not Path(file_path).is_absolute():
                    file_path = str(Path(output_dir) / file_path)
                bam_files.append(file_path)
        
        ref_fasta = self.short_sim_results.get("reference_fasta_used") or self.sr_ref_fasta_entry.get()
        
        if bam_files and ref_fasta:
            self.update_status(f"Loading {len(bam_files)} BAM file(s) in IGV...", clear_first=True)
            self.send_to_igv_desktop(reference_fasta_path=ref_fasta, data_file_paths=bam_files)
        else:
            self.update_status("No BAM files found in simulation results.", is_error=True, clear_first=True)

    def create_spike_variants_tab(self):
        tab = self.spike_tab
        tab.grid_columnconfigure(1, weight=1)
        row_idx = 0

        self.sv_ref_fasta_entry = self._create_path_entry(tab, "Reference FASTA:", row_idx, tkinter.filedialog.askopenfilename)
        row_idx += 1

        self.sv_input_bam_entry = self._create_path_entry(tab, "Input BAM:", row_idx, tkinter.filedialog.askopenfilename)
        row_idx += 1

        customtkinter.CTkLabel(master=tab, text="Manual Variants (VCF format):").grid(row=row_idx, column=0, padx=10, pady=5, sticky="nw")
        self.sv_manual_variants_textbox = customtkinter.CTkTextbox(master=tab, height=80, wrap="word")
        self.sv_manual_variants_textbox.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew", columnspan=2)
        # Add placeholder text
        self.sv_manual_variants_textbox.insert("1.0", "# Enter variants in VCF format (tab-separated):\n# CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO\n# Example:\n# 7\t140453136\t.\tT\tA\t.\t.\t.")
        row_idx += 1

        self.sv_vcf_entry = self._create_path_entry(tab, "Input VCF (optional):", row_idx, tkinter.filedialog.askopenfilename)
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

        customtkinter.CTkLabel(master=tab, text="Strand Bias:").grid(row=row_idx, column=0, padx=10, pady=5, sticky="w")
        self.sv_strand_bias_entry = customtkinter.CTkEntry(master=tab)
        self.sv_strand_bias_entry.insert(0, "0.5") # Default 0.5 = no bias
        self.sv_strand_bias_entry.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew")
        customtkinter.CTkLabel(master=tab, text="(0=all reverse, 0.5=no bias, 1=all forward)", font=("Arial", 10)).grid(row=row_idx, column=2, padx=5, pady=5, sticky="w")
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

        # Get input BAM from new entry widget
        input_bam_path = self.sv_input_bam_entry.get().strip()
        if not input_bam_path:
            self.update_status("Input BAM path is required.", is_error=True, clear_first=True)
            return None
        input_bam_paths_str = [input_bam_path]  # Keep as list for compatibility

        # Get variants from manual entry or VCF file
        manual_variants_text = self.sv_manual_variants_textbox.get("1.0", "end-1c").strip()
        vcf_file_str = self.sv_vcf_entry.get().strip() or None
        
        # Parse manual variants if provided
        temp_vcf_path = None
        if manual_variants_text:
            # Remove comment lines and empty lines
            variant_lines = [line for line in manual_variants_text.splitlines() 
                           if line.strip() and not line.strip().startswith('#')]
            
            if variant_lines:
                # Create temporary VCF file from manual variants
                import tempfile
                with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
                    f.write("##fileformat=VCFv4.2\n")
                    f.write("##source=BaseBuddy_Manual_Entry\n")
                    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                    for line in variant_lines:
                        f.write(line + "\n")
                    temp_vcf_path = f.name

        # Use manual variants VCF if created, otherwise use provided VCF file
        variants_vcf = temp_vcf_path or vcf_file_str
        
        if not variants_vcf:
            self.update_status("Either manual variants or a VCF file must be provided.", is_error=True, clear_first=True)
            return None
        
        # For now, use the same VCF for both SNPs and indels
        # The backend will handle parsing them appropriately
        snp_vcf_str = variants_vcf
        indel_vcf_str = variants_vcf

        picard_jar_str = self.sv_picard_jar_entry.get().strip() or None

        output_bam_prefix = self.sv_output_prefix_entry.get()
        if not output_bam_prefix:
            self.update_status("Output File Prefix is required.", is_error=True, clear_first=True)
            return None

        run_name_str = self.sv_run_name_entry.get() or None

        vaf_val_str = self.sv_vaf_entry.get()
        seed_val_str = self.sv_seed_entry.get()
        strand_bias_str = self.sv_strand_bias_entry.get()
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
        try:
            strand_bias_val = float(strand_bias_str)
            if not (0 <= strand_bias_val <= 1.0):
                self.update_status("Strand bias must be between 0 and 1.0.", is_error=True, clear_first=True); return None
        except ValueError:
            self.update_status(f"Invalid Strand Bias value: '{strand_bias_str}'. Must be a number.", is_error=True, clear_first=True); return None

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
            "strand_bias": strand_bias_val,
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

        # Create button frame for Run and IGV buttons
        button_frame = customtkinter.CTkFrame(master=tab, fg_color="transparent")
        button_frame.grid(row=row_idx, column=0, columnspan=3, padx=10, pady=20)
        
        self.run_long_sim_button = customtkinter.CTkButton(master=button_frame, text="Run Simulation", command=self.run_long_simulation_thread)
        self.run_long_sim_button.pack(side="left", padx=5)
        
        self.igv_button_long = customtkinter.CTkButton(master=button_frame, text="Import into IGV", command=self.on_igv_button_long_click, state="disabled")
        self.igv_button_long.pack(side="left", padx=5)

    def run_long_simulation_thread(self):
        if hasattr(self, 'igv_button_spike'): self.igv_button_spike.configure(state="disabled")
        if hasattr(self, 'igv_button_long'): self.igv_button_long.configure(state="disabled")
        self._start_runner_thread(
            button_to_disable=self.run_long_sim_button,
            runner_func=src_bb_runner.simulate_long,
            param_extractor_func=self._get_long_sim_params,
            result_handler_func=self._handle_long_sim_results
        )
    
    def _handle_long_sim_results(self, result_dict: Dict[str, Any]):
        # Store results for IGV button
        self.long_sim_results = result_dict
        
        # Long reads typically produce FASTQ, but check if any BAM files were created
        if result_dict.get("output_files"):
            bam_files = [f for f in result_dict["output_files"] if f["type"] == "BAM"]
            if bam_files:
                self.igv_button_long.configure(state="normal")
    
    def on_igv_button_long_click(self):
        if not hasattr(self, 'long_sim_results') or not self.long_sim_results:
            self.update_status("No simulation results available to send.", is_error=True, clear_first=True)
            return
        
        # Get BAM files from results and construct full paths
        output_files = self.long_sim_results.get("output_files", [])
        output_dir = self.long_sim_results.get("output_directory", "")
        
        # Construct full paths for BAM files
        bam_files = []
        for f in output_files:
            if f["type"] == "BAM":
                # If path is relative (just filename), join with output directory
                file_path = f["path"]
                if not Path(file_path).is_absolute():
                    file_path = str(Path(output_dir) / file_path)
                bam_files.append(file_path)
        
        ref_fasta = self.long_sim_results.get("reference_fasta_used") or self.lr_ref_fasta_entry.get()
        
        if bam_files and ref_fasta:
            self.update_status(f"Loading {len(bam_files)} BAM file(s) in IGV...", clear_first=True)
            self.send_to_igv_desktop(reference_fasta_path=ref_fasta, data_file_paths=bam_files)
        else:
            self.update_status("No BAM files found in long read simulation results.", is_error=True, clear_first=True)

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
            logger.info(f"Runner function completed successfully")
            if result_handler_func:
                result_handler_func(result)
            logger.info(f"Sending completion message to GUI")
            self.thread_queue.put(("completion", (result, True), button_to_disable))
            logger.info(f"Completion message sent")
        except (BaseBuddyInputError, BaseBuddyToolError, BaseBuddyFileError, BaseBuddyConfigError, NameError) as e:
            self.thread_queue.put(("error", str(e), button_to_disable))
        except Exception as e:
            tb_str = traceback.format_exc()
            self.thread_queue.put(("error", f"An unexpected error occurred: {str(e)}\nTraceback:\n{tb_str}", button_to_disable))

    def create_apply_signature_tab(self):
        tab = self.apply_sig_tab
        tab.grid_columnconfigure(1, weight=1)
        row_idx = 0

        # Reference FASTA
        self.sig_ref_fasta_entry = self._create_path_entry(tab, "Reference FASTA:", row_idx, tkinter.filedialog.askopenfilename)
        row_idx += 1

        # Signature selection frame
        sig_frame = customtkinter.CTkFrame(tab)
        sig_frame.grid(row=row_idx, column=0, columnspan=3, sticky="ew", padx=10, pady=5)
        sig_frame.grid_columnconfigure(1, weight=1)
        row_idx += 1

        # Bundled or custom radio
        self.sig_source_var = tkinter.StringVar(value="bundled")
        
        bundled_radio = customtkinter.CTkRadioButton(sig_frame, text="Bundled Signature", 
                                                    variable=self.sig_source_var, value="bundled",
                                                    command=self._on_sig_source_change)
        bundled_radio.grid(row=0, column=0, padx=5, pady=5, sticky="w")
        
        custom_radio = customtkinter.CTkRadioButton(sig_frame, text="Custom File", 
                                                   variable=self.sig_source_var, value="custom",
                                                   command=self._on_sig_source_change)
        custom_radio.grid(row=0, column=1, padx=5, pady=5, sticky="w")

        # Bundled type selection
        self.sig_type_label = customtkinter.CTkLabel(sig_frame, text="Signature Type:")
        self.sig_type_label.grid(row=1, column=0, padx=5, pady=5, sticky="e")
        
        self.sig_type_var = tkinter.StringVar(value="sbs")
        self.sig_type_menu = customtkinter.CTkOptionMenu(sig_frame, variable=self.sig_type_var,
                                                        values=["sbs", "dbs", "id"])
        self.sig_type_menu.grid(row=1, column=1, padx=5, pady=5, sticky="ew")

        # Signature name
        self.sig_name_label = customtkinter.CTkLabel(sig_frame, text="Signature Name:")
        self.sig_name_label.grid(row=2, column=0, padx=5, pady=5, sticky="e")
        
        self.sig_name_entry = customtkinter.CTkEntry(sig_frame, placeholder_text="e.g., SBS1")
        self.sig_name_entry.grid(row=2, column=1, padx=5, pady=5, sticky="ew")
        self.sig_name_entry.insert(0, "SBS1")

        # Custom file entry (hidden by default)
        self.sig_file_label = customtkinter.CTkLabel(tab, text="Signature File:")
        self.sig_file_entry = self._create_path_entry(tab, "Signature File:", row_idx, tkinter.filedialog.askopenfilename)
        row_idx += 1
        
        # Initially hide custom file selector
        self.sig_file_label.grid_remove()
        self.sig_file_entry.grid_remove()

        # Output file
        self.sig_output_entry = self._create_path_entry(tab, "Output FASTA:", row_idx, tkinter.filedialog.asksaveasfilename)
        row_idx += 1

        # Number of mutations
        self.sig_num_mut_label = customtkinter.CTkLabel(tab, text="Number of Mutations:")
        self.sig_num_mut_label.grid(row=row_idx, column=0, padx=10, pady=5, sticky="e")
        
        self.sig_num_mut_entry = customtkinter.CTkEntry(tab, placeholder_text="1000")
        self.sig_num_mut_entry.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew")
        self.sig_num_mut_entry.insert(0, "1000")
        row_idx += 1

        # Seed
        self.sig_seed_label = customtkinter.CTkLabel(tab, text="Random Seed:")
        self.sig_seed_label.grid(row=row_idx, column=0, padx=10, pady=5, sticky="e")
        
        self.sig_seed_entry = customtkinter.CTkEntry(tab, placeholder_text="42")
        self.sig_seed_entry.grid(row=row_idx, column=1, padx=10, pady=5, sticky="ew")
        self.sig_seed_entry.insert(0, "42")
        row_idx += 1

        # Run button
        # Create button frame for Run and IGV buttons
        button_frame = customtkinter.CTkFrame(master=tab, fg_color="transparent")
        button_frame.grid(row=row_idx, column=0, columnspan=3, padx=10, pady=20)
        
        self.run_apply_sig_button = customtkinter.CTkButton(button_frame, text="Apply Signature", 
                                                           command=self.on_apply_signature_click)
        self.run_apply_sig_button.pack(side="left", padx=5)
        
        self.igv_button_apply_sig = customtkinter.CTkButton(button_frame, text="Import into IGV", 
                                                           command=self.on_igv_button_apply_sig_click, state="disabled")
        self.igv_button_apply_sig.pack(side="left", padx=5)

    def _on_sig_source_change(self):
        """Toggle visibility of bundled vs custom signature options"""
        if self.sig_source_var.get() == "bundled":
            # Show bundled options
            self.sig_type_label.grid()
            self.sig_type_menu.grid()
            self.sig_name_label.grid()
            self.sig_name_entry.grid()
            # Hide custom file
            self.sig_file_label.grid_remove()
            self.sig_file_entry.grid_remove()
        else:
            # Hide bundled options
            self.sig_type_label.grid_remove()
            self.sig_type_menu.grid_remove()
            # Show custom file and name
            self.sig_file_label.grid(row=3, column=0, padx=10, pady=5, sticky="e")
            self.sig_file_entry.grid(row=3, column=1, padx=10, pady=5, sticky="ew")

    def on_apply_signature_click(self):
        if hasattr(self, 'igv_button_apply_sig'): self.igv_button_apply_sig.configure(state="disabled")
        self._start_runner_thread(self.run_apply_sig_button, src_bb_runner.apply_signature_to_fasta,
                                self._get_apply_signature_params, self._handle_apply_sig_results)

    def _get_apply_signature_params(self):
        ref_fasta = self.sig_ref_fasta_entry.get().strip()
        if not ref_fasta:
            raise ValueError("Reference FASTA is required")
        
        output_fasta = self.sig_output_entry.get().strip()
        if not output_fasta:
            output_fasta = str(Path(ref_fasta).parent / "mutated.fa")
        
        num_mutations = int(self.sig_num_mut_entry.get() or "1000")
        seed = int(self.sig_seed_entry.get() or "42")
        
        params = {
            "reference_fasta": ref_fasta,
            "output_fasta": output_fasta,
            "num_mutations": num_mutations,
            "seed": seed
        }
        
        if self.sig_source_var.get() == "bundled":
            sig_type = self.sig_type_var.get()
            sig_name = self.sig_name_entry.get().strip()
            if not sig_name:
                raise ValueError("Signature name is required")
            
            # Get bundled signature file path
            sig_file = sig_utils.get_bundled_signature_path(sig_type)
            params["signature_file"] = str(sig_file)
            params["signature_name"] = sig_name
        else:
            sig_file = self.sig_file_entry.get().strip()
            if not sig_file:
                raise ValueError("Signature file is required")
            sig_name = self.sig_name_entry.get().strip()
            if not sig_name:
                raise ValueError("Signature name is required")
            
            params["signature_file"] = sig_file
            params["signature_name"] = sig_name
        
        return params
    
    def _handle_apply_sig_results(self, result_dict: Dict[str, Any]):
        # Store results for IGV button
        self.apply_sig_results = result_dict
        
        # Enable IGV button if output FASTA was created
        output_fasta = result_dict.get("output_fasta_path")
        if output_fasta and Path(output_fasta).exists():
            self.igv_button_apply_sig.configure(state="normal")
    
    def on_igv_button_apply_sig_click(self):
        if not hasattr(self, 'apply_sig_results') or not self.apply_sig_results:
            self.update_status("No signature results available to send.", is_error=True, clear_first=True)
            return
        
        output_fasta = self.apply_sig_results.get("output_fasta_path")
        
        if output_fasta and Path(output_fasta).exists():
            self.update_status(f"Loading mutated FASTA in IGV...", clear_first=True)
            # For FASTA files, we load them as genome references
            self.send_to_igv_desktop(reference_fasta_path=output_fasta)
        else:
            self.update_status("No output FASTA found in results.", is_error=True, clear_first=True)

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


def run_gui():
    app = BaseBuddyGUI()
    app.mainloop()


if __name__ == "__main__":
    run_gui()
