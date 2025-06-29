import tkinter as tk
from tkinter import filedialog, scrolledtext, messagebox
from gene_finder import run_analysis
from off_target_finder import run_cas_offinder


class GeneFinder:
    # File type constants
    GENOME_FILE = 0
    ANNOTATION_FILE = 1
    GFF_FILE = 2
    SGRNA_FILE = 3
    CAS_OFFINDER_EXE = 4
    TOTAL_FILES = 5

    # Default values
    DEFAULT_GENE = "cesa"
    DEFAULT_GENE_FAMILY = "cellulose synthase subfamily"
    DEFAULT_CONSERVED_REGION = "WPGN"
    DEFAULT_MAX_OFF_TARGETS = "1"

    def __init__(self, root):
        self.root = root
        self.root.title("GOFFinder")
        self.file_paths = [None] * self.TOTAL_FILES
        self.output_dir = ""

        self._initialize_variables()
        self._create_ui()

    def _initialize_variables(self):
        """Initialize StringVar objects with default values."""
        self.gene = tk.StringVar(value=self.DEFAULT_GENE)
        self.gene_family = tk.StringVar(value=self.DEFAULT_GENE_FAMILY)
        self.conserved_region = tk.StringVar(value=self.DEFAULT_CONSERVED_REGION)
        self.max_off_targets = tk.StringVar(value=self.DEFAULT_MAX_OFF_TARGETS)

        self.file_labels = [
            "Genome.fasta File",
            ".Annotation.gz File",
            ".gff.gz File",
            "sgRNA.fasta File",
            "cas-offinder.exe File"
        ]
        self.file_entries = []

    def _create_ui(self):
        """Create the main user interface."""
        self._create_file_upload_section()
        self._create_output_directory_section()
        self._create_gene_input_section()
        self._create_analysis_section()
        self._create_off_target_section()
        self._create_display_section()

    def _create_file_upload_section(self):
        """Create the file upload section for required files."""
        file_frame = tk.LabelFrame(self.root, text="1. Upload Required Files", padx=10, pady=10)
        file_frame.pack(fill="x", padx=10, pady=5)

        for i in range(3):  # First 3 files go in main file frame
            self._create_file_row(file_frame, i)

    def _create_output_directory_section(self):
        """Create the output directory selection section."""
        output_dir_frame = tk.LabelFrame(self.root, text="2. Select Output Directory", padx=10, pady=10)
        output_dir_frame.pack(fill="x", padx=10, pady=5)

        tk.Button(output_dir_frame, text="Choose Output Directory",
                  command=self.set_output_directory).pack(side="left")
        self.output_dir_label = tk.Label(output_dir_frame, text="No directory selected", anchor="w")
        self.output_dir_label.pack(side="left", padx=10)

    def _create_gene_input_section(self):
        """Create the gene input section."""
        output_frame = tk.LabelFrame(self.root, text="3. Input gene, subfamily name and conserved residues",
                                     padx=10, pady=10)
        output_frame.pack(fill="x", padx=10, pady=5)

        tk.Entry(output_frame, textvariable=self.gene).pack(side="left", padx=10)
        tk.Entry(output_frame, textvariable=self.gene_family).pack(side="left", padx=10)
        tk.Entry(output_frame, textvariable=self.conserved_region).pack(side="left", padx=10)

    def _create_analysis_section(self):
        """Create the analysis buttons section."""
        analysis_frame = tk.LabelFrame(self.root, text="4. Run Analysis", padx=10, pady=10)
        analysis_frame.pack(fill="x", padx=10, pady=5)

        tk.Button(analysis_frame, text="Find Sequences", command=self.find_gene_sequences).pack(side="left", padx=10)
        tk.Button(analysis_frame, text="Find Off Target", command=self.find_off_target_sequences).pack(side="left", padx=10)

    def _create_off_target_section(self):
        """Create the off-target finder section."""
        off_target_frame = tk.LabelFrame(self.root, text="5. Off Target Finder", padx=10, pady=10)
        off_target_frame.pack(fill="x", padx=10, pady=5)

        # Add max off targets input
        max_targets_row = tk.Frame(off_target_frame)
        max_targets_row.pack(fill="x", pady=3)
        tk.Label(max_targets_row, text="Max Off Targets", width=20, anchor="w").pack(side="left")
        tk.Entry(max_targets_row, textvariable=self.max_off_targets).pack(side="left", padx=5)

        # Add remaining file inputs
        self._create_file_row(off_target_frame, self.SGRNA_FILE)
        self._create_file_row(off_target_frame, self.CAS_OFFINDER_EXE)

    def _create_file_row(self, parent, file_index):
        """Create a file input row with label, entry, and browse button."""
        row = tk.Frame(parent)
        row.pack(fill="x", pady=3)

        label = tk.Label(row, text=self.file_labels[file_index], width=20, anchor="w")
        label.pack(side="left")

        entry = tk.Entry(row, width=50, fg="gray")
        entry.insert(0, f"Select {self.file_labels[file_index].lower()}...")
        entry.config(state='readonly')
        entry.pack(side="left", padx=5)
        self.file_entries.append(entry)

        button = tk.Button(row, text="Browse", command=lambda idx=file_index: self.browse_file(idx))
        button.pack(side="left")

    def _create_display_section(self):
        """Create the main display area."""
        display_frame = tk.LabelFrame(self.root, text="Main Display", padx=10, pady=10)
        display_frame.pack(fill="both", expand=True, padx=10, pady=5)

        self.main_display = scrolledtext.ScrolledText(display_frame, wrap=tk.WORD, height=10)
        self.main_display.pack(fill="both", expand=True)

    def browse_file(self, index):
        filepath = filedialog.askopenfilename(title=f"Select {self.file_labels[index]}")
        if filepath:
            self.file_paths[index] = filepath
            self.file_entries[index].config(state='normal', fg='black')
            self.file_entries[index].delete(0, tk.END)
            self.file_entries[index].insert(0, filepath)
            self.file_entries[index].config(state='readonly')
            self.main_display.insert(tk.END, f"{self.file_labels[index]} selected: {filepath}\n")

    def set_output_directory(self):
        dir_selected = filedialog.askdirectory(title="Select Output Directory")
        if dir_selected:
            self.output_dir = dir_selected
            self.output_dir_label.config(text=dir_selected)
            self.main_display.insert(tk.END, f"Output directory set to: {dir_selected}\n")

    def validate_inputs(self):
        if None in self.file_paths:
            messagebox.showwarning("Missing Files", "Please upload all required files.")
            return False
        if not self.output_dir:
            messagebox.showwarning("No Output Directory", "Please select an output directory first.")
            return False
        return True

    def find_gene_sequences(self):
        if not self.validate_inputs():
            return
        self.main_display.insert(tk.END, chars=f"Finding {self.gene.get()} sequences...\n")

        run_analysis(
            self.file_paths[self.GENOME_FILE],
            self.file_paths[self.ANNOTATION_FILE],
            self.file_paths[self.GFF_FILE],
            self.output_dir,
            self.gene.get(),
            self.gene_family.get(),
            self.conserved_region.get()
        )
        self.main_display.insert(tk.END, f"Search complete!\n Saved three files to {self.output_dir} directory\n")

    def find_off_target_sequences(self):
        if not self.validate_inputs():
            return

        self.main_display.insert(tk.END, "Finding potential off target hits in the rosemary genome...\n")
        run_cas_offinder(
            self.file_paths[self.CAS_OFFINDER_EXE],
            self.file_paths[self.GENOME_FILE],
            self.file_paths[self.SGRNA_FILE],
            self.output_dir,
            self.gene.get(),
            self.conserved_region.get(),
            max_off_targets=self.max_off_targets.get()
        )
        self.main_display.insert(tk.END, f"Off target analysis complete! Saving results to {self.output_dir} directory\n")


if __name__ == "__main__":
    root = tk.Tk()
    app = GeneFinder(root)
    root.mainloop()
