import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
from tkinter import messagebox
from file_handler import load_fasta_with_orientation, clean_temp_dir
from endonucleases import get_sticky_end_enzymes, get_blunt_end_enzymes, restriction_enzymes
from logic import cut_blunt, cut_sticky
from logic_varpro import run_variant_analysis  # Nová logika s EVA a NCBI (run_NCBI_variant_anaysis)
from logic_varpro import analyze_variant_enzyme_cutting

class RENDO:
    def __init__(self, root):
        self.root = root
        self.root.title("RENDO")
        self.root.geometry("800x800")

        # Notebook se třemi záložkami
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True)

        # Záložka RENDOSim
        self.sim_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.sim_tab, text="RENDOSim")

        # Záložka RENDOVar
        self.var_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.var_tab, text="RENDOVar")

        # Záložka RENDOVarPro
        self.varpro_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.varpro_tab, text="RENDOVarPro")

        # Proměnné pro RENDOSim
        self.file_path = tk.StringVar()
        self.end_type = tk.StringVar(value="sticky")
        self.selected_enzyme = tk.StringVar()
        self.sequences = {}
        self.current_sequence = None
        self.trimmed_sequences = {}
        self.create_sim_tab_widgets()

        # Proměnné pro RENDOVar
        self.var_file_paths = []
        self.variant_sequences = {}
        self.optimal_results = {}
        self.create_var_tab_widgets()

        # Proměnné pro RENDOVarPro
        self.varpro_email = tk.StringVar()
        self.varpro_organism = tk.StringVar()
        self.varpro_chromosome = tk.StringVar()
        self.varpro_amplicon_pos = tk.StringVar()
        self.varpro_snp_position = tk.StringVar()
        self.varpro_assembly = tk.StringVar() 
        self.create_varpro_tab_widgets()

    # RENDOSim
    def create_sim_tab_widgets(self):
        file_frame = ttk.LabelFrame(self.sim_tab, text="File Selection", padding="5")
        file_frame.pack(fill="x", padx=5, pady=5)
        ttk.Entry(file_frame, textvariable=self.file_path, width=50).pack(side="left", padx=5)
        ttk.Button(file_frame, text="Browse", command=self.browse_file).pack(side="left", padx=5)

        enzyme_frame = ttk.LabelFrame(self.sim_tab, text="Enzyme Selection", padding="5")
        enzyme_frame.pack(fill="x", padx=5, pady=5)
        end_type_frame = ttk.Frame(enzyme_frame)
        end_type_frame.pack(fill="x", padx=5, pady=5)

        ttk.Radiobutton(end_type_frame, text="Sticky Ends", variable=self.end_type,
                        value="sticky", command=self.update_enzyme_list).pack(side="left", padx=5)
        ttk.Radiobutton(end_type_frame, text="Blunt Ends", variable=self.end_type,
                        value="blunt", command=self.update_enzyme_list).pack(side="left", padx=5)

        self.enzyme_dropdown = ttk.Combobox(enzyme_frame, textvariable=self.selected_enzyme)
        self.enzyme_dropdown.pack(fill="x", padx=5, pady=5)
        self.update_enzyme_list()

        sequence_frame = ttk.LabelFrame(self.sim_tab, text="Sequence", padding="5")
        sequence_frame.pack(fill="both", expand=True, padx=5, pady=5)
        self.sequence_text = tk.Text(sequence_frame, wrap="word", height=10)
        self.sequence_text.pack(fill="both", expand=True, padx=5, pady=5)

        results_frame = ttk.LabelFrame(self.sim_tab, text="Analysis Results", padding="5")
        results_frame.pack(fill="both", expand=True, padx=5, pady=5)
        self.results_text = tk.Text(results_frame, wrap="word", height=10)
        self.results_text.pack(fill="both", expand=True, padx=5, pady=5)

        button_frame = ttk.Frame(self.sim_tab)
        button_frame.pack(fill="x", padx=5, pady=5)
        ttk.Button(button_frame, text="Load", command=self.load_sequence).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Cut", command=self.analyze_sequence).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Save Cut", command=self.save_trimmed).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Clear", command=self.clear_all).pack(side="left", padx=5)

    def update_enzyme_list(self):
        from endonucleases import get_sticky_end_enzymes, get_blunt_end_enzymes
        if self.end_type.get() == "sticky":
            enzymes = get_sticky_end_enzymes()
        else:
            enzymes = get_blunt_end_enzymes()
        sorted_enzymes = sorted(enzymes.keys(), key=lambda name: name.lower())
        self.enzyme_dropdown['values'] = sorted_enzymes
        if sorted_enzymes:
            self.enzyme_dropdown.set(sorted_enzymes[0])

    def browse_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])
        if file_path:
            self.file_path.set(file_path)

    def load_sequence(self):
        try:
            file_path = self.file_path.get()
            if not file_path:
                messagebox.showwarning("Warning", "Prosím, vyberte soubor.")
                return
            from file_handler import load_fasta_with_orientation
            self.sequences, _ = load_fasta_with_orientation(file_path, "5_3")
            self.current_sequence = list(self.sequences.values())[0]
            self.sequence_text.delete("1.0", tk.END)
            self.sequence_text.insert(tk.END, self.current_sequence)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def analyze_sequence(self):
        if not self.current_sequence:
            messagebox.showwarning("Warning", "Nejdříve načtěte sekvenci.")
            return
        selected_enzyme = self.selected_enzyme.get()
        from endonucleases import restriction_enzymes
        enzyme_info = restriction_enzymes.get(selected_enzyme)
        if not enzyme_info:
            messagebox.showerror("Error", f"Enzym '{selected_enzyme}' nebyl nalezen v databázi.")
            return
        try:
            from logic import cut_blunt, cut_sticky
            overhang_length = enzyme_info.get('overhang', 0)
            self.results_text.delete("1.0", tk.END)
            if overhang_length > 0:
                f_frags, r_frags = cut_sticky(self.current_sequence, selected_enzyme)
                self.trimmed_sequences = {'forward': f_frags, 'reverse': r_frags}
                self.results_text.insert(tk.END, f"Forward strand fragments generated by {selected_enzyme}:\n")
                for i, fragment in enumerate(f_frags, 1):
                    self.results_text.insert(tk.END, f"Fragment {i} ({len(fragment)} bp): {fragment}\n")
                self.results_text.insert(tk.END, "\nReverse strand fragments:\n")
                for i, fragment in enumerate(r_frags, len(f_frags) + 1):
                    self.results_text.insert(tk.END, f"Fragment {i} ({len(fragment)} bp): {fragment}\n")
            else:
                b_frags = cut_blunt(self.current_sequence, selected_enzyme)
                self.trimmed_sequences = b_frags
                self.results_text.insert(tk.END, f"Fragmenty ({selected_enzyme}):\n")
                for i, fragment in enumerate(b_frags, 1):
                    self.results_text.insert(tk.END, f"Fragment {i} ({len(fragment)} bp): {fragment}\n")
        except Exception as e:
            messagebox.showerror("Error", f"Analýza selhala: {e}")

    def save_trimmed(self):
        if not self.trimmed_sequences:
            messagebox.showwarning("Warning", "Nejsou k dispozici žádné ořezané sekvence.")
            return
        from endonucleases import restriction_enzymes
        selected_enzyme = self.selected_enzyme.get()
        enzyme_info = restriction_enzymes.get(selected_enzyme)
        overhang_length = enzyme_info.get('overhang', 0)
        if overhang_length > 0:
            base_save_path = filedialog.asksaveasfilename(defaultextension=".fasta",
                                                          filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])
            if base_save_path:
                dir_name = os.path.dirname(base_save_path)
                base_name = os.path.splitext(os.path.basename(base_save_path))[0]
                forward_file = os.path.join(dir_name, f"{base_name}_F_cut.fasta")
                reverse_file = os.path.join(dir_name, f"{base_name}_R_cut.fasta")
                try:
                    with open(forward_file, "w") as f_fw:
                        for i, fragment in enumerate(self.trimmed_sequences['forward'], 1):
                            f_fw.write(f">Fragment_{i}\n{fragment}\n")
                    with open(reverse_file, "w") as f_rev:
                        for i, fragment in enumerate(self.trimmed_sequences['reverse'], 1):
                            f_rev.write(f">Fragment_{i}\n{fragment}\n")
                    messagebox.showinfo("Success", f"Uloženo:\n{forward_file}\n{reverse_file}")
                except Exception as e:
                    messagebox.showerror("Error", f"Ukládání selhalo: {e}")
        else:
            save_path = filedialog.asksaveasfilename(defaultextension=".fasta",
                                                     filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])
            if save_path:
                try:
                    with open(save_path, "w") as f:
                        for i, fragment in enumerate(self.trimmed_sequences, 1):
                            f.write(f">Fragment_{i}\n{fragment}\n")
                    messagebox.showinfo("Success", "Ořezané sekvence byly uloženy.")
                except Exception as e:
                    messagebox.showerror("Error", f"Ukládání selhalo: {e}")

    def clear_all(self):
        self.file_path.set("")
        self.sequence_text.delete("1.0", tk.END)
        self.results_text.delete("1.0", tk.END)
        self.current_sequence = None
        self.sequences = {}
        self.trimmed_sequences = {}

    # RENDOVar
    def create_var_tab_widgets(self):
        var_file_frame = ttk.LabelFrame(self.var_tab, text="Multi-FASTA File Selection", padding="5")
        var_file_frame.pack(fill="x", padx=5, pady=5)
        ttk.Button(var_file_frame, text="Browse Files", command=self.browse_files).pack(side="left", padx=5)
        self.var_files_label = ttk.Label(var_file_frame, text="No files selected")
        self.var_files_label.pack(side="left", padx=5)

        var_sequence_frame = ttk.LabelFrame(self.var_tab, text="Loaded Sequences", padding="5")
        var_sequence_frame.pack(fill="both", expand=True, padx=5, pady=5)
        self.var_sequence_text = tk.Text(var_sequence_frame, wrap="word", height=10)
        self.var_sequence_text.pack(fill="both", expand=True, padx=5, pady=5)

        var_results_frame = ttk.LabelFrame(self.var_tab, text="Optimal Enzyme Analysis Results", padding="5")
        var_results_frame.pack(fill="both", expand=True, padx=5, pady=5)
        self.var_results_text = tk.Text(var_results_frame, wrap="word", height=10)
        self.var_results_text.pack(fill="both", expand=True, padx=5, pady=5)

        var_button_frame = ttk.Frame(self.var_tab)
        var_button_frame.pack(fill="x", padx=5, pady=5)
        ttk.Button(var_button_frame, text="Load Sequences", command=self.load_variant_sequences).pack(side="left", padx=5)
        ttk.Button(var_button_frame, text="Analyze Variants", command=self.analyze_variants).pack(side="left", padx=5)
        ttk.Button(var_button_frame, text="Clear", command=self.clear_variants).pack(side="left", padx=5)

    def browse_files(self):
        file_paths = filedialog.askopenfilenames(filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])
        if file_paths:
            self.var_file_paths = file_paths
            self.var_files_label.config(text=f"{len(file_paths)} Files selected")

    def load_variant_sequences(self):
        if not self.var_file_paths:
            messagebox.showwarning("Warning", "Choose at least one file.")
            return
        self.variant_sequences.clear()
        self.var_sequence_text.delete("1.0", tk.END)
        try:
            from file_handler import load_fasta_with_orientation
            for file in self.var_file_paths:
                sequences, _ = load_fasta_with_orientation(file, "5_3")
                for header, seq in sequences.items():
                    key = f"{os.path.basename(file)}: {header}"
                    self.variant_sequences[key] = seq
                    self.var_sequence_text.insert("end", f">{key}\n{seq}\n\n")
        except Exception as e:
            messagebox.showerror("Error", f"Loading failed: {e}")

    def analyze_variants(self):
        self.var_results_text.delete("1.0", tk.END)
        try:
            if not self.variant_sequences:
                messagebox.showwarning("Warning", "No sequences loaded.")
                return
            from logic_var import analyze_variant_enzyme_cutting
            all_res, optimal_enz, candidate = analyze_variant_enzyme_cutting(self.variant_sequences)
            result_text = "Analysis results:\n"
            if candidate is not None:
                ratio = f"{candidate}:{candidate+1}"
                result_text += f"Optimal enzymes (fragment ratio = {ratio}):\n"
                for enzyme, ratio_tuple, enzyme_data in optimal_enz:
                    result_text += f"Enzyme: {enzyme} (ratio: {ratio_tuple[0]}:{ratio_tuple[1]})\n"
                    for seq_id, res in enzyme_data.items():
                        for strand, data in res.items():
                            if "fragment_count" in data:
                                result_text += f"  {seq_id} - {strand}: {data['fragment_count']} fragments, length: {data['fragment_lengths']}\n"
                result_text += "\n"
            else:
                result_text += "No enzymes found with ratio n : n+1.\n\n"

            result_text += "All enzymes with different fragments ratio:\n"
            for enzyme, en_res in all_res.items():
                differences = []
                for strand in ["forward", "reverse"]:
                    counts = []
                    for seq_id, rdata in en_res.items():
                        if strand in rdata and "fragment_count" in rdata[strand]:
                            counts.append(rdata[strand]["fragment_count"])
                    if counts:
                        diff = max(counts) - min(counts)
                        if diff > 0:
                            differences.append(diff)
                if differences:
                    display_diff = min(differences)
                    result_text += f"Enzym: {enzyme} (min defference: {display_diff})\n"
                    for seq_id, rdata in en_res.items():
                        for strand, data in rdata.items():
                            if "fragment_count" in data:
                                result_text += f"  {seq_id} - {strand}: {data['fragment_count']} fragments, lenght: {data['fragment_lengths']}\n"
                    result_text += "\n"

            self.var_results_text.insert("end", result_text)

        except Exception as e:
            messagebox.showerror("Error", f"Variant analysis failed: {e}")

    def clear_variants(self):
        self.var_file_paths = []
        self.variant_sequences.clear()
        self.var_files_label.config(text="No files selected")
        self.var_sequence_text.delete("1.0", tk.END)
        self.var_results_text.delete("1.0", tk.END)

    # RENDOVarPro
    def create_varpro_tab_widgets(self):
        # Přístupové údaje
        email_frame = ttk.LabelFrame(self.varpro_tab, text="Acces data", padding="5")
        email_frame.pack(fill="x", padx=5, pady=5)
        ttk.Label(email_frame, text="Email:").pack(side="left", padx=5)
        ttk.Entry(email_frame, textvariable=self.varpro_email, width=50).pack(side="left", padx=5)

        # Info o genu 
        info_frame = ttk.LabelFrame(self.varpro_tab, text="Species information", padding="5")
        info_frame.pack(fill="x", padx=5, pady=5)

        ttk.Label(info_frame, text="Organism:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        ttk.Entry(info_frame, textvariable=self.varpro_organism, width=30).grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(info_frame, text="Chromosome:").grid(row=0, column=2, padx=5, pady=5, sticky="w")
        ttk.Entry(info_frame, textvariable=self.varpro_chromosome, width=10).grid(row=0, column=3, padx=5, pady=5)

        ttk.Label(info_frame, text="Reference assembly:").grid(row=0, column=4, padx=5, pady=5, sticky="w")
        ttk.Entry(info_frame, textvariable=self.varpro_assembly, width=15).grid(row=0, column=5, padx=5, pady=5)

        # Pozice
        pos_frame = ttk.LabelFrame(self.varpro_tab, text="Position", padding="5")
        pos_frame.pack(fill="x", padx=5, pady=5)
        ttk.Label(pos_frame, text="Amplified positions (format: xxx-xxx):").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        ttk.Entry(pos_frame, textvariable=self.varpro_amplicon_pos, width=20).grid(row=0, column=1, padx=5, pady=5)
        ttk.Label(pos_frame, text="Analyzed SNV position:").grid(row=0, column=2, padx=5, pady=5, sticky="w")
        ttk.Entry(pos_frame, textvariable=self.varpro_snp_position, width=10).grid(row=0, column=3, padx=5, pady=5)

        # Zobrazení FASTA
        fasta_frame = ttk.LabelFrame(self.varpro_tab, text="FASTA sequence", padding="5")
        fasta_frame.pack(fill="both", expand=True, padx=5, pady=5)
        self.varpro_fasta_text = tk.Text(fasta_frame, wrap="word", height=5)
        self.varpro_fasta_text.pack(fill="both", expand=True, padx=5, pady=5)

        # Zobrazení VCF
        variant_file_frame = ttk.LabelFrame(self.varpro_tab, text="VCF file", padding="5")
        variant_file_frame.pack(fill="both", expand=True, padx=5, pady=5)
        self.varpro_variant_file_text = tk.Text(variant_file_frame, wrap="word", height=15)
        self.varpro_variant_file_text.pack(fill="both", expand=True, padx=5, pady=5)

        # Výsledky
        results_frame = ttk.LabelFrame(self.varpro_tab, text="Analysis results", padding="5")
        results_frame.pack(fill="both", expand=True, padx=5, pady=5)
        self.varpro_results_text = tk.Text(results_frame, wrap="word", height=10)
        self.varpro_results_text.pack(fill="both", expand=True, padx=5, pady=5)

        # Tlačítka
        button_frame = ttk.Frame(self.varpro_tab)
        button_frame.pack(fill="x", padx=5, pady=10)
        ttk.Button(button_frame, text="Analyze", command=self.analyze_varpro).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Clear", command=self.clear_varpro).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Save Cut", command=self.save_varpro_results).pack(side="left", padx=5)


    def analyze_varpro(self):
        try:
            # 1) Načtení a validace vstupů
            email = self.varpro_email.get().strip()
            organism = self.varpro_organism.get().strip()
            chromosome = self.varpro_chromosome.get().strip()
            amplicon_positions = self.varpro_amplicon_pos.get().strip()
            snp_pos_str = self.varpro_snp_position.get().strip()
            assembly = self.varpro_assembly.get().strip()
            if not all([email, organism, chromosome, amplicon_positions, snp_pos_str, assembly]):
                messagebox.showwarning("Warning", "Vyplňte prosím všechna pole.")
                return

            region_start = int(amplicon_positions.split('-')[0])
            snp_global_position = int(snp_pos_str)

            # 2) Spuštění analýzy (získáme sekvenci, varianty, iupac, temp_vcf, enzyme_results...)
            result = run_variant_analysis(
                email, organism, assembly,
                chromosome, amplicon_positions,
                snp_global_position
            )

            # 3) Zobrazení FASTA
            self.varpro_fasta_text.delete("1.0", tk.END)
            fasta_seq = result.get('sequence', '')
            self.varpro_fasta_text.insert(tk.END, fasta_seq)

            # 4) Zobrazení dočasného VCF
            self.varpro_variant_file_text.delete("1.0", tk.END)
            temp_vcf_path = result.get('temp_vcf')
            try:
                with open(temp_vcf_path, 'r') as vf:
                    vcf_content = vf.read()
            except Exception as e:
                vcf_content = f"Chyba při čtení VCF: {e}"
            self.varpro_variant_file_text.insert(tk.END, vcf_content)

            # 5) Příprava pro enzyme‐analysis
            enzyme_results = result['enzyme_results']
            optimal = result['optimal']
            candidate = result['candidate']
            seqs = result['variant_sequences']

            # 6) Vymazání starého výstupu
            self.varpro_results_text.delete("1.0", tk.END)

            if not optimal:
                self.varpro_results_text.insert(tk.END, "Žádné optimální enzymy.\n")
                return

            self.varpro_results_text.insert(
                tk.END,
                f"Optimální enzymy (poměr {candidate}:{candidate+1}):\n\n"
            )

            # 7) Výpis fragmentů pro každou variantu
            from endonucleases import restriction_enzymes
            from logic import cut_blunt, cut_sticky

            for enz in optimal:
                self.varpro_results_text.insert(tk.END, f"--- Enzyme: {enz} ---\n")
                overhang = restriction_enzymes[enz].get('overhang', 0)

                for var_id, seq in seqs.items():
                    # zjistíme bázi na SNV pozici
                    idx = snp_global_position - region_start
                    allele = seq[idx] if 0 <= idx < len(seq) else '?'
                    self.varpro_results_text.insert(
                        tk.END,
                        f"Variant: {var_id} (base {allele})\n"
                    )

                    # provést řezání
                    if overhang > 0:
                        ffrags, rfrags = cut_sticky(seq, enz)
                        # forward fragments
                        for i, frag in enumerate(ffrags, 1):
                            self.varpro_results_text.insert(
                                tk.END,
                                f"  Forward fragment {i} ({len(frag)} bp): {frag}\n"
                            )
                        # reverse fragments
                        for i, frag in enumerate(rfrags, 1):
                            self.varpro_results_text.insert(
                                tk.END,
                                f"  Reverse fragment {i} ({len(frag)} bp): {frag}\n"
                            )
                    else:
                        bfrags = cut_blunt(seq, enz)
                        for i, frag in enumerate(bfrags, 1):
                            self.varpro_results_text.insert(
                                tk.END,
                                f"  Fragment {i} ({len(frag)} bp): {frag}\n"
                            )

                self.varpro_results_text.insert(tk.END, "\n")

        except Exception as e:
            messagebox.showerror("Error", f"Analýza selhala: {e}")


    def save_varpro_results(self):
        try:
            results = self.varpro_results_text.get("1.0", tk.END)
            if not results.strip():
                messagebox.showwarning("Warning", "Žádné výsledky k uložení.")
                return
            save_path = filedialog.asksaveasfilename(
                defaultextension=".txt",
                filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
            )
            if save_path:
                with open(save_path, "w", encoding="utf-8") as f:
                    f.write(results)
                messagebox.showinfo("Success", "Výsledky byly úspěšně uloženy.")
        except Exception as e:
            messagebox.showerror("Error", f"Ukládání výsledků selhalo: {e}")

    def clear_varpro(self):
        self.varpro_email.set("")
        self.varpro_organism.set("")
        self.varpro_chromosome.set("")
        self.varpro_amplicon_pos.set("")
        self.varpro_snp_position.set("")
        self.varpro_assembly.set("")
        self.varpro_fasta_text.delete("1.0", tk.END)
        self.varpro_variant_file_text.delete("1.0", tk.END)
        self.varpro_results_text.delete("1.0", tk.END)

def main():
    root = tk.Tk()
    app = RENDO(root)
    root.mainloop()

if __name__ == "__main__":
    main()
