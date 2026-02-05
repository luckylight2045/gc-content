# src/gui.py

from analyzer import (
    calculate_gc,
    sliding_window_gc,
    nucleotide_counts,
    reverse_complement,
    validate_sequence,
    translate_dna,
    find_motif,
    gc_skew,
    find_orfs,
    segment_gc_profile,
)

import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog
from analyzer import calculate_gc, sliding_window_gc
from file_reader import read_fasta
import matplotlib.pyplot as plt
from analyzer import nucleotide_counts
import matplotlib.colors as mcolors

import csv

def export_report():
    seq = text_seq.get("1.0", tk.END).strip()
    if not seq:
        messagebox.showwarning("No sequence", "Paste a sequence or load a FASTA file.")
        return
    if not validate_sequence(seq):
        messagebox.showerror("Invalid Sequence", "Sequence contains invalid characters!")
        return

    counts = nucleotide_counts(seq)
    gc = calculate_gc(seq)
    protein = translate_dna(seq)

    path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files","*.csv")])
    if not path:
        return

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Metric","Value"])
        writer.writerow(["Length", counts['Length']])
        writer.writerow(["A", counts['A']])
        writer.writerow(["T", counts['T']])
        writer.writerow(["G", counts['G']])
        writer.writerow(["C", counts['C']])
        writer.writerow(["GC%", gc])
        writer.writerow(["Protein Translation", protein])

    messagebox.showinfo("Export Complete", f"Report saved to {path}")


def search_motif():
    seq = text_seq.get("1.0", tk.END).strip()
    if not seq:
        messagebox.showwarning("No sequence", "Paste a sequence or load a FASTA file.")
        return
    if not validate_sequence(seq):
        messagebox.showerror("Invalid Sequence", "Sequence contains invalid characters!")
        return
    motif = simpledialog.askstring("Input Motif", "Enter motif (use standard IUPAC bases):")
    if not motif:
        return
    positions = find_motif(seq, motif)
    top = tk.Toplevel(root)
    top.title(f"Motif: {motif}")
    text_motif = tk.Text(top, width=80, height=20)
    text_motif.pack(padx=10, pady=10)
    if not positions:
        text_motif.insert(tk.END, f"Motif '{motif}' not found.")
    else:
        text_motif.insert(tk.END, f"Motif '{motif}' found at positions:\n")
        text_motif.insert(tk.END, ", ".join(map(str, positions)))

def plot_gc_skew():
    seq = text_seq.get("1.0", tk.END).strip()
    if not seq:
        messagebox.showwarning("No sequence", "Paste a sequence or load a FASTA file.")
        return
    if not validate_sequence(seq):
        messagebox.showerror("Invalid Sequence", "Sequence contains invalid characters!")
        return
    data = gc_skew(seq)
    if not data:
        messagebox.showwarning("Sequence too short", "Sequence too short for sliding window GC skew.")
        return
    positions = [pos for pos, skew in data]
    skews = [skew for pos, skew in data]
    plt.figure()
    plt.plot(positions, skews, marker="o")
    plt.xlabel("Position (start of window)")
    plt.ylabel("GC Skew (G-C)/(G+C)")
    plt.title("GC Skew Plot")
    plt.show()


def show_orfs():
    seq = text_seq.get("1.0", tk.END).strip()
    if not seq:
        messagebox.showwarning("No sequence", "Paste a sequence or load a FASTA file.")
        return
    if not validate_sequence(seq):
        messagebox.showerror("Invalid Sequence", "Sequence contains invalid characters!")
        return
    orfs = find_orfs(seq)
    top = tk.Toplevel(root)
    top.title("ORF Finder")
    text_orfs = tk.Text(top, width=80, height=20)
    text_orfs.pack(padx=10, pady=10)
    if not orfs:
        text_orfs.insert(tk.END, "No ORFs found.")
    else:
        for idx, orf in enumerate(orfs, start=1):
            text_orfs.insert(tk.END, f"ORF {idx}: Start={orf['start']}, End={orf['end']}, Length={orf['length']}\n")


def plot_gc_heatmap(seq: str, window: int = 20, step: int = 5):
    data = sliding_window_gc(seq, window, step)
    if not data:
        messagebox.showwarning("Sequence too short", "Sequence shorter than window size.")
        return
    gcs = [gc for pos, gc in data]
    colors = ["red" if gc < 40 else "yellow" if gc < 60 else "green" for gc in gcs]
    plt.figure(figsize=(8,1))
    plt.bar(range(len(gcs)), [1]*len(gcs), color=colors, width=1)
    plt.title("GC Heatmap: Red=low, Yellow=medium, Green=high")
    plt.axis("off")
    plt.show()

def show_protein_translation():
    seq = text_seq.get("1.0", tk.END).strip()
    if not seq:
        messagebox.showwarning("No sequence", "Paste a sequence or load a FASTA file.")
        return
    if not validate_sequence(seq):
        messagebox.showerror("Invalid Sequence", "Sequence contains invalid characters! Only A,T,G,C,N allowed.")
        return
    protein = translate_dna(seq)
    top = tk.Toplevel(root)
    top.title("Protein Translation")
    text_protein = tk.Text(top, width=70, height=12)
    text_protein.pack(padx=10, pady=10)
    text_protein.insert(tk.END, protein)



def show_nucleotide_counts():
    seq = text_seq.get("1.0", tk.END).strip()
    if not seq:
        messagebox.showwarning("No sequence", "Paste a sequence or load a FASTA file.")
        return
    counts = nucleotide_counts(seq)
    msg = (
        f"A: {counts['A']}\n"
        f"T: {counts['T']}\n"
        f"G: {counts['G']}\n"
        f"C: {counts['C']}\n"
        f"Total Length: {counts['Length']}"
    )
    messagebox.showinfo("Nucleotide Counts", msg)

def show_reverse_complement():
    seq = text_seq.get("1.0", tk.END).strip()
    if not seq:
        messagebox.showwarning("No sequence", "Paste a sequence or load a FASTA file.")
        return
    rev_comp = reverse_complement(seq)
    # Show in a new window
    top = tk.Toplevel(root)
    top.title("Reverse Complement")
    text_rev = tk.Text(top, width=70, height=12)
    text_rev.pack(padx=10, pady=10)
    text_rev.insert(tk.END, rev_comp)



def load_file():
    path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa"), ("All files","*.*")])
    if not path:
        return
    seq = read_fasta(path)
    text_seq.delete("1.0", tk.END)
    text_seq.insert(tk.END, seq)

def compute_gc_show():
    seq = text_seq.get("1.0", tk.END).strip()
    if not seq:
        messagebox.showwarning("No sequence", "Paste a sequence or load a FASTA file.")
        return
    gc = calculate_gc(seq)
    messagebox.showinfo("GC Content", f"GC Content: {gc}%")

def plot_sliding_window():
    seq = text_seq.get("1.0", tk.END).strip()
    if not seq:
        messagebox.showwarning("No sequence", "Paste a sequence or load a FASTA file.")
        return
    if not validate_sequence(seq):
        messagebox.showerror("Invalid Sequence", "Sequence contains invalid characters! Only A,T,G,C,N allowed.")
        return
    window = 50  # or let user change later
    step = 10
    data = sliding_window_gc(seq, window=window, step=step)
    if not data:
        messagebox.showwarning("Sequence too short", "Sequence shorter than window size.")
        return
    positions = [pos for pos, gc in data]
    gcs = [gc for pos, gc in data]
    plt.figure()
    plt.plot(positions, gcs, marker="o")
    plt.xlabel("Position (start of window)")
    plt.ylabel("GC%")
    plt.title(f"Sliding Window GC% (window={window}, step={step})")
    plt.show()

def plot_gc_segmentation():
    seq = text_seq.get("1.0", tk.END).strip()
    if not seq:
        messagebox.showwarning("No sequence", "Paste a sequence or load a FASTA file.")
        return
    if not validate_sequence(seq):
        messagebox.showerror("Invalid Sequence", "Sequence contains invalid characters!")
        return

    # parameters (you can tune later)
    window = 50
    step = 10
    threshold = 1.0

    segments, breakpoints, gc_profile = segment_gc_profile(
        seq,
        window=window,
        step=step,
        threshold=threshold
    )

    if not gc_profile:
        messagebox.showwarning("Too short", "Sequence too short for segmentation.")
        return

    positions = [pos for pos, gc in gc_profile]
    gcs = [gc for pos, gc in gc_profile]

    plt.figure(figsize=(10, 4))
    plt.plot(positions, gcs, label="GC%")

    # draw breakpoints
    for bp in breakpoints:
        plt.axvline(bp, color="red", linestyle="--", alpha=0.6)

    plt.xlabel("Position (bp)")
    plt.ylabel("GC%")
    plt.title("GC Segmentation (Change-Point Based)")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # ---- show segments table ----
    top = tk.Toplevel(root)
    top.title("GC Segments")

    text_seg = tk.Text(top, width=80, height=20)
    text_seg.pack(padx=10, pady=10)

    if not segments:
        text_seg.insert(tk.END, "No significant segments detected.")
        return

    text_seg.insert(
        tk.END,
        "Start\tEnd\tLength\tMean GC%\n"
        "-------------------------------------------\n"
    )

    for s in segments:
        text_seg.insert(
            tk.END,
            f"{s['start']}\t{s['end']}\t{s['length']}\t{s['mean_gc']:.2f}\n"
        )



root = tk.Tk()
root.title("GC Content Analyzer")

frame = tk.Frame(root)
frame.pack(padx=10, pady=10)

text_seq = tk.Text(frame, width=70, height=12)
text_seq.grid(row=0, column=0, columnspan=3, padx=5, pady=5)

btn_load = tk.Button(frame, text="Load FASTA", command=load_file)
btn_load.grid(row=1, column=0, sticky="ew", padx=5, pady=5)

btn_calc = tk.Button(frame, text="Calculate GC", command=compute_gc_show)
btn_calc.grid(row=1, column=1, sticky="ew", padx=5, pady=5)

btn_plot = tk.Button(frame, text="Plot Sliding Window", command=plot_sliding_window)
btn_plot.grid(row=1, column=2, sticky="ew", padx=5, pady=5)

btn_counts = tk.Button(frame, text="Show Nucleotide Counts", command=show_nucleotide_counts)
btn_counts.grid(row=2, column=0, columnspan=3, sticky="ew", padx=5, pady=5)

btn_revcomp = tk.Button(frame, text="Reverse Complement", command=show_reverse_complement)
btn_revcomp.grid(row=3, column=0, columnspan=3, sticky="ew", padx=5, pady=5)

btn_heatmap = tk.Button(frame, text="GC Heatmap", command=lambda: plot_gc_heatmap(text_seq.get("1.0", tk.END).strip()))
btn_heatmap.grid(row=4, column=0, columnspan=3, sticky="ew", padx=5, pady=5)

btn_protein = tk.Button(frame, text="Translate DNA → Protein", command=show_protein_translation)
btn_protein.grid(row=5, column=0, columnspan=3, sticky="ew", padx=5, pady=5)

btn_orf = tk.Button(frame, text="Find ORFs", command=show_orfs)
btn_orf.grid(row=6, column=0, columnspan=3, sticky="ew", padx=5, pady=5)

btn_skew = tk.Button(frame, text="GC Skew Plot", command=plot_gc_skew)
btn_skew.grid(row=7, column=0, columnspan=3, sticky="ew", padx=5, pady=5)

btn_motif = tk.Button(frame, text="Find Motif", command=search_motif)
btn_motif.grid(row=8, column=0, columnspan=3, sticky="ew", padx=5, pady=5)

btn_export = tk.Button(frame, text="Export Report", command=export_report)
btn_export.grid(row=9, column=0, columnspan=3, sticky="ew", padx=5, pady=5)

btn_segment = tk.Button(
    frame,
    text="GC Segmentation",
    command=plot_gc_segmentation
)
btn_segment.grid(row=10, column=0, columnspan=3, sticky="ew", padx=5, pady=5)


root.mainloop()
