#!/usr/bin/env python3
# =============================================================================
# MK Flu-Pipe — Interface Gráfica
# Version: 1.0
# Developed by: Jean Phellipe Marques do Nascimento
# Genomic Surveillance Laboratory — LACEN/AL
# =============================================================================
#
# CHANGES v1.0 (GUI):
#   • Seletor "Tipo de dado" substituindo a inferência implícita de tecnologia:
#       - Automatic, Short reads paired, Short reads single, Long reads (ONT)
#       - Clarifica o uso para dados BGI, SRA, Ion Torrent sem ambiguidade
#   • Seção "Análises avançadas" com:
#       - Toggle iVar (variant calling short reads) + frequência e profundidade
#       - Toggle Medaka (variant calling ONT)
#       - Toggle análise de resistência a antivirais (FluSurver/WHO)
#       - Toggle análise de virulência H5 (condicional ao subtipo detectado)
#       - Parâmetros de co-infecção: freq. mínima alelo minoritário e % posições
#   • Seção "Submissão GISAID" com:
#       - Campo de localidade (ex: Brazil-AL)
#       - Campo de ano de coleta
#       - Prévia do formato do header: A/localidade/amostra/ano(HxNy)
#   • FLU-utr adicionado à lista de módulos no ComboBox
#   • Persistência de todas as novas configurações em ~/.mkflupipe_config.json
# =============================================================================
#   • Seção "Controle de Qualidade" com novos campos:
#       - Seletor de arquivo FASTA de adaptadores/primers (opcional)
#       - Toggle para ativar/desativar FastQC nos reads brutos
#       - Comprimento mínimo reads Illumina (fastp)
#       - Comprimento mínimo reads long-read/ONT (Filtlong)
#       - Comprimento máximo reads long-read (0 = sem limite)
#       - Qualidade mínima Phred (fastp / Filtlong)
#       - Toggle para depleção de hospedeiro humano
#   • Seção "QC pós-montagem" com campos:
#       - Cobertura mínima por segmento
#       - Percentual máximo de N
#       - Número mínimo de segmentos (PASS)
#   • Botão "Restaurar padrões" nos parâmetros de QC
#   • Persistência de todas as novas configurações em ~/.mkflupipe_config.json
# =============================================================================

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib, Pango

import threading
import subprocess
import os
import datetime
import json

CONFIG_FILE = os.path.join(os.path.expanduser("~"), ".mkflupipe_config.json")

# Valores padrão dos parâmetros de QC
QC_DEFAULTS = {
    "adapter_fasta":    "",
    "run_fastqc":       True,
    "min_len_illumina": "75",
    "min_len_long":     "200",
    "max_len_long":     "0",
    "min_qual":         "20",
    "host_depletion":   False,
    "min_coverage":     "50",
    "max_n_pct":        "10",
    "min_segments":     "4",
    # v1.0
    "seq_type":         "auto",
    "run_ivar":         False,
    "run_medaka":       False,
    "ivar_freq":        "0.03",
    "ivar_depth":       "10",
    "minority_freq":    "0.20",
    "coinfection_pct":  "5.0",
    "run_antiviral":    True,
    "run_h5_virulence": True,
    "run_fullvarcall":  False,
    "gisaid_location":  "",
    "gisaid_year":      str(datetime.date.today().year),
    "medaka_env":       "medaka_env",
}


def load_config() -> dict:
    try:
        with open(CONFIG_FILE, "r") as f:
            return json.load(f)
    except Exception:
        return {}


def save_config(data: dict):
    try:
        with open(CONFIG_FILE, "w") as f:
            json.dump(data, f, indent=2)
    except Exception:
        pass


class JanelaPrincipal(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self, title="MK Flu-Pipe by Jean Nascimento v1.0")
        self.set_border_width(12)
        self.set_default_size(980, 820)
        self.set_resizable(True)
        self.connect("delete-event", self._on_window_delete)

        self._processo = None
        self._executando = False
        self._config = load_config()

        # Layout principal com scroll
        vbox_outer = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=0)
        self.add(vbox_outer)

        scroll = Gtk.ScrolledWindow()
        scroll.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
        scroll.set_vexpand(True)
        vbox_outer.pack_start(scroll, True, True, 0)

        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
        vbox.set_margin_top(8)
        vbox.set_margin_bottom(8)
        vbox.set_margin_start(8)
        vbox.set_margin_end(8)
        scroll.add(vbox)

        # ── Cabeçalho ─────────────────────────────────────────────────────
        header_label = Gtk.Label()
        header_label.set_markup(
            "<span size='large' weight='bold'>MK Flu-Pipe</span>  "
            "<span size='small' foreground='gray'>v1.0</span>"
        )
        header_label.set_halign(Gtk.Align.START)
        vbox.pack_start(header_label, False, False, 0)

        separator = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)
        vbox.pack_start(separator, False, False, 2)

        # ── Seção: Main parameters ──────────────────────────────────
        frame_main = Gtk.Frame(label=" Main parameters ")
        frame_main.set_label_align(0.01, 0.5)
        vbox.pack_start(frame_main, False, False, 0)

        grid = Gtk.Grid(column_spacing=12, row_spacing=10)
        grid.set_margin_top(8)
        grid.set_margin_bottom(8)
        grid.set_margin_start(8)
        grid.set_margin_end(8)
        frame_main.add(grid)

        # Campo: Diretório de entrada
        lbl_in = Gtk.Label(label="Input directory:")
        lbl_in.set_halign(Gtk.Align.END)
        self.fc_input = Gtk.FileChooserButton(
            title="Select input directory",
            action=Gtk.FileChooserAction.SELECT_FOLDER)
        self.fc_input.set_tooltip_text(
            "Folder containing FASTQ files (.fastq.gz)\n"
            "Input mode is auto-detected."
        )
        self.fc_input.set_hexpand(True)
        if self._config.get("input_dir") and os.path.isdir(self._config["input_dir"]):
            self.fc_input.set_filename(self._config["input_dir"])

        # Campo: Diretório de saída
        lbl_out = Gtk.Label(label="Output directory:")
        lbl_out.set_halign(Gtk.Align.END)
        self.fc_output = Gtk.FileChooserButton(
            title="Select output directory",
            action=Gtk.FileChooserAction.SELECT_FOLDER)
        self.fc_output.set_tooltip_text(
            "Folder where results will be saved.\n"
            "Created automatically if it does not exist."
        )
        self.fc_output.set_hexpand(True)
        if self._config.get("output_dir") and os.path.isdir(self._config["output_dir"]):
            self.fc_output.set_filename(self._config["output_dir"])

        # IRMA module field
        lbl_mod = Gtk.Label(label="IRMA module:")
        lbl_mod.set_halign(Gtk.Align.END)
        self.combo_module = Gtk.ComboBoxText()
        IRMA_MODULES = ["FLU", "FLU-utr", "FLU-lowQC", "FLU_AD", "FLU-minion"]
        for mod in IRMA_MODULES:
            self.combo_module.append_text(mod)
        saved_mod = self._config.get("irma_module", "FLU")
        self.combo_module.set_active(IRMA_MODULES.index(saved_mod) if saved_mod in IRMA_MODULES else 0)
        self.combo_module.set_tooltip_text(
            "FLU        → General use (Illumina/BGI A+B)\n"
            "FLU-utr    → Includes UTR regions (more complete)\n"
            "FLU-lowQC  → Low-quality samples\n"
            "FLU_AD     → Also detects C and D\n"
            "FLU-minion → Oxford Nanopore data (ONT)"
        )

        # Data type field — new in v1.0
        lbl_seqtype = Gtk.Label(label="Data type:")
        lbl_seqtype.set_halign(Gtk.Align.END)
        self.combo_seqtype = Gtk.ComboBoxText()
        SEQ_TYPES = ["auto", "short_paired", "short_single", "long"]
        SEQ_TYPE_LABELS = [
            "Automatic",
            "Short reads paired  (Illumina/BGI — SRA _1/_2 or R1/R2)",
            "Short reads single  (Illumina/BGI/Ion Torrent)",
            "Long reads (ONT — use with FLU-minion)",
        ]
        for lbl in SEQ_TYPE_LABELS:
            self.combo_seqtype.append_text(lbl)
        saved_seqtype = self._config.get("seq_type", "auto")
        self.combo_seqtype.set_active(
            SEQ_TYPES.index(saved_seqtype) if saved_seqtype in SEQ_TYPES else 0
        )
        self.combo_seqtype.set_tooltip_text(
            "Automatic       → Auto-detects by file naming pattern\n"
            "short_paired     → Illumina/BGI paired; SRA (_1/_2); BaseSpace (R1/R2)\n"
            "short_single     → Illumina/BGI/Ion Torrent single-end\n"
            "long             → Oxford Nanopore (ONT) — use with FLU-minion module\n\n"
            "⚠ BGI data has the same format as Illumina — select short_paired."
        )
        self._SEQ_TYPES = SEQ_TYPES  # usado em _get_params

        # Campo: Modo de entrada (padrão de arquivo — mantido para retrocompatibilidade)
        lbl_mode = Gtk.Label(label="File mode:")
        lbl_mode.set_halign(Gtk.Align.END)
        self.combo_mode = Gtk.ComboBoxText()
        FILE_MODES = ["Automatic", "illumina_paired", "sra_paired", "generic_paired", "single"]
        for mode in FILE_MODES:
            self.combo_mode.append_text(mode)
        saved_mode = self._config.get("seq_mode", "Automatic")
        self.combo_mode.set_active(FILE_MODES.index(saved_mode) if saved_mode in FILE_MODES else 0)
        self.combo_mode.set_tooltip_text(
            "Automatic       → Auto-detects by file naming pattern (recommended)\n"
            "illumina_paired  → *_L001_R1_001.fastq.gz (BaseSpace)\n"
            "sra_paired       → *_1.fastq.gz + *_2.fastq.gz (SRA)\n"
            "generic_paired   → *_R1_*.fastq.gz (Illumina/BGI generic)\n"
            "single           → *.fastq.gz (ONT, Ion Torrent, PacBio, SRA single)"
        )

        grid.attach(lbl_in,              0, 0, 1, 1)
        grid.attach(self.fc_input,       1, 0, 3, 1)
        grid.attach(lbl_out,             0, 1, 1, 1)
        grid.attach(self.fc_output,      1, 1, 3, 1)
        grid.attach(lbl_mod,             0, 2, 1, 1)
        grid.attach(self.combo_module,   1, 2, 1, 1)
        grid.attach(lbl_seqtype,         2, 2, 1, 1)
        grid.attach(self.combo_seqtype,  3, 2, 1, 1)
        grid.attach(lbl_mode,            0, 3, 1, 1)
        grid.attach(self.combo_mode,     1, 3, 3, 1)

        # ── Seção: Controle de Qualidade pré-IRMA ─────────────────────────
        frame_qc = Gtk.Frame(label=" Pre-IRMA quality control ")
        frame_qc.set_label_align(0.01, 0.5)
        vbox.pack_start(frame_qc, False, False, 0)

        grid_qc = Gtk.Grid(column_spacing=12, row_spacing=8)
        grid_qc.set_margin_top(8)
        grid_qc.set_margin_bottom(8)
        grid_qc.set_margin_start(8)
        grid_qc.set_margin_end(8)
        frame_qc.add(grid_qc)

        # FastQC toggle
        lbl_fqc = Gtk.Label(label="Raw reads FastQC:")
        lbl_fqc.set_halign(Gtk.Align.END)
        self.chk_fastqc = Gtk.CheckButton(label="Enable FastQC diagnostic before trimming")
        self.chk_fastqc.set_active(self._config.get("run_fastqc", QC_DEFAULTS["run_fastqc"]))
        self.chk_fastqc.set_tooltip_text(
            "Generates quality report for raw reads (FastQC + MultiQC)\n"
            "Useful for identifying quality issues before analysis."
        )

        # Arquivo de adaptadores
        lbl_adp = Gtk.Label(label="Adapters/primers:")
        lbl_adp.set_halign(Gtk.Align.END)

        adp_box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        self.fc_adapter = Gtk.FileChooserButton(
            title="Select FASTA file for adapters",
            action=Gtk.FileChooserAction.OPEN)
        fasta_filter = Gtk.FileFilter()
        fasta_filter.set_name("FASTA files (*.fasta, *.fa, *.fna)")
        fasta_filter.add_pattern("*.fasta")
        fasta_filter.add_pattern("*.fa")
        fasta_filter.add_pattern("*.fna")
        self.fc_adapter.add_filter(fasta_filter)
        self.fc_adapter.set_tooltip_text(
             "FASTA file with adapters and primers for removal by fastp.\n"
             "If not provided, fastp auto-detects adapters.\n\n"
            "Seu arquivo inclui: PrefixNX, Opti1, primers universais Flu A/B,\n"
             "12U/13U primers and the full B-Uni primer collection."
        )
        saved_adapter = self._config.get("adapter_fasta", "")
        if saved_adapter and os.path.isfile(saved_adapter):
            self.fc_adapter.set_filename(saved_adapter)
        self.fc_adapter.set_hexpand(True)

        btn_clear_adp = Gtk.Button(label="Clear")
        btn_clear_adp.set_tooltip_text("Remove adapter file (use auto-detection)")
        btn_clear_adp.connect("clicked", self._on_clear_adapter)

        adp_box.pack_start(self.fc_adapter, True, True, 0)
        adp_box.pack_start(btn_clear_adp, False, False, 0)

        # Min len Illumina
        lbl_ml_il = Gtk.Label(label="Min. length Illumina (nt):")
        lbl_ml_il.set_halign(Gtk.Align.END)
        self.spin_min_len_illumina = Gtk.SpinButton()
        adj_ml_il = Gtk.Adjustment(value=float(self._config.get("min_len_illumina", QC_DEFAULTS["min_len_illumina"])),
                                   lower=30, upper=300, step_increment=5, page_increment=25)
        self.spin_min_len_illumina.set_adjustment(adj_ml_il)
        self.spin_min_len_illumina.set_digits(0)
        self.spin_min_len_illumina.set_tooltip_text(
            "Minimum read length for Illumina after trimming (fastp).\n"
            "Shorter reads are discarded. Default: 75 nt."
        )

        # Min len long-read
        lbl_ml_lr = Gtk.Label(label="Min. length long-read (nt):")
        lbl_ml_lr.set_halign(Gtk.Align.END)
        self.spin_min_len_long = Gtk.SpinButton()
        adj_ml_lr = Gtk.Adjustment(value=float(self._config.get("min_len_long", QC_DEFAULTS["min_len_long"])),
                                   lower=50, upper=2000, step_increment=50, page_increment=200)
        self.spin_min_len_long.set_adjustment(adj_ml_lr)
        self.spin_min_len_long.set_digits(0)
        self.spin_min_len_long.set_tooltip_text(
            "Minimum read length for ONT/PacBio reads (Filtlong).\n"
            "Shorter reads are discarded. Default: 200 nt."
        )

        # Max len long-read
        lbl_maxl = Gtk.Label(label="Max. length long-read (nt):")
        lbl_maxl.set_halign(Gtk.Align.END)
        self.spin_max_len_long = Gtk.SpinButton()
        adj_maxl = Gtk.Adjustment(value=float(self._config.get("max_len_long", QC_DEFAULTS["max_len_long"])),
                                  lower=0, upper=100000, step_increment=100, page_increment=1000)
        self.spin_max_len_long.set_adjustment(adj_maxl)
        self.spin_max_len_long.set_digits(0)
        self.spin_max_len_long.set_tooltip_text(
            "Maximum length for long reads (Filtlong).\n"
            "0 = no limit (recommended for influenza)."
        )

        # Min qual
        lbl_mq = Gtk.Label(label="Minimum quality (Phred):")
        lbl_mq.set_halign(Gtk.Align.END)
        self.spin_min_qual = Gtk.SpinButton()
        adj_mq = Gtk.Adjustment(value=float(self._config.get("min_qual", QC_DEFAULTS["min_qual"])),
                                lower=5, upper=40, step_increment=1, page_increment=5)
        self.spin_min_qual.set_adjustment(adj_mq)
        self.spin_min_qual.set_digits(0)
        self.spin_min_qual.set_tooltip_text(
            "Minimum Phred quality for base trimming at read ends.\n"
            "Q20 = 99% accuracy (default). Q30 = 99.9% (more stringent)."
        )

        # Host depletion toggle
        lbl_dep = Gtk.Label(label="Host depletion:")
        lbl_dep.set_halign(Gtk.Align.END)
        self.chk_host_dep = Gtk.CheckButton(
            label="Remover reads humanos antes do IRMA (Bowtie2 / minimap2)")
        self.chk_host_dep.set_active(self._config.get("host_depletion", QC_DEFAULTS["host_depletion"]))
        self.chk_host_dep.set_tooltip_text(
            "Remove reads do hospedeiro humano (GRCh38) antes da montagem.\n"
            "Recommended for clinical samples (nasopharyngeal swab).\n"
            "Illumina/Ion Torrent: usa Bowtie2.\n"
            "ONT/PacBio (FLU-minion): usa minimap2.\n\n"
            "⚠ On first run, the human genome (~3 GB) will be\n"
            "downloaded and indexed automatically (~30-60 min)."
        )

        # Adiciona campos ao grid_qc
        grid_qc.attach(lbl_fqc,                    0, 0, 1, 1)
        grid_qc.attach(self.chk_fastqc,             1, 0, 3, 1)
        grid_qc.attach(lbl_adp,                     0, 1, 1, 1)
        grid_qc.attach(adp_box,                     1, 1, 3, 1)
        grid_qc.attach(lbl_ml_il,                   0, 2, 1, 1)
        grid_qc.attach(self.spin_min_len_illumina,  1, 2, 1, 1)
        grid_qc.attach(lbl_ml_lr,                   2, 2, 1, 1)
        grid_qc.attach(self.spin_min_len_long,      3, 2, 1, 1)
        grid_qc.attach(lbl_maxl,                    0, 3, 1, 1)
        grid_qc.attach(self.spin_max_len_long,      1, 3, 1, 1)
        grid_qc.attach(lbl_mq,                      2, 3, 1, 1)
        grid_qc.attach(self.spin_min_qual,          3, 3, 1, 1)
        grid_qc.attach(lbl_dep,                     0, 4, 1, 1)
        grid_qc.attach(self.chk_host_dep,           1, 4, 3, 1)

        # ── Seção: QC pós-montagem ─────────────────────────────────────────
        frame_postqc = Gtk.Frame(label=" Post-assembly QC (acceptance filters) ")
        frame_postqc.set_label_align(0.01, 0.5)
        vbox.pack_start(frame_postqc, False, False, 0)

        grid_postqc = Gtk.Grid(column_spacing=12, row_spacing=8)
        grid_postqc.set_margin_top(8)
        grid_postqc.set_margin_bottom(8)
        grid_postqc.set_margin_start(8)
        grid_postqc.set_margin_end(8)
        frame_postqc.add(grid_postqc)

        # Cobertura mínima
        lbl_cov = Gtk.Label(label="Minimum coverage (×):")
        lbl_cov.set_halign(Gtk.Align.END)
        self.spin_min_cov = Gtk.SpinButton()
        adj_cov = Gtk.Adjustment(value=float(self._config.get("min_coverage", QC_DEFAULTS["min_coverage"])),
                                 lower=1, upper=10000, step_increment=10, page_increment=100)
        self.spin_min_cov.set_adjustment(adj_cov)
        self.spin_min_cov.set_digits(0)
        self.spin_min_cov.set_tooltip_text(
            "Minimum mean coverage per assembled segment.\n"
            "Segmentos abaixo deste valor recebem alerta WARN.\n"
            "Default: 50×."
        )

        # % max N
        lbl_npct = Gtk.Label(label="Maximum % N (%):")
        lbl_npct.set_halign(Gtk.Align.END)
        self.spin_max_n = Gtk.SpinButton()
        adj_npct = Gtk.Adjustment(value=float(self._config.get("max_n_pct", QC_DEFAULTS["max_n_pct"])),
                                  lower=0, upper=100, step_increment=1, page_increment=5)
        self.spin_max_n.set_adjustment(adj_npct)
        self.spin_max_n.set_digits(1)
        self.spin_max_n.set_tooltip_text(
            "Maximum percentage of N bases allowed in an assembled segment.\n"
            "Segmentos acima deste valor recebem alerta WARN.\n"
            "Default: 10%."
        )

        # Min segmentos
        lbl_minseg = Gtk.Label(label="Minimum segments (PASS):")
        lbl_minseg.set_halign(Gtk.Align.END)
        self.spin_min_segs = Gtk.SpinButton()
        adj_minseg = Gtk.Adjustment(value=float(self._config.get("min_segments", QC_DEFAULTS["min_segments"])),
                                    lower=1, upper=8, step_increment=1, page_increment=1)
        self.spin_min_segs.set_adjustment(adj_minseg)
        self.spin_min_segs.set_digits(0)
        self.spin_min_segs.set_tooltip_text(
            "Minimum number of assembled segments for PASS status.\n"
            "Amostras com menos segmentos recebem status FAIL.\n"
            "Default: 4 (half-genome). 8 = full genome required."
        )

        # Botão restaurar padrões
        btn_reset = Gtk.Button(label="↺ Restore QC defaults")
        btn_reset.set_tooltip_text("Restores all QC parameters to default values")
        btn_reset.connect("clicked", self._on_reset_qc_defaults)

        grid_postqc.attach(lbl_cov,           0, 0, 1, 1)
        grid_postqc.attach(self.spin_min_cov, 1, 0, 1, 1)
        grid_postqc.attach(lbl_npct,          2, 0, 1, 1)
        grid_postqc.attach(self.spin_max_n,   3, 0, 1, 1)
        grid_postqc.attach(lbl_minseg,        0, 1, 1, 1)
        grid_postqc.attach(self.spin_min_segs,1, 1, 1, 1)
        grid_postqc.attach(btn_reset,         3, 1, 1, 1)

        # ── Seção: Advanced analyses (v1.0) ────────────────────────────
        frame_adv = Gtk.Frame(label=" Advanced analyses ")
        frame_adv.set_label_align(0.01, 0.5)
        vbox.pack_start(frame_adv, False, False, 0)

        grid_adv = Gtk.Grid(column_spacing=12, row_spacing=8)
        grid_adv.set_margin_top(8)
        grid_adv.set_margin_bottom(8)
        grid_adv.set_margin_start(8)
        grid_adv.set_margin_end(8)
        frame_adv.add(grid_adv)

        # iVar variant calling
        lbl_ivar = Gtk.Label(label="Variant calling (short):")
        lbl_ivar.set_halign(Gtk.Align.END)
        self.chk_ivar = Gtk.CheckButton(label="iVar — variant calling in short reads (Illumina/BGI)")
        self.chk_ivar.set_active(self._config.get("run_ivar", QC_DEFAULTS["run_ivar"]))
        self.chk_ivar.set_tooltip_text(
            "Uses BAMs generated by IRMA to call variants with iVar.\n"
            "Complements IRMA internal variant calling with configurable\n"
            "frequency and depth thresholds.\n"
            "Requer: iVar instalado no ambiente mk_flu."
        )

        lbl_ivar_freq = Gtk.Label(label="Min. freq. iVar:")
        lbl_ivar_freq.set_halign(Gtk.Align.END)
        self.spin_ivar_freq = Gtk.SpinButton()
        adj_ivf = Gtk.Adjustment(
            value=float(self._config.get("ivar_freq", QC_DEFAULTS["ivar_freq"])),
            lower=0.01, upper=1.0, step_increment=0.01, page_increment=0.05)
        self.spin_ivar_freq.set_adjustment(adj_ivf)
        self.spin_ivar_freq.set_digits(2)
        self.spin_ivar_freq.set_tooltip_text("Minimum alternative allele frequency (default: 0.03 = 3%)")

        lbl_ivar_depth = Gtk.Label(label="Min. depth iVar:")
        lbl_ivar_depth.set_halign(Gtk.Align.END)
        self.spin_ivar_depth = Gtk.SpinButton()
        adj_ivd = Gtk.Adjustment(
            value=float(self._config.get("ivar_depth", QC_DEFAULTS["ivar_depth"])),
            lower=1, upper=10000, step_increment=5, page_increment=50)
        self.spin_ivar_depth.set_adjustment(adj_ivd)
        self.spin_ivar_depth.set_digits(0)
        self.spin_ivar_depth.set_tooltip_text("Minimum coverage depth to call variant (default: 10)")

        # Medaka variant calling
        lbl_medaka = Gtk.Label(label="Variant calling (ONT):")
        lbl_medaka.set_halign(Gtk.Align.END)
        self.chk_medaka = Gtk.CheckButton(label="Medaka — variant calling in long reads (ONT)")
        self.chk_medaka.set_active(self._config.get("run_medaka", QC_DEFAULTS["run_medaka"]))
        self.chk_medaka.set_tooltip_text(
            "Usa os BAMs do IRMA para polishing e variant calling com Medaka.\n"
            "Applicable only when using FLU-minion (ONT long reads).\n"
            "Requer: medaka instalado no ambiente conda separado abaixo."
        )

        lbl_medaka_env = Gtk.Label(label="Medaka conda environment:")
        lbl_medaka_env.set_halign(Gtk.Align.END)
        self.entry_medaka_env = Gtk.Entry()
        self.entry_medaka_env.set_text(
            self._config.get("medaka_env", QC_DEFAULTS["medaka_env"]))
        self.entry_medaka_env.set_hexpand(True)
        self.entry_medaka_env.set_tooltip_text(
            "Name of the conda environment where Medaka is installed.\n"
            "Default: medaka_env\n\n"
            "Medaka has dependencies incompatible with mk_flu,\n"
            "so it requires a separate environment.\n\n"
            "Para criar: conda create -n medaka_env -c conda-forge -c bioconda medaka"
        )

        # Co-infecção
        lbl_coinf = Gtk.Label(label="Co-infection detection:")
        lbl_coinf.set_halign(Gtk.Align.END)
        lbl_coinf_info = Gtk.Label(
            label="Min. minority allele freq.:")
        lbl_coinf_info.set_halign(Gtk.Align.END)
        self.spin_minor_freq = Gtk.SpinButton()
        adj_mf = Gtk.Adjustment(
            value=float(self._config.get("minority_freq", QC_DEFAULTS["minority_freq"])),
            lower=0.05, upper=0.50, step_increment=0.01, page_increment=0.05)
        self.spin_minor_freq.set_adjustment(adj_mf)
        self.spin_minor_freq.set_digits(2)
        self.spin_minor_freq.set_tooltip_text(
             "Minimum minority allele frequency to consider\n"
             "a position as bimodal (default: 0.20 = 20%).\n"
            "Analysis uses IRMA allAlleles.txt."
        )

        lbl_coinf_pct = Gtk.Label(label="% bimodal positions for alert:")
        lbl_coinf_pct.set_halign(Gtk.Align.END)
        self.spin_coinf_pct = Gtk.SpinButton()
        adj_cp = Gtk.Adjustment(
            value=float(self._config.get("coinfection_pct", QC_DEFAULTS["coinfection_pct"])),
            lower=0.1, upper=100.0, step_increment=0.5, page_increment=5.0)
        self.spin_coinf_pct.set_adjustment(adj_cp)
        self.spin_coinf_pct.set_digits(1)
        self.spin_coinf_pct.set_tooltip_text(
            "Minimum percentage of bimodal positions in a segment\n"
            "to trigger a possible co-infection alert (default: 5.0%)."
        )

        # Resistência a antivirais
        lbl_antiv = Gtk.Label(label="Antiviral resistance:")
        lbl_antiv.set_halign(Gtk.Align.END)
        self.chk_antiviral = Gtk.CheckButton(
            label="Analyze resistance mutations (FluSurver/WHO: oseltamivir, baloxavir…)")
        self.chk_antiviral.set_active(self._config.get("run_antiviral", QC_DEFAULTS["run_antiviral"]))
        self.chk_antiviral.set_tooltip_text(
            "Uses samtools mpileup on IRMA BAMs to check\n"
            "resistance positions curated by FluSurver/WHO:\n"
            "  • NA H275Y, E119V, R292K, N294S (oseltamivir)\n"
            "  • PA I38T, I38M (baloxavir)\n"
            "  • M2 S31N, A30T (amantadinas)\n"
            "Database created automatically on first run."
        )

        # Virulência H5
        lbl_h5 = Gtk.Label(label="H5 virulence:")
        lbl_h5.set_halign(Gtk.Align.END)
        self.chk_h5 = Gtk.CheckButton(
            label="Analyze H5 virulence markers (enabled only if H5 detected)")
        self.chk_h5.set_active(self._config.get("run_h5_virulence", QC_DEFAULTS["run_h5_virulence"]))
        self.chk_h5.set_tooltip_text(
            "If BLAST/Nextclade identifies H5, automatically checks:\n"
            "  • PB2 E627K and D701N (mammalian adaptation)\n"
            "  • PA I97V (virulence)\n"
            "  • PB1-F2 N66S (pathogenicity)\n"
            "  • NS1 D92E (interferon resistance)\n"
            "Analysis is conditional — does not run if H5 is absent."
        )

        # Full genome variant calling (RefSeq + GFF3)
        lbl_fvc = Gtk.Label(label="Full genome varcall:")
        lbl_fvc.set_halign(Gtk.Align.END)
        self.chk_fullvarcall = Gtk.CheckButton(
            label="All 8 segments — NCBI RefSeq per segment + GFF3 protein annotation")
        self.chk_fullvarcall.set_active(
            self._config.get("run_fullvarcall", QC_DEFAULTS["run_fullvarcall"]))
        self.chk_fullvarcall.set_tooltip_text(
            "For each assembled segment, downloads the best NCBI RefSeq + GFF3,\n"
            "re-aligns reads with minimap2, calls variants with iVar, and\n"
            "translates nt variants into protein mutations (e.g. K189Q).\n\n"
            "Works for Flu A (H1N1, H3N2, H5N1) and Flu B (Victoria/Yamagata).\n"
            "Compatible with short reads (Illumina/BGI) and long reads (ONT).\n\n"
            "Output:\n"
            "  full_variant_calls/<sample>/*_protein_mutations.tsv\n"
            "  full_variant_calls/all_samples_protein_mutations.tsv\n"
            "  Included in surveillance_report.html (Protein mutations tab)\n\n"
            "⚠ Requires internet access on first run (downloads RefSeq + GFF3).\n"
            "   Cached locally after first download.\n"
            "⚠ Adds ~5–15 min per sample depending on coverage depth.\n"
            "   Requires: iVar + minimap2 + samtools in the mk_flu environment."
        )

        # Layout da seção avançada
        grid_adv.attach(lbl_ivar,             0, 0, 1, 1)
        grid_adv.attach(self.chk_ivar,        1, 0, 3, 1)
        grid_adv.attach(lbl_ivar_freq,        0, 1, 1, 1)
        grid_adv.attach(self.spin_ivar_freq,  1, 1, 1, 1)
        grid_adv.attach(lbl_ivar_depth,       2, 1, 1, 1)
        grid_adv.attach(self.spin_ivar_depth, 3, 1, 1, 1)
        grid_adv.attach(lbl_medaka,           0, 2, 1, 1)
        grid_adv.attach(self.chk_medaka,      1, 2, 3, 1)
        grid_adv.attach(lbl_medaka_env,       0, 3, 1, 1)
        grid_adv.attach(self.entry_medaka_env,1, 3, 3, 1)
        grid_adv.attach(lbl_coinf,            0, 4, 1, 1)
        grid_adv.attach(lbl_coinf_info,       0, 5, 1, 1)
        grid_adv.attach(self.spin_minor_freq, 1, 5, 1, 1)
        grid_adv.attach(lbl_coinf_pct,        2, 5, 1, 1)
        grid_adv.attach(self.spin_coinf_pct,  3, 5, 1, 1)
        grid_adv.attach(lbl_antiv,            0, 6, 1, 1)
        grid_adv.attach(self.chk_antiviral,   1, 6, 3, 1)
        grid_adv.attach(lbl_h5,               0, 7, 1, 1)
        grid_adv.attach(self.chk_h5,          1, 7, 3, 1)
        grid_adv.attach(lbl_fvc,              0, 8, 1, 1)
        grid_adv.attach(self.chk_fullvarcall, 1, 8, 3, 1)

        # ── Seção: Submissão GISAID (v1.0) ──────────────────────────────
        frame_gisaid = Gtk.Frame(label=" GISAID EpiFlu submission ")
        frame_gisaid.set_label_align(0.01, 0.5)
        vbox.pack_start(frame_gisaid, False, False, 0)

        grid_gisaid = Gtk.Grid(column_spacing=12, row_spacing=8)
        grid_gisaid.set_margin_top(8)
        grid_gisaid.set_margin_bottom(8)
        grid_gisaid.set_margin_start(8)
        grid_gisaid.set_margin_end(8)
        frame_gisaid.add(grid_gisaid)

        lbl_gloc = Gtk.Label(label="Location:")
        lbl_gloc.set_halign(Gtk.Align.END)
        self.entry_gisaid_loc = Gtk.Entry()
        self.entry_gisaid_loc.set_placeholder_text("ex: Brazil-AL")
        self.entry_gisaid_loc.set_text(self._config.get("gisaid_location", ""))
        self.entry_gisaid_loc.set_hexpand(True)
        self.entry_gisaid_loc.set_tooltip_text(
            "Location for the GISAID header. Examples:\n"
            "  Brazil-AL  → A/Brazil-AL/AMOSTRA/2025(H3N2)\n"
            "  Alagoas    → A/Alagoas/AMOSTRA/2025(H1N1)\n"
            "If empty, GISAID-ready generation will be skipped."
        )
        self.entry_gisaid_loc.connect("changed", self._on_gisaid_preview_changed)

        lbl_gyear = Gtk.Label(label="Collection year:")
        lbl_gyear.set_halign(Gtk.Align.END)
        self.spin_gisaid_year = Gtk.SpinButton()
        adj_gy = Gtk.Adjustment(
            value=float(self._config.get("gisaid_year", datetime.date.today().year)),
            lower=2000, upper=2100, step_increment=1, page_increment=1)
        self.spin_gisaid_year.set_adjustment(adj_gy)
        self.spin_gisaid_year.set_digits(0)
        self.spin_gisaid_year.set_tooltip_text("Sample collection year for the GISAID header")
        self.spin_gisaid_year.connect("value-changed", self._on_gisaid_preview_changed)

        lbl_prev = Gtk.Label(label="Header preview:")
        lbl_prev.set_halign(Gtk.Align.END)
        self.lbl_gisaid_preview = Gtk.Label(label="A/Brazil-AL/AMOSTRA/2025(H3N2)")
        self.lbl_gisaid_preview.set_halign(Gtk.Align.START)
        self.lbl_gisaid_preview.set_markup(
            "<span font_family='monospace' foreground='#666'>A/location/SAMPLE/year(HxNy)</span>"
        )

        grid_gisaid.attach(lbl_gloc,                   0, 0, 1, 1)
        grid_gisaid.attach(self.entry_gisaid_loc,       1, 0, 1, 1)
        grid_gisaid.attach(lbl_gyear,                  2, 0, 1, 1)
        grid_gisaid.attach(self.spin_gisaid_year,       3, 0, 1, 1)
        grid_gisaid.attach(lbl_prev,                   0, 1, 1, 1)
        grid_gisaid.attach(self.lbl_gisaid_preview,    1, 1, 3, 1)

        # Botão restaurar padrões
        btn_reset = Gtk.Button(label="↺ Restore all QC and analysis defaults")
        btn_reset.set_tooltip_text("Restores all parameters to default values")
        btn_reset.connect("clicked", self._on_reset_qc_defaults)
        vbox.pack_start(btn_reset, False, False, 2)

        sep2 = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)
        vbox.pack_start(sep2, False, False, 2)

        btn_box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
        vbox.pack_start(btn_box, False, False, 0)

        self.btn_run = Gtk.Button(label="▶ Run Pipeline")
        self.btn_run.set_tooltip_text("Starts full pipeline execution")
        self.btn_run.connect("clicked", self._on_run_clicked)

        self.btn_stop = Gtk.Button(label="⏹ Stop")
        self.btn_stop.set_tooltip_text("Interrupts execution (can be resumed later)")
        self.btn_stop.set_sensitive(False)
        self.btn_stop.connect("clicked", self._on_stop_clicked)

        self.btn_save_log = Gtk.Button(label="💾 Save Log")
        self.btn_save_log.set_tooltip_text("Saves the execution log to a .txt file")
        self.btn_save_log.connect("clicked", self._on_save_log)

        self.btn_clear_log = Gtk.Button(label="🗑 Clear Log")
        self.btn_clear_log.connect("clicked", self._on_clear_log)

        btn_box.pack_start(self.btn_run,      False, False, 0)
        btn_box.pack_start(self.btn_stop,     False, False, 0)
        btn_box.pack_start(self.btn_save_log, False, False, 0)
        btn_box.pack_start(self.btn_clear_log, False, False, 0)

        # ── Área de log ───────────────────────────────────────────────────
        log_frame = Gtk.Frame(label=" Execution log ")
        log_frame.set_label_align(0.01, 0.5)
        vbox_outer.pack_start(log_frame, True, True, 4)

        log_scroll = Gtk.ScrolledWindow()
        log_scroll.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        log_scroll.set_min_content_height(220)
        log_frame.add(log_scroll)

        self.txt_log = Gtk.TextView()
        self.txt_log.set_editable(False)
        self.txt_log.set_monospace(True)
        self.txt_log.set_wrap_mode(Gtk.WrapMode.WORD_CHAR)
        mono_font = Pango.FontDescription("Monospace 9")
        self.txt_log.override_font(mono_font)
        log_scroll.add(self.txt_log)

        self.show_all()

    # ── Helpers ───────────────────────────────────────────────────────────

    def _log(self, text: str):
        """Appends text to the log area (thread-safe)."""
        def _append():
            buf = self.txt_log.get_buffer()
            buf.insert(buf.get_end_iter(), text + "\n")
            adj = self.txt_log.get_parent().get_vadjustment()
            adj.set_value(adj.get_upper() - adj.get_page_size())
            return False
        GLib.idle_add(_append)

    def _set_ui_running(self, running: bool):
        def _set():
            self.btn_run.set_sensitive(not running)
            self.btn_stop.set_sensitive(running)
            return False
        GLib.idle_add(_set)

    def _on_clear_adapter(self, _btn):
        """Clears the adapter file selection."""
        self.fc_adapter.unselect_all()

    def _on_gisaid_preview_changed(self, _widget):
        """Updates the GISAID header preview in real time."""
        loc = self.entry_gisaid_loc.get_text().strip() or "localidade"
        year = int(self.spin_gisaid_year.get_value())
        preview = f"A/{loc}/AMOSTRA/{year}(H3N2)   |   B/{loc}/AMOSTRA/{year}"
        self.lbl_gisaid_preview.set_markup(
            f"<span font_family='monospace' foreground='#448'>{preview}</span>"
        )

    def _on_reset_qc_defaults(self, _btn):
        """Restores all QC and analysis parameters to default values."""
        self.chk_fastqc.set_active(QC_DEFAULTS["run_fastqc"])
        self.fc_adapter.unselect_all()
        self.spin_min_len_illumina.set_value(float(QC_DEFAULTS["min_len_illumina"]))
        self.spin_min_len_long.set_value(float(QC_DEFAULTS["min_len_long"]))
        self.spin_max_len_long.set_value(float(QC_DEFAULTS["max_len_long"]))
        self.spin_min_qual.set_value(float(QC_DEFAULTS["min_qual"]))
        self.chk_host_dep.set_active(QC_DEFAULTS["host_depletion"])
        self.spin_min_cov.set_value(float(QC_DEFAULTS["min_coverage"]))
        self.spin_max_n.set_value(float(QC_DEFAULTS["max_n_pct"]))
        self.spin_min_segs.set_value(float(QC_DEFAULTS["min_segments"]))
        # v1.0
        self.chk_ivar.set_active(QC_DEFAULTS["run_ivar"])
        self.chk_medaka.set_active(QC_DEFAULTS["run_medaka"])
        self.entry_medaka_env.set_text(QC_DEFAULTS["medaka_env"])
        self.spin_ivar_freq.set_value(float(QC_DEFAULTS["ivar_freq"]))
        self.spin_ivar_depth.set_value(float(QC_DEFAULTS["ivar_depth"]))
        self.spin_minor_freq.set_value(float(QC_DEFAULTS["minority_freq"]))
        self.spin_coinf_pct.set_value(float(QC_DEFAULTS["coinfection_pct"]))
        self.chk_antiviral.set_active(QC_DEFAULTS["run_antiviral"])
        self.chk_h5.set_active(QC_DEFAULTS["run_h5_virulence"])
        self.chk_fullvarcall.set_active(QC_DEFAULTS["run_fullvarcall"])
        self.entry_gisaid_loc.set_text(QC_DEFAULTS["gisaid_location"])
        self.spin_gisaid_year.set_value(float(QC_DEFAULTS["gisaid_year"]))
        self._log("[GUI] All parameters restored to default values.")

    # ── Leitura de parâmetros da GUI ──────────────────────────────────────

    def _get_params(self):
        """Collects all parameters from the interface."""
        input_dir  = self.fc_input.get_filename() or ""
        output_dir = self.fc_output.get_filename() or ""
        irma_mod   = self.combo_module.get_active_text() or "FLU"

        # Data type v1.0
        seqtype_idx = self.combo_seqtype.get_active()
        seq_type = self._SEQ_TYPES[seqtype_idx] if 0 <= seqtype_idx < len(self._SEQ_TYPES) else "auto"

        seq_mode_txt = self.combo_mode.get_active_text() or "Automatic"
        seq_mode   = "" if seq_mode_txt == "Automatic" else seq_mode_txt

        adapter_fasta  = self.fc_adapter.get_filename() or ""
        run_fastqc     = "yes" if self.chk_fastqc.get_active() else "no"
        min_len_il     = int(self.spin_min_len_illumina.get_value())
        min_len_lr     = int(self.spin_min_len_long.get_value())
        max_len_lr     = int(self.spin_max_len_long.get_value())
        min_qual       = int(self.spin_min_qual.get_value())
        host_dep       = "yes" if self.chk_host_dep.get_active() else "no"
        min_cov        = int(self.spin_min_cov.get_value())
        max_n_pct      = self.spin_max_n.get_value()
        min_segs       = int(self.spin_min_segs.get_value())

        # v1.0
        run_ivar       = "yes" if self.chk_ivar.get_active() else "no"
        run_medaka     = "yes" if self.chk_medaka.get_active() else "no"
        medaka_env     = self.entry_medaka_env.get_text().strip() or "medaka_env"
        ivar_freq      = round(self.spin_ivar_freq.get_value(), 2)
        ivar_depth     = int(self.spin_ivar_depth.get_value())
        minority_freq  = round(self.spin_minor_freq.get_value(), 2)
        coinf_pct      = round(self.spin_coinf_pct.get_value(), 1)
        run_antiviral  = "yes" if self.chk_antiviral.get_active() else "no"
        run_h5         = "yes" if self.chk_h5.get_active() else "no"
        run_fullvarcall= "yes" if self.chk_fullvarcall.get_active() else "no"
        gisaid_loc     = self.entry_gisaid_loc.get_text().strip()
        gisaid_year    = int(self.spin_gisaid_year.get_value())

        return {
            "input_dir":       input_dir,
            "output_dir":      output_dir,
            "irma_module":     irma_mod,
            "seq_type":        seq_type,
            "seq_mode":        seq_mode,
            "adapter_fasta":   adapter_fasta,
            "run_fastqc":      run_fastqc,
            "min_len_illumina": min_len_il,
            "min_len_long":    min_len_lr,
            "max_len_long":    max_len_lr,
            "min_qual":        min_qual,
            "host_depletion":  host_dep,
            "min_coverage":    min_cov,
            "max_n_pct":       max_n_pct,
            "min_segments":    min_segs,
            # v1.0
            "run_ivar":        run_ivar,
            "run_medaka":      run_medaka,
            "medaka_env":      medaka_env,
            "ivar_freq":       ivar_freq,
            "ivar_depth":      ivar_depth,
            "minority_freq":   minority_freq,
            "coinfection_pct": coinf_pct,
            "run_antiviral":   run_antiviral,
            "run_h5_virulence": run_h5,
            "run_fullvarcall": run_fullvarcall,
            "gisaid_location": gisaid_loc,
            "gisaid_year":     gisaid_year,
        }

    def _build_command(self, params: dict) -> list:
        """Builds the argument list for the bash script."""
        script_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "script_influenza_gui.sh"
        )

        cmd = [
            "bash", script_path,
            params["input_dir"],
            params["output_dir"],
            params["irma_module"],
        ]

        if params["seq_mode"]:
            cmd.append(params["seq_mode"])
        else:
            cmd.append("")  # placeholder to force correct parsing

        # Flags de QC existentes
        cmd += ["--run_fastqc",    params["run_fastqc"]]
        cmd += ["--min_len_short", str(params["min_len_illumina"])]
        cmd += ["--min_len_long",  str(params["min_len_long"])]
        cmd += ["--max_len_long",  str(params["max_len_long"])]
        cmd += ["--min_qual",      str(params["min_qual"])]
        cmd += ["--host_depletion", params["host_depletion"]]
        cmd += ["--min_coverage",  str(params["min_coverage"])]
        cmd += ["--max_n_pct",     str(params["max_n_pct"])]
        cmd += ["--min_segments",  str(params["min_segments"])]

        if params["adapter_fasta"]:
            cmd += ["--adapter_fasta", params["adapter_fasta"]]

        # Flags v1.0
        cmd += ["--seq_type",        params["seq_type"]]
        cmd += ["--ivar",            params["run_ivar"]]
        cmd += ["--medaka",          params["run_medaka"]]
        cmd += ["--medaka_env",      params["medaka_env"]]
        cmd += ["--ivar_freq",       str(params["ivar_freq"])]
        cmd += ["--ivar_depth",      str(params["ivar_depth"])]
        cmd += ["--minority_freq",   str(params["minority_freq"])]
        cmd += ["--coinfection_pct", str(params["coinfection_pct"])]
        cmd += ["--antiviral",       params["run_antiviral"]]
        cmd += ["--h5_virulence",    params["run_h5_virulence"]]
        cmd += ["--fullvarcall",     params["run_fullvarcall"]]

        if params["gisaid_location"]:
            cmd += ["--gisaid_location", params["gisaid_location"]]
        cmd += ["--gisaid_year", str(params["gisaid_year"])]

        return cmd

    # ── Execução do pipeline ──────────────────────────────────────────────

    def _on_run_clicked(self, _btn):
        params = self._get_params()

        # Validações básicas
        if not params["input_dir"] or not os.path.isdir(params["input_dir"]):
            self._show_error("Please select a valid input directory.")
            return
        if not params["output_dir"]:
            self._show_error("Please select an output directory.")
            return
        if not params["irma_module"]:
            self._show_error("Please select an IRMA module.")
            return

        # Salva configurações
        save_config({
            "input_dir":       params["input_dir"],
            "output_dir":      params["output_dir"],
            "irma_module":     params["irma_module"],
            "seq_type":        params["seq_type"],
            "seq_mode":        self.combo_mode.get_active_text(),
            "adapter_fasta":   params["adapter_fasta"],
            "run_fastqc":      self.chk_fastqc.get_active(),
            "min_len_illumina": str(params["min_len_illumina"]),
            "min_len_long":    str(params["min_len_long"]),
            "max_len_long":    str(params["max_len_long"]),
            "min_qual":        str(params["min_qual"]),
            "host_depletion":  self.chk_host_dep.get_active(),
            "min_coverage":    str(params["min_coverage"]),
            "max_n_pct":       str(params["max_n_pct"]),
            "min_segments":    str(params["min_segments"]),
            # v1.0
            "run_ivar":        self.chk_ivar.get_active(),
            "run_medaka":      self.chk_medaka.get_active(),
            "medaka_env":      params["medaka_env"],
            "ivar_freq":       str(params["ivar_freq"]),
            "ivar_depth":      str(params["ivar_depth"]),
            "minority_freq":   str(params["minority_freq"]),
            "coinfection_pct": str(params["coinfection_pct"]),
            "run_antiviral":   self.chk_antiviral.get_active(),
            "run_h5_virulence": self.chk_h5.get_active(),
            "run_fullvarcall": self.chk_fullvarcall.get_active(),
            "gisaid_location": params["gisaid_location"],
            "gisaid_year":     str(params["gisaid_year"]),
        })

        cmd = self._build_command(params)

        # Log de início
        self._log("=" * 70)
        self._log(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Starting MK Flu-Pipe v1.0...")
        self._log(f"  Input        : {params['input_dir']}")
        self._log(f"  Output       : {params['output_dir']}")
        self._log(f"  IRMA Module  : {params['irma_module']}")
        self._log(f"  Data type    : {params['seq_type']}")
        self._log(f"  File mode    : {params['seq_mode'] or 'auto'}")
        self._log(f"  FastQC       : {params['run_fastqc']}")
        self._log(f"  Adapter      : {params['adapter_fasta'] or 'auto-detect'}")
        self._log(f"  Min len I    : {params['min_len_illumina']} nt")
        self._log(f"  Min len L    : {params['min_len_long']} nt")
        self._log(f"  Min qual     : Q{params['min_qual']}")
        self._log(f"  Host dep     : {params['host_depletion']}")
        self._log(f"  Min cov      : {params['min_coverage']}×")
        self._log(f"  Max N%       : {params['max_n_pct']}%")
        self._log(f"  Min segs     : {params['min_segments']}")
        self._log(f"  iVar         : {params['run_ivar']} (freq≥{params['ivar_freq']}, depth≥{params['ivar_depth']})")
        self._log(f"  Medaka       : {params['run_medaka']} (env: {params['medaka_env']})")
        self._log(f"  Co-infection : minor≥{params['minority_freq']}, alert≥{params['coinfection_pct']}% positions")
        self._log(f"  Antivirals   : {params['run_antiviral']}")
        self._log(f"  H5 virulence: {params['run_h5_virulence']} (condicional)")
        self._log(f"  Full varcall : {params['run_fullvarcall']} (all segments RefSeq+GFF3)")
        self._log(f"  GISAID loc   : {params['gisaid_location'] or 'not set'}")
        self._log(f"  GISAID year  : {params['gisaid_year']}")
        self._log("=" * 70)

        self._executando = True
        self._set_ui_running(True)

        def _run():
            try:
                self._processo = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )
                for line in iter(self._processo.stdout.readline, ""):
                    self._log(line.rstrip())
                self._processo.wait()
                rc = self._processo.returncode
                if rc == 0:
                    self._log("\n✅ Pipeline completed successfully!")
                    self._show_finished_dialog(params["output_dir"])
                else:
                    self._log(f"\n⚠ Pipeline exited with code {rc}.")
                    self._show_pipeline_error_dialog(
                        f"Pipeline exited with return code {rc}.\n"
                        f"Check execution log at:\n{params['output_dir']}/run_log.txt"
                    )
            except Exception as e:
                self._log(f"\n❌ Error running pipeline: {e}")
                self._show_pipeline_error_dialog(str(e))
            finally:
                self._executando = False
                self._set_ui_running(False)

        threading.Thread(target=_run, daemon=True).start()

    def _on_stop_clicked(self, _btn):
        if self._processo and self._processo.poll() is None:
            self._processo.terminate()
            self._log("\n⏹ Pipeline stopped by user.")
            self._executando = False
            self._set_ui_running(False)

    def _on_save_log(self, _btn):
        buf = self.txt_log.get_buffer()
        text = buf.get_text(buf.get_start_iter(), buf.get_end_iter(), True)
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        log_path = os.path.join(os.path.expanduser("~"), f"mkflupipe_log_{ts}.txt")
        try:
            with open(log_path, "w") as f:
                f.write(text)
            self._log(f"\n💾 Log saved at: {log_path}")
        except Exception as e:
            self._log(f"\n❌ Error saving log: {e}")

    def _on_clear_log(self, _btn):
        self.txt_log.get_buffer().set_text("")

    def _on_window_delete(self, _win, _event):
        if self._executando:
            dialog = Gtk.MessageDialog(
                transient_for=self,
                flags=0,
                message_type=Gtk.MessageType.QUESTION,
                buttons=Gtk.ButtonsType.YES_NO,
                text="The pipeline is running.",
            )
            dialog.format_secondary_text(
                "Do you want to stop execution and close the window?"
            )
            resp = dialog.run()
            dialog.destroy()
            if resp == Gtk.ResponseType.YES:
                self._on_stop_clicked(None)
                Gtk.main_quit()
            return True
        Gtk.main_quit()
        return False

    def _show_error(self, msg: str):
        dialog = Gtk.MessageDialog(
            transient_for=self,
            flags=0,
            message_type=Gtk.MessageType.ERROR,
            buttons=Gtk.ButtonsType.OK,
            text=msg,
        )
        dialog.run()
        dialog.destroy()

    def _show_finished_dialog(self, output_dir: str):
        """Exibe aviso sonoro e visual ao concluir o pipeline com sucesso."""
        def _show():
            self.get_display().beep()
            dialog = Gtk.MessageDialog(
                transient_for=self,
                flags=0,
                message_type=Gtk.MessageType.INFO,
                buttons=Gtk.ButtonsType.OK,
                text="✅  Analysis Complete!",
            )
            dialog.format_secondary_text(
                f"Check results at:\n{output_dir}"
            )
            dialog.run()
            dialog.destroy()
            return False
        GLib.idle_add(_show)

    def _show_pipeline_error_dialog(self, msg: str):
        """Exibe aviso sonoro e visual quando o pipeline encerra com erro."""
        def _show():
            self.get_display().beep()
            dialog = Gtk.MessageDialog(
                transient_for=self,
                flags=0,
                message_type=Gtk.MessageType.ERROR,
                buttons=Gtk.ButtonsType.OK,
                text="❌  Pipeline Execution Error",
            )
            dialog.format_secondary_text(msg)
            dialog.run()
            dialog.destroy()
            return False
        GLib.idle_add(_show)


def main():
    win = JanelaPrincipal()
    win.show_all()
    Gtk.main()


if __name__ == "__main__":
    main()
