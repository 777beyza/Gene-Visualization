import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pandas as pd
import numpy as np
from visualizer import manhattan_plot, venn_diagram, heatmap, volcano_plot


class GeneVizApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GeneViz: Genomik Veri Görselleştirme")
        self.root.geometry("800x600")
        self.data = None

        self.create_menu()
        self.setup_ui()

    def create_menu(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Yardım", menu=help_menu)
        help_menu.add_command(
            label="Manhattan Plot", command=lambda: self.show_help("Manhattan Plot")
        )
        help_menu.add_command(
            label="Venn Diyagramı", command=lambda: self.show_help("Venn Diyagramı")
        )
        help_menu.add_command(
            label="Isı Haritası", command=lambda: self.show_help("Isı Haritası")
        )
        help_menu.add_command(
            label="Volcano Plot", command=lambda: self.show_help("Volcano Plot")
        )
        help_menu.add_separator()
        help_menu.add_command(label="Genel Bilgi", command=self.show_general_help)

    def show_help(self, plot_type):
        help_texts = {
            "Manhattan Plot": {
                "title": "Manhattan Plot Yardım",
                "text": "Manhattan plot, genom çapında ilişkilendirme çalışmaları (GWAS) sonuçlarını görselleştirir.\n\n"
                "Gerekli CSV Sütunları:\n"
                "- 'chromosome': Kromozom numarası\n"
                "- 'position': Kromozom üzerindeki pozisyon\n"
                "- 'p_value': İstatistiksel önem derecesi",
            },
            "Venn Diyagramı": {
                "title": "Venn Diyagramı Yardım",
                "text": "Venn diyagramı, gen setleri arasındaki kesişim ve farkları gösterir. Özel giriş alanlarını kullanarak gen setlerini virgülle ayırarak girin.",
            },
            "Isı Haritası": {
                "title": "Isı Haritası Yardım",
                "text": "Isı haritası, gen ekspresyon verilerini renk yoğunluğu ile görselleştirir.\n\n"
                "Gerekli CSV Sütunları:\n"
                "- 'Gene': Gen adı\n"
                "- 'Condition': Örnek koşulu\n"
                "- 'Expression': Ekspresyon seviyesi",
            },
            "Volcano Plot": {
                "title": "Volcano Plot Yardım",
                "text": "Volcano plot, iki koşul arasındaki gen ekspresyon değişimlerini gösterir. 'Anlamlı' genler kırmızı renkte vurgulanır.\n\n"
                "Gerekli CSV Sütunları:\n"
                "- 'log2_fold_change': İki koşul arasındaki kat değişiminin log2 değeri\n"
                "- 'p_value': İstatistiksel önem derecesi",
            },
        }
        info = help_texts.get(plot_type)
        if info:
            messagebox.showinfo(info["title"], info["text"])

    def show_general_help(self):
        messagebox.showinfo(
            "Genel Bilgi",
            "GeneViz uygulaması, genomik verileri farklı biyoinformatik grafik türleriyle görselleştirmenizi sağlar.\n\n"
            "1. 'Dosya Seç' ile bir CSV dosyası yükleyin.\n"
            "2. 'Grafik Tipi Seç' bölümünden istediğiniz grafiği seçin.\n"
            "3. Varsa ilgili ayarları yapın.\n"
            "4. 'Grafik Çiz' butonuna basarak görselleştirmeyi oluşturun.",
        )

    def setup_ui(self):
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(pady=10, fill=tk.X)
        file_frame = ttk.LabelFrame(control_frame, text="1. Veri Dosyası Seç")
        file_frame.pack(side=tk.LEFT, padx=5, pady=5, fill=tk.X, expand=True)
        self.file_path_label = ttk.Label(file_frame, text="Henüz dosya seçilmedi.")
        self.file_path_label.pack(side=tk.LEFT, padx=5, pady=5)
        select_btn = ttk.Button(file_frame, text="Dosya Seç", command=self.load_data)
        select_btn.pack(side=tk.RIGHT, padx=5, pady=5)

        graph_frame = ttk.LabelFrame(control_frame, text="2. Grafik Tipi Seç")
        graph_frame.pack(side=tk.LEFT, padx=5, pady=5, fill=tk.X, expand=True)
        self.graph_type = tk.StringVar(value="Manhattan Plot")
        radio_frame = ttk.Frame(graph_frame)
        radio_frame.pack(pady=5, padx=5, fill=tk.X)
        ttk.Radiobutton(
            radio_frame,
            text="Manhattan Plot",
            variable=self.graph_type,
            value="Manhattan Plot",
            command=self.toggle_settings,
        ).pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(
            radio_frame,
            text="Venn Diyagramı",
            variable=self.graph_type,
            value="Venn Diyagramı",
            command=self.toggle_settings,
        ).pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(
            radio_frame,
            text="Isı Haritası",
            variable=self.graph_type,
            value="Isı Haritası",
            command=self.toggle_settings,
        ).pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(
            radio_frame,
            text="Volcano Plot",
            variable=self.graph_type,
            value="Volcano Plot",
            command=self.toggle_settings,
        ).pack(side=tk.LEFT, padx=5)

        settings_and_draw_frame = ttk.Frame(main_frame)
        settings_and_draw_frame.pack(fill=tk.BOTH, expand=True)
        self.common_settings_frame = ttk.LabelFrame(
            settings_and_draw_frame, text="Genel Ayarlar"
        )
        self.common_settings_frame.pack(pady=10, fill=tk.X)
        ttk.Label(self.common_settings_frame, text="Başlık:").pack(
            side=tk.LEFT, padx=5, pady=5
        )
        self.title_entry = ttk.Entry(self.common_settings_frame, width=40)
        self.title_entry.pack(side=tk.LEFT, padx=5, pady=5)
        self.title_entry.insert(0, "Grafik Başlığı")
        self.manhattan_settings_frame = ttk.LabelFrame(
            settings_and_draw_frame, text="Manhattan Ayarları"
        )
        self.venn_settings_frame = ttk.LabelFrame(
            settings_and_draw_frame, text="Venn Diyagramı Ayarları"
        )
        self.heatmap_settings_frame = ttk.LabelFrame(
            settings_and_draw_frame, text="Isı Haritası Ayarları"
        )
        self.volcano_settings_frame = ttk.LabelFrame(
            settings_and_draw_frame, text="Volcano Plot Ayarları"
        )
        ttk.Label(self.volcano_settings_frame, text="log2 Kat Değişimi Eşiği:").pack(
            side=tk.LEFT, padx=5, pady=5
        )
        self.fold_change_entry = ttk.Entry(self.volcano_settings_frame, width=10)
        self.fold_change_entry.pack(side=tk.LEFT, padx=5, pady=5)
        self.fold_change_entry.insert(0, "1.0")

        ttk.Label(self.volcano_settings_frame, text="p-Değeri Eşiği:").pack(
            side=tk.LEFT, padx=5, pady=5
        )
        self.p_value_entry = ttk.Entry(self.volcano_settings_frame, width=10)
        self.p_value_entry.pack(side=tk.LEFT, padx=5, pady=5)
        self.p_value_entry.insert(0, "0.05")
        ttk.Label(self.venn_settings_frame, text="Set 1 (Etiket):").pack(
            side=tk.TOP, anchor=tk.W, padx=5, pady=2
        )
        self.venn_label1_entry = ttk.Entry(self.venn_settings_frame, width=20)
        self.venn_label1_entry.pack(side=tk.TOP, anchor=tk.W, padx=5, pady=2)
        self.venn_label1_entry.insert(0, "Set 1")

        ttk.Label(self.venn_settings_frame, text="Genler (Virgülle Ayırın):").pack(
            side=tk.TOP, anchor=tk.W, padx=5, pady=2
        )
        self.venn_set1_entry = tk.Text(self.venn_settings_frame, height=2, width=40)
        self.venn_set1_entry.pack(side=tk.TOP, fill=tk.X, padx=5, pady=2)

        ttk.Label(self.venn_settings_frame, text="Set 2 (Etiket):").pack(
            side=tk.TOP, anchor=tk.W, padx=5, pady=2
        )
        self.venn_label2_entry = ttk.Entry(self.venn_settings_frame, width=20)
        self.venn_label2_entry.pack(side=tk.TOP, anchor=tk.W, padx=5, pady=2)
        self.venn_label2_entry.insert(0, "Set 2")

        ttk.Label(self.venn_settings_frame, text="Genler (Virgülle Ayırın):").pack(
            side=tk.TOP, anchor=tk.W, padx=5, pady=2
        )
        self.venn_set2_entry = tk.Text(self.venn_settings_frame, height=2, width=40)
        self.venn_set2_entry.pack(side=tk.TOP, fill=tk.X, padx=5, pady=2)

        ttk.Label(self.venn_settings_frame, text="Set 3 (Opsiyonel):").pack(
            side=tk.TOP, anchor=tk.W, padx=5, pady=2
        )
        self.venn_set3_entry = tk.Text(self.venn_settings_frame, height=2, width=40)
        self.venn_set3_entry.pack(side=tk.TOP, fill=tk.X, padx=5, pady=2)
        draw_btn = ttk.Button(main_frame, text="Grafik Çiz", command=self.draw_plot)
        draw_btn.pack(pady=10)
        self.status_label = ttk.Label(main_frame, text="", foreground="blue")
        self.status_label.pack(pady=10)
        self.toggle_settings()

    def toggle_settings(self):
        self.manhattan_settings_frame.pack_forget()
        self.venn_settings_frame.pack_forget()
        self.heatmap_settings_frame.pack_forget()
        self.volcano_settings_frame.pack_forget()
        selected_plot = self.graph_type.get()
        if selected_plot == "Manhattan Plot":
            self.manhattan_settings_frame.pack(pady=10, fill=tk.X)
        elif selected_plot == "Venn Diyagramı":
            self.venn_settings_frame.pack(pady=10, fill=tk.X)
        elif selected_plot == "Isı Haritası":
            self.heatmap_settings_frame.pack(pady=10, fill=tk.X)
        elif selected_plot == "Volcano Plot":
            self.volcano_settings_frame.pack(pady=10, fill=tk.X)

    def load_data(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
        if file_path:
            try:
                self.data = pd.read_csv(file_path)
                self.file_path_label.config(text=f"Dosya: {file_path.split('/')[-1]}")
                self.status_label.config(text="Dosya başarıyla yüklendi!")
            except Exception as e:
                messagebox.showerror("Hata", f"Dosya yüklenirken bir hata oluştu:\n{e}")
                self.data = None
                self.file_path_label.config(text="Henüz dosya seçilmedi.")

    def draw_plot(self):
        if self.data is None and self.graph_type.get() != "Venn Diyagramı":
            messagebox.showwarning("Uyarı", "Lütfen önce bir veri dosyası seçin!")
            return

        plot_title = (
            self.title_entry.get() if self.title_entry.get() else "Grafik Başlığı"
        )
        selected_plot = self.graph_type.get()

        try:
            if selected_plot == "Manhattan Plot":
                manhattan_plot(self.data.copy(), title=plot_title)
            elif selected_plot == "Isı Haritası":
                heatmap(self.data.copy(), title=plot_title)
            elif selected_plot == "Volcano Plot":
                fold_change = float(self.fold_change_entry.get())
                p_value = float(self.p_value_entry.get())
                volcano_plot(
                    self.data.copy(),
                    fold_change_threshold=fold_change,
                    p_value_threshold=p_value,
                    title=plot_title,
                )
            elif selected_plot == "Venn Diyagramı":
                self.draw_venn_diagram(plot_title)
        except ValueError:
            messagebox.showerror(
                "Hata", "Lütfen eşik değerleri için geçerli sayılar girin."
            )
        except Exception as e:
            messagebox.showerror(
                "Hata",
                f"Grafik çizilirken bir hata oluştu:\n{e}\n\n"
                f"Lütfen verilerinizin '{selected_plot}' grafiği için gerekli sütunları içerdiğinden emin olun.",
            )

    def draw_venn_diagram(self, title):
        try:
            set1_genes = {
                gene.strip()
                for gene in self.venn_set1_entry.get("1.0", tk.END).split(",")
                if gene.strip()
            }
            set2_genes = {
                gene.strip()
                for gene in self.venn_set2_entry.get("1.0", tk.END).split(",")
                if gene.strip()
            }
            set3_genes_text = self.venn_set3_entry.get("1.0", tk.END).strip()

            sets = [set1_genes, set2_genes]
            labels = [self.venn_label1_entry.get(), self.venn_label2_entry.get()]

            if set3_genes_text:
                set3_genes = {
                    gene.strip() for gene in set3_genes_text.split(",") if gene.strip()
                }
                sets.append(set3_genes)
                labels.append(
                    self.venn_label3_entry.get()
                )  

            if not set1_genes or not set2_genes:
                messagebox.showwarning(
                    "Uyarı", "Lütfen en az iki gen seti için genleri girin."
                )
                return

            venn_diagram(sets, labels, title=title)
        except Exception as e:
            messagebox.showerror(
                "Hata", f"Venn diyagramı çizilirken bir hata oluştu:\n{e}"
            )


if __name__ == "__main__":
    root = tk.Tk()
    app = GeneVizApp(root)
    root.mainloop()
