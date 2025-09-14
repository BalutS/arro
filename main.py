import tkinter as tk
from tkinter import filedialog, messagebox
import logging

from gcode_reader import GcodeReader, GcodeType, CONFIG
from gcode_visualizer import GcodeVisualizer
from ui import InterfaceA, InterfaceB

class MainApplication:
    def __init__(self, root):
        self.root = root
        self.root.title("G-code Reader Optimizado")
        self.root.geometry("1000x800")

        self.gcode_reader = None
        self.visualizer = None
        self.min_cut_info = {}

        self.current_ui_type = None
        self.ui_frame = None

        self.switch_ui('A') # Start with Interface A

    def switch_ui(self, ui_type):
        if self.ui_frame:
            self.ui_frame.destroy()

        if ui_type == 'A':
            self.ui_frame = InterfaceA(self.root, self)
            self.current_ui_type = 'A'
        elif ui_type == 'B':
            self.ui_frame = InterfaceB(self.root, self)
            self.current_ui_type = 'B'
        else:
            raise ValueError("Unknown UI type")

    def load_gcode_file(self):
        filename = filedialog.askopenfilename(
            title="Seleccione un archivo .gcode",
            filetypes=[("G-code files", "*.gcode"), ("All files", "*.*")]
        )
        if not filename:
            return

        try:
            self.gcode_reader = GcodeReader(filename, GcodeType.FDM_REGULAR)
            self.visualizer = GcodeVisualizer(self.gcode_reader)

            # Perform analysis once
            self._analyze_gcode()

            # Automatically show the first layer plot
            self.show_first_layer()

        except FileNotFoundError as e:
            messagebox.showerror("Error", str(e))
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo procesar el archivo: {e}")
            logging.error(f"Failed to process file: {e}", exc_info=True)

    def _analyze_gcode(self):
        if not self.gcode_reader:
            return

        # This can be slow, might be worth running in a separate thread in a real app
        areaCortes, minP, minArea, minCut_solidArea, minCut_points = self.gcode_reader.search_minorArea(
            CONFIG["delta"], CONFIG["step"], CONFIG["ejeMenor"], CONFIG["ejeMayor"]
        )
        self.min_cut_info = {
            "minP": minP,
            "minArea": minArea,
            "minCut_solidArea": minCut_solidArea,
            "minCut_points": minCut_points
        }

    def _check_if_loaded(self):
        if not self.gcode_reader:
            messagebox.showinfo("Información", "Por favor, cargue primero un archivo G-code.")
            return False
        return True

    def show_first_layer(self):
        if self._check_if_loaded():
            minP = self.min_cut_info.get("minP", 0)
            self.ui_frame._draw_plot(self.visualizer.plot_first_layer, minP)

    def show_cross_section(self):
        if self._check_if_loaded():
            cut_points = self.min_cut_info.get("minCut_points", [])
            self.ui_frame._draw_plot(self.visualizer.plot_cross_section, cut_points)

    def show_3d_view(self):
        if self._check_if_loaded():
            minP = self.min_cut_info.get("minP", 0)
            self.ui_frame._draw_plot(self.visualizer.plot_3d_view, minP)

    def animate_layer(self, layer_index):
        if self._check_if_loaded():
            if self.gcode_reader.n_layers <= layer_index:
                messagebox.showwarning("Advertencia", f"La capa {layer_index} no existe. La pieza solo tiene {self.gcode_reader.n_layers} capas.")
                return
            # Note: this will block the GUI
            self.ui_frame._draw_plot(self.visualizer.animate_layer_step_by_step, layer_index)

    def show_full_2d_simulation(self):
        if self._check_if_loaded():
            minP = self.min_cut_info.get("minP", 0)
            # Note: this will block the GUI
            self.ui_frame._draw_plot(self.visualizer.plot_full_2d_simulation, minP)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    logging.info("Iniciando Gcode Reader Optimizado...")

    root = tk.Tk()
    app = MainApplication(root)
    root.mainloop()
