import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

class BaseInterface(tk.Frame):
    def __init__(self, parent, app_controller):
        super().__init__(parent)
        self.app = app_controller
        self.pack(fill=tk.BOTH, expand=True)
        self._create_widgets()

    def _create_widgets(self):
        raise NotImplementedError

    def _create_matplotlib_canvas(self):
        self.figure = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.figure, self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(self.canvas, self)
        toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def _clear_canvas(self):
        self.figure.clear()
        self.canvas.draw()

    def _draw_plot(self, plot_function, *args):
        self._clear_canvas()
        plot_function(self.figure, *args)
        self.canvas.draw()


class InterfaceA(BaseInterface):
    def _create_widgets(self):
        # Top frame for controls
        control_frame = tk.Frame(self)
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

        # Matplotlib Canvas
        self._create_matplotlib_canvas()

        # --- Control Buttons ---

        # Load file button
        load_button = ttk.Button(control_frame, text="Cargar G-code", command=self.app.load_gcode_file)
        load_button.pack(side=tk.LEFT, padx=5)

        # Separator
        ttk.Separator(control_frame, orient=tk.VERTICAL).pack(side=tk.LEFT, padx=5, fill='y')

        # Visualization buttons
        vis_label = ttk.Label(control_frame, text="Visualizaciones:")
        vis_label.pack(side=tk.LEFT, padx=(10, 2))

        vis_buttons = {
            "Primera Capa": self.app.show_first_layer,
            "Corte Transversal": self.app.show_cross_section,
            "Vista 3D": self.app.show_3d_view,
            "Animar Capa 6": lambda: self.app.animate_layer(6),
            "Simulación 2D Completa": self.app.show_full_2d_simulation,
        }

        for text, command in vis_buttons.items():
            button = ttk.Button(control_frame, text=text, command=command)
            button.pack(side=tk.LEFT, padx=2)

        # Switch UI button on the far right
        switch_button = ttk.Button(control_frame, text="Cambiar a UI B", command=lambda: self.app.switch_ui('B'))
        switch_button.pack(side=tk.RIGHT, padx=5)


class InterfaceB(BaseInterface):
    def _create_widgets(self):
        # Left frame for controls
        control_frame = tk.Frame(self, bg="lightgrey")
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)

        # Matplotlib Canvas
        self._create_matplotlib_canvas()

        # --- Control Buttons ---

        # Load file button
        load_button = tk.Button(control_frame, text="Cargar Archivo", command=self.app.load_gcode_file, bg="#cce", relief=tk.FLAT)
        load_button.pack(pady=10, padx=10, fill=tk.X)

        # Separator
        ttk.Separator(control_frame, orient=tk.HORIZONTAL).pack(pady=5, fill='x')

        # Visualization buttons
        vis_label = tk.Label(control_frame, text="Vistas:", bg="lightgrey", font=("Helvetica", 12, "bold"))
        vis_label.pack(pady=(10, 5))

        vis_buttons = {
            "Capa Inicial": self.app.show_first_layer,
            "Sección de Corte": self.app.show_cross_section,
            "Modelo 3D": self.app.show_3d_view,
            "Animar Capa (6)": lambda: self.app.animate_layer(6),
            "Simulación 2D": self.app.show_full_2d_simulation,
        }

        for text, command in vis_buttons.items():
            button = tk.Button(control_frame, text=text, command=command, anchor="w", bg="lightblue", relief=tk.GROOVE)
            button.pack(pady=2, padx=10, fill=tk.X)

        # Switch UI button at the bottom
        switch_button = tk.Button(control_frame, text="Cambiar a UI A", command=lambda: self.app.switch_ui('A'), bg="pink")
        switch_button.pack(side=tk.BOTTOM, pady=20, padx=10)
