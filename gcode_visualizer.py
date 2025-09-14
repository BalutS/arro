import logging
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

# Configuración global
CONFIG = {
    "margin_ratio": 0.2
}

def add_margin_to_axis_limits(min_v, max_v, margin_ratio=CONFIG["margin_ratio"]):
    if min_v == max_v:
        return min_v - 0.5, max_v + 0.5
    dv = (max_v - min_v) * margin_ratio
    return (min_v - dv, max_v + dv)

class GcodeVisualizer:
    def __init__(self, gcode_reader):
        self.reader = gcode_reader
        plt.style.use('seaborn-v0_8-darkgrid')

    def _create_axis(self, fig, projection='2d', title=""):
        projection = projection.lower()
        if projection == '2d':
            ax = fig.add_subplot(111)
        elif projection == '3d':
            ax = fig.add_subplot(111, projection='3d')
        else:
            raise ValueError("Proyección no válida, use '2d' o '3d'")

        ax.set_title(title)
        return ax

    def plot_first_layer(self, fig, minP):
        """Visualización de la primera capa en 2D."""
        ax = self._create_axis(fig, projection='2d', title="Primera Capa")

        if self.reader.n_layers == 0:
            ax.text(0.5, 0.5, "No hay capas para mostrar", ha='center', va='center', transform=ax.transAxes)
            return

        try:
            z0 = float(self.reader.seg_index[0])
            layer0 = self.reader.get_layerSegs(z0, z0)

            if not layer0:
                ax.text(0.5, 0.5, "La primera capa no contiene segmentos", ha='center', va='center', transform=ax.transAxes)
                return

            for x0, y0, x1, y1, z in layer0:
                ax.plot([x0, x1], [y0, y1], 'k-', linewidth=0.5)

            # Línea de corte transversal
            ax.axvline(x=minP, color='g', linestyle='--', linewidth=1)
            ax.set_aspect('equal', adjustable='datalim')
            ax.set_title(f'Primera capa z={z0:.3f}')

            xmin, xmax, ymin, ymax, _, _ = self.reader.xyzlimits
            ax.set_xlim(add_margin_to_axis_limits(xmin, xmax))
            ax.set_ylim(add_margin_to_axis_limits(ymin, ymax))

        except Exception as e:
            logging.warning(f"No se pudo graficar la primera capa: {e}")
            ax.text(0.5, 0.5, f"Error al graficar la primera capa:\n{e}", ha='center', va='center', transform=ax.transAxes)

    def animate_layer_step_by_step(self, fig, layer_index, animation_time=5):
        """Animación de una capa, segmento por segmento."""
        ax = self._create_axis(fig, projection='2d', title=f"Animación de la capa {layer_index}")
        xmin, xmax, ymin, ymax, _, _ = self.reader.xyzlimits
        ax.set_xlim(add_margin_to_axis_limits(xmin, xmax))
        ax.set_ylim(add_margin_to_axis_limits(ymin, ymax))
        ax.set_xlabel('x')
        ax.set_ylabel('y')

        try:
            layer_z = self.reader.seg_index[layer_index]
            temp = self.reader.get_layerSegs(layer_z, layer_z)

            if not temp:
                logging.warning(f"No hay segmentos en la capa {layer_index}")
                ax.text(0.5, 0.5, f"No hay segmentos en la capa {layer_index}", ha='center', va='center', transform=ax.transAxes)
                return

            # This is a blocking animation, not ideal for tkinter, but let's keep it for now
            # A better implementation would use FuncAnimation
            for i, (x0, y0, x1, y1, z) in enumerate(temp):
                ax.plot([x0, x1], [y0, y1], 'y-')
                fig.canvas.draw()
                fig.canvas.flush_events() # Update the plot
                plt.pause(0.01) # Small pause

            ax.set_title(f"Animación de la capa {layer_index} (Completa)")

        except IndexError:
            logging.error(f"El índice de capa {layer_index} está fuera de rango.")
            ax.text(0.5, 0.5, f"Índice de capa {layer_index} fuera de rango", ha='center', va='center', transform=ax.transAxes)

    def plot_3d_view(self, fig, minP):
        """Visualización 3D de las trayectorias."""
        ax = self._create_axis(fig, projection='3d', title="Vista 3D")

        if self.reader.segs is None or len(self.reader.segs) == 0:
            ax.text2D(0.5, 0.5, "No hay segmentos para mostrar", ha='center', va='center', transform=ax.transAxes)
            return

        try:
            # Limitar a 10 capas for performance
            if self.reader.n_layers > 10:
                z_max_layer = self.reader.seg_index[9]
                segs_to_plot = self.reader.get_layerSegs(self.reader.seg_index[0], z_max_layer)
                title = 'Trayectorias 3D (Primeras 10 capas)'
            else:
                segs_to_plot = self.reader.segs
                title = 'Trayectorias 3D (G-code)'

            ax.set_title(title)

            if segs_to_plot:
                arr = np.array(segs_to_plot)
                for seg in arr:
                    x0, y0, x1, y1, z = seg
                    ax3d.plot([x0, x1], [y0, y1], [z, z], 'c-', linewidth=0.3)

                # Plano de corte transversal
                _, _, ymin, ymax, _, _ = self.reader.xyzlimits
                zmin_plot, zmax_plot = np.min(arr[:, 4]), np.max(arr[:, 4])
                Y = np.array([[ymin, ymax], [ymin, ymax]])
                Z = np.array([[zmin_plot, zmin_plot], [zmax_plot, zmax_plot]])
                X = np.full_like(Y, minP)
                ax.plot_surface(X, Y, Z, color='g', alpha=0.3, shade=False)

        except Exception as e:
            logging.warning(f"No se pudo graficar en 3D: {e}")
            ax.text2D(0.5, 0.5, f"Error al graficar en 3D:\n{e}", ha='center', va='center', transform=ax.transAxes)

    def plot_cross_section(self, fig, cut_points):
        """Dibuja los puntos de corte y un rectángulo que los contiene."""
        ax = self._create_axis(fig, projection='2d', title="Vista de Corte Transversal")
        ax.axis('on') # Override 'off' from original

        if not cut_points:
            logging.warning("No hay puntos de corte para graficar.")
            ax.text(0.5, 0.5, "No hay puntos de corte para graficar", ha='center', va='center', transform=ax.transAxes)
            return

        y_coords = [p[1] for p in cut_points]
        z_coords = [p[2] for p in cut_points]

        miny, maxy = min(y_coords), max(y_coords)
        minz, maxz = min(z_coords), max(z_coords)

        ax.add_patch(Rectangle((miny, minz), (maxy - miny), (maxz - minz), facecolor='k', fill=False, lw=2))
        ax.scatter(y_coords, z_coords, color=['red'], s=10)
        ax.set_xlabel("Eje Y (mm)")
        ax.set_ylabel("Eje Z (mm)")
        ax.set_aspect('equal', adjustable='datalim')

        xlim = add_margin_to_axis_limits(miny, maxy)
        ylim = add_margin_to_axis_limits(minz, maxz)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    def plot_full_2d_simulation(self, fig, minP):
        """Animación de la impresión completa en 2D."""
        ax = self._create_axis(fig, projection='2d', title="Simulación de impresión (2D)")
        ax.set_aspect('equal', adjustable='datalim')

        if self.reader.segs is None or len(self.reader.segs) == 0:
            ax.text(0.5, 0.5, "No hay segmentos para mostrar", ha='center', va='center', transform=ax.transAxes)
            return

        xmin, xmax, ymin, ymax, _, _ = self.reader.xyzlimits
        ax.set_xlim(add_margin_to_axis_limits(xmin, xmax))
        ax.set_ylim(add_margin_to_axis_limits(ymin, ymax))

        # Línea de corte transversal
        ax.axvline(x=minP, color='g', linestyle='--', linewidth=1)

        # This is also blocking. For a real GUI, this should be handled with a non-blocking animation.
        for idx, z in enumerate(self.reader.seg_index):
            layer = self.reader.get_layerSegs(z, z)
            for x0, y0, x1, y1, _ in layer:
                ax.plot([x0, x1], [y0, y1], 'k-', linewidth=0.4)
            ax.set_title(f"Simulación de impresión (2D) - Capa {idx+1}/{self.reader.n_layers}")
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(0.05)

        ax.set_title("Impresión completa (2D)")
