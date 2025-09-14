import collections
from enum import Enum
import math
import os.path
import sys
import logging

import numpy as np

# Configuración global
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

CONFIG = {
    "delta": 7.62,
    "step": 0.1,
    "ejeMenor": 0.119 * 2,
    "ejeMayor": 0.191 * 2,
    "minDist_y": 0.338,
    "margin_ratio": 0.2
}

# Definiciones base
Element = collections.namedtuple('Element', ['x0', 'y0', 'x1', 'y1', 'z'])

class GcodeType(Enum):
    """Tipos de G-code soportados."""
    FDM_REGULAR = 1
    FDM_STRATASYS = 2
    LPBF_REGULAR = 3
    LPBF_SCODE = 4

    @classmethod
    def has_value(cls, value: int) -> bool:
        return any(value == item.value for item in cls)

class GcodeReader:
    """Clase lectora y analítica de archivos G-code."""

    def __init__(self, filename: str, filetype: GcodeType = GcodeType.FDM_REGULAR):
        if not os.path.exists(filename):
            logging.error(f"{filename} no existe!")
            # Instead of exiting, raise an exception
            raise FileNotFoundError(f"{filename} not found!")

        self.filename = filename
        self.filetype = filetype

        self.n_segs = 0
        self.segs = None
        self.n_layers = 0
        self.seg_index = []
        self.xyzlimits = None
        self.minDimensions = None

        # Leer archivo
        self.segs = self._read(filename)
        if self.segs.size == 0:
            logging.warning("No segments found in the gcode file.")
            self.xyzlimits = (0,0,0,0,0,0)
            self.minDimensions = [0,0,0,0]
            self.n_layers = 0
            self.seg_index = []
        else:
            self.xyzlimits = self._compute_xyzlimits(self.segs)
            self.minDimensions = self.get_specimenDimensions()
            self.remove_skirt() # Remove skirt by default

            # Recalculate after skirt removal
            self.xyzlimits = self._compute_xyzlimits(self.segs)
            self.minDimensions = self.get_specimenDimensions()


        logging.info(f"Dimensiones mínimas: {self.minDimensions}")

    def _read(self, filename: str):
        if self.filetype == GcodeType.FDM_REGULAR:
            segs = self._read_fdm_regular(filename)
        else:
            logging.error("Tipo de archivo no soportado")
            # Return empty array instead of exiting
            return np.array([])
        return segs

    def _read_fdm_regular(self, filename: str):
        """Lee un archivo de G-code FDM estándar y extrae segmentos."""
        segs = []
        temp = -float('inf')
        gxyzef = [temp] * 7
        d = dict(zip(['G', 'X', 'Y', 'Z', 'E', 'F', 'S'], range(7)))

        x0 = y0 = temp
        z = -math.inf

        with open(filename, 'r') as infile:
            for raw_line in infile:
                line = raw_line.strip()
                if not line or not line.startswith('G'):
                    continue
                if ';' in line:
                    line = line.split(';')[0]

                for token in line.split():
                    if token[0] in d:
                        try:
                            gxyzef[d[token[0]]] = float(token[1:])
                        except (ValueError, IndexError):
                            # Ignore malformed tokens
                            pass

                if gxyzef[0] == 1:  # G1 (movimiento con extrusión)
                    if np.isfinite(gxyzef[3]):
                        z = gxyzef[3]
                    if np.isfinite(gxyzef[1]) and np.isfinite(gxyzef[2]) and not np.isfinite(gxyzef[4]):
                        x0, y0 = gxyzef[1], gxyzef[2]
                    elif np.isfinite(gxyzef[1]) and np.isfinite(gxyzef[2]) and (gxyzef[4] > 0):
                        if np.isfinite(x0) and np.isfinite(y0):
                            segs.append((x0, y0, gxyzef[1], gxyzef[2], z))
                        x0, y0 = gxyzef[1], gxyzef[2]

                # Reset for next line
                gxyzef[1:3] = [temp, temp]
                gxyzef[4:] = [temp, temp, temp]


        segs = np.array(segs)
        if segs.size > 0:
            self.n_segs = len(segs)
            self.seg_index = np.unique(segs[:, 4])
            self.n_layers = len(self.seg_index)
        else:
            self.n_segs = 0
            self.seg_index = []
            self.n_layers = 0


        logging.info(f"Número de segmentos: {self.n_segs}")
        logging.info(f"Número de capas: {self.n_layers - 1 if self.n_layers > 0 else 0}")

        return segs

    def _compute_xyzlimits(self, seg_list: np.ndarray):
        """Calcula los límites XYZ de todos los segmentos."""
        if seg_list.size == 0:
            return (0, 0, 0, 0, 0, 0)
        arr = np.array(seg_list)
        xmin, xmax = np.min(arr[:, [0, 2]]), np.max(arr[:, [0, 2]])
        ymin, ymax = np.min(arr[:, [1, 3]]), np.max(arr[:, [1, 3]])
        zmin, zmax = np.min(arr[:, 4]), np.max(arr[:, 4])
        return xmin, xmax, ymin, ymax, zmin, zmax

    def get_specimenDimensions(self) -> list[float]:
        """Obtiene las dimensiones mínimas de la probeta en XY."""
        if self.n_layers == 0:
            return [0, 0, 0, 0]

        n_zCoords = len(self.seg_index)
        mz = int(n_zCoords / 2)
        mz_idx = self.seg_index[mz]
        mz_layerSegs = self.get_layerSegs(mz_idx, mz_idx)
        if not mz_layerSegs:
            return [0,0,0,0]
        arr = np.array(mz_layerSegs)
        minx, miny = np.min(arr[:, [0, 2]]), np.min(arr[:, [1, 3]])
        maxx, maxy = np.max(arr[:, [0, 2]]), np.max(arr[:, [1, 3]])
        return [minx, miny, maxx, maxy]

    def get_layerSegs(self, min_layer: float, max_layer: float):
        """Devuelve los segmentos de capa entre min_layer y max_layer."""
        return [(x0, y0, x1, y1, z) for (x0, y0, x1, y1, z) in self.segs if min_layer <= z <= max_layer]

    def remove_skirt(self):
        """Elimina líneas externas (skirt) fuera de la pieza."""
        if self.segs is None or len(self.segs) == 0:
            return
        minx, miny, maxx, maxy = self.minDimensions
        # Check if dimensions are valid
        if minx == maxx or miny == maxy:
            logging.warning("Cannot remove skirt, invalid specimen dimensions.")
            return

        original_count = len(self.segs)
        new_segs = [seg for seg in self.segs if not self.is_skirt(seg)]
        self.segs = np.array(new_segs)
        self.n_segs = len(self.segs)
        if self.n_segs > 0:
            self.seg_index = np.unique(self.segs[:, 4])
            self.n_layers = len(self.seg_index)
        else:
            self.seg_index = []
            self.n_layers = 0

        logging.info(f"Skirt removal: {original_count - self.n_segs} segments removed.")

    def is_skirt(self, seg: tuple) -> bool:
        """Verifica si un segmento está fuera de las dimensiones mínimas."""
        minx, miny, maxx, maxy = self.minDimensions
        # A simple check: if both points of the segment are outside the bounding box
        return (seg[0] < minx and seg[2] < minx) or \
               (seg[1] < miny and seg[3] < miny) or \
               (seg[0] > maxx and seg[2] > maxx) or \
               (seg[1] > maxy and seg[3] > maxy)

    def search_minorArea(self, delta, step, ejeMenor, ejeMayor):
        """Busca la menor área transversal en la pieza."""
        if self.segs is None or len(self.segs) == 0:
            return [], 0, np.inf, 0, []

        minx, maxx = self.minDimensions[0], self.minDimensions[2]
        middleP = minx + ((maxx - minx) / 2)
        limInf, limSup = middleP - delta / 2, middleP + delta / 2

        ptoCortes = np.arange(limInf, limSup, step)
        if len(ptoCortes) == 0:
            return [], 0, np.inf, 0, []

        minArea, minP = np.inf, 0
        minCut_solidArea = 0
        minCut_points = []
        areaCortes = []

        for p in ptoCortes:
            areaP, nCutPoints, areaSolida, cutPoints = self.apply_cutPoint(p, ejeMenor, ejeMayor)
            areaCortes.append((p, areaP, nCutPoints, areaSolida))
            if areaP < minArea:
                minArea, minP = areaP, p
                minCut_solidArea = areaSolida
                minCut_points = cutPoints

        logging.info(f"Menor área encontrada: {minArea:.4f} en corte {minP:.4f}")
        return areaCortes, minP, minArea, minCut_solidArea, minCut_points

    def apply_cutPoint(self, xcorte, ejeMenor, ejeMayor, verbose=False):
        """Aplica un corte vertical y calcula su área."""
        if self.segs is None or len(self.segs) == 0:
            return 0, 0, 0, []

        miny, maxy = self.minDimensions[1], self.minDimensions[3]
        cutSeg = [xcorte, miny, maxy]
        cutPoints = self.apply_cutSeg(cutSeg)

        if not cutPoints:
            return 0, 0, 0, []

        extremePoints = self.elispse_extremePoints(cutPoints, ejeMenor, ejeMayor)
        area_totalSolida = self.calcular_areaTotal_solida(extremePoints)

        minDist_y = CONFIG["minDist_y"]
        areaP = self.estimate_proportionalArea(cutPoints, area_totalSolida, minDist_y)

        if verbose:
            logging.info(f"Área proporcional del corte: {areaP:.4f}")

        return areaP, len(cutPoints), area_totalSolida, cutPoints

    def apply_cutSeg(self, cutSeg):
        """Calcula los puntos de intersección de un corte con los segmentos."""
        cutPoints = []
        if self.segs is None:
            return cutPoints

        for (x0, y0, x1, y1, z) in self.segs:
            if x0 == x1: # Vertical segment
                continue
            # Check if the cut line is within the x-range of the segment
            if min(x0, x1) <= cutSeg[0] <= max(x0, x1):
                mseg = (y1 - y0) / (x1 - x0)
                y = mseg * (cutSeg[0] - x0) + y0
                # Check if the intersection point is within the y-range of the cut segment (and the segment itself)
                if min(y0, y1) <= y <= max(y0, y1):
                     cutPoints.append([cutSeg[0], y, z])

        return cutPoints

    def estimate_proportionalArea(self, cutPoints, areaSolida, minDist_y):
        """Estima el área proporcional de un corte."""
        if not cutPoints:
            return 0
        y_coords = [p[1] for p in cutPoints]
        z_coords = [p[2] for p in cutPoints]

        miny, maxy = min(y_coords), max(y_coords)
        nPoints_y = round((maxy - miny) / minDist_y) if minDist_y > 0 else 0
        nPoints_z = len(np.unique(z_coords))
        nCutPoints = len(cutPoints)

        nGridPoints = nPoints_y * nPoints_z
        areaEstimada = (areaSolida * nCutPoints) / max(1, nGridPoints)
        return min(areaEstimada, areaSolida) if areaSolida > 0 else 0

    def elispse_extremePoints(self, cutPoints, ejeMenor, ejeMayor):
        """Genera puntos extremos simulando elipse de filamento."""
        extremePoints = []
        for x, y, z in cutPoints:
            extremePoints.extend([
                [x, y, z + ejeMenor],
                [x, y, z - ejeMenor],
                [x, y + ejeMayor, z],
                [x, y - ejeMayor, z]
            ])
        return extremePoints

    def calcular_areaTotal_solida(self, extremePoints):
        """Calcula área sólida total a partir de puntos extremos."""
        if not extremePoints:
            return 0
        y_coords = [p[1] for p in extremePoints]
        z_coords = [p[2] for p in extremePoints]
        miny, maxy = min(y_coords), max(y_coords)
        minz, maxz = min(z_coords), max(z_coords)
        a, b = maxy - miny, maxz - minz
        return a * b
