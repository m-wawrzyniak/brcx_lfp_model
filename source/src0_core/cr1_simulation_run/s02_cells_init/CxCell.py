from neuron import h
import math

import config_templates.conf01_simulation_parameters as conf01


class CxCell:
    """
    CX cell templates used in (Markram, 2015) do not have any stochastic mechanisms or random assignments.
    Templates initialized as they are, only define channels/conductances etc. only based on the morphology. And so,
    USING THE SAME MORPHOLOGY CREATES THE SAME CELL. To verify.

    Globals:
        NONE

    Attributes:
        cell_id (str): e.g. <0>
        cell_name (str): Shorten me-type e.g. 'L4_SS_cAD'
        h_cell (h.HocObject): Initialized from h.Template, instance of cell.
        center (tuple): (x, y, z) of soma[0]
        rotation (tuple): Rotation in (x, y, z) axis. In radians.
        spike_times (h.Vector): Times of spikes during simulation.
        spike_detector (h.NetCon): NetCon object used for spike detection.

    """

    def __init__(self, cell_id: str, cell_name: str, cell_temp, offset: tuple, rotation: tuple):
        self.cell_id = cell_id
        self.cell_name = cell_name
        self.h_cell = cell_temp(0)

        self._set_initial_potential()

        self.center = offset
        self.rotation = rotation

        self._rotate_cell()
        self._move_cell()

        self.spike_times = h.Vector()
        self.spike_detector = None

    def __repr__(self):
        return f"{self.cell_id}_{self.cell_name}"

    def _move_cell(self):
        """
        Sets all cell h.Sections by an offset provided during init.
        """
        for sec in self.h_cell.all:  # or cell.soma, cell.dend, etc. if section lists are available
            npt = int(h.n3d(sec=sec))  # number of 3D points
            if npt == 0:
                continue  # skip if no 3D points

            for i in range(npt):
                x = h.x3d(i, sec=sec) + self.center[0]
                y = h.y3d(i, sec=sec) + self.center[1]
                z = h.z3d(i, sec=sec) + self.center[2]
                diam = h.diam3d(i, sec=sec)
                h.pt3dchange(i, x, y, z, diam, sec=sec)

    def _rotate_cell(self):
        """
        Rotates all cell h.Sections along all three axles, provided during init.
        Has to be done prior to _move_cell(), as the rotation occurs with respect to (0, 0, 0)
        """
        rx, ry, rz = self.rotation

        for sec in self.h_cell.all:
            npt = int(h.n3d(sec=sec))
            if npt == 0:
                continue
            for i in range(npt):
                x = h.x3d(i, sec=sec)
                y = h.y3d(i, sec=sec)
                z = h.z3d(i, sec=sec)
                diam = h.diam3d(i, sec=sec)

                # X rotation
                y1 = y * math.cos(rx) - z * math.sin(rx)
                z1 = y * math.sin(rx) + z * math.cos(rx)
                x1 = x

                # Y rotation
                x2 = x1 * math.cos(ry) + z1 * math.sin(ry)
                z2 = -x1 * math.sin(ry) + z1 * math.cos(ry)
                y2 = y1

                # Z rotation
                x3 = x2 * math.cos(rz) - y2 * math.sin(rz)
                y3 = x2 * math.sin(rz) + y2 * math.cos(rz)
                z3 = z2

                h.pt3dchange(i, x3, y3, z3, diam, sec=sec)

    def _set_initial_potential(self):
        """
        Set the initial membrane potential for all sections and segments of a NEURON cell.

        Parameters:
        cell : h.Section or list of Sections
            The neuron section(s) to initialize.
        v_init : float
            Desired initial potential in mV.
        """

        rest_pot = conf01.RESTING_POTENTIALS[self.cell_name]
        for sec in self.h_cell.all:
            for seg in sec:
                seg.v = rest_pot

    def setup_spike_detector(self, threshold=0.0):
        """
        Sets up a NetCon to detect APs in the soma.

        Args:
            threshold (float): Potential threshold for detecting spike.
        """
        # You can customize this depending on the morphology (e.g., h_cell.soma[0](0.5))
        soma_seg = self.h_cell.soma[0](0.5)
        self.spike_detector = h.NetCon(soma_seg._ref_v, None, sec=self.h_cell.soma[0])
        self.spike_detector.threshold = threshold
        self.spike_detector.record(self.spike_times)

    def get_spike_times(self):
        """
        Returns spiking times after simulation as list.

        Returns:
            list: List of spikes. Probably sorted.
        """
        return list(self.spike_times)

    def start_recording_imem(self):
        """
        Record transmembrane current density at every segment (i_membrane_) and time.
        Also keep segment geometry for later LFP reconstruction.
        """
        # 1) Enable fast i_mem (gives seg._ref_i_membrane_)
        h.cvode.use_fast_imem(1)

        self.tvec = h.Vector()
        self.tvec.record(h._ref_t)

        self.imem_vecs = []  # list of h.Vector, one per segment (records mA/cm2)
        self.segs_meta = []  # list of dicts with geometry per segment

        for sec in self.h_cell.all:
            # cache 3D points for interpolation of segment midpoints
            n3d = int(h.n3d(sec=sec))
            pts = [(h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec),
                    h.arc3d(i, sec=sec)) for i in range(n3d)]

            for seg in sec:  # iterate segments (0..1 normalized)
                # record i_membrane_ (mA/cm2)
                v = h.Vector()
                v.record(seg._ref_i_membrane_)
                self.imem_vecs.append(v)

                # segment area in cm^2 (NEURON’s .area is um^2)
                area_um2 = h.area(seg.x, sec=sec)
                area_cm2 = area_um2 * 1e-8

                # midpoint xyz (interpolate along pt3d by arc length)
                x, y, z = self._seg_xyz_from_pt3d(sec, seg.x, pts)

                # approximate segment length in μm from area and diameter
                # (ok as metadata; line-source kernels can use true ends if needed)
                diam_um = seg.diam
                # length_um ≈ area / (π * diam)
                length_um = area_um2 / (math.pi * max(diam_um, 1e-9))

                self.segs_meta.append({
                    "sec": sec,
                    "x": float(seg.x),
                    "xyz_um": (float(x), float(y), float(z)),
                    "area_cm2": float(area_cm2),
                    "length_um": float(length_um),
                    "stype": sec.name().split(".")[1].split("[")[0]  # e.g. "soma", "dend", "axon", "apic"
                })

    def _seg_xyz_from_pt3d(self, sec, seg_x, pts):
        """Interpolate the 3D midpoint (μm) of a segment using pt3d and arc length."""
        if len(pts) == 0:
            # fallback: put at soma center if no 3D data
            return (0.0, 0.0, 0.0)
        L = sec.L
        target_arc = seg_x * L
        # find bracketing points by arc
        for i in range(len(pts) - 1):
            if pts[i][3] <= target_arc <= pts[i + 1][3]:
                a0, a1 = pts[i][3], pts[i + 1][3]
                t = 0.0 if a1 == a0 else (target_arc - a0) / (a1 - a0)
                x = pts[i][0] + t * (pts[i + 1][0] - pts[i][0])
                y = pts[i][1] + t * (pts[i + 1][1] - pts[i][1])
                z = pts[i][2] + t * (pts[i + 1][2] - pts[i][2])
                return x, y, z
        # if target beyond last arc (numerical edge), return last point
        return pts[-1][0], pts[-1][1], pts[-1][2]

    def get_imem_array(self):
        """
        Return numpy array of transmembrane current per segment in **Amperes**,
        shape: (n_segments, n_time). Converts from mA/cm2 * area(cm2).
        Also returns segment metadata and time vector (ms).
        """
        import numpy as np
        assert hasattr(self, "imem_vecs"), "Call start_recording_imem() before sim."
        nseg = len(self.imem_vecs)
        nt = len(self.tvec)

        imem_A = np.zeros((nseg, nt), dtype=np.float64)
        for k, vec in enumerate(self.imem_vecs):
            area = self.segs_meta[k]["area_cm2"]
            imem_A[k, :] = np.array(vec.as_numpy(), dtype=np.float64) * area * 1e-3

        t_ms = np.array(self.tvec.as_numpy(), dtype=np.float64)
        return imem_A, t_ms, self.segs_meta

    def populate_segs_meta_for_synapses(self):
        """Populate segment metadata needed for synapse placement without enabling recording."""
        self.segs_meta = []

        for sec in self.h_cell.all:
            sec_name = sec.name().split(".")[1].split("[")[0]  # stype: 'apic', 'dend', etc.
            for seg in sec:
                # Get 3D coordinates safely
                try:
                    x, y, z = self._seg_xyz_from_pt3d(sec, seg.x)
                except Exception:
                    # If no 3D points are available, fallback to (0,0,0)
                    x = y = z = 0.0

                diam_um = seg.diam
                area_um2 = h.area(seg.x, sec=sec)
                area_cm2 = area_um2 * 1e-8
                length_um = area_um2 / (math.pi * max(diam_um, 1e-9))

                self.segs_meta.append({
                    "sec": sec,
                    "x": seg.x,
                    "xyz_um": (x, y, z),
                    "area_cm2": area_cm2,
                    "length_um": length_um,
                    "stype": sec_name
                })