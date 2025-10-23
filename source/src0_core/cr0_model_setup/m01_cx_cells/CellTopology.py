import os
import math
import numpy as np
from neuron import h
import sys
import contextlib
h.load_file("import3d.hoc")

@contextlib.contextmanager
def suppress_stdout_stderr():
    with open(os.devnull, 'w') as fnull:
        fd_stdout = sys.__stdout__.fileno()
        fd_stderr = sys.__stderr__.fileno()

        def _redirect_fd(fd, target_fd):
            os.dup2(target_fd.fileno(), fd)

        orig_stdout = os.dup(fd_stdout)
        orig_stderr = os.dup(fd_stderr)

        try:
            _redirect_fd(fd_stdout, fnull)
            _redirect_fd(fd_stderr, fnull)
            yield
        finally:
            os.dup2(orig_stdout, fd_stdout)
            os.dup2(orig_stderr, fd_stderr)
            os.close(orig_stdout)
            os.close(orig_stderr)


class CellTopology:
    def __init__(self, morph_path: str, offset:tuple=(0,0,0), rotation:tuple=(0,0,0)):
        self.morph_file = morph_path
        self.offset = np.array(offset)
        self.rotation = np.array(rotation)

        self.axon = []
        self.dendrite = []
        self.soma = []

        self.axon_pts_with_sec = []
        self.dendrite_pts_with_sec = []
        self.soma_pts_with_sec = []

        self._load_and_transform()
        self.axon_len = self._compute_axon_length()

    def _load_and_transform(self):
        importer = h.Import3d_Neurolucida3()
        importer.quiet = 1
        importer.input(self.morph_file)

        gui = h.Import3d_GUI(importer, 0)
        gui.instantiate(None)

        section_types = set()
        for sec in h.allsec():
            sec_name = sec.name()
            stype = sec_name.split('[')[0]
            section_types.add(stype)

            n3d = int(h.n3d(sec=sec))
            if n3d <= 0:
                continue

            try:
                pts = np.array([[h.x3d(i, sec=sec),
                                 h.y3d(i, sec=sec),
                                 h.z3d(i, sec=sec)] for i in range(n3d)])
            except Exception as e:
                print(f"err reading pts in sec {sec_name}: {e}")
                continue

            pts = self._apply_rotation(pts)
            pts = self._apply_offset(pts)

            if stype == "axon":
                self.axon.extend(pts.tolist())
                self.axon_pts_with_sec.extend([(p[0], p[1], p[2], sec_name) for p in pts])
            elif stype == "soma":
                self.soma.extend(pts.tolist())
                self.soma_pts_with_sec.extend([(p[0], p[1], p[2], sec_name) for p in pts])
            elif stype in ("dend", "apic"):
                self.dendrite.extend(pts.tolist())
                self.dendrite_pts_with_sec.extend([(p[0], p[1], p[2], sec_name) for p in pts])


    def _apply_offset(self, pts:np.ndarray):
        return pts + self.offset

    def _apply_rotation(self, pts:np.ndarray):
        rx, ry, rz = self.rotation
        Rx = np.array([
            [1, 0, 0],
            [0, math.cos(rx), -math.sin(rx)],
            [0, math.sin(rx), math.cos(rx)]
        ])
        Ry = np.array([
            [math.cos(ry), 0, math.sin(ry)],
            [0, 1, 0],
            [-math.sin(ry), 0, math.cos(ry)]
        ])
        Rz = np.array([
            [math.cos(rz), -math.sin(rz), 0],
            [math.sin(rz), math.cos(rz), 0],
            [0, 0, 1]
        ])
        R = Rz @ Ry @ Rx
        return pts @ R.T

    def get_axon_points(self) -> np.ndarray:
        return np.array(self.axon)

    def get_dendrite_points(self) -> np.ndarray:
        return np.array(self.dendrite)

    def _compute_axon_length(self) -> float:
        if not self.axon or len(self.axon) < 2:
            return 0.0
        axon_array = np.array(self.axon)
        diffs = axon_array[1:] - axon_array[:-1]
        segment_lengths = np.linalg.norm(diffs, axis=1)
        return np.sum(segment_lengths)