# -*- coding: utf-8 -*-
"""
Agent-based model for voxel-level alpha-syn spreading on mouse connectome.

By: Vikram Nathan, 03/24/2026
Modified from original model used in Rahayel et al., 2022. to remove dependencies on voxel size 
and instead use voxel-level data.

OPTIMIZED for time complexity.

"""

import numpy as np

class AgentBasedModel:
    _HISTORY_INIT_CAP = 256  # initial row capacity for history buffers; doubles on overflow

    def __init__(
            self, weights, distance,
            sources, targets, dt=0.1
    ):
        """
        construct the object. Note all the params are those that are fixed
        we don't need to fit them in model-fitting stage
        :param weights: ndarray-like, row-source column-target
        :param distance: ndarray-like, same as above; the index must match
        :param sources: list, labels of source voxels
        :param targets: list, labels of target voxels
        :param dt: float
        """
        self.adj = np.where(weights != 0, 1, 0)
        self.weights = np.copy(weights)
        self.distance = np.copy(distance)

        # voxel counts
        self.sources = sources
        self.n_sources = len(sources)
        self.targets = targets
        self.n_targets = len(targets)
        
        # spread rates (load pre-normalized data)
        self.spread_weights = self.weights / \
            self.weights.sum(axis=1)[:, np.newaxis]
        # index of diagonal matrix -- may be useful
        self.diagonal = np.eye(self.n_sources, self.n_targets)
        self.voxel_to_edge_weights = np.copy(self.spread_weights)
        self.voxel_to_edge_weights[self.diagonal == 1] = 0

        self.adj_dist = self.adj * self.distance
        self.dist_inv = np.zeros(self.distance.shape)
        self.dist_inv[self.adj_dist != 0] = 1 / \
            self.adj_dist[self.adj_dist != 0]
        self.edge_to_voxel_weights = np.copy(self.dist_inv)

        self.dt = dt

        # rates
        self.growth_rate = 0
        self.clearance_rate = 0.5 
        self.trans_rate = 0.005 ###similar to the mean transition rate in regional SIR model of 0.002209

        # initialize population
        n = max(self.n_targets, self.n_sources)
        self.s_voxel = np.zeros(n)
        self.i_voxel = np.zeros(n)
        self.s_edge = np.zeros(self.adj.shape)
        self.i_edge = np.zeros(self.adj.shape)

        # Pre-allocated reusable buffer for spread steps — avoids per-step
        # np.append(..., np.zeros(...)) allocation. Trailing elements beyond
        # n_sources stay 0 permanently (set at init, never written).
        self._v2e_row_buf = np.zeros(n)

        # Pre-allocated history buffers — grow by doubling instead of np.append
        cap = self._HISTORY_INIT_CAP
        self._s_voxel_hist_buf = np.empty((cap, n))
        self._i_voxel_hist_buf = np.empty((cap, n))
        self._hist_len = 0

    # ------------------------------------------------------------------
    # History properties — expose only the recorded portion, preserving
    # the same array interface as the original (read and write).
    # ------------------------------------------------------------------

    @property
    def s_voxel_history(self):
        return self._s_voxel_hist_buf[:self._hist_len]

    @property
    def i_voxel_history(self):
        return self._i_voxel_hist_buf[:self._hist_len]

    @s_voxel_history.setter
    def s_voxel_history(self, value):
        value = np.atleast_2d(value)
        n = self.s_voxel.shape[0]
        cap = max(value.shape[0], self._HISTORY_INIT_CAP)
        self._s_voxel_hist_buf = np.empty((cap, n))
        self._s_voxel_hist_buf[:value.shape[0]] = value
        self._hist_len = value.shape[0]

    @i_voxel_history.setter
    def i_voxel_history(self, value):
        value = np.atleast_2d(value)
        n = self.i_voxel.shape[0]
        cap = max(value.shape[0], self._HISTORY_INIT_CAP)
        self._i_voxel_hist_buf = np.empty((cap, n))
        self._i_voxel_hist_buf[:value.shape[0]] = value
        self._hist_len = value.shape[0]

    def _grow_history_buffers(self):
        n = self.s_voxel.shape[0]
        new_cap = self._s_voxel_hist_buf.shape[0] * 2
        for attr in ('_s_voxel_hist_buf', '_i_voxel_hist_buf'):
            old = getattr(self, attr)
            new = np.empty((new_cap, n))
            new[:old.shape[0]] = old
            setattr(self, attr, new)

    # ------------------------------------------------------------------

    def set_spread_process(self, v):
        """v spread"""
        self.edge_to_voxel_weights *= v

    def set_growth_process(self, growth_rate):
        """growth_rate involves params to be  fitted"""
        self.growth_rate = np.array(growth_rate)

    def set_clearance_process(self, clearance_rate):
        """clearance_rate involves params to be fitted"""
        self.clearance_rate = np.array(clearance_rate)

    def set_trans_process(self, trans_rate):
        """trans_rate involves params to be fitted"""
        self.trans_rate *= trans_rate

    def update_spread_process(self, v_scale=1, spread_scale=1): 
        """
        slow_down, exit_down involve params to be fitted
        :param v_scale: float, by which the spreading in edges is discounted
        :param spread_scale: float, by which the probability of leaving voxels
            is discounted
        """
        self.edge_to_voxel_weights *= v_scale
        self.voxel_to_edge_weights *= spread_scale

    def s_spread_step(self):
        """spread step in each time step"""
        voxel_to_edge = self.voxel_to_edge_weights * self.dt * \
            self.s_voxel[:self.n_sources][:, np.newaxis]

        edge_to_voxel = self.edge_to_voxel_weights * \
            self.s_edge * self.dt

        # update edges and voxels
        self.s_edge += voxel_to_edge - edge_to_voxel

        # Reuse pre-allocated buffer instead of np.append(..., np.zeros(...))
        self._v2e_row_buf[:self.n_sources] = voxel_to_edge.sum(axis=1)
        self.s_voxel += edge_to_voxel.sum(axis=0) - self._v2e_row_buf

    def i_spread_step(self): #same as s_spread_step but after injection 
        """spread step in each time step"""
        voxel_to_edge = self.voxel_to_edge_weights * self.dt * \
            self.i_voxel[:self.n_sources][:, np.newaxis]

        edge_to_voxel = self.edge_to_voxel_weights * \
            self.i_edge * self.dt

        # update edges and voxels
        self.i_edge += voxel_to_edge - edge_to_voxel

        # Reuse pre-allocated buffer instead of np.append(..., np.zeros(...))
        self._v2e_row_buf[:self.n_sources] = voxel_to_edge.sum(axis=1)
        self.i_voxel += edge_to_voxel.sum(axis=0) - self._v2e_row_buf

    def growth_step(self):
        self.s_voxel += self.growth_rate * self.dt

    def clearance_step(self):
        """clearance step"""
        factor = np.exp(-self.clearance_rate * self.dt)
        self.s_voxel *= factor
        self.i_voxel *= factor

    def injection(self, seed, amount=1):
        """inject infected agents into seed"""
        self.i_voxel[seed] = amount

    def trans_step(self):
        infected = self.s_voxel * (
            1 - np.exp(-self.trans_rate * self.dt * self.i_voxel)
        )

        self.s_voxel -= infected
        self.i_voxel += infected

    def record_history_voxel(self):
        if self._hist_len >= self._s_voxel_hist_buf.shape[0]:
            self._grow_history_buffers()
        self._s_voxel_hist_buf[self._hist_len] = self.s_voxel
        self._i_voxel_hist_buf[self._hist_len] = self.i_voxel
        self._hist_len += 1