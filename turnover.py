from __future__ import annotations

import os
from typing import Optional, Dict, Any, cast

import numpy as np
import pandas as pd

from cobaya.likelihood import Likelihood
from cobaya.conventions import Const
from cobaya.log import LoggedError


class TurnOverLikelihood(Likelihood):
    """
    Base Cobaya likelihood for turnover observable: DV(z) * k_TO(z).

    Subclasses or YAML-configured classes should supply:
      - path: directory containing measurement/covariance files
      - measurements_file: (columns: z value [observable optional])
      - cov_file: covariance matrix file
    """

    path: Optional[str] = None
    measurements_file: Optional[str] = None
    cov_file: Optional[str] = None
    aliases: Optional[list[str]] = None

    # Attribute types for linters/type-checkers
    data: pd.DataFrame
    cov: np.ndarray
    invcov: np.ndarray
    k: np.ndarray

    def initialize(self, k: Optional[np.ndarray] = None) -> None:
        if k is None:
            k = np.logspace(-3.0, 1.0, 400)
        self.k = np.array(k, dtype="f8")

        data_dir = self._resolve_data_dir()
        if not self.measurements_file or not self.cov_file:
            raise LoggedError(self.log, "Both 'measurements_file' and 'cov_file' must be provided (via YAML)")

        self.data = self._load_measurements(data_dir, self.measurements_file)
        self.cov = self._load_covariance(data_dir, self.cov_file)

        self.invcov = np.linalg.inv(np.atleast_2d(self.cov))
        self.logpdf = lambda _x: (lambda x_: -0.5 * x_.dot(self.invcov).dot(x_))(
            _x - self.data["value"].values
        )
        self.log.info("Initialized turnover likelihood with data dir=%s", data_dir)

    def _resolve_data_dir(self) -> str:
        if not getattr(self, "path", None):
            raise LoggedError(self.log, "No 'path' given to turnover data. Set the likelihood property in YAML.")
        return cast(str, self.path)

    def _load_measurements(self, directory: str, filename: str) -> pd.DataFrame:
        file_path = os.path.join(directory, filename)
        try:
            # Load the first two numeric columns (z, value)
            numeric = np.loadtxt(file_path, comments="#", usecols=(0, 1), dtype=float)
            if numeric.ndim == 1:
                numeric = numeric.reshape(1, -1)
            # Try to load an optional third column (observable) as strings
            observable: Optional[np.ndarray]
            try:
                observable = np.loadtxt(file_path, comments="#", usecols=(2,), dtype=str)
                if observable.ndim == 0:
                    observable = observable.reshape(1)
            except Exception:
                observable = None
            if observable is None:
                observable = np.full(numeric.shape[0], "DV_times_kT0", dtype=object)
            df = pd.DataFrame({
                "z": numeric[:, 0].astype(float),
                "value": numeric[:, 1].astype(float),
                "observable": observable,
            })
        except OSError:
            raise LoggedError(
                self.log,
                f"Couldn't find measurements file '{filename}' in folder '{directory}'. Check your paths.",
            )
        return df

    def _load_covariance(self, directory: str, filename: str) -> np.ndarray:
        file_path = os.path.join(directory, filename)
        try:
            cov = np.loadtxt(file_path)
        except OSError:
            raise LoggedError(
                self.log,
                f"Couldn't find covariance file '{filename}' in folder '{directory}'. Check your paths.",
            )
        return np.atleast_2d(cov)

    @property
    def rs(self):
        return self.provider.get_param("rdrag")

    def get_requirements(self):
        zs = {obs: self.data.loc[self.data["observable"] == obs, "z"].values for obs in self.data["observable"].unique()}
        
        theory_reqs = {
            "DV_times_kT0": {
                "angular_diameter_distance": {"z": zs.get("DV_times_kT0", None)},
                "Hubble": {"z": zs.get("DV_times_kT0", None)},
                # 'fourier': {'z': zs.get("DV_times_kT0", None),  'of': [('delta_cb', 'delta_cb')]}},
                # 'Pk_interpolator': {'z': zs.get("DV_times_kT0", None), 'k_max': self.k.max(), 'vars_pairs': [('delta_cb', 'delta_cb')],"nonlinear": True}},
                'Pk_interpolator': {'z': zs.get("DV_times_kT0", None), 'k_max': self.k.max(), 'vars_pairs': [('delta_tot', 'delta_tot')],
                                    # "nonlinear": False,
                                    }},
        }
        requisites = {}
        # if self.has_type:
        for observable in self.data["observable"].unique():
            for req, req_values in theory_reqs[observable].items():
                if req not in requisites:
                    requisites[req] = req_values
                else:
                    if isinstance(req_values, dict):
                        for k, v in req_values.items():
                            if v is not None:
                                requisites[req][k] = np.unique(
                                    np.concatenate((requisites[req][k], v)))
        return requisites

    def find_turn_over(self, pk_interpolator, z: float):
        # Find turnover via local parabolic interpolation around the maximum
        power_values = pk_interpolator.P(z, pk_interpolator.k)
        index_max = np.argmax(power_values)
        logk = np.log10(pk_interpolator.k[index_max - 1 : index_max + 2])
        logpk = np.log10(power_values[index_max - 1 : index_max + 2])
        c0 = logpk[0] / (logk[0] - logk[1]) / (logk[0] - logk[2])
        c1 = logpk[1] / (logk[1] - logk[0]) / (logk[1] - logk[2])
        c2 = logpk[2] / (logk[2] - logk[0]) / (logk[2] - logk[1])
        a = c0 + c1 + c2
        logk0 = (c0 * (logk[1] + logk[2]) + c1 * (logk[0] + logk[2]) + c2 * (logk[0] + logk[1])) / (2 * a)
        if not (a <= 0.0 and logk[0] <= logk0 <= logk[2]):
            raise LoggedError(self.log, "Parabolic interpolation around turnover failed.")
        k0 = 10 ** logk0
        return k0, pk_interpolator.P(z, k0)

    def theory_fun(self, z: float, observable: str) -> float:
        if observable == "DV_times_kT0":
            DV = np.cbrt(
                ((1 + z) * self.provider.get_angular_diameter_distance(z)) ** 2
                * Const.c_km_s
                * z
                / self.provider.get_Hubble(z, units="km/s/Mpc")
            )
            pk_interp = self.provider.get_Pk_interpolator()
            kTO, _ = self.find_turn_over(pk_interp, z)
            DV_times_kTO = DV * kTO
            return float(DV_times_kTO)
        raise LoggedError(self.log, f"Unsupported observable '{observable}'")

    def logp(self, **params_values):
        theory_values = np.array(
            [
                self.theory_fun(z, obs)
                for z, obs in zip(
                    self.data["z"].values, self.data["observable"].values
                )
            ]
        )
        if self.is_debug():
            for i, (z, obs, theo) in enumerate(
                zip(
                    self.data["z"].values,
                    self.data["observable"].values,
                    theory_values,
                )
            ):
                self.log.debug(
                    "%s at z=%g : %g (theo) ; %g (data)", obs, z, theo, self.data.iloc[i, 1]
                )
        return self.logpdf(theory_values) 