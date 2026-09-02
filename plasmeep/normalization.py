"""
Stored incident-coefficient normalization for modal S-parameters.

StudentGuide / PMM-Design philosophy: obtain c_i^+ from a reference
run and store it. Later physical runs compute

    S_ji = c_j^-(out) / c_i^+(norm)

Never recompute the denominator from the physical simulation.
Keep complex amplitude. Do not convert to power or dB.

c_i^+ is the toward-aperture coefficient at the driven port.
c_j^- is the away-from-aperture coefficient at a receive port.
S_ii uses the away-from-aperture (reflected) coefficient at the
source port over the same stored c_i^+.
"""

import numpy as np


NEAR_ZERO = 1e-18


def _as_complex_1d(values, name):
    arr = np.asarray(values, dtype=complex).reshape(-1)
    if arr.size == 0:
        raise ValueError("{} must be a nonempty complex vector.".format(name))
    if np.any(~np.isfinite(arr.real)) or np.any(~np.isfinite(arr.imag)):
        raise ValueError("{} must be finite.".format(name))
    return arr


def _as_freq_1d(values):
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size == 0 or np.any(~np.isfinite(arr)) or np.any(arr <= 0.0):
        raise ValueError(
            "frequencies must be a nonempty list of finite "
            "MEEP cyclic frequencies > 0."
        )
    return arr


def _require_int_band(eig_band):
    try:
        eig_band = int(eig_band)
    except (TypeError, ValueError) as exc:
        raise ValueError("eig_band must be an integer >= 1.") from exc
    if eig_band < 1:
        raise ValueError("eig_band must be an integer >= 1.")
    return eig_band


class IncidentNormalization:
    """
    Cached c_i^+(norm) for one source port and one eigenmode setting.

    Attributes are the identity of the stored reference. apply / normalize
    refuse silently mixing a different source, frequency list, band, or
    parity.
    """

    def __init__(
        self,
        coefficients,
        frequencies,
        source_id,
        eig_band=1,
        eig_parity=None,
        direction="toward_aperture",
    ):
        if source_id is None or source_id == "":
            raise ValueError("source_id is required to identify the driven port.")
        self.source_id = source_id
        self.coefficients = _as_complex_1d(coefficients, "coefficients")
        self.frequencies = _as_freq_1d(frequencies)
        if self.coefficients.size != self.frequencies.size:
            raise ValueError(
                "coefficients and frequencies must have the same length "
                "(got {} coeffs, {} frequencies).".format(
                    self.coefficients.size, self.frequencies.size
                )
            )
        if np.any(np.abs(self.coefficients) < NEAR_ZERO):
            raise ValueError(
                "stored incident coefficient is zero or nearly zero at "
                "one or more frequencies; cannot use it as a denominator."
            )
        self.eig_band = _require_int_band(eig_band)
        if eig_parity is None:
            raise ValueError(
                "eig_parity is required so a later apply() cannot mix modes."
            )
        try:
            self.eig_parity = int(eig_parity)
        except (TypeError, ValueError) as exc:
            raise ValueError("eig_parity must be a MEEP parity flag.") from exc
        self.direction = direction

    @classmethod
    def from_port(cls, port, coefficients, source_id, frequencies=None):
        """Build a store from an Add_Port descriptor and measured c_i^+."""
        if not isinstance(port, dict):
            raise ValueError("port must be a dict returned by Add_Port.")
        if frequencies is None:
            if port.get("frequencies") is not None:
                frequencies = port["frequencies"]
            elif port.get("frequency") is not None:
                frequencies = [port["frequency"]]
            else:
                raise ValueError(
                    "frequencies must be provided or stored on the port."
                )
        return cls(
            coefficients,
            frequencies,
            source_id=source_id,
            eig_band=port.get("eig_band", 1),
            eig_parity=port.get("eig_parity"),
            direction=port.get("direction", "toward_aperture"),
        )

    def _check_compatible(self, source_id, frequencies, eig_band, eig_parity):
        if source_id != self.source_id:
            raise ValueError(
                "normalization source_id {!r} does not match "
                "requested source_id {!r}.".format(self.source_id, source_id)
            )
        freqs = _as_freq_1d(frequencies)
        if freqs.size != self.frequencies.size or not np.allclose(
            freqs, self.frequencies, rtol=0.0, atol=1e-12
        ):
            raise ValueError(
                "frequency list does not match the stored normalization "
                "run; refusing to mix incompatible settings."
            )
        eig_band = _require_int_band(eig_band)
        if eig_band != self.eig_band:
            raise ValueError(
                "eig_band {} does not match stored normalization "
                "eig_band {}.".format(eig_band, self.eig_band)
            )
        try:
            eig_parity = int(eig_parity)
        except (TypeError, ValueError) as exc:
            raise ValueError("eig_parity must be a MEEP parity flag.") from exc
        if eig_parity != self.eig_parity:
            raise ValueError(
                "eig_parity does not match the stored normalization run."
            )
        return freqs

    def normalize(
        self,
        outgoing,
        source_id,
        frequencies=None,
        eig_band=None,
        eig_parity=None,
    ):
        """
        Complex S = outgoing / c_i^+(norm).

        outgoing is c_j^- (transmission) or the source-port backward
        coefficient (reflection). Same stored denominator in both cases.
        """
        if frequencies is None:
            frequencies = self.frequencies
        if eig_band is None:
            eig_band = self.eig_band
        if eig_parity is None:
            eig_parity = self.eig_parity
        self._check_compatible(source_id, frequencies, eig_band, eig_parity)
        outgoing = _as_complex_1d(outgoing, "outgoing")
        if outgoing.size != self.coefficients.size:
            raise ValueError(
                "outgoing coefficient length {} does not match stored "
                "incident length {}.".format(
                    outgoing.size, self.coefficients.size
                )
            )
        return outgoing / self.coefficients

    def s_parameter(self, outgoing, **kwargs):
        """Alias for normalize(); S_ji or S_ii from stored c_i^+."""
        return self.normalize(outgoing, **kwargs)


class NormalizationTable:
    """One stored IncidentNormalization per source port."""

    def __init__(self):
        self._by_source = {}

    def add(self, norm):
        if not isinstance(norm, IncidentNormalization):
            raise TypeError("norm must be an IncidentNormalization.")
        if norm.source_id in self._by_source:
            raise ValueError(
                "normalization for source_id {!r} is already stored; "
                "refusing to overwrite.".format(norm.source_id)
            )
        self._by_source[norm.source_id] = norm
        return norm

    def get(self, source_id):
        try:
            return self._by_source[source_id]
        except KeyError as exc:
            raise ValueError(
                "no stored normalization for source_id {!r}.".format(source_id)
            ) from exc

    def normalize(self, outgoing, source_id, **kwargs):
        return self.get(source_id).normalize(
            outgoing, source_id=source_id, **kwargs
        )

    def __contains__(self, source_id):
        return source_id in self._by_source

    def __len__(self):
        return len(self._by_source)
