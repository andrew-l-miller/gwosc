from __future__ import annotations

from dataclasses import dataclass
from typing import BinaryIO, Dict, Any, Tuple, Optional

import numpy as np


@dataclass
class SFDBheader:
    # matches MATLAB piahead.* fields (plus eof)
    eof: int = 0
    def __repr__(self):
        lines = ["SFDBheader:"]
        for k, v in self.__dict__.items():
            lines.append(f"  {k:<15} = {v}")
        return "\n".join(lines)


def _read_exact(f: BinaryIO, nbytes: int) -> bytes:
    b = f.read(nbytes)
    if len(b) != nbytes:
        raise EOFError("Unexpected EOF while reading sfdb09 block")
    return b


def _read_scalar(f: BinaryIO, dtype: np.dtype) -> Any:
    """Read one scalar of dtype from file (little-endian assumed)."""
    dt = np.dtype(dtype)
    b = f.read(dt.itemsize)
    if len(b) == 0:
        return None  # signals EOF at block start (matches MATLAB count==0)
    if len(b) != dt.itemsize:
        raise EOFError("Unexpected EOF while reading scalar")
    return np.frombuffer(b, dtype=dt, count=1)[0]


def _read_array(f: BinaryIO, dtype: np.dtype, count: int) -> np.ndarray:
    """Read `count` items of dtype from file into a numpy array."""
    dt = np.dtype(dtype)
    b = _read_exact(f, dt.itemsize * int(count))
    return np.frombuffer(b, dtype=dt, count=int(count))


def sfdb_read_an_FFT(
    fid: BinaryIO,
    det: int = 0,
    want_sft: bool = True,
    endian: str = "<",
) -> Tuple[SFDBheader | int, np.ndarray | int, np.ndarray | int, Optional[np.ndarray] | int]:
    """
    Python translation of MATLAB sfdb_read_an_FFT(fid, det).

    Parameters
    ----------
    fid : BinaryIO
        Open file handle in binary mode ('rb').
    det : int
        If 1, check detector type and use "bar old" extra fields when detector==0.
        Default 0.
    want_sft : bool
        If False, skip reading SFT data (seek past it) and return sft=0.
        This mirrors MATLAB's `if nargout < 4` branch.
    endian : str
        '<' little-endian (default) or '>' big-endian. MATLAB code reads an
        endian marker but doesn't use it; choose the correct one here.

    Returns
    -------
    sfdb_header : SFDBheader or 0
        Header object; if EOF at block start, returns 0 (to match MATLAB behavior).
    tps : np.ndarray or 0
    sps : np.ndarray or 0
    sft : np.ndarray (complex64) or 0
    """
    h = SFDBheader()
    # define dtypes with explicit endianness
    f64 = np.dtype(endian + "f8")
    f32 = np.dtype(endian + "f4")
    i32 = np.dtype(endian + "i4")

    # --- read endian marker (double) ---
    endian_val = _read_scalar(fid, f64)
    if endian_val is None:
        # end of file at block start
        h.eof = 1
        return 0, 0, 0, 0

    # populate header fields
    hd: Dict[str, Any] = {"eof": 0, "endian": float(endian_val)}

    hd["detector"] = int(_read_scalar(fid, i32))
    hd["gps_sec"] = int(_read_scalar(fid, i32))
    hd["gps_nsec"] = int(_read_scalar(fid, i32))
    hd["tbase"] = float(_read_scalar(fid, f64))
    hd["firstfrind"] = int(_read_scalar(fid, i32))
    hd["nsamples"] = int(_read_scalar(fid, i32))

    hd["red"] = int(_read_scalar(fid, i32))
    hd["typ"] = int(_read_scalar(fid, i32))
    hd["n_flag"] = float(_read_scalar(fid, f32))
    hd["einstein"] = float(_read_scalar(fid, f32))
    hd["mjdtime"] = float(_read_scalar(fid, f64))
    hd["nfft"] = int(_read_scalar(fid, i32))
    hd["wink"] = int(_read_scalar(fid, i32))
    hd["normd"] = float(_read_scalar(fid, f32))
    hd["normw"] = float(_read_scalar(fid, f32))
    hd["frinit"] = float(_read_scalar(fid, f64))
    hd["tsamplu"] = float(_read_scalar(fid, f64))
    hd["deltanu"] = float(_read_scalar(fid, f64))

    # conditional extra fields
    if hd["detector"] == 0 and det == 1:
        # bar old
        hd["frcal"] = float(_read_scalar(fid, f64))
        hd["freqm"] = float(_read_scalar(fid, f64))
        hd["freqp"] = float(_read_scalar(fid, f64))
        hd["taum"] = float(_read_scalar(fid, f64))
        hd["taup"] = float(_read_scalar(fid, f64))
        # NOTE: MATLAB doesn't read lavesp in this branch, but later uses sfdb_header.lavesp.
        # In practice, many files likely follow the "else" layout. If you truly have
        # bar-old blocks, you'll need the exact format spec for where lavesp lives.
        hd["lavesp"] = 0
    else:
        hd["vx_eq"] = float(_read_scalar(fid, f64))
        hd["vy_eq"] = float(_read_scalar(fid, f64))
        hd["vz_eq"] = float(_read_scalar(fid, f64))
        hd["px_eq"] = float(_read_scalar(fid, f64))
        hd["py_eq"] = float(_read_scalar(fid, f64))
        hd["pz_eq"] = float(_read_scalar(fid, f64))
        hd["n_zeroes"] = int(_read_scalar(fid, i32))
        hd["sat_howmany"] = float(_read_scalar(fid, f64))
        hd["spare1"] = float(_read_scalar(fid, f64))
        hd["spare2"] = float(_read_scalar(fid, f64))
        hd["spare3"] = float(_read_scalar(fid, f64))
        hd["spare4"] = float(_read_scalar(fid, f32))
        hd["spare5"] = float(_read_scalar(fid, f32))
        hd["spare6"] = float(_read_scalar(fid, f32))
        hd["lavesp"] = int(_read_scalar(fid, i32))
        hd["spare8"] = int(_read_scalar(fid, i32))
        hd["spare9"] = int(_read_scalar(fid, i32))

    # --- tps / sps ---
    ltps = int(hd["red"])
    if int(hd.get("lavesp", 0)) > 0:
        lsps = int(hd["lavesp"])
        tps = _read_array(fid, f32, lsps).astype(np.float32, copy=False)
    else:
        tps = _read_array(fid, f32, ltps).astype(np.float32, copy=False)
        lsps = int(hd["nsamples"] // ltps)

    sps = _read_array(fid, f32, lsps).astype(np.float32, copy=False)

    # --- sft (complex) or skip ---
    n = int(hd["nsamples"])
    n2 = 2 * n  # interleaved re/im float32 pairs

    if not want_sft:
        # skip 4*n2 bytes (float32)
        fid.seek(4 * n2, 1)  # relative seek
        sft = 0
    else:
        cdat = _read_array(fid, f32, n2).astype(np.float32, copy=False)
        # interleaved -> complex
        sft = cdat[0:n2:2] + 1j * cdat[1:n2:2]
        sft = sft.astype(np.complex64, copy=False)

    # package header into dataclass-like object (keeping MATLAB field names)
    sfdb_header = SFDBheader()
    for k, v in hd.items():
        setattr(sfdb_header, k, v)

    return sfdb_header, tps, sps, sft


# --- Example usage ---
# with open("your.sfdb09", "rb") as f:
#     while True:
#         sfdb_header, tps, sps, sft = sfdb_read_an_FFT(f, det=0, want_sft=True, endian="<")
#         if sfdb_header == 0:
#             break
#         # do something with sfdb_header, tps, sps, sft
