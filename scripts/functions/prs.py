# Documentation: scripts/README.md.
import numpy as np


def _orient_effect_allele(row):
    """orient effect allele."""
    if row["effect_allele"] == row["A1"]:
        return "A1"
    if row["effect_allele"] == row["A2"]:
        return "A2"
    return "MISMATCH"


def _decode_bed_block(raw_bytes, N):
    byte_arr = np.frombuffer(raw_bytes, dtype=np.uint8)
    codes = np.zeros((len(byte_arr), 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k] = (byte_arr >> (2 * k)) & 0b11
    codes = codes.reshape(-1)[:N]

    dosage_A1 = np.select([codes == 0, codes == 2, codes == 3], [2.0, 1.0, 0.0], default=np.nan)
    return dosage_A1


def _encode_bed_bytes(dosage_mat):
    """encode bed bytes."""
    n_var, n_samp = dosage_mat.shape
    codes = np.ones(dosage_mat.shape, dtype=np.uint8)
    codes[dosage_mat == 2] = 0
    codes[dosage_mat == 1] = 2
    codes[dosage_mat == 0] = 3
    record_size = -(-n_samp // 4)
    pad = (-n_samp) % 4
    if pad:
        codes = np.concatenate([codes, np.zeros((n_var, pad), dtype=np.uint8)], axis=1)
    codes = codes.reshape(n_var, -1, 4)
    byte_vals = (
        codes[:, :, 0] | (codes[:, :, 1] << 2) | (codes[:, :, 2] << 4) | (codes[:, :, 3] << 6)
    ).astype(np.uint8)
    out = bytearray(3 + n_var * record_size)
    out[0:3] = bytes([0x6C, 0x1B, 0x01])
    out[3:] = byte_vals.tobytes()
    return bytes(out), record_size
