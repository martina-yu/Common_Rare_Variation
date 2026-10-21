"""Input checks without printing identifiers, counts or values."""
import numpy as np
import pandas as pd


def unique_ids(frame, column='person_id'):
    if frame[column].isna().any() or frame[column].duplicated().any():
        raise ValueError('Missing or duplicate identifiers')


def same_membership(left, right, column='person_id'):
    unique_ids(left, column)
    unique_ids(right, column)
    if set(left[column]) != set(right[column]):
        raise ValueError('Frozen membership changed')


def aligned(frame, index):
    """No implicit dropping; duplicate identifiers or mismatched membership fail."""
    if frame.index.has_duplicates or index.has_duplicates or frame.index.isna().any():
        raise ValueError('Invalid index')
    if set(frame.index) != set(index):
        raise ValueError('Input membership differs')
    return frame.loc[index].copy()


def finite(values):
    x = np.asarray(values, dtype=float)
    if not np.isfinite(x).all():
        raise ValueError('Required numeric input is missing or nonfinite')
    return x


def standardize(values, ddof=1, reference=None):
    x = finite(values)
    r = x if reference is None else finite(reference)
    mean, sd = r.mean(axis=0), r.std(axis=0, ddof=ddof)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError('Cannot standardize constant input')
    return (x - mean) / sd, {'mean': mean, 'sd': sd, 'ddof': ddof}
