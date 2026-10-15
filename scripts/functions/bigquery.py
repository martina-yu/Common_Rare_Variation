"""Explicit-client query helpers; importing this module performs no query."""

import time


def timed_query(client, label, query):
    started = time.time()
    frame = client.query(query).to_dataframe()
    print(f'[{label}] elapsed={time.time()-started:.1f}s  n_rows={len(frame)}')
    return frame
