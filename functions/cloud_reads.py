"""Read-only, in-memory GCS access. Configure complete URIs privately."""
import subprocess


def read_bytes(uri, billing_project, byte_start=None, byte_end=None):
    if not uri.startswith('gs://') or not billing_project:
        raise ValueError('A configured GCS URI and billing project are required')
    command = ['gsutil', '-u', billing_project, 'cat']
    if byte_start is not None:
        if byte_end is None or byte_start < 0 or byte_end < byte_start:
            raise ValueError('Invalid inclusive byte range')
        command += ['-r', f'{byte_start}-{byte_end}']
    process = subprocess.run(command + [uri], capture_output=True)
    if process.returncode:
        raise RuntimeError('Controlled read failed; inspect within the authorized environment')
    if byte_start is not None and len(process.stdout) != byte_end - byte_start + 1:
        raise RuntimeError('Incomplete controlled byte-range read')
    return process.stdout


def bed_record(uri, billing_project, sample_count, variant_index):
    if sample_count <= 0 or variant_index < 0:
        raise ValueError('Invalid BED dimensions')
    size = (sample_count + 3) // 4
    start = 3 + variant_index * size
    return read_bytes(uri, billing_project, start, start + size - 1)


def tabix_region(configured_file, region, executable='tabix'):
    """Returns text in memory; paths and stderr are not printed."""
    process = subprocess.run([executable, str(configured_file), region], capture_output=True, text=True)
    if process.returncode:
        raise RuntimeError('Controlled region query failed')
    return process.stdout
