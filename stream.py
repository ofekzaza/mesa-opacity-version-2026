import os
import mesa_reader as mr
import utils


def stream_profiles(logs_path: str, collectors: list) -> int:
    """Single streaming pass through all profile files.

    Calls collector.consume(prof) on each profile in order (number ordering,
    with ``profile_final.data`` appended last when it exists).
    Each MesaData object is discarded after the call — only the
    numpy arrays extracted by collectors remain in memory.

    Returns the number of profiles processed.
    """
    count = 0
    for f_path in utils.list_profiles(logs_path):
        try:
            prof = mr.MesaData(f_path)
        except Exception as e:
            print(f"    stream: could not read {f_path}: {e}")
            continue
        for collector in collectors:
            collector.consume(prof)
        count += 1
        # prof goes out of scope here — MesaData is discarded
    return count
