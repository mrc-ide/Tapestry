import os
from unravel.sample.main import sample


def population(input_dir: str) -> None:
    sample_dirs = sorted([
        f"{input_dir}/{d}"
        for d in os.listdir(input_dir)
        if os.path.isdir(f"{input_dir}/{d}")
    ])
    print(f"Found {len(sample_dirs)} samples to process.")
    for sample_dir in sample_dirs:
        sample(sample_dir)
    print("Done.")

