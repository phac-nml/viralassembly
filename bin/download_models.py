#!/usr/bin/env python3

"""
Rewritten version of artic's get_models.py - https://github.com/artic-network/fieldbioinformatics/blob/master/artic/get_models.py

Changes are to allow pulling one model at a time along with some small simplifications to the structure of the MODELs
"""
import argparse
import logging
import os
import requests
import tarfile
import shutil
import sys

from pathlib import Path

MODELS = [
    "r1041_e82_260bps_fast_g632", "r1041_e82_260bps_hac_g632",
    "r1041_e82_260bps_hac_v400", "r1041_e82_260bps_hac_v410",
    "r1041_e82_260bps_sup_g632", "r1041_e82_260bps_sup_v400",
    "r1041_e82_260bps_sup_v410", "r1041_e82_400bps_fast_g615",
    "r1041_e82_400bps_fast_g632", "r1041_e82_400bps_hac_g615",
    "r1041_e82_400bps_hac_g632", "r1041_e82_400bps_hac_v400",
    "r1041_e82_400bps_hac_v410", "r1041_e82_400bps_hac_v420",
    "r1041_e82_400bps_hac_v430", "r1041_e82_400bps_hac_v500",
    "r1041_e82_400bps_hac_v520", "r1041_e82_400bps_sup_g615",
    "r1041_e82_400bps_sup_v400", "r1041_e82_400bps_sup_v410",
    "r1041_e82_400bps_sup_v420", "r1041_e82_400bps_sup_v430",
    "r1041_e82_400bps_sup_v500", "r1041_e82_400bps_sup_v520",
    "r104_e81_hac_g5015", "r104_e81_sup_g5015"
]

# Really basic logging setup
logging.basicConfig(
    format="[%(asctime)s] %(levelname)s: %(message)s",
    level=logging.INFO,
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def download_file(url: str, local_path: Path) -> None:
    """
    Purpose
    -------
    Download given URL to specific local path

    Parameters
    ----------
    url: str
        The URL to stream data from
    local_path: Path
        The path to save the streamed data
    """
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)


def download_model(download_dir: Path, model: str) -> None:
    """
    Purpose
    -------
    Download model to the given directory based on the oxfordnanoportal model formatting

    Parameters
    ----------
    download_dir: path
        The path to download the model to
    model: str
        The name of the model to download
    """
    # All models follow the same format for the moment
    #  From https://github.com/nanoporetech/rerio/tree/master/clair3_models
    model_fname = f"{model}.tar.gz"
    model_url = f"https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/{model_fname}"

    # Check if downloaded already and if not, setup, download, and extract
    final_model_path = Path(download_dir, model)
    if (final_model_path.exists()) and (len(list(final_model_path.iterdir())) > 0):
        logger.info(f"Model {model} already downloaded, skipping")
        return

    model_download_path = Path(download_dir, model_fname)
    download_dir.mkdir(exist_ok=True)
    download_file(model_url, model_download_path)

    # Extract
    with tarfile.open(model_download_path, "r") as tar:
        paths = [Path(x) for x in tar.getnames()]
        root_paths = [str(x.parent) for x in paths if str(x.parent) != "."]
        if len(set(root_paths)) != 1:
            raise ValueError(
                f"The Clair3 model tarfile {model_fname} contains multiple root directories (there can only be one), please check the tar file."
            )

        tar.extractall(download_dir)

        if root_paths[0] != model:
            shutil.move(Path(download_dir, root_paths[0]), final_model_path)
    model_download_path.unlink()

    logger.info(f"Downloaded model: {model} to {final_model_path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=f"{os.getenv('CONDA_PREFIX')}/bin/models/",
        help="Directory to download the model to, default is: %(default)s",
    )
    parser.add_argument(
        "--model",
        type=str,
        default='',
        help="Download given specific model only instead of all",
    )
    args = parser.parse_args()

    if not os.getenv("CONDA_PREFIX"):
        logger.warning(
            "CONDA_PREFIX is not set, this probably means you are not running this inside a conda environment, "
            "if you have not provided a model path argument '--model-dir' the models might be downloaded "
            "somewhere you don't want them to be."
        )

    # Download specific model
    if args.model:
        if args.model in MODELS:
            download_model(
                args.model_dir,
                args.model
            )
        else:
            logger.error(f"Given model {args.model} does not exist or has yet to be added to the list of available models")
            sys.exit(1)

    # Download all
    else:
        for model in MODELS:
            download_model(
                args.model_dir,
                model
            )

if __name__ == "__main__":
    main()
