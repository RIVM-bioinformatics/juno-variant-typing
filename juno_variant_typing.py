"""
Juno variant typing
Authors: Boas van der Putten
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), IBR     
Date: 14-06-2023   
"""

from pathlib import Path
import pathlib
import yaml
import argparse
import sys
from dataclasses import dataclass, field
from juno_library import Pipeline
from typing import Optional, Union, Callable
from version import __package_name__, __version__, __description__


def main() -> None:
    juno_variant_typing = JunoVariantTyping()
    juno_variant_typing.run()


def check_number_within_range(
    minimum: float = 0, maximum: float = 1
) -> Union[Callable[[str], str], argparse.FileType]:
    """
    Creates a function to check whether a numeric value is within a range, inclusive.

    The generated function can be used by the `type` parameter in argparse.ArgumentParser.
    See https://stackoverflow.com/a/53446178.

    Args:
        value: the numeric value to check.
        minimum: minimum of allowed range, inclusive.
        maximum: maximum of allowed range, inclusive.

    Returns:
        A function which takes a single argument and checks this against the range.

    Raises:
        argparse.ArgumentTypeError: if the value is outside the range.
        ValueError: if the value cannot be converted to float.
    """

    def generated_func_check_range(value: str) -> str:
        value_f = float(value)
        if (value_f < minimum) or (value_f > maximum):
            raise argparse.ArgumentTypeError(
                f"Supplied value {value} is not within expected range {minimum} to {maximum}."
            )
        return str(value)

    return generated_func_check_range


@dataclass
class JunoVariantTyping(Pipeline):
    pipeline_name: str = __package_name__
    pipeline_version: str = __version__
    input_type: str = "bam_and_vcf"

    def _add_args_to_parser(self) -> None:
        super()._add_args_to_parser()

        self.parser.description = (
            "Juno variant typing pipeline for interpretation of genomic variants"
        )

        self.add_argument(
            "-m",
            "--metadata",
            type=Path,
            default=None,
            required=False,
            metavar="FILE",
            help="Relative or absolute path to the metadata csv file. If "
            "provided, it must contain at least one column named 'sample' "
            "with the name of the sample (same than file name but removing "
            "the suffix _R1.fastq.gz), a column called "
            "'genus' and a column called 'species'. The genus and species "
            "provided will be used to choose the serotyper and the MLST schema(s)."
            "If a metadata file is provided, it will overwrite the --species "
            "argument for the samples present in the metadata file.",
        )
        self.add_argument(
            "-s",
            "--species",
            type=lambda s: s.strip().lower(),
            nargs=2,
            default=["mycobacterium", "tuberculosis"],
            required=False,
            metavar=("GENUS", "SPECIES"),
            help="Species name (any species in the metadata file will overwrite"
            " this argument). It should be given as two words (e.g. --species "
            "Mycobacterium tuberculosis)",
        )
        self.add_argument(
            "-d",
            "--db_dir",
            type=Path,
            required=False,
            metavar="DIR",
            default="/mnt/db/juno/variant-typing",
            help="Relative or absolute path to the directory that contains the"
            " databases for all the tools used in this pipeline or where they"
            " should be downloaded. Default is: /mnt/db/juno/variant-typing",
        )
        self.add_argument(
            "--presets-path",
            type=Path,
            required=False,
            metavar="PATH",
            help="Relative or absolute path to custom presets.yaml to use. If"
            " none is provided, the default (config/presets.yaml) is used.",
        )

    def _parse_args(self) -> argparse.Namespace:
        args = super()._parse_args()

        # Optional arguments are loaded into self here
        self.db_dir: Path = args.db_dir.resolve()

        self.genus: Optional[str]
        self.species: Optional[str]
        self.genus, self.species = args.species
        self.metadata_file: Path = args.metadata
        self.presets_path: Optional[Path] = args.presets_path
        # self.single_copy_bed: Optional[Path] = args.single_copy_bed

        return args

    # # Extra class methods for this pipeline can be defined here
    # def example_class_method(self):
    #     print(f"example option is set to {self.example}")

    def setup(self) -> None:
        super().setup()
        self.update_sample_dict_with_metadata()
        self.set_presets()

        if self.snakemake_args["use_singularity"]:
            self.snakemake_args["singularity_args"] = " ".join(
                [
                    self.snakemake_args["singularity_args"],
                    f"--bind {self.db_dir}:{self.db_dir}",
                ]  # paths that singularity should be able to read from can be bound by adding to the above list
            )

        # # Extra class methods for this pipeline can be invoked here
        # if self.example:
        #     self.example_class_method()

        with open(
            Path(__file__).parent.joinpath("config/pipeline_parameters.yaml")
        ) as f:
            parameters_dict = yaml.safe_load(f)
        self.snakemake_config.update(parameters_dict)

        self.user_parameters = {
            "input_dir": str(self.input_dir),
            "output_dir": str(self.output_dir),
            "exclusion_file": str(self.exclusion_file),
            "custom_presets_file": str(self.presets_path),
            "use_singularity": str(self.snakemake_args["use_singularity"]),
            # "example": str(self.example), # other user parameters can be included in user_parameters.yaml here
        }

    def update_sample_dict_with_metadata(self) -> None:
        self.get_metadata_from_csv_file(
            filepath=self.metadata_file,
            expected_colnames=["sample", "genus", "species"],
        )
        # Add metadata
        for sample in self.sample_dict:
            if self.genus is not None and self.species is not None:
                self.sample_dict[sample]["genus"] = self.genus
                self.sample_dict[sample]["species"] = self.species
            else:
                try:
                    self.sample_dict[sample].update(self.juno_metadata[sample])
                except (KeyError, TypeError):
                    raise ValueError(
                        f"One of your samples is not in the metadata file "
                        f"({self.metadata_file}). Please ensure that all "
                        "samples are present in the metadata file or provide "
                        "a --species argument."
                    )
                self.sample_dict[sample]["genus"] = (
                    self.sample_dict[sample]["genus"].strip().lower()
                )
                self.sample_dict[sample]["species"] = (
                    self.sample_dict[sample]["species"].strip().lower()
                )

    def set_presets(self) -> None:
        if self.presets_path is None:
            self.presets_path = Path(__file__).parent.joinpath("config/presets.yaml")

        with open(self.presets_path) as f:
            presets_dict = yaml.safe_load(f)

        for sample in self.sample_dict:
            complete_species_name = "_".join(
                [self.sample_dict[sample]["genus"], self.sample_dict[sample]["species"]]
            )

            if complete_species_name in presets_dict.keys():
                for key, value in presets_dict[complete_species_name].items():
                    self.sample_dict[sample][key] = value


if __name__ == "__main__":
    main()
