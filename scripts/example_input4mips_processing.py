from pathlib import Path

from climateset.processing.abstract_processor_step import process_steps
from climateset.processing.raw.input4mips.input4mips_processor import (
    CorrectUnitsStep,
    CreateAnthroFireFilesStep,
    Input4MipsEmissionProcessor,
    RenameBiomassBurningFilesStep,
    ReorderSSPCO2DimensionsStep,
)

#
# Simple example with Input4MipsEmissionProcessor
#

# Stand in for arguments
processor_arguments = {"input_directory": "<SOME_PATH_HERE>"}

# Create processor for a specific directory
processor = Input4MipsEmissionProcessor(**processor_arguments)

# Add processing steps - These will be ordered by Input4MipsEmissionProcessor
processor.add_correct_calendar_step()
processor.add_reorder_ssp_co2_dimensions_step()
processor.add_correct_names_step()
processor.add_correct_units_step()
processor.add_create_fire_files()

# Execute steps
processor.process_directory()

#
# More manual example with manual setting of steps
#

# Stand ins for arguments
input_directory = Path("SOME_PATH_HERE/")
anthrofire_arguments = {"metadata_directory": "<SOME_PATH_HERE>", "cdo_operators": []}
correct_units_arguments = {"desired_units": "<DESIRED_UNITS>"}

# Create steps
step_1 = RenameBiomassBurningFilesStep()
step_2 = ReorderSSPCO2DimensionsStep()
step_3 = CorrectUnitsStep(**correct_units_arguments)
step_4 = CreateAnthroFireFilesStep(**anthrofire_arguments)

process_list = [step_1, step_2, step_3, step_4]

# Execute steps
output_path = process_steps(input_directory=input_directory, list_of_steps=process_list)
