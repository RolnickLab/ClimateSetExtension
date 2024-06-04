# Configuration

This folder is for configurations, like:

- Application configurations
- Training parameters
- Script parameters
- etc.

The following configs need to be adapted by the user:
- main_config.yaml
- downloader 
- cmip6 processor
- input4mips processor 

If you are adding a new variables, please edit:
- units 

If you are using a new output resolution (i.e. you are not using NorESM as "role model" for the output resolution), please edit:
- resolutions
- output_rolemodel 

## Output Role Model
This is the only config directory that has not a yaml file, but a netcdf file. Please store here the role model of your desired output data. The default is NorESM, i.e. all files will be aligned to have the same resolution format like NorESM (in terms of where longitude / latitude points start, step size, etc.).


