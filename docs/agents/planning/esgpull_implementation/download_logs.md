 $ python scripts/download_example.py download-basic
[2026-03-03 16:18:24] INFO       [MainThread][climateset.utils] Yaml config file [/home/francispelletier/projects/ClimateSetExtension/configs/downloader/constants/cmip6.yaml] found.
[2026-03-03 16:18:24] INFO       [MainThread][climateset.utils] Loading YAML config file [/home/francispelletier/projects/ClimateSetExtension/configs/downloader/constants/cmip6.yaml].
[2026-03-03 16:18:24] INFO       [MainThread][climateset.utils] Yaml config file [/home/francispelletier/projects/ClimateSetExtension/configs/downloader/constants/cmip6plus.yaml] found.
[2026-03-03 16:18:24] INFO       [MainThread][climateset.utils] Loading YAML config file [/home/francispelletier/projects/ClimateSetExtension/configs/downloader/constants/cmip6plus.yaml].
[2026-03-03 16:18:24] INFO       [MainThread][climateset.utils] Yaml config file [/home/francispelletier/projects/ClimateSetExtension/configs/downloader/constants/imput4MIPs.yaml] found.
[2026-03-03 16:18:24] INFO       [MainThread][climateset.utils] Loading YAML config file [/home/francispelletier/projects/ClimateSetExtension/configs/downloader/constants/imput4MIPs.yaml].
[2026-03-03 16:18:24] INFO       [MainThread][climateset.utils] Yaml config file [/home/francispelletier/projects/ClimateSetExtension/configs/micro_dataset.yaml] found.
[2026-03-03 16:18:24] INFO       [MainThread][climateset.utils] Loading YAML config file [/home/francispelletier/projects/ClimateSetExtension/configs/micro_dataset.yaml].
[2026-03-03 16:18:24] INFO       [MainThread][climateset.download.downloader_config] Cleaned variables : ['CO2']
[2026-03-03 16:18:24] INFO       [MainThread][climateset.download.downloader_config] Emission variables to download: ['CO2_em_anthro', 'CO2_em_AIR_anthro']
[2026-03-03 16:18:24] INFO       [MainThread][climateset.download.downloader_config] Biomass burning vars to download: ['CO2']
[2026-03-03 16:18:24] INFO       [MainThread][climateset.download.downloader_config] Meta emission vars to download:
        []
        []
[2026-03-03 16:18:24] INFO       [MainThread][climateset.utils] Yaml config file [/home/francispelletier/projects/ClimateSetExtension/configs/micro_dataset.yaml] found.
[2026-03-03 16:18:24] INFO       [MainThread][climateset.utils] Loading YAML config file [/home/francispelletier/projects/ClimateSetExtension/configs/micro_dataset.yaml].
[2026-03-03 16:18:24] INFO       [MainThread][climateset.download.input4mips_downloader] Downloading data for variable: CO2_em_anthro
[2026-03-03 16:18:24] INFO       [MainThread][climateset.download.input4mips_downloader] Using download_raw_input_single_var() function
[2026-03-03 16:18:24] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf-node.llnl.gov/esg-search
[2026-03-03 16:18:25] WARNING    [MainThread][climateset.download.client] Error fetching facets from https://esgf-node.llnl.gov/esg-search: 422 Client Error: Unprocessable Content for url: https://esgf-node.ornl.gov/esgf-1-5-bridge?format=application%2Fsolr%2Bjson&limit=0&distrib=false&type=Dataset&project=input4MIPs&variable=CO2_em_anthro&institution_id=PNNL-JGCRI&facets=%2A
[2026-03-03 16:18:25] INFO       [MainThread][climateset.download.client] Rotating to next ESGF node...
[2026-03-03 16:18:25] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf.ceda.ac.uk/esg-search
[2026-03-03 16:18:28] INFO       [MainThread][climateset.download.utils] Available grid labels : ['gn']
[2026-03-03 16:18:28] INFO       [MainThread][climateset.download.utils] Choosing grid : gn
[2026-03-03 16:18:31] INFO       [MainThread][climateset.download.utils] Available nominal resolution : ['50 km']
[2026-03-03 16:18:31] INFO       [MainThread][climateset.download.utils] Choosing nominal resolution : 50 km
[2026-03-03 16:18:34] INFO       [MainThread][climateset.download.utils] Available frequencies : ['mon']
[2026-03-03 16:18:34] INFO       [MainThread][climateset.download.utils] Choosing default frequency : mon
[2026-03-03 16:18:37] INFO       [MainThread][climateset.download.utils] Available target mips: ['CMIP']
[2026-03-03 16:18:37] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf-node.llnl.gov/esg-search
[2026-03-03 16:18:38] WARNING    [MainThread][climateset.download.client] Error fetching facets from https://esgf-node.llnl.gov/esg-search: 422 Client Error: Unprocessable Content for url: https://esgf-node.ornl.gov/esgf-1-5-bridge?format=application%2Fsolr%2Bjson&limit=0&distrib=false&type=Dataset&project=input4MIPs&variable=CO2_em_anthro&institution_id=PNNL-JGCRI&grid_label=gn&nominal_resolution=50+km&frequency=mon&target_mip=CMIP&facets=%2A
[2026-03-03 16:18:38] INFO       [MainThread][climateset.download.client] Rotating to next ESGF node...
[2026-03-03 16:18:38] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf.ceda.ac.uk/esg-search
[2026-03-03 16:18:41] INFO       [MainThread][climateset.download.utils] Available versions : ['20250421', '20250325', '20241203']
[2026-03-03 16:18:41] INFO       [MainThread][climateset.download.utils] Choosing latest version: 20250421
[2026-03-03 16:18:45] INFO       [MainThread][climateset.download.utils] Result len for target CMIP: 1
********************************************************************************
*                                                                              *
* Note that new functionality to allow authentication without the need for     *
* certificates is available with this version of the wget script.  To enable,  *
* use the "-H" option and enter your OpenID and password when prompted:        *
*                                                                              *
* $ download -H [options...]                                     *
*                                                                              *
* For a full description of the available options use the help option:         *
*                                                                              *
* $ download -h                                                  *
*                                                                              *
********************************************************************************
Running download version: 1.3.2
Use download -h for help.

Script created for 6 file(s)
(The count won't match if you manually edit this file!)



CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-179912.nc ...Already downloaded and verified
CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_180001-184912.nc ...Already downloaded and verified
CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_185001-189912.nc ...Already downloaded and verified
CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_190001-194912.nc ...Already downloaded and verified
CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_195001-199912.nc ...Already downloaded and verified
CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_200001-202312.nc ...Already downloaded and verified
done
[2026-03-03 16:18:52] INFO       [MainThread][climateset.download.input4mips_downloader] Download results: [<pyesgf.search.results.ResultSet object at 0x7edb0f197d10>]
[2026-03-03 16:18:52] INFO       [MainThread][climateset.download.input4mips_downloader] Downloading data for variable: CO2_em_AIR_anthro
[2026-03-03 16:18:52] INFO       [MainThread][climateset.download.input4mips_downloader] Using download_raw_input_single_var() function
[2026-03-03 16:18:52] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf-node.llnl.gov/esg-search
[2026-03-03 16:18:53] WARNING    [MainThread][climateset.download.client] Error fetching facets from https://esgf-node.llnl.gov/esg-search: 422 Client Error: Unprocessable Content for url: https://esgf-node.ornl.gov/esgf-1-5-bridge?format=application%2Fsolr%2Bjson&limit=0&distrib=false&type=Dataset&project=input4MIPs&variable=CO2_em_AIR_anthro&institution_id=PNNL-JGCRI&facets=%2A
[2026-03-03 16:18:53] INFO       [MainThread][climateset.download.client] Rotating to next ESGF node...
[2026-03-03 16:18:53] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf.ceda.ac.uk/esg-search
[2026-03-03 16:18:56] INFO       [MainThread][climateset.download.utils] Available grid labels : ['gn']
[2026-03-03 16:18:56] INFO       [MainThread][climateset.download.utils] Choosing grid : gn
[2026-03-03 16:18:59] INFO       [MainThread][climateset.download.utils] Available nominal resolution : ['50 km']
[2026-03-03 16:18:59] INFO       [MainThread][climateset.download.utils] Choosing nominal resolution : 50 km
[2026-03-03 16:19:03] INFO       [MainThread][climateset.download.utils] Available frequencies : ['mon']
[2026-03-03 16:19:03] INFO       [MainThread][climateset.download.utils] Choosing default frequency : mon
[2026-03-03 16:19:06] INFO       [MainThread][climateset.download.utils] Available target mips: ['CMIP']
[2026-03-03 16:19:06] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf-node.llnl.gov/esg-search
[2026-03-03 16:19:06] WARNING    [MainThread][climateset.download.client] Error fetching facets from https://esgf-node.llnl.gov/esg-search: 422 Client Error: Unprocessable Content for url: https://esgf-node.ornl.gov/esgf-1-5-bridge?format=application%2Fsolr%2Bjson&limit=0&distrib=false&type=Dataset&project=input4MIPs&variable=CO2_em_AIR_anthro&institution_id=PNNL-JGCRI&grid_label=gn&nominal_resolution=50+km&frequency=mon&target_mip=CMIP&facets=%2A
[2026-03-03 16:19:06] INFO       [MainThread][climateset.download.client] Rotating to next ESGF node...
[2026-03-03 16:19:06] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf.ceda.ac.uk/esg-search
[2026-03-03 16:19:09] INFO       [MainThread][climateset.download.utils] Available versions : ['20250421', '20250325', '20241109']
[2026-03-03 16:19:09] INFO       [MainThread][climateset.download.utils] Choosing latest version: 20250421
[2026-03-03 16:19:13] INFO       [MainThread][climateset.download.utils] Result len for target CMIP: 1
********************************************************************************
*                                                                              *
* Note that new functionality to allow authentication without the need for     *
* certificates is available with this version of the wget script.  To enable,  *
* use the "-H" option and enter your OpenID and password when prompted:        *
*                                                                              *
* $ download -H [options...]                                     *
*                                                                              *
* For a full description of the available options use the help option:         *
*                                                                              *
* $ download -h                                                  *
*                                                                              *
********************************************************************************
Running download version: 1.3.2
Use download -h for help.

Script created for 6 file(s)
(The count won't match if you manually edit this file!)



CO2-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_175001-179912.nc ...Already downloaded and verified
CO2-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_180001-184912.nc ...Already downloaded and verified
CO2-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_185001-189912.nc ...Already downloaded and verified
CO2-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_190001-194912.nc ...Already downloaded and verified
CO2-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_195001-199912.nc ...Already downloaded and verified
CO2-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_200001-202312.nc ...Already downloaded and verified
done
[2026-03-03 16:19:20] INFO       [MainThread][climateset.download.input4mips_downloader] Download results: [<pyesgf.search.results.ResultSet object at 0x7edae27abc10>]
[2026-03-03 16:19:20] INFO       [MainThread][climateset.download.input4mips_downloader] Downloading biomassburing data for variable: CO2
[2026-03-03 16:19:20] INFO       [MainThread][climateset.download.input4mips_downloader] Using download_raw_input_single_var() function
[2026-03-03 16:19:20] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf-node.llnl.gov/esg-search
[2026-03-03 16:19:20] WARNING    [MainThread][climateset.download.client] Error fetching facets from https://esgf-node.llnl.gov/esg-search: 422 Client Error: Unprocessable Content for url: https://esgf-node.ornl.gov/esgf-1-5-bridge?format=application%2Fsolr%2Bjson&limit=0&distrib=false&type=Dataset&project=input4MIPs&variable=CO2&institution_id=VUA&facets=%2A
[2026-03-03 16:19:20] INFO       [MainThread][climateset.download.client] Rotating to next ESGF node...
[2026-03-03 16:19:20] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf.ceda.ac.uk/esg-search
[2026-03-03 16:19:23] INFO       [MainThread][climateset.download.utils] Available grid labels : ['gn']
[2026-03-03 16:19:23] INFO       [MainThread][climateset.download.utils] Choosing grid : gn
[2026-03-03 16:19:27] INFO       [MainThread][climateset.download.utils] Available nominal resolution : ['25 km']
[2026-03-03 16:19:27] INFO       [MainThread][climateset.download.utils] Choosing nominal resolution : 25 km
[2026-03-03 16:19:30] INFO       [MainThread][climateset.download.utils] Available frequencies : ['mon']
[2026-03-03 16:19:30] INFO       [MainThread][climateset.download.utils] Choosing default frequency : mon
[2026-03-03 16:19:33] INFO       [MainThread][climateset.download.utils] Available target mips: ['CMIP']
[2026-03-03 16:19:33] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf-node.llnl.gov/esg-search
[2026-03-03 16:19:33] WARNING    [MainThread][climateset.download.client] Error fetching facets from https://esgf-node.llnl.gov/esg-search: 422 Client Error: Unprocessable Content for url: https://esgf-node.ornl.gov/esgf-1-5-bridge?format=application%2Fsolr%2Bjson&limit=0&distrib=false&type=Dataset&project=input4MIPs&variable=CO2&institution_id=VUA&grid_label=gn&nominal_resolution=25+km&frequency=mon&target_mip=CMIP&facets=%2A
[2026-03-03 16:19:33] INFO       [MainThread][climateset.download.client] Rotating to next ESGF node...
[2026-03-03 16:19:33] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf.ceda.ac.uk/esg-search
[2026-03-03 16:19:36] INFO       [MainThread][climateset.download.utils] Available versions : ['20161002', '20160705']
[2026-03-03 16:19:36] INFO       [MainThread][climateset.download.utils] Choosing latest version: 20161002
[2026-03-03 16:19:40] INFO       [MainThread][climateset.download.utils] Result len for target CMIP: 1
********************************************************************************
*                                                                              *
* Note that new functionality to allow authentication without the need for     *
* certificates is available with this version of the wget script.  To enable,  *
* use the "-H" option and enter your OpenID and password when prompted:        *
*                                                                              *
* $ download -H [options...]                                     *
*                                                                              *
* For a full description of the available options use the help option:         *
*                                                                              *
* $ download -h                                                  *
*                                                                              *
********************************************************************************
Running download version: 1.3.2
Use download -h for help.

Script created for 2 file(s)
(The count won't match if you manually edit this file!)



CO2-em-biomassburning_input4MIPs_emissions_CMIP_VUA-CMIP-BB4CMIP6-1-1_gn_175001-184912.nc ...Already downloaded and verified
CO2-em-biomassburning_input4MIPs_emissions_CMIP_VUA-CMIP-BB4CMIP6-1-1_gn_185001-201512.nc ...Already downloaded and verified
done
[2026-03-03 16:19:47] INFO       [MainThread][climateset.download.input4mips_downloader] Download results: [<pyesgf.search.results.ResultSet object at 0x7edb0c944490>]
[2026-03-03 16:19:47] INFO       [MainThread][climateset.download.cmip6_downloader] Downloading data for model: [NorESM2-LM]
[2026-03-03 16:19:47] INFO       [MainThread][climateset.download.cmip6_downloader] Downloading data for variable: [tas]
[2026-03-03 16:19:47] INFO       [MainThread][climateset.download.cmip6_downloader] Downloading data for experiment: [historical]
[2026-03-03 16:19:47] INFO       [MainThread][climateset.download.utils] Using download_from_model_single_var() function
[2026-03-03 16:19:47] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf-node.llnl.gov/esg-search

-------------------------------------------------------------------------------
Warning - defaulting to search with facets=*

This behavior is kept for backward-compatibility, but ESGF indexes might not
successfully perform a distributed search when this option is used, so some
results may be missing.  For full results, it is recommended to pass a list of
facets of interest when instantiating a context object.  For example,

      ctx = conn.new_context(facets='project,experiment_id')

Only the facets that you specify will be present in the facets_counts dictionary.

This warning is displayed when a distributed search is performed while using the
facets=* default, a maximum of once per context object.  To suppress this warning,
set the environment variable ESGF_PYCLIENT_NO_FACETS_STAR_WARNING to any value
or explicitly use  conn.new_context(facets='*')

-------------------------------------------------------------------------------
[2026-03-03 16:19:47] WARNING    [MainThread][climateset.download.client] Error fetching facets from https://esgf-node.llnl.gov/esg-search: 422 Client Error: Unprocessable Content for url: https://esgf-node.ornl.gov/esgf-1-5-bridge?format=application%2Fsolr%2Bjson&limit=0&distrib=true&type=Dataset&project=CMIP6&variable=tas&experiment_id=historical&source_id=NorESM2-LM&facets=%2A
[2026-03-03 16:19:47] INFO       [MainThread][climateset.download.client] Rotating to next ESGF node...
[2026-03-03 16:19:47] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf.ceda.ac.uk/esg-search

-------------------------------------------------------------------------------
Warning - defaulting to search with facets=*

This behavior is kept for backward-compatibility, but ESGF indexes might not
successfully perform a distributed search when this option is used, so some
results may be missing.  For full results, it is recommended to pass a list of
facets of interest when instantiating a context object.  For example,

      ctx = conn.new_context(facets='project,experiment_id')

Only the facets that you specify will be present in the facets_counts dictionary.

This warning is displayed when a distributed search is performed while using the
facets=* default, a maximum of once per context object.  To suppress this warning,
set the environment variable ESGF_PYCLIENT_NO_FACETS_STAR_WARNING to any value
or explicitly use  conn.new_context(facets='*')

-------------------------------------------------------------------------------
[2026-03-03 16:19:56] INFO       [MainThread][climateset.download.utils] Available frequencies : ['mon', 'day', '6hrPt', '6hr']
[2026-03-03 16:19:56] INFO       [MainThread][climateset.download.utils] Choosing default frequency : mon
[2026-03-03 16:20:04] INFO       [MainThread][climateset.download.utils] Available grid labels : ['gn']
[2026-03-03 16:20:04] INFO       [MainThread][climateset.download.utils] Choosing grid : gn
[2026-03-03 16:20:12] INFO       [MainThread][climateset.download.utils] Available variants : ['r9i1p1f1', 'r8i1p1f1', 'r7i1p1f1', 'r6i1p1f1', 'r5i1p1f1', 'r4i1p1f1', 'r43i1p1f1', 'r42i1p1f1', 'r41i1p1f1', 'r40i1p1f1', 'r3i1p1f1', 'r39i1p1f1', 'r38i1p1f1', 'r37i1p1f1', 'r36i1p1f1', 'r35i1p1f1', 'r34i1p1f1', 'r33i1p1f1', 'r32i1p1f1', 'r31i1p1f1', 'r30i1p1f1', 'r2i1p1f1', 'r29i1p1f1', 'r28i1p1f1', 'r27i1p1f1', 'r26i1p1f1', 'r25i1p1f1', 'r24i1p1f1', 'r23i1p1f1', 'r22i1p1f1', 'r21i1p1f1', 'r20i1p1f1', 'r1i1p4f1', 'r1i1p1f1', 'r19i1p1f1', 'r18i1p1f1', 'r17i1p1f1', 'r16i1p1f1', 'r15i1p1f1', 'r14i1p1f1', 'r13i1p1f1', 'r12i1p1f1', 'r11i1p1f1', 'r10i1p1f1']

[2026-03-03 16:20:12] INFO       [MainThread][climateset.download.utils] Length : 44
[2026-03-03 16:20:12] INFO       [MainThread][climateset.download.utils] Desired list of ensemble members given: ['r2i1p1f1']
[2026-03-03 16:20:12] INFO       [MainThread][climateset.download.utils] Ensembles member: r2i1p1f1
[2026-03-03 16:20:12] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf-node.llnl.gov/esg-search

-------------------------------------------------------------------------------
Warning - defaulting to search with facets=*

This behavior is kept for backward-compatibility, but ESGF indexes might not
successfully perform a distributed search when this option is used, so some
results may be missing.  For full results, it is recommended to pass a list of
facets of interest when instantiating a context object.  For example,

      ctx = conn.new_context(facets='project,experiment_id')

Only the facets that you specify will be present in the facets_counts dictionary.

This warning is displayed when a distributed search is performed while using the
facets=* default, a maximum of once per context object.  To suppress this warning,
set the environment variable ESGF_PYCLIENT_NO_FACETS_STAR_WARNING to any value
or explicitly use  conn.new_context(facets='*')

-------------------------------------------------------------------------------
[2026-03-03 16:20:13] WARNING    [MainThread][climateset.download.client] Error fetching facets from https://esgf-node.llnl.gov/esg-search: 422 Client Error: Unprocessable Content for url: https://esgf-node.ornl.gov/esgf-1-5-bridge?format=application%2Fsolr%2Bjson&limit=0&distrib=true&type=Dataset&project=CMIP6&variable=tas&experiment_id=historical&source_id=NorESM2-LM&frequency=mon&grid_label=gn&variant_label=r2i1p1f1&facets=%2A
[2026-03-03 16:20:13] INFO       [MainThread][climateset.download.client] Rotating to next ESGF node...
[2026-03-03 16:20:13] INFO       [MainThread][climateset.download.client] Connecting to ESGF node: https://esgf.ceda.ac.uk/esg-search

-------------------------------------------------------------------------------
Warning - defaulting to search with facets=*

This behavior is kept for backward-compatibility, but ESGF indexes might not
successfully perform a distributed search when this option is used, so some
results may be missing.  For full results, it is recommended to pass a list of
facets of interest when instantiating a context object.  For example,

      ctx = conn.new_context(facets='project,experiment_id')

Only the facets that you specify will be present in the facets_counts dictionary.

This warning is displayed when a distributed search is performed while using the
facets=* default, a maximum of once per context object.  To suppress this warning,
set the environment variable ESGF_PYCLIENT_NO_FACETS_STAR_WARNING to any value
or explicitly use  conn.new_context(facets='*')

-------------------------------------------------------------------------------
[2026-03-03 16:20:21] INFO       [MainThread][climateset.download.utils] Available versions : ['20190920']
[2026-03-03 16:20:21] INFO       [MainThread][climateset.download.utils] Choosing latest version: 20190920
[2026-03-03 16:20:31] INFO       [MainThread][climateset.download.utils] Result len 3
[2026-03-03 16:20:31] INFO       [MainThread][climateset.download.utils] [<pyesgf.search.results.ResultSet object at 0x7edaf66cffd0>]

-------------------------------------------------------------------------------
Warning - defaulting to search with facets=*

This behavior is kept for backward-compatibility, but ESGF indexes might not
successfully perform a distributed search when this option is used, so some
results may be missing.  For full results, it is recommended to pass a list of
facets of interest when instantiating a context object.  For example,

      ctx = conn.new_context(facets='project,experiment_id')

Only the facets that you specify will be present in the facets_counts dictionary.

This warning is displayed when a distributed search is performed while using the
facets=* default, a maximum of once per context object.  To suppress this warning,
set the environment variable ESGF_PYCLIENT_NO_FACETS_STAR_WARNING to any value
or explicitly use  conn.new_context(facets='*')

-------------------------------------------------------------------------------
********************************************************************************
*                                                                              *
* Note that new functionality to allow authentication without the need for     *
* certificates is available with this version of the wget script.  To enable,  *
* use the "-H" option and enter your OpenID and password when prompted:        *
*                                                                              *
* $ download -H [options...]                                     *
*                                                                              *
* For a full description of the available options use the help option:         *
*                                                                              *
* $ download -h                                                  *
*                                                                              *
********************************************************************************
Running download version: 1.3.2
Use download -h for help.

Script created for 17 file(s)
(The count won't match if you manually edit this file!)



tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_185001-185912.nc ...Downloading
--2026-03-03 16:20:38--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_185001-185912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3921956 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_185001-185912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_185001-185912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M   925KB/s    in 4,1s    

2026-03-03 16:20:43 (925 KB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_185001-185912.nc’ saved [3921956/3921956]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_186001-186912.nc ...Downloading
--2026-03-03 16:20:43--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_186001-186912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3921991 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_186001-186912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_186001-186912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M  1,04MB/s    in 3,6s    

2026-03-03 16:20:47 (1,04 MB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_186001-186912.nc’ saved [3921991/3921991]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_187001-187912.nc ...Downloading
--2026-03-03 16:20:47--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_187001-187912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3922810 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_187001-187912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_187001-187912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M   956KB/s    in 4,1s    

2026-03-03 16:20:52 (926 KB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_187001-187912.nc’ saved [3922810/3922810]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_188001-188912.nc ...Downloading
--2026-03-03 16:20:52--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_188001-188912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3920723 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_188001-188912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_188001-188912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M   953KB/s    in 4,0s    

2026-03-03 16:20:56 (953 KB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_188001-188912.nc’ saved [3920723/3920723]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_189001-189912.nc ...Downloading
--2026-03-03 16:20:56--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_189001-189912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3920346 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_189001-189912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_189001-189912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M   958KB/s    in 4,2s    

2026-03-03 16:21:01 (905 KB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_189001-189912.nc’ saved [3920346/3920346]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_190001-190912.nc ...Downloading
--2026-03-03 16:21:01--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_190001-190912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3922637 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_190001-190912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_190001-190912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M  1019KB/s    in 3,8s    

2026-03-03 16:21:06 (1019 KB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_190001-190912.nc’ saved [3922637/3922637]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_191001-191912.nc ...Downloading
--2026-03-03 16:21:06--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_191001-191912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3921878 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_191001-191912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_191001-191912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M   966KB/s    in 4,1s    

2026-03-03 16:21:10 (933 KB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_191001-191912.nc’ saved [3921878/3921878]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_192001-192912.nc ...Downloading
--2026-03-03 16:21:10--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_192001-192912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3918591 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_192001-192912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_192001-192912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M   968KB/s    in 4,1s    

2026-03-03 16:21:15 (928 KB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_192001-192912.nc’ saved [3918591/3918591]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_193001-193912.nc ...Downloading
--2026-03-03 16:21:15--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_193001-193912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3920337 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_193001-193912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_193001-193912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M  1,02MB/s    in 3,7s    

2026-03-03 16:21:19 (1,02 MB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_193001-193912.nc’ saved [3920337/3920337]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_194001-194912.nc ...Downloading
--2026-03-03 16:21:19--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_194001-194912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3921350 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_194001-194912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_194001-194912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M  1011KB/s    in 3,8s    

2026-03-03 16:21:24 (1011 KB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_194001-194912.nc’ saved [3921350/3921350]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_195001-195912.nc ...Downloading
--2026-03-03 16:21:24--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_195001-195912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3920415 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_195001-195912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_195001-195912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M  2,85MB/s    in 1,3s    

2026-03-03 16:21:26 (2,85 MB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_195001-195912.nc’ saved [3920415/3920415]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_196001-196912.nc ...Downloading
--2026-03-03 16:21:26--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_196001-196912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3920774 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_196001-196912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_196001-196912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M  1,06MB/s    in 3,5s    

2026-03-03 16:21:30 (1,06 MB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_196001-196912.nc’ saved [3920774/3920774]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_197001-197912.nc ...Downloading
--2026-03-03 16:21:30--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_197001-197912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3920052 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_197001-197912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_197001-197912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M  1,03MB/s    in 3,6s    

2026-03-03 16:21:34 (1,03 MB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_197001-197912.nc’ saved [3920052/3920052]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_198001-198912.nc ...Downloading
--2026-03-03 16:21:34--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_198001-198912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3919588 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_198001-198912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_198001-198912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M   969KB/s    in 3,9s    

2026-03-03 16:21:39 (969 KB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_198001-198912.nc’ saved [3919588/3919588]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_199001-199912.nc ...Downloading
--2026-03-03 16:21:39--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_199001-199912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3919135 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_199001-199912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_199001-199912.nc                            100%[=========================================================================================================================================================================================================================>]   3,74M   821KB/s    in 5,1s    

2026-03-03 16:21:44 (748 KB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_199001-199912.nc’ saved [3919135/3919135]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_200001-200912.nc ...Downloading
--2026-03-03 16:21:45--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_200001-200912.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3915211 (3,7M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_200001-200912.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_200001-200912.nc                            100%[=========================================================================================================================================================================================================================>]   3,73M  3,63MB/s    in 1,0s    

2026-03-03 16:21:46 (3,63 MB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_200001-200912.nc’ saved [3915211/3915211]

  sha256 ok. done!
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_201001-201412.nc ...Downloading
--2026-03-03 16:21:46--  https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/NCC/NorESM2-LM/historical/r2i1p1f1/Amon/tas/gn/v20190920/tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_201001-201412.nc
Resolving esgf.ceda.ac.uk (esgf.ceda.ac.uk)... 130.246.128.97
Connecting to esgf.ceda.ac.uk (esgf.ceda.ac.uk)|130.246.128.97|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 1978507 (1,9M) [application/octet-stream]
Saving to: ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_201001-201412.nc’

tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_201001-201412.nc                            100%[=========================================================================================================================================================================================================================>]   1,89M  2,38MB/s    in 0,8s    

2026-03-03 16:21:47 (2,38 MB/s) - ‘tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_201001-201412.nc’ saved [1978507/1978507]

  sha256 ok. done!
done

-------------------------------------------------------------------------------
Warning - defaulting to search with facets=*

This behavior is kept for backward-compatibility, but ESGF indexes might not
successfully perform a distributed search when this option is used, so some
results may be missing.  For full results, it is recommended to pass a list of
facets of interest when instantiating a context object.  For example,

      ctx = conn.new_context(facets='project,experiment_id')

Only the facets that you specify will be present in the facets_counts dictionary.

This warning is displayed when a distributed search is performed while using the
facets=* default, a maximum of once per context object.  To suppress this warning,
set the environment variable ESGF_PYCLIENT_NO_FACETS_STAR_WARNING to any value
or explicitly use  conn.new_context(facets='*')

-------------------------------------------------------------------------------
********************************************************************************
*                                                                              *
* Note that new functionality to allow authentication without the need for     *
* certificates is available with this version of the wget script.  To enable,  *
* use the "-H" option and enter your OpenID and password when prompted:        *
*                                                                              *
* $ download -H [options...]                                     *
*                                                                              *
* For a full description of the available options use the help option:         *
*                                                                              *
* $ download -h                                                  *
*                                                                              *
********************************************************************************
Running download version: 1.3.2
Use download -h for help.

Script created for 17 file(s)
(The count won't match if you manually edit this file!)



tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_185001-185912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_186001-186912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_187001-187912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_188001-188912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_189001-189912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_190001-190912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_191001-191912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_192001-192912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_193001-193912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_194001-194912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_195001-195912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_196001-196912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_197001-197912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_198001-198912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_199001-199912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_200001-200912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_201001-201412.nc ...Already downloaded and verified
done

-------------------------------------------------------------------------------
Warning - defaulting to search with facets=*

This behavior is kept for backward-compatibility, but ESGF indexes might not
successfully perform a distributed search when this option is used, so some
results may be missing.  For full results, it is recommended to pass a list of
facets of interest when instantiating a context object.  For example,

      ctx = conn.new_context(facets='project,experiment_id')

Only the facets that you specify will be present in the facets_counts dictionary.

This warning is displayed when a distributed search is performed while using the
facets=* default, a maximum of once per context object.  To suppress this warning,
set the environment variable ESGF_PYCLIENT_NO_FACETS_STAR_WARNING to any value
or explicitly use  conn.new_context(facets='*')

-------------------------------------------------------------------------------
********************************************************************************
*                                                                              *
* Note that new functionality to allow authentication without the need for     *
* certificates is available with this version of the wget script.  To enable,  *
* use the "-H" option and enter your OpenID and password when prompted:        *
*                                                                              *
* $ download -H [options...]                                     *
*                                                                              *
* For a full description of the available options use the help option:         *
*                                                                              *
* $ download -h                                                  *
*                                                                              *
********************************************************************************
Running download version: 1.3.2
Use download -h for help.

Script created for 17 file(s)
(The count won't match if you manually edit this file!)



tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_185001-185912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_186001-186912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_187001-187912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_188001-188912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_189001-189912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_190001-190912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_191001-191912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_192001-192912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_193001-193912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_194001-194912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_195001-195912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_196001-196912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_197001-197912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_198001-198912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_199001-199912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_200001-200912.nc ...Already downloaded and verified
tas_Amon_NorESM2-LM_historical_r2i1p1f1_gn_201001-201412.nc ...Already downloaded and verified
done