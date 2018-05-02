# Coorga
## Description
Collection of tools to perform object-based analysis with special emphasis on metric of convection organization.

## Workflow
The typical workflow involves

* __time stacking__: slot files are packed into daily data stacks (also variable transformations and secondary variable calculations are included, here)
* __segmentation__: each daily data stack is segmented into categorized fields (involves to define a threshold) and also daily cluster files are saved
 * e.g. script `save_cluster_properties.py`
* __analysis stack__: at the next step, a analysis period (e.g. a month) is defined and selected properties from the cluster files are stacked into a larger cluster properties file
 * e.g. script `save_interesting_cluster_props.py`
* __average_statistics__: average statistics like number densities are calculated 
 * e.g. notebook `18-Average_NARVAL-SMF_NumberDensites.ipynb`
* __addition of fields__: at this point, several fields are added to the daily cluster files 
 * _pair-correlation function_: e.g. script `sort_pcf2cluster.py`
 * environmental variables