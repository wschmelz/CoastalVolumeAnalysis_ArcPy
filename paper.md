title: "A Geospatial Workflow for Quantifying Coastal Geomorphological Change"
tags:
  - coastal geomorphology
  - GIS
  - 3D volumetric analysis
  - shoreline monitoring
  - ArcGIS workflow
authors:
  - name: "FirstName LastName"
    affiliation: 1
  - name: "SecondName LastName"
    affiliation: 2
affiliations:
  - index: 1
    name: "Affiliation 1, Some Institution, Country"
  - index: 2
    name: "Affiliation 2, Another Institution, Country"
date: 9 December 2024
bibliography: paper.bib
---

## Summary

Shorelines and coastal landforms are dynamic systems subject to continuous changes driven by natural processes such as storms, sea-level rise, and sediment supply, as well as by human interventions (Carter, 1988; Psuty & Ofiara, 2002; Nicholls et al., 2007). Understanding the spatial and temporal patterns of coastal geomorphological change is critical for effective coastal management and conservation. Accurate quantification of erosion, deposition, and net volumetric change provides coastal scientists and resource managers with insights into sediment budgets, landform evolution, and habitat dynamics.

We present a standardized geospatial workflow that leverages high-density coastal topographical surveys, geodetic control, and robust geospatial analysis tools to derive three-dimensional volumetric changes in coastal landscapes. Implemented primarily using ESRI ArcGIS software and survey-grade GNSS data, this methodology transforms raw topographic survey points into raster-based Digital Elevation Models (DEMs), enabling the quantification of volumetric changes through space and time. The workflow includes: (1) generation of DEMs from survey data using Delaunay triangulation, (2) creation of a standardized base raster for spatial consistency, (3) spatial segmentation of the alongshore domain into compartments for localized analyses, (4) inter-comparison of DEMs from multiple survey dates to calculate volumetric erosion, deposition, and net change, and (5) reporting and visualization of results via graphical and statistical summaries.

This framework adheres to scientific principles of repeatability and reproducibility, facilitating integration with existing coastal monitoring protocols and datasets. It can be readily applied to multiple sites and temporal scales, supporting short-term event-based assessments (e.g., storm impacts) as well as long-term trend analyses (Roman & Barrett, 1999; Psuty & Silveira, 2010). The methodology provides coastal managers with essential metrics that inform resource protection strategies, infrastructure planning, and ecological restoration efforts.

## Statement of need

As coastal parks and protected areas face increasing pressures from climate change, sea-level rise, and storm events, park managers and scientists require reliable, spatially explicit, and temporally resolved information to guide decision-making (National Research Council, 1995; Psuty & Pace, 2009; Rodriguez, 2004). Traditional coastal monitoring methods that rely solely on shoreline-position metrics or sparse topographic data are often inadequate to capture the complexities of sediment movement and morphological change.

The proposed workflow addresses these needs by integrating:  
- State-of-the-art GNSS surveying to ensure high-accuracy vertical and horizontal positioning,  
- ArcGIS-based tools for efficient generation of DEMs and spatial analyses, and  
- A standardized process for segmenting the coast into compartments, enabling spatially focused analyses of sediment budgets and geomorphological evolution.

This approach is particularly valuable for national parks, coastal reserves, and other coastal management units seeking to establish long-term monitoring protocols. It complements and enhances existing datasets, including historical aerial imagery, LiDAR, and other remote sensing products. The open, step-by-step nature of the methodology ensures flexibility, allowing practitioners to incorporate additional datasets, adapt to different software packages, and refine the parameters as needed.

## Acknowledgements

Development of this workflow was informed by ongoing coastal geomorphological research programs and guidance from local and federal resource management agencies. We acknowledge the contributions of field teams collecting survey data, GIS specialists who tested the procedures, and the support of stakeholders who provided feedback on data interpretation and reporting formats.

## References