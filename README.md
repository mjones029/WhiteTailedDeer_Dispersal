# WhiteTailedDeer_Dispersal
Public repository for code identifying and analyzing dispersal events from white tailed deer movement data


This code corresponds to analyses from Gilbertson et al 2022, "Agricultural land use shapes dispersal in white-tailed deer (*Odocoileus virginianus*)," with is available from the journal [*Movement Ecology*](https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-022-00342-5).

This code is also archived with Zenodo: [![DOI](https://zenodo.org/badge/551513859.svg)](https://zenodo.org/badge/latestdoi/551513859)


## Dispersal detection
Movement data used in Gilbertson et al 2022 is not publicly available but interested parties can contact the Wisconsin Department of Natural Resources for more information.

To demonstrate the functionality of the dispersal detection algorithm used in our publication, we've included code for simulating simple movement paths. Simulated movement data and metadata is included in this repository, as well as the code to simulate this data. 

### Workflow

1. Simulate simple movement data with the script **Simulate_data.R**
2. Detect dispersal events among trajectories using the script **HR_dispersal.R**. This script also estimates the timing of dispersal events. 
3. Calculate "traversability" metrics using the script **Traversability.R**.
4. Calculate Proximity scores using the script **ProxRate.R**
5. Calculate proportion of different habitats in natal ranges using the script **InRange_cover.R**

*Additional code associated with this publication is forthcoming.*

