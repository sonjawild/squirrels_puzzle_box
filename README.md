# Code and Data accompanying 'Behavioral type predicts three stages of anthropogenic resource discovery and exploitation in California ground squirrels'

## R Code
Analysis, tables and supplementary figure can be fully reproduced with the available code. 

## Data
### latency.data.discovery: 
Data detailing latency to discovery of the resource in the two populations (more disturbed = Crow, less disturbed = Paradise)
- id: PIT tag identifier
- timestamp: Time and date of discovery (NA if never discovered)
- box: Box1-3 (refers to the three experimental puzzle boxes simultanously deployed)
- time: Discovery time in experimental hours (max time if never discovered)
- t_end: time at the end of the experiment
- other_present: 1 if another squirrel was present at the same box at time of discovery (or had departed within the 10s before), 0 i no conspecific present
- eigenvector: eigenvector centrlaity based on co-feeding networks at RFID feeders
- mobility: average number of antenna visited during social network mapping days (called 'mobility' in MS)
- n_days_net: # of network days during which the squirrel was registered (not used)
- log_num_trapped_per_day: the log of the average trapping numbers per trapping day (part of boldness calculation)
- any_beh_prop: proportion of trapping events during which squirrel showed fear responses (part of boldness calculation)
- age: A for adult, P for juvenile
- sex: M for male, F for female
- colony: population where squirrels was most often trapped
- prop_days_trapped: proportion of days a squirrel was trapped during trapping efforts (not used)
- bold: boldness score from PCA

### latency.data.solve: 
Data detailing latency to solve of the tax in the two populations (more disturbed = Crow, less disturbed = Paradise)
- id: PIT tag identifier
- timestamp: Time and date of third solve (NA if never discovered)
- time: Solving time in experimental hours (max time if never solved)
- t_end: time at the end of the experiment
- training.solves: # of solves during the training phase (when additional peanut butter was supplied on levers). Serves as a proxy for exposure during training and opportunities for tria-and-error learning.
- solve.rate: # of solves after the third solves divided by the number of arrivals to the box
- solve.rate.day: # solves after the third solve divided by the number of unique days discovered at the puzzle box
- total.solves: # of solves after the third solve until the end of the experiment
- total.arrivals: # of arrivals to the box after the third solve (used to calculate solving rate)
- n_observation_opportunities: # of observations of another squirrel performing a solve. Considered to have had opportunity if a solve occurred within 10s of arrival or while the focal was already present at the box
- eigenvector: eigenvector centrlaity based on co-feeding networks at RFID feeders
- mobility: average number of antenna visited during social network mapping days (called 'mobility' in MS)
- n_days_net: # of network days during which the squirrel was registered (not used)
- log_num_trapped_per_day: the log of the average trapping numbers per trapping day (part of boldness calculation)
- any_beh_prop: proportion of trapping events during which squirrel showed fear responses (part of boldness calculation)
- age: A for adult, P for juvenile
- sex: M for male, F for female
- colony: population where squirrels was most often trapped
- prop_days_trapped: proportion of days a squirrel was trapped during trapping efforts (not used)
- bold: boldness score from PCA

### Crow.net.data & Paradise.net.data
Co-feeding data from grids of RFID antennae used for mapping of social networks and space use.
- Date: mm/dd/yyyy
- Time: HH:MM:SS
- Antenna_ID: which of the three antenna the squirrel was registered on
- PIT: PIT tag identifier
- logger: which RFID logger
- date_time: yyyymmddHHMMSS
- date_logger: yyyymmdd_logger#_antenna#
- Network_location: L1-L6 (see Figure 1b in the main manuscript for map of study area)
- week: which week of network data collection (1-5)

### Crow.discovery.net & Paradise.discovery.net
R Data objects containing the gmm objects for computing social networks. These models are computationally intense, so we provide the output. These objects contain three slots:
- gmm: gmm object resulting from the gmm_events function (gaussian mixture model)
- net: the social network based on the simple ratio index
- eigen: eigenvector centralities for all individuals in the network (this value corresponds to the column 'eigen' in the data frames loaded above)



