# Arimanator Integration
- Final file: R:\regulatory\Stress Testing 2019\MRM\Model Documents\template
- Typical path:
	load ST2019_.R in R:\Projects\ArimanatorRunner\Stress Test 2019\CHIMPS ->
	parse model data ->
 	runArimanator.R (arimanatorEngage()) to run model, export RDS ->
	convertArimanatorRDS.R (Not used?) ->
	modelDoctext.R (Not used?) ->
	renderModelDoc.R (calls ModelDoc.Rmd)

- modelDoc.Rmd
-- rds$id
-- rds$
- Installed packages: XLConnectJars, rJava, XLConnect