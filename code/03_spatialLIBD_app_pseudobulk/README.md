
## Habenula_Visium

### Internal

JHPCE location: `/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium`

<br>

#### To deploy a pseudobulk `Habenula` spatialLIBD shiny app follow next steps:

<br>

1. Be sure to create a smaller spe (RDS, RData, etc) object for the shiny app: For Ex. spe_harmony_shiny.rds or spe_pseudobulk_shiny.rds

``run .../01_make_spe_shiny_app.R``

2. Create data (model) and prepare soft-links pointing to the specific model (ex. sce_pseudo_k09.rds), as well as other directories, Readme and documentation for the spatialLIBD app

``Run .../select_BayesSpaceK_brainModel.R``

3. Check and test the app before deploy

``Run .../app.R``

4. Deploy if ready

``Run .../deploy.R``

<br>




