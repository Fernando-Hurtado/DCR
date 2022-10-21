# Index:

0. General information
1. Dataset_Indicators' description
2. Dataset_Traits' description
3. Script description

# 0 - General information

The datasets and code were compiled by Fernando Hurtado Bocanegra in colaboration with Belén Estébanez, Pedro Aragón, Joaquín Hortal, Manuel Molina-Bustamante and Nagore G. Medina.

Both datasets and code were use in the article: "*Moss establishment success is determined by the interaction between propagule size and species identity*".

Link of the article: https://assets.researchsquare.com/files/rs-1358759/v1/6d414dbb-521d-4b6e-bbb0-c29072bb5461.pdf?c=1647365290.

The datasets and code are published under Creative Commons BY-SA 4.0 International license. 
When using the dataset, please cite the main manuscript.

If you need further assistance with the datasets or code do not hesitate to contact me at: ferhurboc@gmail.com

# 1 - Dataset_Indicators' description

Data gathered from a factorial design culture experiment of 3 fragments classes of 6 moss species with 12 replicates for each one. The experiment tests the establishment success of the mosses and use some indicators based on the growth in terms of biomass and surface colonized and characterize the propagules of the mosses to better understand their establishment.

The dataset contains 217 rows and 15 columns of data. 
The first row of the dataset contains the name of the columns, and the rest correspond to the 216 samples. 
First 4 columns are dedicated to the identification of the sample.
The rest of columns correspond to measurements or obtained through calculations with the collected data.

Column names and descriptions are:
 - (1) Block -> The number of the block in the chamber experiment corresponding to the sample.
 - (2) Block_Row -> The number of the row within the block corresponding to the sample.
 - (3) Block_Column -> The number of the column within the block corresponding to the sample.
 - (4) Sample -> The name of the species and the fragment size corresponding to the sample.
 - (5) Shoots_established -> The number of shoots established (result obtained after finishing the culture period).
 - (6) S_shoots -> The surface occupied by the shoots (in mm²).
 - (7) S_total -> The surface occupied by the total biomass in the sample (in mm²). Includes viable and non-viable biomass.
 - (8) S_substratum -> The surface occupied by the substratum (in mm²).
 - (9) Sp. -> The species name, one of the following: D. scoparium, H. aureum, H. cupressiforme , P. capillare, S. ruralis or T. squarrosa.
 - (10) Size -> The size class, one of the following: Small (Ø 0.25 mm – 0.16 mm), Medium (Ø 0.45 mm – 0.25 mm) or Large (Ø 0.75 mm – 0.63 mm). 
 - (11) Biomass_V -> The viable moss biomass (in g). Biomass that seem alive/viable (result obtained after finishing the culture period).
 - (12) Biomass_NV -> The non-viable moss biomass (in g). Biomass that seemed death/non-viable (result obtained after finishing the culture period).
 - (13) Biom_total -> The total biomass viable+non-viable (result obtained after finishing the culture period).
 - (14) per_BiomV -> The percentage of viable biomass. Viable biomass divided by the cultured biomass (0.004 g) multiplied by 100.
 - (15) RGR -> Relative Growth Rate (g/month). Calculated using (log(Biom_total) - log(cultured biomass))/2. Cultured biomass was 0.004 g.

# 2 - Dataset_Traits' description

The dataset is a compilation of measurements of 30 fragments of 3 size classes of 6 moss species done with the microscope and image analyses (ImageJ) in the stereomicroscope.
For further information of the image analyses please look at the section Analyze Particles of ImageJ: https://imagej.nih.gov/ij/docs/menus/analyze.html#ap

The dataset contains 541 rows and 73 columns.
The first row contains the name of the columns, and the rest correspond to the 540 samples.
First 3 columns and column 6 are dedicated to the general identification of the samples. Columns 4 and 5 specifies the type and morfology.
Columns 8 to 37 are image analysis measurements in wet propagules, and 39 to 68 are image analysis measurements in the dry propagules. Note that the propagules are the same but with different hydrated status.
Columns 7 and 38 are for the identification of the propagules when subsetting the propagules in wet and dry conditions, respectively.
The rest of columns (69-73) correspond to measurements or obtained through calculations with the collected data.

Column names and descriptions are:
 - (1) Sp. -> The species name, one of the following: D. scoparium, H. aureum, H. cupressiforme , P. capillare, S. ruralis or T. squarrosa.
 - (2) Size -> The size class, one of the following: Small (Ø 0.25 mm – 0.16 mm), Medium (Ø 0.45 mm – 0.25 mm) or Large (Ø 0.75 mm – 0.63 mm).
 - (3) Number -> The number of the fragment being analysed. Note that is by Specie and Size, from 1 to 30.
 - (4) Type -> The type of fragment. One of the following: S (for shoot type) or L (for leaf type).
 - (5) Morfology -> The morfology of the fragment in relation to the origin part of the moss. Can be Ápice (Apex of the moss), F.Apical (from the apical region), F. Intermedio (Not basal, not apical, medium part) or Macho Enano (stands for small males, dimorphism present in some species).
 - (6) Sample -> The name of the species and the fragment size corresponding to the sample.
 - (7) Label_wet -> The label used in the image analysis with ImageJ.
 
 - (8-37) Image analysis with ImageJ in wet fragments. Important Columns used directly in the paper were:
   - (8) Area_wet -> Area of the fragment.
   - (9) Circ_wet -> Circularity of the fragment. 4π*area/perimeter^2. A value of 1.0 indicates a perfect circle. As the value approaches 0.0, it indicates an increasingly elongated shape.
   - (26) Feret_wet -> Length of the fragment.
 - (38) Label_dry -> The label used in the image analysis with ImageJ.
 - (39-68) Image analysis with ImageJ in dry fragments. Important Columns used directly in the paper were:
   - (39) Area_dry -> Area of the fragment.
   - (57) Feret_dry -> Length of the fragment.
 - (69) Index of viability -> Apparent viability of the propagules, analysed at the microscope. The index and values are fully described in the article.
 - (70) Number of shoots -> Number of meristems/shoots present in the propagule. 1R means the propagule had one branch apart from the principal (only in shoot type propagules, not leaf ones).
 - (71-72) Index1 and Index2 -> Indexes calculated using the index of viability and the contamination index.
 - (73) Index of contamination -> Index in which it was evaluated how contaminated was the propagule, between 1 (low contamination) and 3 (highly invaded).

# 3 - Script description

The script has all the analysis and figures related to the article "*Moss establishment success is determined by the interaction between propagule size and species identity*".
It also presents an Index and each section commented, dedicated to the analysis, tables and figures of the article. Further information is in the code.
