* Encoding: UTF-8.
*LEIDA analyse

*Baseline characteristics, to check for table. 
*file used: LEIDAcharac.sav

*05-09-2018: fout! Nieuwe file gemaakt: LEiDAcharacnew.sav

MATCH FILES /FILE= 'DataSet1'
  /TABLE='DataSet2'
  /BY ppnr.
EXECUTE.

CROSSTABS
  /TABLES=patient BY geslacht
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ PHI
  /CELLS=COUNT ROW COLUMN TOTAL
  /COUNT ROUND CELL.

**niet significant verschil in hele groep, niet sign in LEIDA groep** 

*** calculate differences in age ***

T-TEST GROUPS=patient(0 1)
  /MISSING=ANALYSIS
  /VARIABLES=leeftijd
  /CRITERIA=CI(.95).

**niet significant verschil in hele groep, niet sign in LEIDA groep** 
*** calculate differences in education distribution ***

CROSSTABS
  /TABLES=patient BY opleiding
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ PHI
  /CELLS=COUNT ROW COLUMN TOTAL
  /COUNT ROUND CELL.

**niet significant verschil in hele groep, niet sign in LEIDA groep** 
*** calculate differences in IQ ***

MEANS TABLES= IQ_schatting BY patient
  /CELLS=MEAN COUNT STDDEV MIN MAX.

NPAR TESTS
  /M-W=  IQ_schatting BY patient(1 0)
  /MISSING ANALYSIS.

**niet significant verschil in hele groep, niet sign in LEIDA groep** 
*** calculate differences in living situation distribution ***

CROSSTABS
  /TABLES=patient BY woonsituatie
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ PHI
  /CELLS=COUNT ROW COLUMN TOTAL
  /COUNT ROUND CELL.

**niet significant verschil in hele groep, niet sign in LEIDA groep** 
*** calculate differences in working class distribution ***

CROSSTABS
  /TABLES=patient BY beroepskls
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ PHI
  /CELLS=COUNT ROW COLUMN TOTAL
  /COUNT ROUND CELL.

**niet significant verschil in hele groep, niet sign in LEIDA groep** 
*** calculate differences in handedness distribution ***

CROSSTABS
  /TABLES=patient BY handvoorkeur
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ PHI
  /CELLS=COUNT ROW COLUMN TOTAL
  /COUNT ROUND CELL.

*** age of onset for patients ***

COMPUTE gbjaar=XDATE.YEAR(geboortedt).
EXECUTE.

COMPUTE age_onset=Bjrederp-gbjaar.
EXECUTE.
DESCRIPTIVES VARIABLES=age_onset
  /STATISTICS=MEAN STDDEV MIN MAX.


*** calculate mean previous episodes for patients ***

MEANS TABLES=deprepleven BY patient
  /CELLS=MEAN COUNT STDDEV.

*** calculate differences in HDRS scores ***

MEANS TABLES=HDRS_T0_combined BY patient
  /CELLS=MEAN COUNT STDDEV MIN MAX.

NPAR TESTS
  /M-W= HDRS_totaal_MRI BY patient(1 0)
  /MISSING ANALYSIS.

NPAR TESTS
  /M-W= HDRS_T0_combined BY patient(1 0)
  /MISSING ANALYSIS.

* nog uit CRF halen




FREQUENCIES VARIABLES=HDRS_totaal_MRI HDRS_T0_combined deprepleven
  /NTILES=4
  /STATISTICS=MEDIAN
  /ORDER=ANALYSIS.




SORT CASES  BY patient.
SPLIT FILE LAYERED BY patient.

DATASET ACTIVATE DataSet3.
FREQUENCIES VARIABLES=HDRS_T0_combined
  /FORMAT=NOTABLE
  /NTILES=4
  /STATISTICS=MINIMUM MAXIMUM MEDIAN
  /ORDER=ANALYSIS.

DATASET ACTIVATE DataSet3.
FREQUENCIES VARIABLES=HDRS_totaal_MRI
  /FORMAT=NOTABLE
  /NTILES=4
  /STATISTICS=MINIMUM MAXIMUM MEDIAN
  /ORDER=ANALYSIS.

COMPUTE intake_MRI_tijd=dtMRI - dtintake.
EXECUTE.

COMPUTE intake_MRI_tijd=TIME.DAYS(dtMRI - dtintake).
EXECUTE.

 COMPUTE days = CTIME.DAYS(dtMRI - dtintake) .
  EXECUTE .

*Mood ratings:

Mood repeated measures*

DATASET ACTIVATE DataSet2.

COMPUTE DeltaneutralMood= MRI_cijfer_st2 - MRI_cijfer_st1 .
EXECUTE.

COMPUTE DeltasadMood= MRI_cijfer_st4 - MRI_cijfer_st3.
EXECUTE. 

SORT CASES  BY patient.
SPLIT FILE LAYERED BY patient.

FREQUENCIES VARIABLES=MRI_cijfer_st1 MRI_cijfer_st2 MRI_cijfer_st3 MRI_cijfer_st4 DeltaneutralMood DeltasadMood MRI_cijfer_st5
  /FORMAT=NOTABLE
  /STATISTICS=STDDEV MEAN
  /ORDER=ANALYSIS.

SPLIT FILE OFF.

GLM MRI_cijfer_st2 MRI_cijfer_st1 BY patient
  /WSFACTOR=mood 2 Polynomial 
  /METHOD=SSTYPE(3)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=mood 
  /DESIGN=patient.

GLM MRI_cijfer_st3 MRI_cijfer_st4 BY patient
  /WSFACTOR=mood 2 Polynomial 
  /METHOD=SSTYPE(3)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=mood 
  /DESIGN=patient.

GLM MRI_cijfer_st1 MRI_cijfer_st2 MRI_cijfer_st3 MRI_cijfer_nst4 BY patient
  /WSFACTOR=mood 4 Polynomial 
  /METHOD=SSTYPE(3)
  /POSTHOC=patient(BONFERRONI) 
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=mood 
  /DESIGN=patient.

SPLIT FILE OFF.

NPAR TESTS
  /M-W= MRI_cijfer_st1 BY patient(0 1)
  /MISSING ANALYSIS.

NPAR TESTS
  /M-W= MRI_cijfer_st2 BY patient(0 1)
  /MISSING ANALYSIS.

NPAR TESTS
  /M-W= MRI_cijfer_st3 BY patient(0 1)
  /MISSING ANALYSIS.

NPAR TESTS
  /M-W= MRI_cijfer_st4 BY patient(0 1)
  /MISSING ANALYSIS.

NPAR TESTS
  /M-W= MRI_cijfer_st5 BY patient(0 1)
  /MISSING ANALYSIS.

*after revision: check if correlations with HDRS:

*check interactions



DATASET ACTIVATE DataSet4.
GLM ifetime10states ifetime10states_sad BY Patient
  /WSFACTOR=mood 2 Polynomial 
  /METHOD=SSTYPE(3) 
  /PLOT=PROFILE(mood*Patient) TYPE=LINE ERRORBAR=NO MEANREFERENCE=NO YAXIS=AUTO
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=mood 
  /DESIGN=Patient.





GLM Probability_10states Probability_10states_sad BY Patient
  /WSFACTOR=mood 2 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(mood*Patient) TYPE=LINE ERRORBAR=NO MEANREFERENCE=NO YAXIS=AUTO
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=mood 
  /DESIGN=Patient.

*regression analysis HDRS influence


DATASET ACTIVATE DataSet1.
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT Patient
  /METHOD=ENTER HDRS Probability_10states.


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT Patient
  /METHOD=ENTER HDRS ifetime10states.

*number of psychiatric diagnoses correlated to lifetime or probability?

SORT CASES  BY Patient.
SPLIT FILE LAYERED BY Patient.

.
NONPAR CORR
  /VARIABLES=Probability_10states psychcm_ndsm4 ifetime10states
  /PRINT=SPEARMAN TWOTAIL NOSIG
  /MISSING=PAIRWISE.

GRAPH
  /HISTOGRAM=psychcm_ndsm4
  /PANEL ROWVAR=Patient ROWOP=CROSS.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT Patient
  /METHOD=ENTER Probability_10states psychcm_ndsm4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT Patient
  /METHOD=ENTER psychcm_ndsm4 ifetime10states.

*effect sizes:

SORT CASES  BY Patient.
SPLIT FILE LAYERED BY Patient.


FREQUENCIES VARIABLES=ifetime10states Probability_10states Probability_10states_sad 
    ifetime10states_sad
  /FORMAT=NOTABLE
  /STATISTICS=STDDEV MEAN
  /ORDER=ANALYSIS.
