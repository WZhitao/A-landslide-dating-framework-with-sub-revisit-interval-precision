**A landslide dating framework with sub-revisit interval precision through integration of remote sensing imagery and natural triggers**
This project includes the executable program for landslide dating associated with the article A landslide dating framework with sub-revisit interval precision through integration of remote sensing imagery and natural triggers and the test data available for use. Once the manuscript is accepted, we will release the source code here.
The approach begins by constructing a dense vegetation index (VI) series through the fusion of multi source remote sensing imagery. From this sequence, we extract signature anomalies that quantify the intensity of surface disturbances induced by landslides and provide timing information at the revisit interval. This preliminary date is further refined by integrating anomalies from natural triggers, which enables precision beyond the satellite revisit cycle. Additional methodological details are available in our paper.

**How to Use**
**Step 1:** Download time series data from GEE
In our GEE App (Fig. 1), download the signature and trigger time series as follows:
Select the method for inputting landslide boundaries. Two modes are supported: Draw Polygon and Upload Vector, the latter allowing batch processing of multiple polygons in one vector file.
Select the time range for exporting the series. Very long or very short series are not recommended. Long series significantly increase GEE export time, while short series (e.g., less than 3 years) may reduce dating accuracy in Step 2. 
Choose the vegetation index based on the study area characteristics. NDVI is recommended for its stable and reliable performance.
Select the trigger types as needed. If prior knowledge of landslide triggers exists, you may exclude unrelated trigger types to save computation time.
Click Start Export to begin the task.
After completion, you will obtain the landslide time series data (Fig. 2) for use in the next step.

**Step 2:** Determine the landslide date based on time series
Before this step, please download our executable program (Fig. 3).
Choose the input table. Both single file and batch processing are supported. For batch processing, store multiple CSV files in a single folder and select that folder.
Choose the output file (optional). The program can save the dating process as images and summarize the landslide dates in a single table. Leave this field empty if you do not need to save the results.
Click Start Analysis to begin the calculation. Once completed, you can view the dating results directly in the app (Fig. 4) or in the saved output files (Fig. 5).

**Cite Us**
If this project is helpful to your research, please cite our article:
A landslide dating framework with sub revisit interval precision through integration of remote sensing imagery and natural triggers.
<img width="976" height="619" alt="image" src="https://github.com/user-attachments/assets/ab959bf5-591e-4e5f-82ba-9ec0234bf893" />
