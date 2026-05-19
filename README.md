*A landslide dating framework with sub-revisit interval precision through integration of remote sensing imagery and natural triggers*

This project includes the executable program for landslide dating associated with the article A landslide dating framework with sub-revisit interval precision through integration of remote sensing imagery and natural triggers and the test data available for use. Once the manuscript is accepted, we will release the source code here.

The approach begins by constructing a dense vegetation index (VI) series through the fusion of multi source remote sensing imagery. From this sequence, we extract signature anomalies that quantify the intensity of surface disturbances induced by landslides and provide timing information at the revisit interval. This preliminary date is further refined by integrating anomalies from natural triggers, which enables precision beyond the satellite revisit cycle. Additional methodological details are available in our paper.

**How to Use**

**Step 1:** Download time series data from GEE

In our GEE App (Fig. 1, https://idyllic-nova-460013-e0.projects.earthengine.app/view/landslide-time-series-downloader), download the signature and trigger time series as follows:

1. Select the method for inputting landslide boundaries. Two modes are supported: Draw Polygon and Upload Vector, the latter allowing batch processing of multiple polygons in one vector file.

2. Select the time range for exporting the series. Very long or very short series are not recommended. Long series significantly increase GEE export time, while short series (e.g., less than 3 years) may reduce dating accuracy in Step 2.

3. Choose the vegetation index based on the study area characteristics. NDVI is recommended for its stable and reliable performance.
  
4. Select the trigger types as needed. If prior knowledge of landslide triggers exists, you may exclude unrelated trigger types to save computation time.
   
5. Click Start Export to begin the task.
    
After completion, you will obtain the landslide time series data (Fig. 2) for use in the next step.



**Step 2:** Determine the landslide date based on time series

Before this step, please download our executable program (Fig. 3).

1. Choose the input table. Both single file and batch processing are supported. For batch processing, store multiple CSV files in a single folder and select that folder.

2. Choose the output file (optional). The program can save the dating process as images and summarize the landslide dates in a single table. Leave this field empty if you do not need to save the results.

3. Click Start Analysis to begin the calculation. Once completed, you can view the dating results directly in the app (Fig. 4) or in the saved output files (Fig. 5).


**Cite Us**

If this project is helpful to your research, please cite our article:

A landslide dating framework with sub revisit interval precision through integration of remote sensing imagery and natural triggers.




<img width="974" height="620" alt="image" src="https://github.com/user-attachments/assets/04da5420-4a3e-4fc0-ae95-49c3f17f24bd" />

Fig. 1 Description of the Landslide Time Series Download Tool.

 <img width="984" height="692" alt="image" src="https://github.com/user-attachments/assets/a9074515-6ca3-4b17-8d80-3f44e714a3b0" />

Fig. 2 Landslide time series table obtained from Step 1.

 <img width="980" height="511" alt="image" src="https://github.com/user-attachments/assets/47e50b2b-5756-40f9-8ba4-cbf5bded1b30" />

Fig. 3 Local landslide dating executable program.

 <img width="917" height="1368" alt="image" src="https://github.com/user-attachments/assets/e3e74731-6c08-4e85-b499-552a190646d2" />

Fig. 4 Visualization interface of the dating result.  

<img width="976" height="619" alt="image" src="https://github.com/user-attachments/assets/de831df7-0b3d-4df9-9802-5b243b651efa" />

Fig. 5 Example of saved output.
