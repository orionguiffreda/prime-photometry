# Automated Detection Pipeline PRIME

A data processing and machine learning pipeline for automatically detecting transients in PRIME observations. 

Pass in a list of externally detected transients, and it will identify and reduce PRIME candidate observation fields, subtract science and reference epochs, build a training dataset of images and relevant analysis features (metadata), train and test random forrest and CNN classifiers on said data, and simulate using these models on unseen observations. 


### Pipeline Architecture
![Architecture Diagram](https://github.com/orionguiffreda/prime-photometry/blob/ml_analysis/automated_detection/images/architecture_chart.png)

### Resources

Tutorial notebook: automated_detection/AutomatedDetection.ipynb
Methods paper: 

