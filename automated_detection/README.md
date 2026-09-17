# Automated Detection Pipeline for PRIME

A data processing and machine learning pipeline for automatically detecting transients in PRIME observations by Fiona McDaniel 

### Resources

Tutorial notebook: automated_detection/AutomatedDetection.ipynb

Methods paper: docs/AutomatedDetection.pdf (will update)


### Pipeline Architecture

Input: Pass in a list of externally detected transients
1. Reduce PRIME candidate observation fields
2. Subtract science and reference epochs
3. Build a training dataset of images and relevant analysis features (metadata)
4. Train and test random forrest and CNN classifiers on image and metadata
5. Simulate using these models on unseen observations

![Architecture Diagram](https://github.com/orionguiffreda/prime-photometry/blob/ml_analysis/automated_detection/images/architecture_chart.png)

