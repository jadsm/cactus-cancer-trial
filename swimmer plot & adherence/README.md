# Swimmer plot
 
 This project is a swimmerplot, which is a type of graph used to visualize the survival data of patients in a clinical study. 
 It displays the time to an event (such as death or disease progression) on the x-axis and the status of the event (alive or dead) on the y-axis.
 Each patient is represented by a line segment that starts at the time of entry into the study and ends at the time of the event or censoring.
 The swimmerplot is useful for comparing the survival outcomes of different groups or treatments in the study.

## Installation
- Install python 3.11 or above
- Install pip 
- Install all required packages:
    ```bash
    pip install -r requirements.txt
    ```

## Run the code
- run the code
```bash
    python swimmerplot.py
    ```

## View the chart
- start a simple server. There are may ways, one way is with python: 
```bash
    python -m http.server
    ```
- navigate in your favourite browser to: `http://localhost:8000/timeline_swimmer.html`
