import os
from flask import Flask, render_template
import matplotlib.pyplot as plt
import networkx as nx


network_image_path = r"D:\OneDrive - University of Toronto\19 mTb project\mTb\web\static\network_graph.png"


app = Flask(__name__)

@app.route("/")
def index():
    # The directory for the SAFE results
    safe_output_dir = "data_analysis/results/safe"

    # Bring the file list in there.
    results_files = os.listdir(safe_output_dir)

    return render_template('index.html', results_files=results_files)




print("app.py is running")

if __name__ == "__main__":
    app.run(debug=True)