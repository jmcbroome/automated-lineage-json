# automated-lineage-json
Streamlit webapp for the annotation of Nextstrain Auspice v2 JSON with 

## Summary

This repository contain both a CLI implementation and a Streamlit-based GUI for the application of a lineage nomenclature based on the [genotype representation heuristic](https://github.com/jmcbroome/automate-lineages-prototype).

A public version of the app can be found at: https://jmcbroome-automated-lineage-json-streamlit-app-3adskh.streamlit.app/

## How To Use

The app takes a set of parameters, mostly numeric, and an input [Auspice v2 JSON](https://docs.nextstrain.org/projects/auspice/en/stable/releases/v2.html). The user sets the parameters in the indicated fields and uploads the JSON from local storage, then clicks the run button.
After completion, an annotated JSON and a simple table of sample-lineage associations become available to download in a zip archive. Once downloaded, the JSON can be drag-and-dropped into the Auspice view directly below the main app, or into a separate tab of auspice.us.

### Explanation of Method and Parameters

This method works by computing the genotype representation index, explained in detail [here](https://github.com/jmcbroome/automate-lineages-prototype#mathematical-underpinnings). Essentially, this metric is high for nodes that represent clusters of samples that are distinct from their parent lineage and relatively uniform within the group. It is computed for every node on the input tree, then the node with the highest value is chosen as a new lineage root. This is repeated many times. 

The resulting nomenclature is genotype-based and hierarchical, with a simplified Pango-style naming schema. For example, the lineage A.1.1 is a sublineage of A.1, which in turn is a sublineage of group A. Each of the these are considered a 'level' of annotation- in this example, A.1.1 is in the 3rd level of annotation on the tree.

The nomenclature is generated iteratively; each level is generated as a series of mutually exclusive lineage labels (A,B,C...). After the minimum proportion of samples are labeled with mutually exclusive lineages, each of the resulting labels is independently subdivided by the same process (e.g. A is divided into A.1, A.2, A.3... until the minimum proportion of A samples are labeled wih a child lineage). 

## Run Locally

### Installation

This app has minimal dependencies, relying only on streamlit and Python's standard libraries.

```
pip install streamlit
git clone https://github.com/jmcbroome/automated-lineage-json
cd automated-lineage-json
```

### Command Line Interface

We provide a Python script which has the same parameters and behavior as the app proper.

```
python3 annotate_json.py --help
```

### App

You can spin up a local instance of the Streamlit GUI.

```
streamlit run streamlit_app.py
```