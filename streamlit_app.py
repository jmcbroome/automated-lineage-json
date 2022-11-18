import streamlit as st
import os
import sys
import json
from annotate_json import *

from streamlit.runtime.scriptrunner.script_run_context import get_script_run_ctx
import streamlit.components.v1 as components

def _get_session():
    ctx = get_script_run_ctx()
    if ctx is None:
        raise Exception("Failed to get the thread context")            
    return ctx.session_id
    
st.set_page_config(layout='wide')
with st.form(key="autolin"):
    st.markdown("# AUTOLIN")
    st.markdown("This app is a tool that uses the genotype representation score heuristic to add lineage nomenclature labels to a Nextstrain Auspice JSON.")
    st.markdown("The Nextstrain JSON files produced by this tool can be uploaded to [Auspice](https://auspice.us/) for viewing.")
    size = st.number_input("Minimum number of samples to define a lineage.",min_value=1)
    distinction = st.number_input("Minimum number of distinguishing mutations to define a lineage.",min_value=1)
    cutoff = st.number_input("Proportion of samples that should be covered at each level of lineage annotation.",min_value=0.0,max_value=1.0,value=0.95)
    levels = st.number_input("Maximum number of levels to generate. Set to 0 to generate as many as possible.",min_value=0)
    floor = st.number_input("Minimum genotype representation score to annotate a lineage.",min_value=0)
    missense = st.checkbox("Consider amino-acid altering mutations across the genome only.")
    gene = st.text_input("Limit considered mutations to amino-acid altering mutations in a specific gene. Leave blank to consider mutations in any gene.",value="")
    uploaded_file = st.file_uploader("Choose a JSON to generate lineage labels from.")
    runbutton = st.form_submit_button(label='Generate the labeled JSON and table.')
    st.markdown("Once downloaded, you can drag and drop the results file into the view below, or to a [separate tab.](https://auspice.us/)")
    st.markdown("You have to download the results file first as Auspice is rendered client-side.")

pref = _get_session()
if runbutton:
    if uploaded_file == None:
        st.write("ERROR: Upload a file first!")
    else:
        ijd = json.load(uploaded_file)
        if gene == "":
            genearg = None
        else:
            genearg = gene
        pipeline(ijd,pref+"_annotated.json",floor,size,distinction,cutoff,missense,genearg,levels,pref+"_labels.tsv")
        with open(pref+'_annotated.json', 'r') as f:
            with open(pref+"_labels.tsv", 'r') as f2:
                db = st.download_button(label="Download Annotated JSON", file_name="annotated.json", data=f.read())
                db2 = st.download_button(label="Download Sample Lineage Association Table", file_name="labels.tsv", data=f2.read())
                if db:
                    seshfile = pref+"_annotated.json"
                    if os.path.exists(seshfile):
                        print("Clearing temporary file: " + seshfile,file=sys.stderr)
                        os.remove(seshfile)
                if db2:
                    seshfile = pref+"_labels.tsv"
                    if os.path.exists(seshfile):
                        print("Clearing temporary file: " + seshfile,file=sys.stderr)
                        os.remove(seshfile)

components.iframe("https://auspice.us/", height=1000, scrolling=True)