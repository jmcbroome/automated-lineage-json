import streamlit as st
import os
import sys
import json
from annotate_json import *
import zipfile

from streamlit.scriptrunner.scriptrunner_utils.script_run_context import get_script_run_ctx
import streamlit.components.v1 as components

def _get_session():
    ctx = get_script_run_ctx()
    if ctx is None:
        raise Exception("Failed to get the thread context")            
    return ctx.session_id
    
st.set_page_config(layout='wide')
with st.form(key="autolin"):
    st.markdown("# AUTOLIN")
    st.markdown("This app is a tool that uses the [genotype representation index heuristic](https://github.com/jmcbroome/automate-lineages-prototype#mathematical-underpinnings) to add lineage nomenclature labels to a Nextstrain Auspice JSON.")
    st.markdown(" ".join([
        "The generated nomenclature is genotype-based and hierarchical, with a simplified Pango-style naming schema.",
        "For example, the lineage A.1.1 is a sublineage of A.1, which in turn is a sublineage of group A.",
        "Each of the these would be considered a 'level' of annotation.",
        "The nomenclature is generated iteratively; each 'level' is generated as a series of mutually exclusive lineage labels (A,B,C...).",
        "After the minimum proportion of samples are labeled with mutually exclusive lineages, each of the resulting labels is independently",
        "subdivided by the same process (e.g. A is divided into A.1, A.2, A.3... until the minimum proportion of A samples are labeled with an A.X lineage).",
        "Lineage label generation ceases when no candidate lineage roots fulfill conditions set by the user or the maximum number of levels have been generated."
    ]))
    st.markdown("This tool takes specifically Auspice v2 format JSON that include mutation annotations, ideally including at least one set of amino acid change translations. One example is Michael Wolfinger's excellent CHIKV Nextstrain build, which can be found [here](https://github.com/ViennaRNA/CHIKV/blob/master/auspice/CHIKV.json). Numerous others can be found under Nextstrains community builds on Github, built from raw read data with the Augur pipeline, or exported from a [MAT](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/) with [matUtils](https://github.com/yatisht/usher).")
    st.markdown("The Nextstrain JSON files produced by this tool can be uploaded to [Auspice](https://auspice.us/) for viewing. For convenience, a view of auspice.us is embedded below.")
    size = st.number_input("Output lineages will contain at least this many samples.",min_value=1)
    distinction = st.number_input("Output lineages will have at least this many mutations distinguishing them from their parent lineage or the tree root.",min_value=1)
    cutoff = st.number_input("Proportion of samples that should be covered at each level of lineage annotation.",min_value=0.0,max_value=1.0,value=0.90)
    levels = st.number_input("Maximum number of levels to generate. Set to 0 to generate as many as possible.",min_value=0)
    floor = st.number_input("Minimum genotype representation index to annotate a lineage. This value considers both the number and distinction of descendent samples- a value of 1 means a lineage that represents an average of 1 mutation for a randomly chosen sample from the tree. Set to higher values to exclude small, marginal lineages.",min_value=0)
    missense = st.checkbox("Consider only amino-acid altering mutations across the genome.")
    gene = st.text_input("Limit considered mutations to amino-acid altering mutations in one or more specific genes, comma delinated, named here. Leave blank to consider mutations in any gene. Ensure that the genes are present in your input JSON!",value="")
    uploaded_file = st.file_uploader("Upload a JSON to generate lineage labels from.")
    runbutton = st.form_submit_button(label='Generate the labeled JSON and tables.')
    st.markdown("Once downloaded, you can drag and drop the annotated JSON into the view below, or to a [separate tab.](https://auspice.us/)")
    st.markdown("You have to download the results file first as Auspice is rendered client-side.")
    st.markdown("Once uploaded to Auspice, the different lineage label levels can be viewed using the 'Color By' dropdown menu, as the 'GRI Lineage Level X' labels.")
    st.markdown("Any questions, problems, or suggestions can be posted on the [Github repo!](https://github.com/jmcbroome/automated-lineage-json/issues)")

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
        pipeline(ijd,"annotated.json",floor,size,distinction,cutoff,missense,genearg,levels,"labels.tsv","report.tsv")
        with zipfile.ZipFile(pref+'_results.zip','w') as zipf:
            zipf.write("annotated.json")
            zipf.write("labels.tsv")
            zipf.write("report.tsv")
        with open(pref+"_results.zip","rb") as f:
            db = st.download_button(label="Download Annotated JSON and Table in ZIP Format", file_name="results.zip", data=f.read())
            if db != None:
                print("Attempting to clear temp files.")
                for seshfile in ["annotated.json","labels.tsv","report.tsv",pref+"_results.zip"]:
                    if os.path.exists(seshfile):
                        print("Clearing temporary file: " + seshfile,file=sys.stderr)
                        os.remove(seshfile)
components.iframe("https://auspice.us/", height=1000, scrolling=True)
