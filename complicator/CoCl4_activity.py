import pandas as pd
import math
from .autopp93s4 import complicate
import plotly.express as px
import plotly.graph_objects as go
from ipywidgets import interactive, interact_manual
interact_manual.opts['manual_name'] = 'Calculate'

import pychnosz
_ = pychnosz.thermo("WORM", messages=False)

@interact_manual(log_beta=(-5, 5, 0.01), S=(-100, 100, 1), Cp=(-100, 100, 1), V=(-100, 150, 1))
def CoCl4_2_activity(log_beta=0, S=0, Cp=0, V=0):

    df_input = pd.DataFrame({
            "Metal":["Co+2"],
            "Ligand":["Cl-"],
            "BETA_1":[0],
            "BETA_2":[0],
            "BETA_3":[0],
            "BETA_4":[log_beta],
            "S_4":[S],
            "Cp_4":[Cp],
            "V_4":[V],
            })
    
    df_out, _, _, _ = complicate(df_in=df_input, print_warnings=False)

    G = df_out[df_out["name"] == "Co(Cl)4-2"]["G"].values[0]
    H = df_out[df_out["name"] == "Co(Cl)4-2"]["H"].values[0]
    
    _ = pychnosz.add_OBIGT(df_out, force=True, messages=False)
    T = [0.01, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350]
    out_estimate_psat = pychnosz.subcrt(["Co+2", "Cl-", "Co(Cl)4-2"], [-1, -4, 1], T=T, show=False, messages=False).out
    out_estimate_600bars = pychnosz.subcrt(["Co+2", "Cl-", "Co(Cl)4-2"], [-1, -4, 1], T=T, P=600, show=False, messages=False).out

    out_estimate_psat["type"] = "Psat (model)"
    out_estimate_600bars["type"] = "600 bars (model)"
    df_plot = pd.concat([out_estimate_psat, out_estimate_600bars], axis=0)
    
    fig = px.line(df_plot, x="T", y="logK", color="type", template="simple_white")

    
    fig.add_trace(go.Scatter(x=[25, 50, 90],
                             y=[-2.09, -1.72, -1.41],
                             mode='markers', name="Psat (Pan and Susak, 1990)",
                             marker=dict(
                                line_color='#1F77B4',
                                color='#17BECF',
                                size=10,
                                symbol="circle",
                                line_width=2,
                                ),
                             error_y=dict(
                                type='data',
                                array=[0.05, 0.05, 0.05],
                                color='#1F77B4',
                                visible=True)
                             )
                            
                 )
    
    fig.add_trace(go.Scatter(x=[150, 250, 350],
                             y=[-0.69, 2.16, 5.68],
                             mode='markers', name="600 bars (Liu et al., 2011)",
                             marker=dict(
                                line_color='#FF7F0E',
                                color='#FECB52',
                                size=10,
                                symbol="circle",
                                line_width=2,
                                ),
                             error_y=dict(
                                type='data',
                                array=[0.2, 0.15, 0.6],
                                color='#FF7F0E',
                                visible=True)
                             )
                 )
    


    a1 = df_out[df_out["name"]=="Co(Cl)4-2"]["a1.a"].values[0]
    a2 = df_out[df_out["name"]=="Co(Cl)4-2"]["a2.b"].values[0]
    a3 = df_out[df_out["name"]=="Co(Cl)4-2"]["a3.c"].values[0]
    a4 = df_out[df_out["name"]=="Co(Cl)4-2"]["a4.d"].values[0]
    c1 = df_out[df_out["name"]=="Co(Cl)4-2"]["c1.e"].values[0]
    c2 = df_out[df_out["name"]=="Co(Cl)4-2"]["c2.f"].values[0]
    omega = df_out[df_out["name"]=="Co(Cl)4-2"]["omega.lambda"].values[0]
    
    annotation_text = ("Co<sup>2+</sup> + 4Cl<sup>-</sup> = CoCl<sub>4</sub><sup>2-</sup>"
                       "<br>logβ<sub>4</sub> = {}".format(round(log_beta, 2))+""
                       "<br>∆<sub>f</sub>G° = {} cal mol<sup>-1</sup>".format(int(G))+""
                       "<br>∆<sub>f</sub>H° = {} cal mol<sup>-1</sup>".format(int(H))+""
                       "<br>S° = {} cal mol<sup>-1</sup>K<sup>-1</sup>".format(S)+""
                       "<br>Cp° = {} cal mol<sup>-1</sup>K<sup>-1</sup>".format(Cp)+""
                       "<br>V° = {} cm<sup>3</sup> mol<sup>-1</sup>".format(V)+""
                       "<br>a<sub>1</sub> × 10 = {} cal mol<sup>-1</sup> bar<sup>-1</sup>".format(a1)+""
                       "<br>a<sub>2</sub> × 10<sup>-2</sup> = {} cal mol<sup>-1</sup>".format(a2)+""
                       "<br>a<sub>3</sub> = {} cal K mol<sup>-1</sup> bar<sup>-1</sup>".format(a3)+""
                       "<br>a<sub>4</sub> × 10<sup>-4</sup> = {} cal K mol<sup>-1</sup>".format(a4)+""
                       "<br>c<sub>1</sub> = {} cal mol<sup>-1</sup>".format(c1)+""
                       "<br>c<sub>2</sub> × 10<sup>-4</sup> = {} cal K mol<sup>-1</sup>".format(c2)+""
                       "<br>ω × 10<sup>-5</sup> = {} cal K mol<sup>-1</sup>".format(omega)+""
                       )
    
    fig.add_annotation(
        x=0.05, y=1,
        xref='paper', yref='paper',
        showarrow=False,
        text=annotation_text,
        align="left",
        font = dict(size = 10),
    )
    

    fig.update_layout(
        xaxis_title="T, °C",
        yaxis_title="logβ<sub>4</sub>",
        title=None,
        width=600,
        height=350,
        legend_title=None,
        font=dict(size=12),
        title_x=0.5,
        xaxis_range=[0, 360],
        yaxis_range=[-4, 8],
        margin=dict(l=50,r=50,b=0,t=1),
    )
    
    config = {'toImageButtonOptions': {'format': "png", # one of png, svg, jpeg, webp
                                      },
             }
    
    
    fig.show(config=config)