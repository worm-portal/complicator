import pandas as pd
import math
from scipy.optimize import minimize
from .autopp93s4 import complicate
import plotly.express as px
import plotly.graph_objects as go
import plotly.colors
from wormutils import format_equation
from ipywidgets import interactive, interact_manual
interact_manual.opts['manual_name'] = 'Calculate'
import pychnosz
_ = pychnosz.thermo("WORM", messages=False)

@interact_manual(n_complex=(1, 4, 1), metal="Co+2", ligand="Cl-", logK_data="demo.csv", data_path="WORM", log_beta=(-5, 5, 0.001), S=(-1000, 1000, 0.1), Cp=(-1000, 1000, 0.1), V=(-1000, 1000, 0.1), auto_fit=False)
def complex_fit(n_complex=4, metal="Co+2", ligand="Cl-", logK_data="demo.csv", data_path="WORM", log_beta=0, S=0, Cp=0, V=0, auto_fit=False):

    if data_path == "WORM":
        data_path = None
    else:
        _ = pychnosz.add_OBIGT(data_path, force=True, messages=False)
    
    df_logK = pd.read_csv(logK_data)
    df_logK["legendgroup"] = df_logK.apply(lambda x: str(x["P"]) + " bars (" + str(x["ref"]) + ")" if str(x["P"]) != "Psat" else "Psat (" + str(x["ref"]) + ")", axis=1)

    if auto_fit:
        def _get_predicted_logK(lb, s, cp, v):
            """Run complicate + subcrt and return list of predicted logK values, or None on failure."""
            d = {
                "Metal": [metal],
                "Ligand": [ligand],
                "BETA_1": [0], "BETA_2": [0], "BETA_3": [0], "BETA_4": [0],
                "S_"+str(n_complex): [s],
                "Cp_"+str(n_complex): [cp],
                "V_"+str(n_complex): [v],
            }
            d["BETA_"+str(n_complex)] = [lb]
            df_in = pd.DataFrame(d)
            df_o, _, _, _ = complicate(df_in=df_in, data_path=data_path, print_warnings=False)
            if df_o is None or df_o.empty or len(df_o) < n_complex:
                return None
            cname = df_o["name"].iloc[n_complex-1]
            _ = pychnosz.add_OBIGT(df_o, force=True, messages=False)
            predictions = []
            for i in range(len(df_logK)):
                P_val = df_logK["P"].iloc[i]
                if P_val != "Psat":
                    P_val = float(P_val)
                try:
                    pred = pychnosz.subcrt(
                        [metal, ligand, cname], [-1, -n_complex, 1],
                        T=df_logK["T,C"].iloc[i], P=P_val,
                        show=False, messages=False).out
                    if "logK" in pred.columns:
                        predictions.append(pred["logK"].values[0])
                    else:
                        predictions.append(float('nan'))
                except:
                    predictions.append(float('nan'))
            return predictions

        # Fix log_beta from 25°C Psat data point if available
        ref_mask = (df_logK["T,C"] == 25) & (df_logK["P"].astype(str) == "Psat")
        if ref_mask.any():
            log_beta = df_logK.loc[ref_mask, "logK"].iloc[0]
        elif log_beta == 0:
            idx_lowT = df_logK["T,C"].values.argmin()
            log_beta = df_logK["logK"].iloc[idx_lowT]

        # Test that the model works with the fixed log_beta
        test_pred = _get_predicted_logK(log_beta, S, Cp, V)
        if test_pred is None:
            print("Warning: model cannot be evaluated. Check metal/ligand/data_path settings.")
        else:
            def _objective(params):
                s, cp, v = params
                preds = _get_predicted_logK(log_beta, s, cp, v)
                if preds is None:
                    return 1e10
                residuals = 0.0
                n_valid = 0
                for i in range(len(df_logK)):
                    if preds[i] is None or math.isnan(preds[i]):
                        continue
                    err_val = df_logK["err"].iloc[i]
                    if pd.isna(err_val) or err_val <= 0:
                        err_val = 1.0
                    weight = 1.0 / err_val
                    residuals += (weight * (df_logK["logK"].iloc[i] - preds[i])) ** 2
                    n_valid += 1
                if n_valid == 0:
                    return 1e10
                return residuals / n_valid

            x0 = [S, Cp, V]
            print("Fixing log_beta={:.4f}. Optimizing S, Cp, V...".format(log_beta), end="", flush=True)

            _iter_count = [0]
            def _callback(xk):
                _iter_count[0] += 1
                if _iter_count[0] % 5 == 0:
                    print(".", end="", flush=True)

            result = minimize(_objective, x0, method='Powell',
                              callback=_callback,
                              options={'xtol': 0.001, 'ftol': 0.01, 'maxiter': 2000})
            print(" done.")
            S, Cp, V = result.x
            S = max(-1000, min(1000, S))
            Cp = max(-1000, min(1000, Cp))
            V = max(-1000, min(1000, V))
            print("Auto-fit result: log_beta={:.4f}, S={:.2f}, Cp={:.2f}, V={:.2f}".format(
                log_beta, S, Cp, V))

    dict_input = {
            "Metal":[metal],
            "Ligand":[ligand],
            "BETA_1":[0],
            "BETA_2":[0],
            "BETA_3":[0],
            "BETA_4":[0],
            "S_"+str(n_complex):[S],
            "Cp_"+str(n_complex):[Cp],
            "V_"+str(n_complex):[V],
            }
    dict_input["BETA_"+str(n_complex)] = [log_beta]
    
    df_input = pd.DataFrame(dict_input)
    
    df_out, _, _, _ = complicate(df_in=df_input, data_path=data_path, print_warnings=False)

    if df_out.empty:
        print("The metal or ligand was not found in the data path.")
        return
    
    complex_name = df_out["name"][n_complex-1]

    G = df_out[df_out["name"] == complex_name]["G"].values[0]
    H = df_out[df_out["name"] == complex_name]["H"].values[0]
    
    _ = pychnosz.add_OBIGT(df_out, force=True, messages=False)

    for i,P in enumerate(list(dict.fromkeys(df_logK["P"]))):
        if P != "Psat":
            Pfloat = float(P)
            est_name = "{} bar (model)".format(str(P))
        else:
            Pfloat = P
            est_name = "{} (model)".format(str(P))
        T_min, T_max, T_res = min(df_logK["T,C"]), max(df_logK["T,C"]), 15
        T = [T_min + x*(T_max-T_min)/T_res for x in range(T_res)]+[T_max]
        out_estimate = pychnosz.subcrt([metal, ligand, complex_name], [-1, -n_complex, 1],
                                       T=T, P=Pfloat,
                                       show=False, messages=False).out
        out_estimate["est_name"] = est_name
        if i == 0:
            df_plot = out_estimate
        else:
            df_plot = pd.concat([df_plot, out_estimate], axis=0)

    fig = px.line(df_plot, x="T", y="logK", color="est_name", template="simple_white")

    legendgroups_added = []
    for i,P in enumerate(df_logK["P"]):
        if P != "Psat":
            Pfloat = float(P)
        else:
            Pfloat = P
        
        logK_pred = pychnosz.subcrt([metal, ligand, complex_name], [-1, -n_complex, 1],
                                    T=df_logK["T,C"][i], P=Pfloat,
                                    show=False, messages=False).out
        logK_pred = logK_pred["logK"].values[0]

        if df_logK["legendgroup"][i] in legendgroups_added:
            showlegend = False
        else:
            showlegend = True
            legendgroups_added.append(df_logK["legendgroup"][i])
        
        line_col = plotly.colors.qualitative.D3[len(legendgroups_added)-1]
        err_val = df_logK["err"][i]
        has_err = not pd.isna(err_val)
        if math.isclose(df_logK["logK"][i], logK_pred, rel_tol=0.1) or (has_err and math.isclose(df_logK["logK"][i], logK_pred, abs_tol=err_val)):
            # filled circle
            fill_col = line_col
        else:
            # open circle
            fill_col = 'rgba(255, 0, 0, 0.0)'

        error_y_dict = dict(
                        type='data',
                        array=[err_val if has_err else 0],
                        color=line_col,
                        visible=has_err)

        fig.add_trace(go.Scatter(
            x=[df_logK['T,C'][i]],
            y=[df_logK['logK'][i]],
            legendgroup=df_logK["legendgroup"][i],
            name=df_logK["legendgroup"][i],
            showlegend=showlegend,
            mode="markers",
            marker=dict(
                        size=10,
                        symbol='circle',
                        line_width=2,
                        line_color=line_col,
                        color=fill_col),
                error_y=error_y_dict,
                )
        )
    
    #fig.update_layout(legend=dict(groupclick="toggleitem"))

    a1 = df_out[df_out["name"]==complex_name]["a1.a"].values[0]
    a2 = df_out[df_out["name"]==complex_name]["a2.b"].values[0]
    a3 = df_out[df_out["name"]==complex_name]["a3.c"].values[0]
    a4 = df_out[df_out["name"]==complex_name]["a4.d"].values[0]
    c1 = df_out[df_out["name"]==complex_name]["c1.e"].values[0]
    c2 = df_out[df_out["name"]==complex_name]["c2.f"].values[0]
    omega = df_out[df_out["name"]==complex_name]["omega.lambda"].values[0]

    formatted_rxn = format_equation(species=[metal, ligand, complex_name],
                                       stoich=[-1, -n_complex, 1], sign="equilibrium")
    
    annotation_text = (formatted_rxn+""
                       "<br>logβ<sub>{}</sub> = {}".format(n_complex, round(log_beta, 2))+""
                       "<br>∆<sub>f</sub>G° = {} cal mol<sup>-1</sup>".format(int(G))+""
                       "<br>∆<sub>f</sub>H° = {} cal mol<sup>-1</sup>".format(int(H))+""
                       "<br>S° = {} cal mol<sup>-1</sup>K<sup>-1</sup>".format(round(S, 3))+""
                       "<br>Cp° = {} cal mol<sup>-1</sup>K<sup>-1</sup>".format(round(Cp, 3))+""
                       "<br>V° = {} cm<sup>3</sup> mol<sup>-1</sup>".format(round(V, 3))+""
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
        yaxis_title="logβ<sub>{}</sub>".format(n_complex),
        title=None,
        width=800,
        height=350,
        legend_title=None,
        font=dict(size=12),
        title_x=0.5,
        # xaxis_range=[0, 360],
        # yaxis_range=[-4, 8],
        margin=dict(l=50,r=50,b=0,t=1),
    )
    
    config = {'toImageButtonOptions': {'format': "png", # one of png, svg, jpeg, webp
                                      },
             }
    
    
    fig.show(config=config)