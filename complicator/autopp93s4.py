# load libraries
import sys, os
import re
from datetime import datetime
import pandas as pd
import decimal
import math
import chemparse
import wormutils


def sf(val, sf):
    val = "{0:.{prec}g}".format(val, prec=sf)
    if "e" in val and abs(float(val)) >= 1:
        val = str(int(float(val)))
    elif "e" not in val and abs(float(val)) < 1:
        val = val + "0"*(sf - find_sigfigs(val))
    elif "e" in val and abs(float(val)) < 1:
        val = "{0:.{prec}f}".format(float(val), prec=20) # add arbitrarily high trailing zeros
        val = re.sub("0+$", "", val) # trim off all trailing zeros
        val = val + "0"*(sf - find_sigfigs(val)) # add back correct number of trailing zeros
    return val


def __get_formula_ox_dict(name, df):
    
    entry = df[df["name"]==name]
    formula_ox = list(entry["formula_ox"])[0]
    formula_ox_split = formula_ox.split(" ")

    # remove any "" (empty strings) from list
    formula_ox_split = list(filter(("").__ne__, formula_ox_split))
    
    formula_ox_dict = {}
    for item in formula_ox_split:
        coeff = re.search(r'^\D*(\d+(?:\.\d+)?)', item)
        if coeff != None:
            coeff = coeff.group()
        if coeff == item or coeff == None:
            coeff = 1
            elem = item
        else:
            elem = item.split(coeff)[1]
        formula_ox_dict[elem] = float(coeff)
        
    return formula_ox_dict


def __write_output(filename, cation, ligand, nth_complex, G, H, S, CP, V,
                   a1, a2, a3, a4, c1, c2, wcon, Z, azero,
                   cation_dissrxn_dict, ligand_dissrxn_dict, thermo_data,
                   cation_formula_ox_dict, ligand_formula_ox_dict,
                   skip_duplicates, ligand_name_abbrv_pairs, organic):
    """
    Write output to a CSV.
    """

    duplicate_list = []
    cat_nocharge = re.sub(r"\+.*$", "", cation)
    lig_nocharge = re.sub(r"\-.*$", "", ligand)

    cat_formula = thermo_data[thermo_data["name"]==cation]["formula"].values[0]
    lig_formula = thermo_data[thermo_data["name"]==ligand]["formula"].values[0]
    cat_formula_nocharge = re.sub(r"\+.*$", "", cat_formula)
    lig_formula_nocharge = re.sub(r"\-.*$", "", lig_formula)
    
    if Z > 1:
        this_charge = "+" + str(Z)
    elif Z == 1:
        this_charge = "+"
    elif Z == 0:
        this_charge = ""
    elif Z == -1:
        this_charge = "-"
    else:
        this_charge = str(Z)
    
    azero_val = 4
    if Z == 2:
        azero_val = 6
    elif Z == 3:
        azero_val = 9
    elif Z == 4:
        azero_val = 11
    
    if nth_complex == 1:
        ligand_subscript = ""
    else:
        ligand_subscript = str(nth_complex)
    
    this_date = datetime.today().strftime('%Y%m%d') 

    complex_name = cat_nocharge + "(" + lig_nocharge + ")" + ligand_subscript + this_charge
    complex_formula_ox = ""
    
    if ligand == "OH-":

        if "(OH)2" in complex_name:
            complex_name = complex_name.replace("(OH)2", "(O)")
            complex_formula_ox = cat_formula + " O-2"
            ligand_dissrxn_dict = {"H+": -2, "H2O": 1}
        elif "(OH)3" in complex_name:
            complex_name = "H" + cat_nocharge + "(O)2" + this_charge
            complex_formula_ox = cat_formula + " H+ 2O-2"
            ligand_dissrxn_dict = {"H+": -3, "H2O": 2}
        elif "(OH)4" in complex_name:
            complex_name = cat_nocharge + "(O)2" + this_charge
            complex_formula_ox = cat_formula + " 2O-2"
            ligand_dissrxn_dict = {"H+": -4, "H2O": 2}

        complex_abbrv = complex_name
        complex_formula = complex_name
    else:
        complex_abbrv = cat_nocharge + "(" + lig_nocharge + ")" + ligand_subscript + this_charge
        complex_formula = cat_formula_nocharge + "(" + lig_formula_nocharge + ")" + ligand_subscript + this_charge

    
    cation_dissrxn_dict_coeff = {name:coeff*1 for name, coeff in cation_dissrxn_dict.items()} # future: "1" can be replaced with number of cations
    
    if ligand != "OH-":
        ligand_dissrxn_dict_coeff = {name:coeff*nth_complex for name, coeff in ligand_dissrxn_dict.items()}
    else:
        ligand_dissrxn_dict_coeff = ligand_dissrxn_dict
    
    cation_ligand_dissrxn_dict = {k: cation_dissrxn_dict_coeff.get(k, 0) + ligand_dissrxn_dict_coeff.get(k, 0) for k in set(cation_dissrxn_dict_coeff) | set(ligand_dissrxn_dict_coeff)}
    cation_ligand_dissrxn_formatted = ["{:.4f}".format(coeff) + " " + name for name, coeff in cation_ligand_dissrxn_dict.items()]
    
    cation_ligand_dissrxn_formatted = " ".join(cation_ligand_dissrxn_formatted)
    
    complex_dissrxn = "-1.0000 " + complex_name + " " + cation_ligand_dissrxn_formatted

    if complex_formula_ox == "":
        cation_formula_ox_dict_coeff = {name:coeff*1 for name, coeff in cation_formula_ox_dict.items()} # future: "1" can be replaced with number of cations
        ligand_formula_ox_dict_coeff = {name:coeff*nth_complex for name, coeff in ligand_formula_ox_dict.items()}
        complex_formula_ox_dict = {k: cation_formula_ox_dict_coeff.get(k, 0) + ligand_formula_ox_dict_coeff.get(k, 0) for k in set(cation_formula_ox_dict_coeff) | set(ligand_formula_ox_dict_coeff)}

        complex_formula_ox = []
        for elem,coeff in complex_formula_ox_dict.items():
            if coeff == 1:
                complex_formula_ox.append(str(elem))
            elif coeff == int(coeff):
                complex_formula_ox.append(str(int(coeff))+str(elem))
            else:
                complex_formula_ox.append(str(coeff)+str(elem))
                
        complex_formula_ox = " ".join(complex_formula_ox)

    if organic:
        cat_1 = "organic_aq"
    else:
        cat_1 = "inorganic_aq"
    
    data={
        "name": [complex_name],
        "abbrv": [complex_abbrv],
        "formula": [complex_formula],
        "state": ["aq"],
        "ref1": ["WORM Complicator"],
        "ref2": ["NA"],
        "date": [this_date],
        "model": ["HKF"],
        "E_units": ["cal"],
        "G": [G],
        "H": [H],
        "S": [S],
        "Cp": [CP],
        "V": [V],
        "a1.a": [a1],
        "a2.b": [a2],
        "a3.c": [a3],
        "a4.d": [a4],
        "c1.e": [c1],
        "c2.f": [c2],
        "omega.lambda": [wcon],
        "z.T":[Z],
        "azero":[azero_val],
        "neutral_ion_type":[0],
        "dissrxn":[complex_dissrxn],
        "tag":[""],
        "formula_ox":[complex_formula_ox],
        "category_1":[cat_1],
        "category_2":[""],
    }
    df = pd.DataFrame(data)

    if isinstance(filename, str):
        
        file_exists = os.path.isfile(filename+'.csv')
        
        if file_exists:
            with open(filename+'.csv', 'a') as f:
    
                f_df = pd.read_csv(filename+'.csv')
                
                # check that a monoligand complex without parentheses isn't already
                # in the thermodynamic database. Warn or skip.
                keep = True
                complex_name_no_parentheses = complex_name.replace("(", "").replace(")", "")
    
                if complex_name_no_parentheses in list(thermo_data["name"]):
                    if skip_duplicates:
                        keep = False
                    duplicate_list.append([complex_name, complex_name_no_parentheses])
                
                complex_name_organic_abbv = ""
                for name,abbv in ligand_name_abbrv_pairs:

                    if name in complex_name:
                        # replace "Mg(acetate)+" with "Mg(Ac)+", etc. for the sake of
                        # checking for duplicate entries
                        complex_name_organic_abbv = complex_name.replace(name, abbv)
                
                    if complex_name_organic_abbv in list(thermo_data["name"]):
                        if skip_duplicates:
                            keep = False
                        duplicate_list.append([complex_name, complex_name_organic_abbv])

                        break
                
                # if keep and replace and complex_name in list(f_df["name"]):
                #     cols = list(f_df.columns) 
                #     f_df.iloc[f_df["name"]==complex_name, :] = df
    
                #     f_df.to_csv(filename+'.csv', header=True, index=False)
    
                if keep: # was elif prior to the block above being commented out
                    df.to_csv(f, header=False, index=False)
    
        else:
            df.to_csv(filename+'.csv', index=False)
            
    return df, duplicate_list


def find_sigfigs(number_string):
    """
    Count significant figures (cannot handle scientific notation)
    """
    
    if isinstance(number_string, float):
        return ""

    if "-" in number_string:
        # delete negative symbol
        number_string = number_string.replace("-", "")
    
    if "." in number_string:
        # delete decimal point
        number_string = number_string.replace(".", "")
        # delete leading zeros
        number_string = number_string.lstrip("0")
    else:
        # delete leading and trailing zeros
        number_string = number_string.strip("0")

    return len(number_string)



def calc_params(Z, G, H, S, CP, V,
                a1, a2, a3, a4, c1, c2, omega,
                water_model="SUPCRT92", organic=False):

    kwargs = {"V": V, "Cp": CP, "Gf": G, "Hf": H, "Saq": S,
              "a1":a1, "a2":a2, "a3":a3, "a4":a4, "c1":c1, "c2":c2,
              "aq_complex": True, "HKF_scale":False, "Z":Z, "organic":organic}

    if water_model == "DEW":
        kwargs["DEW"] = True

    if isinstance(omega, float):
        kwargs["omega"] = omega
    
    kwargs["phase_TrPr"] = "cr"

    out_dict, eqn_list = wormutils.find_HKF(**kwargs)
    
    G = out_dict["G"]
    H = out_dict["H"]
    S = out_dict["S"]
    a1 = out_dict["a1"]
    a2 = out_dict["a2"]
    a3 = out_dict["a3"]
    a4 = out_dict["a4"]
    c1 = out_dict["c1"]
    c2 = out_dict["c2"]
    wcon = out_dict["omega"]

    return G, H, S, a1, a2, a3, a4, c1, c2, wcon, eqn_list


def complicate(cation=None, ligand=None, beta=None, sass=None, out_name=None,
               G=None, H=None, S=None, V=None, Cp=None, a1=None, a2=None, a3=None, a4=None,
               c1=None, c2=None, omega=None, azero=4, rt=3, sigfigs=False, df_in=None,
               data_path=None, ligand_abbrv_data_path=None, water_model="SUPCRT92",
               correct_basis=True, skip_duplicates=True, print_warnings=True,
               SSH97_cp_eqns=False, focus_on=None, metal_organic_complex=False):
    
    """
    Estimate the thermodynamic properties and Helgeson-Kirkham-Flowers (HKF)
    equation of state parameters for aqueous inorganic complexes.
    
    Parameters
    ----------
    cation : str
        Name of the cation in the aqueous complex.
        
    ligand : str
        Name of the ligand in the aqueous complex.
    
    beta : list of float
        List of four association equilibrium constants for the first, second,
        third, and fourth association, respectively. The first value should
        correspond to the equilibrium constant of a reaction like this:
        Mg+2 + Cl- = MgCl+
        While the second should represent a reaction written like this:
        Mg+2 + 2Cl- = MgCl2
        And so on, up to the fourth association.
       
    sass : list of float, optional
        List of four entropies of association (in cal mol-1 K-1) for the first,
        second, third, and fourth aqueous complex, respectively. If one or more
        values are not provided, they will be estimated from the properties of
        the cation, ligand, and association constants. The first entropy of
        association would be for a reaction like this:
        Mg+2 + Cl- = MgCl+
        While the second would be for a reaction written like this:
        Mg+2 + 2Cl- = MgCl2
        And so on, up to the fourth association.

    out_name : str, optional
        Name of the CSV file in which results should be exported. If out_name is
        not specified, a file will not be written

    G : list of numeric, optional
        List of four partial molal Gibbs free energies of formation (in cal
        mol-1) corresponding to the aqueous complex formed by the first, second,
        third, and fourth association. If values are not provided, they will be
        estimated from the properties of the cation, ligand, and association
        constants.

    H : list of numeric, optional
        List of four partial molal enthalpies of formation (in cal
        mol-1) corresponding to the aqueous complex formed by the first, second,
        third, and fourth association. If values are not provided, they will be
        estimated from the properties of the cation, ligand, and association
        constants.

    Cp : list of float, optional
        List of four partial molal isobaric heat capacities (in cal mol-1 K-1)
        corresponding to the aqueous complex formed by the first, second, third,
        and fourth association. If values are not provided, they will be
        estimated from the properties of the cation and ligand.

    V : list of float, optional
        List of four partial molal volumes (in cm3 mol-1) corresponding to the
        aqueous complex formed by the first, second, third, and fourth
        association. If values are not provided, they will be estimated from the
        properties of the cation and ligand.

    a1, a2, a3, a4, c1, c2, omega : list of float, optional
        Parameters for the revised Helgeson Kirkam Flowers
        equation of state. The user may provide a list of four values
        corresponding to the aqueous complex formed by the first, second, third,
        and fourth association. If values are not provided, they will be
        estimated using published correlation methods.
    
    azero : float, default 4
        Aqueous species hard core diameter, a parameter in the Debeye-Hückel
        B-dot equation (Helgeson 1969) used to calculate activity coefficients.
        Units are in angstroms.
    
    rt : int, default 3
        Round output values to how many decimal places? Ignored if `sigfigs` is
        True.
    
    sigfigs : bool, default False
        Experimental. Round output values to a number of significant figures
        determined by number of significant figures in values used to make the
        estimation?

    df_in : pandas.DataFrame or filename str, optional
        A dataframe containing columns for "Metal", "Ligand", "BETA_1", "BETA_2",
        "BETA_3", "BETA_4", "S_1", "S_2", "S_3", "S_4". If a dataframe is
        provided to this parameter, thermodynamic estimates will be created for
        the metal-ligand complex given on each row.

        Alternatively, the name of a CSV file can be specified that will be
        loaded as a dataframe.
    
    data_path : str, optional
        File path and name of a custom thermodynamic database CSV. If undefined,
        the default WORM database will be used.

    water_model : str, default "SUPCRT92"
        Can be either "SUPCRT92" or "DEW" (Deep Earth Water) model. Influences
        the compatibility of estimated aqueous complex equation of state
        parameters for the chosen water model. If "DEW" is selected, the a1
        parameter is estimated according to Eq 16 in Sverjensky et al 2014
        (https://doi.org/10.1016/j.gca.2013.12.019)
    
    correct_basis : bool, default True
        If a cation or ligand is not a strict or auxiliary basis species in the
        thermodynamic database, ensure the dissociation reaction of the
        resulting complex includes only basis species? The default is True
        because complexes that are sent to aqueous equilibration package AqEquil
        require dissociation reactions into strict or auxiliary basis species.
    
    skip_duplicates : bool, default True
        Detect whether a complex with a monovalent ligand is already in the
        thermodynamic database without parentheses. E.g., if "Ca(HCO3)+" is being
        added, check that a species named "CaHCO3+" isn't already in the database.
        If set to False, the species will be added anyway, but the user will
        be warned. If set to True, the species will not be added.

    print_warnings : bool, default True
        Print Complicator warnings?
    
    SSH97_cp_eqns : bool, default False
        Use equations given in Sverjensky et al. 1997 for estimating heat
        capacities of chloride complexes, or of acetate complexes with divalent
        cations? If False, the Complicator will use a built-in equations to
        estimate complex heat capacities. These equations can be viewed in the
        "equations" output of the complicator.
    
    focus_on : list of str, optional
        Allows the user to provide a list of metal or ligand names. The
        Complicator will only make estimates for complexes that contain those
        metals or ligands.
    
    metal_organic_complex : bool, default False
        Is this a metal-organic complex? If True, estimation assumptions from
        Prapaipong and Shock (2001) will be used.
    
    Returns
    -------
    df_out : pandas dataframe
        A dataframe with estimated properties and HKF parameters of the aqueous
        complex.

    equations : dict
        A dictionary organizing equations that were used to estimate the
        thermodynamic properties and parameters of aqueous complexes.
    
    warnings : list
        A list of warnings that were encountered.
    
    duplicates : list
        A list of estimated complexes that appear to match existing complexes in
        the thermodynamic database, reported as duplicate pairs, e.g.,
        [["Mg(OH)+", "MgOH+"], ...].
    """

    if print_warnings:
        print("The WORM Complicator is in currently in beta testing. "
              "Estimates should be treated with caution until the full release.")

    if isinstance(df_in, str):
        df_in = pd.read_csv(df_in)
    
    if isinstance(df_in, pd.DataFrame):

        def complicator_check_assign(prop, row):
            if prop in row.index:
                if not math.isnan(row[prop]):
                    return row[prop]
                else:
                    return float('NaN')
            else:
                return float('NaN')

        equations = {}
        warnings = []
        duplicates = []
        df_out_list = []
        for index, row in df_in.iterrows():

            if isinstance(focus_on, list):
                if row['Metal'] not in focus_on:
                    continue
                
                if row['Ligand'] not in focus_on:
                    continue
            
            if 'BETA_1' in row.index:
                BETA1 = row['BETA_1']
            else:
                BETA1 = float('NaN')
            if 'BETA_2' in row.index:
                BETA2 = row['BETA_2']
            else:
                BETA2 = float('NaN')
            if 'BETA_3' in row.index:
                BETA3 = row['BETA_3']
            else:
                BETA3 = float('NaN')
            if 'BETA_4' in row.index:
                BETA4 = row['BETA_4']
            else:
                BETA4 = float('NaN')
            if math.isnan(BETA1) and math.isnan(BETA2) and math.isnan(BETA3) and math.isnan(BETA4):
                continue

            if 'Sass_1' in row.index:
                sass1 = row['Sass_1']
            else:
                sass1 = float('NaN')
            if 'Sass_2' in row.index:
                sass2 = row['Sass_2']
            else:
                sass2 = float('NaN')
            if 'Sass_3' in row.index:
                sass3 = row['Sass_3']
            else:
                sass3 = float('NaN')
            if 'Sass_4' in row.index:
                sass4 = row['Sass_4']
            else:
                sass4 = float('NaN')

            G_1 = complicator_check_assign("G_1", row)
            H_1 = complicator_check_assign("H_1", row)
            S_1 = complicator_check_assign("S_1", row)
            V_1 = complicator_check_assign("V_1", row)
            Cp_1 = complicator_check_assign("Cp_1", row)
            a1_1 = complicator_check_assign("a1_1", row)
            a2_1 = complicator_check_assign("a2_1", row)
            a3_1 = complicator_check_assign("a3_1", row)
            a4_1 = complicator_check_assign("a4_1", row)
            c1_1 = complicator_check_assign("c1_1", row)
            c2_1 = complicator_check_assign("c2_1", row)
            omega_1 = complicator_check_assign("omega_1", row)
            G_2 = complicator_check_assign("G_2", row)
            H_2 = complicator_check_assign("H_2", row)
            S_2 = complicator_check_assign("S_2", row)
            V_2 = complicator_check_assign("V_2", row)
            Cp_2 = complicator_check_assign("Cp_2", row)
            a1_2 = complicator_check_assign("a1_2", row)
            a2_2 = complicator_check_assign("a2_2", row)
            a3_2 = complicator_check_assign("a3_2", row)
            a4_2 = complicator_check_assign("a4_2", row)
            c1_2 = complicator_check_assign("c1_2", row)
            c2_2 = complicator_check_assign("c2_2", row)
            omega_2 = complicator_check_assign("omega_2", row)
            G_3 = complicator_check_assign("G_3", row)
            H_3 = complicator_check_assign("H_3", row)
            S_3 = complicator_check_assign("S_3", row)
            V_3 = complicator_check_assign("V_3", row)
            Cp_3 = complicator_check_assign("Cp_3", row)
            a1_3 = complicator_check_assign("a1_3", row)
            a2_3 = complicator_check_assign("a2_3", row)
            a3_3 = complicator_check_assign("a3_3", row)
            a4_3 = complicator_check_assign("a4_3", row)
            c1_3 = complicator_check_assign("c1_3", row)
            c2_3 = complicator_check_assign("c2_3", row)
            omega_3 = complicator_check_assign("omega_3", row)
            G_4 = complicator_check_assign("G_4", row)
            H_4 = complicator_check_assign("H_4", row)
            S_4 = complicator_check_assign("S_4", row)
            V_4 = complicator_check_assign("V_4", row)
            Cp_4 = complicator_check_assign("Cp_4", row)
            a1_4 = complicator_check_assign("a1_4", row)
            a2_4 = complicator_check_assign("a2_4", row)
            a3_4 = complicator_check_assign("a3_4", row)
            a4_4 = complicator_check_assign("a4_4", row)
            c1_4 = complicator_check_assign("c1_4", row)
            c2_4 = complicator_check_assign("c2_4", row)
            omega_4 = complicator_check_assign("omega_4", row)

            a1_1 = a1_1/10
            a2_1 = a2_1/10**-2
            a3_1 = a3_1
            a4_1 = a4_1/10**-4
            c1_1 = c1_1
            c2_1 = c2_1/10**-4
            omega_1 = omega_1/10**-5
            a1_2 = a1_2/10
            a2_2 = a2_2/10**-2
            a3_2 = a3_2
            a4_2 = a4_2/10**-4
            c1_2 = c1_2
            c2_2 = c2_2/10**-4
            omega_2 = omega_2/10**-5
            a1_3 = a1_3/10
            a2_3 = a2_3/10**-2
            a3_3 = a3_3
            a4_3 = a4_3/10**-4
            c1_3 = c1_3
            c2_3 = c2_3/10**-4
            omega_3 = omega_3/10**-5
            a1_4 = a1_4/10
            a2_4 = a2_4/10**-2
            a3_4 = a3_4
            a4_4 = a4_4/10**-4
            c1_4 = c1_4
            c2_4 = c2_4/10**-4
            omega_4 = omega_4/10**-5

            cation = str(row['Metal'])
            ligand = str(row['Ligand'])
            beta = [BETA1, BETA2, BETA3, BETA4]
            G = [G_1, G_2, G_3, G_4]
            H = [H_1, H_2, H_3, H_4]
            S = [S_1, S_2, S_3, S_4]
            sass = [sass1, sass2, sass3, sass4]
            V = [V_1, V_2, V_3, V_4]
            Cp = [Cp_1, Cp_2, Cp_3, Cp_4]
            a1 = [a1_1, a1_2, a1_3, a1_4]
            a2 = [a2_1, a2_2, a2_3, a2_4]
            a3 = [a3_1, a3_2, a3_3, a3_4]
            a4 = [a4_1, a4_2, a4_3, a4_4]
            c1 = [c1_1, c1_2, c1_3, c1_4]
            c2 = [c2_1, c2_2, c2_3, c2_4]
            omega = [omega_1, omega_2, omega_3, omega_4]

            df_out_row, eqn_row, warn_row, dupe_row, = complicate(cation, ligand, beta, sass,
                    out_name, df_in=None, sigfigs=sigfigs, data_path=data_path,
                    skip_duplicates=skip_duplicates,
                    water_model=water_model, SSH97_cp_eqns=SSH97_cp_eqns,
                    G=G, H=H, S=S, V=V, Cp=Cp, a1=a1, a2=a2, a3=a3, a4=a4, c1=c1, c2=c2,
                    omega=omega, ligand_abbrv_data_path=ligand_abbrv_data_path,
                    metal_organic_complex=metal_organic_complex,
                    print_warnings=False)
            
            df_out_list.append(df_out_row)

            # counts first instance of an equation but removes duplicates
            counted_eq = []
            for eqn in eqn_row:
                if eqn not in counted_eq:
                    counted_eq.append(eqn)
            
            equations[cation + " " + ligand] = counted_eq
            warnings += warn_row
            duplicates += dupe_row

        if len(df_out_list) > 0:
            df_out = pd.concat(df_out_list, ignore_index=True)
        else:
            df_out = None

        warnings = list(dict.fromkeys(warnings))
        if print_warnings:
            for warning in warnings:
                print(warning)
        
        return df_out, equations, warnings, duplicates

    warning_list = []
    BETA1, BETA2, BETA3, BETA4 = beta

    # convert association constants (input) into dissociation constants used by
    # the autopp93s4 code.
    BETA1 = -BETA1
    BETA2 = -BETA2
    BETA3 = -BETA3
    BETA4 = -BETA4
    
    sass1, sass2, sass3, sass4 = sass

    G_1, G_2, G_3, G_4 = G
    H_1, H_2, H_3, H_4 = H
    S_1, S_2, S_3, S_4 = S
    Cp_1, Cp_2, Cp_3, Cp_4 = Cp
    V_1, V_2, V_3, V_4 = V
    a1_1, a1_2, a1_3, a1_4 = a1
    a2_1, a2_2, a2_3, a2_4 = a2
    a3_1, a3_2, a3_3, a3_4 = a3
    a4_1, a4_2, a4_3, a4_4 = a4
    c1_1, c1_2, c1_3, c1_4 = c1
    c2_1, c2_2, c2_3, c2_4 = c2
    omega_1, omega_2, omega_3, omega_4 = omega 

    eqs_used = [] # store list of strings describing which equations were used in the estimation
    duplicate_list = [] # store a list of complexes that appear to be duplicates of complexes that already exist in the thermodynamic database
    duplicate_list_1 = []
    duplicate_list_2 = []
    duplicate_list_3 = []
    duplicate_list_4 = []

    df_out_1 = pd.DataFrame()
    df_out_2 = pd.DataFrame()
    df_out_3 = pd.DataFrame()
    df_out_4 = pd.DataFrame()

    if data_path == None:
        thermo_data = wormutils.import_package_file(__name__, "wrm_data.csv")
    else:
        thermo_data = pd.read_csv(data_path)

    if ligand_abbrv_data_path == None:
        ligand_abbrv = wormutils.import_package_file(__name__, "wrm_ligand_abbrv.csv")
    elif isinstance(ligand_abbrv_data_path, str):
        ligand_abbrv = pd.read_csv(ligand_abbrv_data_path)
        
    # a list of ligand names and abbreviations,
    # e.g., [["formate", "For"], ["acetate", "Ac"]], for the sake of
    # checking for duplicate entries in the thermodynamic database
    ligand_name_abbrv_pairs = [[name,abbv] for name,abbv in zip(ligand_abbrv["ligand"], ligand_abbrv["abbrv"])]
    
    cation_entry = thermo_data[thermo_data["name"]==cation]
    ligand_entry = thermo_data[thermo_data["name"]==ligand]
    
    if cation_entry.empty:
        msg = "Could not find an entry for " + cation + " in reference sheet."
        warning_list.append(msg)
        
    if ligand_entry.empty:
        msg = "Could not find an entry for " + ligand + " in reference sheet."
        warning_list.append(msg)

    if cation_entry.shape[0] == 0 or ligand_entry.shape[0] == 0:
        return pd.DataFrame(), eqs_used, warning_list, duplicate_list

    V_or_Cp = True
    if math.isnan(float(cation_entry["V"].values[0])):
        msg = "Cation " + cation + " does not have a partial molal volume (V) in the database. Skipping estimation for " + cation + " " + ligand + " complexes."
        warning_list.append(msg)
        V_or_Cp = False
    if math.isnan(float(ligand_entry["V"].values[0])):
        msg = "Ligand " + ligand + " does not have a partial molal volume (V) in the database. Skipping estimation for " + cation + " " + ligand + " complexes."
        warning_list.append(msg)
        V_or_Cp = False
    if math.isnan(float(cation_entry["Cp"].values[0])):
        msg = "Cation " + cation + " does not have a partial molal isobaric heat capacity (Cp) in the database. Skipping estimation for " + cation + " " + ligand + " complexes."
        warning_list.append(msg)
        V_or_Cp = False
    if math.isnan(float(ligand_entry["Cp"].values[0])):
        msg = "Ligand " + ligand + " does not have a partial molal isobaric heat capacity (Cp) in the database. Skipping estimation for " + cation + " " + ligand + " complexes."
        warning_list.append(msg)
        V_or_Cp = False
    
    if not V_or_Cp:
        return pd.DataFrame(), eqs_used, warning_list, duplicate_list

    if "tag" in list(cation_entry.keys()):
        if cation_entry["tag"].values[0] != "basis" and cation_entry["tag"].values[0] != "aux" and correct_basis:
            cation_dissrxn = cation_entry["dissrxn"].values[0]
            cation_dissrxn_split = cation_dissrxn.split(" ")
            cation_dissrxn_names = cation_dissrxn_split[1::2]
            cation_dissrxn_coeffs = cation_dissrxn_split[::2]
            
            cation_dissrxn_dict_prediv = {name:float(coeff) for name,coeff in zip(cation_dissrxn_names, cation_dissrxn_coeffs)}
            div_val = abs(cation_dissrxn_dict_prediv[cation])
            cation_dissrxn_dict = {name:coeff/div_val for name,coeff in cation_dissrxn_dict_prediv.items()}
            del cation_dissrxn_dict[cation]
            
            cation_dissrxn_dict = cation_dissrxn_dict/abs(cation_dissrxn_dict[cation])
    
            msg = (cation + " is not a strict or auxiliary basis species. "
                  "Dissociation reactions of the complex will assume the "
                  "cation dissociates according to: " + cation_dissrxn)
            warning_list.append(msg)
        else:
            cation_dissrxn_dict = {cation:1}


        if ligand_entry["tag"].values[0] != "basis" and ligand_entry["tag"].values[0] != "aux" and correct_basis:
            ligand_dissrxn = ligand_entry["dissrxn"].values[0]
            ligand_dissrxn_split = ligand_dissrxn.split(" ")
            ligand_dissrxn_names = ligand_dissrxn_split[1::2]
            ligand_dissrxn_coeffs = ligand_dissrxn_split[::2]
            
            ligand_dissrxn_dict_prediv = {name:float(coeff) for name,coeff in zip(ligand_dissrxn_names, ligand_dissrxn_coeffs)}
            div_val = abs(ligand_dissrxn_dict_prediv[ligand])
            ligand_dissrxn_dict = {name:coeff/div_val for name,coeff in ligand_dissrxn_dict_prediv.items()}
            del ligand_dissrxn_dict[ligand]
    
            msg = (ligand + " is not a strict or auxiliary basis species. "
                  "Dissociation reactions of the complex will assume the "
                  "ligand dissociates according to: " + ligand_dissrxn)
            warning_list.append(msg)
        else:
            ligand_dissrxn_dict = {ligand:1}

        # handle complex formula_ox
        cation_formula_ox_dict = __get_formula_ox_dict(cation, df=thermo_data)
        ligand_formula_ox_dict = __get_formula_ox_dict(ligand, df=thermo_data)
    
    else:
        cation_dissrxn_dict = {cation:1}
        ligand_dissrxn_dict = {ligand:1}
        cation_formula_ox_dict = {}
        ligand_formula_ox_dict = {}
    

    
    ZC = str(cation_entry["z.T"].values[0])
    GC = str(cation_entry["G"].values[0])
    HC = str(cation_entry["H"].values[0])
    SC = str(cation_entry["S"].values[0])
    CPC = str(cation_entry["Cp"].values[0])
    VC = str(cation_entry["V"].values[0])

    if cation_entry["E_units"].values[0] == "J":
        GC = GC/4.184
        HC = HC/4.184
        SC = SC/4.184
        CPC = CPC/4.18
    elif cation_entry["E_units"].values[0] == "cal":
        pass
    else:
        msg = "Cation " + cation + " is not identified as being Joule or calorie based in the thermodynamic database (E_units). Skipping this complex..."
        warning_list.append(msg)
        return pd.DataFrame(), eqs_used, warning_list, duplicate_list
    
    if sigfigs:
        GC_sf = find_sigfigs(GC)
        HC_sf = find_sigfigs(HC)
        SC_sf = find_sigfigs(SC)
        CPC_sf = find_sigfigs(CPC)
        VC_sf = find_sigfigs(VC)
    else:
        GC_sf = GC
        HC_sf = HC
        SC_sf = SC
        CPC_sf = CPC
        VC_sf = VC
    
    ZL = str(ligand_entry["z.T"].values[0])
    GL = str(ligand_entry["G"].values[0])
    HL = str(ligand_entry["H"].values[0])
    SL = str(ligand_entry["S"].values[0])
    CPL = str(ligand_entry["Cp"].values[0])
    VL = str(ligand_entry["V"].values[0])

    if ligand_entry["E_units"].values[0] == "J":
        GL = GL/4.184
        HL = HL/4.184
        SL = SL/4.184
        CPL = CPL/4.18
    elif ligand_entry["E_units"].values[0] == "cal":
        pass
    else:
        msg = "Ligand " + ligand + " is not identified as Joule or calorie based in the thermodynamic database (E_units). Skipping this complex..."
        warning_list.append(msg)
        return pd.DataFrame(), eqs_used, warning_list, duplicate_list
    
    GL_sf = find_sigfigs(GL)
    HL_sf = find_sigfigs(HL)
    SL_sf = find_sigfigs(SL)
    CPL_sf = find_sigfigs(CPL)
    VL_sf = find_sigfigs(VL)

    ZC, GC, HC, SC, CPC, VC = int(float(ZC)), float(GC), float(HC), float(SC), float(CPC), float(VC)
    ZL, GL, HL, SL, CPL, VL = int(float(ZL)), float(GL), float(HL), float(SL), float(CPL), float(VL)

    G_water_298K = -56688 # cal/mol, standard Gibbs free energy of formation of water at 298.15K
    H_water_298K = -68317 # cal/mol, standard enthalpy of formation of water at 298.15K
    S_water_298K = 16.712 # cal/K/mol, standard entropy of formation of water at 298.15K
    V_water_298K = 18.07 # cm3/mol, standard volume of formation of water at 298.15K

    first_row_transition_metals = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu"]
    transition_metals = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zr", "Nb", "Mo",
                         "Tc", "Ru", "Rh", "Pd", "Ag", "Hf", "Ta", "W", "Re", "Os", "Ir",
                         "Pt", "Au", "Hg"]

    alkaline_earths = ["Be", "Mg", "Ca", "Sr", "Ba", "Ra"]

    if not math.isnan(BETA1):
        
        ###### Calculations for the first complex
        eqs_used.append("Beginning calculations for the first complex...")
        
        Z = ZC + ZL
        Z1 = Z
        
        if not math.isnan(G_1):
            G1 = G_1
            DELGR1 = -G1 + GC + GL
            eqs_used.append("G1 = {} cal/mol, standard Gibbs free energy of formation of the first complex".format("{0:.5g}".format(G1)))
            eqs_used.append("DELGR1 = {} cal/mol = -G1 + GC + GL, Gibbs free energy of association of the first complex".format("{0:.5g}".format(DELGR1)))
        else:
            DELGR1 = (2.30259)*(1.98719)*(298.15)*BETA1
            G1 = DELGR1 + GC + GL
            eqs_used.append("DELGR1 = {} cal/mol = (2.30259)*(1.98719)*(298.15)*BETA1, Gibbs free energy of association of the first complex".format("{0:.5g}".format(DELGR1)))
            eqs_used.append("Standard Gibbs free energy of formation of the cation GC = {} cal/mol, and the ligand GL = {} cal/mol".format(GC, GL))
            eqs_used.append("G1 = {} cal/mol = DELGR1 + GC + GL, standard state partial molal Gibbs free energy of formation of the first complex".format("{0:.5g}".format(G1)))
 
        # retrieve the cation's element (e.g., "Fe" for "Fe+2")
        cation_formula_dict = chemparse.parse_formula(cation_entry["formula"].values[0])
        cation_element = [k for k in cation_formula_dict.keys() if k != "+"][0]
        
        # Now the entropy predictor starts for the first complex
        # if sass1 is available, use it. If not, predict it.
        if not math.isnan(sass1) or not math.isnan(S_1):
            # if user provides an entropy of association of the first complex...
            if math.isnan(S_1):
                # if user does not provide a standard entropy of the first complex...
                DELSR1 = sass1
                DELS1 = DELSR1
                S1 = SC + SL + DELS1
                eqs_used.append("Entropy of association Sass_1 = {}, third law entropy of the cation SC = {}, and ligand SL = {}, cal/mol/K".format(sass1, SC, SL))
                eqs_used.append("DELS1 = DELSR1 = {} cal/mol/K, entropy of association of the first complex".format("{0:.5g}".format(DELS1)))
                eqs_used.append("S1 = {} cal/mol/K = SC + SL + DELS1, standard state partial molal third law entropy of the first complex".format("{0:.5g}".format(S1)))
            else:
                # if user provides a third law entropy
                S1 = S_1
                DELSR1 = S1 - SC - SL
                DELS1 = DELSR1
                eqs_used.append("Third law entropy of the complex S1 = {}, cation SC = {}, and ligand SL = {}, cal/mol/K".format(S_1, SC, SL))
                eqs_used.append("DELS1 = DELSR1 = {} cal/mol/K = S1 - SC - SL, entropy of association of the first complex".format("{0:.5g}".format(DELS1))) 
        else:
            # if user provides neither a standard entropy of the first complex, nor an entropy of association...
            if ligand == "OH-":
                # correlation from Shock et al., 1997 for hydroxide complexes (table 9)
                eqs_used.append("Beginning standard entropy estimation for the first hydroxide complex according to Shock et al., 1997...")
                if ZC == 1:
                    S1 = 1.32*SC - 6
                    eqs_used.append("S1 = {} cal/K/mol = 1.32*SC + 6, Eqn for monovalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S1)))
                elif ZC == 2:
                    S1 = 1.32*SC + 24.5
                    eqs_used.append("S1 = {} cal/K/mol = 1.32*SC + 24.5, Eqn for divalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S1)))
                elif ZC == 3:
                    if cation_element in ["Ga", "In", "Tl", "Bi"]:
                        S1 = 1.32*SC + 37
                        eqs_used.append("S1 = {} cal/K/mol = 1.32*SC + 37, Eqn for Ga+3, In+3, Tl+3, Bi+3 from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S1)))
                    else:
                        S1 = 1.32*SC + 62
                        eqs_used.append("S1 = {} cal/K/mol = 1.32*SC + 62, Eqn for trivalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S1)))
                elif ZC == 4:
                    S1 = 1.32*SC + 37
                    eqs_used.append("S1 = {} cal/K/mol = 1.32*SC + 37, Eqn for tetravalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S1)))
                else:
                    # cation charges of 5+ are not supported
                    msg = "Estimations for hydroxide complexes for cations with 5+ charge are not supported. Skipping estimation for " + cation + " complexation with OH-..."
                    warning_list.append(msg)
                    return pd.DataFrame(), eqs_used, warning_list, duplicate_list
                DELSR1 = S1 - SC - SL
                DELS1 = DELSR1
                eqs_used.append("DELS1 = DELSR1 = {} cal/mol/K = S1 - SC - SL, entropy of association of the first complex".format("{0:.5g}".format(DELS1)))
            
            else:
                # Eq 74-77 in Sverjensky et al 1997
                AZ= 0.016241*Z - 0.000479
                AZP= -0.36097*Z + 0.3209
                BZ= 0.32102*Z  - 0.05996
                BZP= 8.2198*Z - 1.557
                ALPHA = AZ*(SL+(-5.0*ZL)) + AZP
                BETA = BZ*(SL+(-5.0*ZL)) + BZP
                DELS1 = ALPHA*(SC+(-5.0*ZC)) + BETA
                DELSR1= DELS1
                S1 = DELSR1 + SC + SL
                eqs_used.append("Calculating parameters to estimate the standard third law entropy of the first complex...")
                eqs_used.append("AZ = {} = 0.016241*Z - 0.000479, Eqn 74 in Sverjensky et al. 1997".format("{0:.5g}".format(AZ)))
                eqs_used.append("AZP = {} = -0.36097*Z + 0.3209, Eqn 75 in Sverjensky et al. 1997".format("{0:.5g}".format(AZP)))
                eqs_used.append("BZ = {} = 0.32102*Z  - 0.05996, Eqn 76 in Sverjensky et al. 1997".format("{0:.5g}".format(BZ)))
                eqs_used.append("BZP = {} = 8.2198*Z - 1.557, Eqn 77 in Sverjensky et al. 1997".format("{0:.5g}".format(BZP)))
                eqs_used.append("ALPHA = {} = AZ*(SL+(-5.0*ZL)) + AZP, Eqn pertaining to Sverjensky et al. 1997".format("{0:.5g}".format(ALPHA)))
                eqs_used.append("BETA = {} = BZ*(SL+(-5.0*ZL)) + BZP, Eqn pertaining to Sverjensky et al. 1997".format("{0:.5g}".format(BETA)))
                eqs_used.append("DELS1 = {} = ALPHA*(SC+(-5.0*ZC)) + BETA, Eqn pertaining to Sverjensky et al. 1997".format("{0:.5g}".format(DELS1)))
                eqs_used.append("DELSR1 = DELS1 = {} cal/mol/K = S1 - SC - SL, entropy of association of the first complex".format("{0:.5g}".format(DELS1)))

        if not math.isnan(H_1):
            H1 = H_1
            eqs_used.append("H1 = {} cal/mol, standard enthalpy of the first complex".format("{0:.5g}".format(H1)))
        else:
            DELHR1 = DELGR1 + (298.15*DELSR1)
            H1 = DELHR1 + HC + (1.0*HL)
            eqs_used.append("DELHR1 = {} cal/mol = DELGR1 + (298.15*DELSR1), enthalpy of the first association".format("{0:.5g}".format(DELHR1)))
            eqs_used.append("H1 = {} cal/mol = DELHR1 + HC + (1.0*HL), standard enthalpy of the first complex".format("{0:.5g}".format(H1)))

        if not math.isnan(Cp_1):
            # if user provides a heat capacity
            CP1 = Cp_1
            DELCP1 = CP1 - CPC - CPL
            eqs_used.append("Isobaric heat capacity of the complex CP1 = {}, cation CPC = {}, and ligand CPL = {}, cal/mol/K".format(Cp_1, CPC, CPL))
            eqs_used.append("DELCP1 = {} cal/mol/K = CP1 - CPC - CPL, heat capacity of association of the first complex".format("{0:.5g}".format(DELCP1)))
        else:
            if ligand == "OH-":
                if ZC in [1, 2]:
                    CP1 = -1.14*S1 + 9
                    eqs_used.append("CP1 = {} cal/mol/K = -1.14*S1 + 9, Eqn for monovalent and divalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(CP1)))
                elif ZC in [3, 4]:
                    CP1 = -1.14*S1 - 37.2
                    eqs_used.append("CP1 = {} cal/mol/K = -1.14*S1 - 37.2, Eqn for trivalent and tetravalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(CP1)))
                else:
                    # cations with >4 charge have been weeded out in entropy calculation already
                    CP1 = float('NaN')
                DELCP1 = CP1 - CPC - CPL
            elif metal_organic_complex:
                DELCP1 = 80 # Prapaipong and Shock 2001
                CP1 = DELCP1 + CPC + CPL
                eqs_used.append("DELCP1 = {} cal/mol/K = 80, assumption for metal organic complexes made in Prapaipong and Shock 2001".format("{0:.5g}".format(DELCP1)))
                eqs_used.append("CP1 = {} cal/mol/K = DELCP1 + CPC + CPL, standard isobaric heat capacity of the first complex".format("{0:.5g}".format(CP1)))
            else:
                if SSH97_cp_eqns and ligand == "Cl-":
                    # GB modified 2024 to add equations as they appear in Sverjensky et al., 1997:
                    DELCP1 = 1.25*CPC + 45.3*(Z1+1) - 27.3
                    CP1 = DELCP1 + CPC + CPL
                    eqs_used.append("DELCP1 = {} cal/mol/K = 1.25*CPC + 45.3*(Z1+1) - 27.3, Eqn 53 in Sverjensky et al. 1997".format("{0:.5g}".format(DELCP1)))
                    eqs_used.append("CP1 = {} cal/mol/K = DELCP1 + CPC + CPL, standard isobaric heat capacity of the first complex".format("{0:.5g}".format(CP1)))
                elif SSH97_cp_eqns and ligand =="acetate" and ZC == 2:
                    # GB modified 2024 to add equations as they appear in Sverjensky et al., 1997:
                    DELCP1 = 1.25*CPC + +93.8
                    CP1 = DELCP1 + CPC + CPL
                    eqs_used.append("DELCP1 = {} cal/mol/K = 1.25*CPC + 93.8, Eqn 54 in Sverjensky et al. 1997".format("{0:.5g}".format(DELCP1)))
                    eqs_used.append("CP1 = {} cal/mol/K = DELCP1 + CPC + CPL, standard isobaric heat capacity of the first complex".format("{0:.5g}".format(CP1)))
                else:
                    # ES notes from original fortran code: Here the Cp predictor starts for the first complex.
                    # Modified 8 Aug 1992 after latest revision Modified 17 March 1992 to
                    # incorporate changes from last July 1991
                    dz = 0.856*CPL - 2.1 + 45.3*ZC
                    DELCP1 = 1.25*CPC + dz
                    CP1 = DELCP1 + CPC + CPL
                    eqs_used.append("dz = {} = 0.856*CPL - 2.1 + 45.3*ZC, unpublished equation pertaining to Haas et al 1995, Shock et al. 1995, and Sverjensky et al. 1997".format("{0:.5g}".format(dz)))
                    eqs_used.append("DELCP1 = {} cal/mol/K = 1.25*CPC + dz, unpublished equation pertaining to Haas et al 1995, Shock et al. 1995, and Sverjensky et al. 1997".format("{0:.5g}".format(DELCP1)))
                    eqs_used.append("CP1 = {} cal/mol/K = DELCP1 + CPC + CPL, standard isobaric heat capacity of the first complex".format("{0:.5g}".format(CP1)))

        if not math.isnan(V_1):
            # if user provides a volume
            V1 = V_1
            DELVR1 = V1 - VC - VL
            eqs_used.append("Volume of the complex V1 = {}, cation VC = {}, and ligand VL = {}, cal/mol/K".format(V1, VC, VL))
            eqs_used.append("DELVR1 = {} cm3/mol = V1 - VC - VL, volume of association of the first complex".format("{0:.5g}".format(DELVR1)))
        else:
            use_sverjensky_V_methods = True
            if ligand == "OH-":
                use_sverjensky_V_methods = False
                if ZC == 1:
                    use_sverjensky_V_methods = True
                elif ZC == 2:
                    if cation_element in transition_metals:
                        V1 = 0.45*S1 - 12
                        eqs_used.append("V1 = {} cm3/mol = 0.45*S1 - 12, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V1)))
                    elif cation_element in alkaline_earths:
                        V1 = 0.16*S1 + 4.9
                        eqs_used.append("V1 = {} cm3/mol = 0.16*S1 + 4.9, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V1)))
                    else:
                        V1 = float('NaN')
                elif Z == 3:
                    if cation_element in first_row_transition_metals:
                        V1 = 0.45*S1 - 12
                        eqs_used.append("V1 = {} cm3/mol = 0.45*S1 - 12, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V1)))
                    else:
                        V1 = 0.16*S1 + 4.9
                        eqs_used.append("V1 = {} cm3/mol = 0.16*S1 + 4.9, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V1)))
                elif Z == 4:
                    V1 = 0.16*S1 + 4.9
                    eqs_used.append("V1 = {} cm3/mol = 0.16*S1 + 4.9, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V1)))
                else:
                    # this cation element is not supported by existing methods
                    V1 = float("NaN")
    
            if ligand != "OH-" or use_sverjensky_V_methods:
                # Here the Volume predictor starts for the first complex
                DELVR1  = 0.11419*VC + 8.9432 # Eq. 62 in Sverjensky et al 1997
                V1 = DELVR1 + VC + VL
                eqs_used.append("DELVR1 = {} cm3/mol = 0.11419*VC + 8.9432, Eq. 62 in Sverjensky et al., 1997".format("{0:.5g}".format(DELVR1)))
                eqs_used.append("V1 = {} cm3/mol = DELVR1 + VC + VL, standard volume of the first complex".format("{0:.5g}".format(V1)))

    if not math.isnan(BETA2):
        ###### Calculations for the second complex
        eqs_used.append("Beginning calculations for the second complex...")
        Z = Z + ZL
        Z2 = Z


        if not math.isnan(G_2):
            G2 = G_2
            DELGR2 = -G2 + GC + (2.0*GL)
            eqs_used.append("G2 = {} cal/mol, standard Gibbs free energy of formation of the second complex".format("{0:.5g}".format(G2)))
            eqs_used.append("DELGR2 = {} cal/mol = -G2 + GC + (2.0*GL), Gibbs free energy of association of the second complex".format("{0:.5g}".format(DELGR2)))
        else:
            DELGR2 = 2.30259*1.98719*298.15*BETA2
            G2 = DELGR2 + GC + (2.0*GL)
            eqs_used.append("DELGR2 = {} cal/mol = 2.30259*1.98719*298.15*BETA2, Gibbs free energy of association of the second complex".format("{0:.5g}".format(DELGR2)))
            eqs_used.append("G2 = {} cal/mol = DELGR2 + GC + (2.0*GL), standard state partial molal Gibbs free energy of formation of the second complex".format("{0:.5g}".format(G2)))
    
        if ligand == "OH-":
            # subtract the Gibbs free energy of 1 water
            G2 = G2 - G_water_298K
            DELGR2 = G2 - GC - (2.0*GL) + G_water_298K
            eqs_used.append("G2 = {} cal/mol = G2 - G_water_298K, subtract the standard Gibbs free energy of H2O at 298.15K as per the convention outlined in Shock et al., 1997".format("{0:.5g}".format(G2)))
    
        # if sass2 is available, use it. If not, predict it.
        if not math.isnan(sass2) or not math.isnan(S_2):
            if math.isnan(S_2):
                DELSR2 = sass2
                DELS2 = DELSR2 - DELS1
                S2 = SC + 2.0*SL + DELS1 + DELS2
                eqs_used.append("Sass_2 = {} cal/mol/K, entropy of association of the second complex".format(sass2))
                eqs_used.append("DELS2 = {} cal/mol/K = DELSR2 - DELS1, entropy of association of the second complex from the first complex".format("{0:.5g}".format(DELS2)))
                eqs_used.append("S2 = {} cal/mol/K = SC + 2.0*SL + DELS1 + DELS2, standard state partial molal third law entropy of the second complex".format("{0:.5g}".format(S2)))
            else:
                S2 = S_2
                DELSR2 = S2 - SC - (2.0*SL)
                DELS2 = DELSR2 - DELS1
                eqs_used.append("Third law entropy of the second complex S2 = {}, cation SC = {}, and ligand SL = {}, cal/mol/K".format(S_2, SC, SL))
                eqs_used.append("DELSR2 = {} cal/mol/K = S2 - SC - (2.0*SL), entropy of association of the second complex".format("{0:.5g}".format(DELSR2)))
                eqs_used.append("DELS2 = {} cal/mol/K = DELSR2 - DELS1, entropy of association of the second complex from the first complex".format("{0:.5g}".format(DELS2)))
        else:
            if ligand == "OH-":
                # correlation from Shock et al., 1997 for second hydroxide complexes (table 9)
                eqs_used.append("Beginning standard entropy estimation for the second hydroxide complex according to Shock et al., 1997...")
                if ZC == 1:
                    S2 = 1.42*SC - 11
                    eqs_used.append("S2 = {} cal/K/mol = 1.42*SC - 11, Eqn for monovalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S2)))
                elif ZC == 2:
                    S2 = 1.42*SC + 20.5
                    eqs_used.append("S2 = {} cal/K/mol = 1.42*SC + 20.5, Eqn for divalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S2)))
                elif ZC == 3:
                    S2 = 1.42*SC + 83
                    eqs_used.append("S2 = {} cal/K/mol = 1.42*SC + 83, Eqn for trivalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S2)))
                elif ZC == 4:
                    S2 = 1.42*SC + 108
                    eqs_used.append("S2 = {} cal/K/mol = 1.42*SC + 108, Eqn for tetravalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S2)))
                else:
                    # cation charges of 5+ are not supported but have already been weeded out...
                    S2 = float('NaN')
                
                DELSR2 = S2 - SC - (2.0*SL) + S_water_298K
                DELS2 = DELSR2 - DELS1
                eqs_used.append("DELSR2 = {} cal/mol/K = S2 - SC - (2.0*SL) + S_water_298K, entropy of association of the second complex".format("{0:.5g}".format(DELSR2)))
                eqs_used.append("DELS2 = {} cal/mol/K = DELSR2 - DELS1, entropy of association of the second complex from the first complex".format("{0:.5g}".format(DELS2)))
            else:
                # Eq 74-77 in Sverjensky et al 1997
                AZ = 0.016241*Z - 0.000479
                AZP = -0.36097*Z + 0.3209
                BZ = 0.32102*Z  - 0.05996
                BZP = 8.2198*Z - 1.557
                ALPHA = AZ*(SL+(-5.0*ZL)) + AZP
                BETA = BZ*(SL+(-5.0*ZL)) + BZP
                DELS2 = ALPHA*(S1+(-5.0*(ZC+ZL))) + BETA
                DELSR2= DELS1 + DELS2
                S2 = DELSR2 + SC + (2.0*SL)
                eqs_used.append("Calculating parameters to estimate the standard third law entropy of the second complex...")
                eqs_used.append("AZ = {} = 0.016241*Z - 0.000479, Eqn 74 in Sverjensky et al. 1997".format("{0:.5g}".format(AZ)))
                eqs_used.append("AZP = {} = -0.36097*Z + 0.3209, Eqn 75 in Sverjensky et al. 1997".format("{0:.5g}".format(AZP)))
                eqs_used.append("BZ = {} = 0.32102*Z  - 0.05996, Eqn 76 in Sverjensky et al. 1997".format("{0:.5g}".format(BZ)))
                eqs_used.append("BZP = {} = 8.2198*Z - 1.557, Eqn 77 in Sverjensky et al. 1997".format("{0:.5g}".format(BZP)))
                eqs_used.append("ALPHA = {} = AZ*(SL+(-5.0*ZL)) + AZP, Eqn pertaining to Sverjensky et al. 1997".format("{0:.5g}".format(ALPHA)))
                eqs_used.append("BETA = {} = BZ*(SL+(-5.0*ZL)) + BZP, Eqn pertaining to Sverjensky et al. 1997".format("{0:.5g}".format(BETA)))
                eqs_used.append("DELS2 = {} cal/mol/K = ALPHA*(S1+(-5.0*(ZC+ZL))) + BETA, Eqn pertaining to Sverjensky et al. 1997".format("{0:.5g}".format(DELS2)))
                eqs_used.append("DELSR2 = {} cal/mol/K = DELS1 + DELS2, entropy of association for the second complex (cal/mol/K)".format("{0:.5g}".format(DELSR2)))
                eqs_used.append("S2 = {} cal/mol/K = DELSR2 + SC + (2.0*SL), standard third law entropy of the second complex".format("{0:.5g}".format(S2)))

        if not math.isnan(H_2):
            H2 = H_2
            eqs_used.append("H2 = {} cal/mol, standard enthalpy of the second complex".format("{0:.5g}".format(H2)))
        else:
            DELHR2 = DELGR2 + (298.15*DELSR2)
            H2 = DELHR2 + HC + (2.0*HL)
            eqs_used.append("DELHR2 = {} cal/mol = DELGR2 + (298.15*DELSR2), enthalpy of the second association".format("{0:.5g}".format(DELHR2)))
            eqs_used.append("H2 = {} cal/mol = DELHR2 + HC + (2.0*HL), standard enthalpy of the second complex".format("{0:.5g}".format(H2)))
        
        if ligand == "OH-":
            H2 = H2 - H_water_298K
            eqs_used.append("H2 = {} cal/mol = H2 - H_water_298K, subtract the standard enthalpy of H2O at 298.15K as per the convention outlined in Shock et al., 1997".format("{0:.5g}".format(G2)))
        
        if not math.isnan(Cp_2):
            # if user provides a heat capacity
            CP2 = Cp_2
            DELCP2 = CP2 - CPC - CPL
            eqs_used.append("Isobaric heat capacity of the complex CP2 = {}, cation CPC = {}, and ligand CPL = {}, cal/mol/K".format(Cp_2, CPC, CPL))
            eqs_used.append("DELCP2 = {} cal/mol/K = CP2 - CPC - CPL, heat capacity of association of the second complex".format("{0:.5g}".format(DELCP2)))
        else:
            if ligand == "OH-":
                if ZC in [1, 2]:
                    CP2 = -1.14*S2 - 15.5
                    eqs_used.append("CP2 = {} cal/mol/K = -1.14*S2 - 15.5, Eqn for monovalent and divalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(CP2)))
                elif ZC in [3, 4]:
                    CP2 = -1.14*S2 - 60.8
                    eqs_used.append("CP1 = {} cal/mol/K = -1.14*S2 - 60.8, Eqn for trivalent and tetravalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(CP2)))
                else:
                    # cations with >4 charge have been weeded out in entropy calculation already
                    CP2 = float('NaN')
                DELCPR2 = CP2 - CPC - (2.0*CPL)
                DELCP2 = DELCPR2 - DELCP1
            else:
                if SSH97_cp_eqns and ligand == "Cl-":
                    # GB added 2024: Cp equations as they appear in Sverjensky et al., 1997
                    DELCP2 = (0.89*CPC - 4.9)*1 + DELCP1
                    CP2 = DELCP2 + CP1 + CPL
                    eqs_used.append("DELCP2 = {} cal/mol/K = (0.89*CPC - 4.9)*1 + DELCP1, Eqn 58 in Sverjensky et al. 1997".format("{0:.5g}".format(DELCP2)))
                    eqs_used.append("CP2 = {} cal/mol/K = DELCP1 + CPC + CPL, standard isobaric heat capacity of the second complex".format("{0:.5g}".format(CP2)))
                elif SSH97_cp_eqns and ligand == "acetate" and ZC == 2:
                    # GB modified 2024 to add equations as they appear in Sverjensky et al., 1997:
                    DELCP2 = (0.89*CPC + 20.6)*1 + DELCP1
                    CP2 = DELCP2 + CP1 + CPL
                    eqs_used.append("DELCP2 = {} cal/mol/K = (0.89*CPC + 20.6)*1 + DELCP1, Eqn 59 in Sverjensky et al. 1997".format("{0:.5g}".format(DELCP2)))
                    eqs_used.append("CP2 = {} cal/mol/K = DELCP1 + CPC + CPL, standard isobaric heat capacity of the second complex".format("{0:.5g}".format(CP2)))
                else:
                    # ES note from original fortran code: CHANGED 17 MARCH 1992!!!
                    gz = 0.89*CPC + 0.72*CPL + 16.3
                    DELCP2 = DELCP1 + gz
                    CP2 = DELCP2 + CP1 + CPL
                    eqs_used.append("gz = {} = 0.89*CPC + 0.72*CPL + 16.3, unpublished equation pertaining to Haas et al 1995, Shock et al. 1995, and Sverjensky et al. 1997".format("{0:.5g}".format(gz)))
                    eqs_used.append("DELCP2 = {} cal/mol/K = DELCP1 + gz, unpublished equation pertaining to Haas et al 1995, Shock et al. 1995, and Sverjensky et al. 1997".format("{0:.5g}".format(DELCP2)))
                    eqs_used.append("CP2 = {} cal/mol/K = DELCP2 + CPC + CPL, standard isobaric heat capacity of the second complex".format("{0:.5g}".format(CP2)))
                
        if not math.isnan(V_2):
            # if user provides a volume
            V2 = V_2
            DELVR2 = V2 - VC - VL
            eqs_used.append("Volume of the complex V2 = {}, cation VC = {}, and ligand VL = {}, cal/mol/K".format(V2, VC, VL))
            eqs_used.append("DELVR2 = {} cm3/mol = V2 - VC - VL, volume of association of the second complex".format("{0:.5g}".format(DELVR2)))
        else:
            if ligand == "OH-":
                if ZC == 1:
                    V2 = 0.45*S2 - 12
                    eqs_used.append("V2 = {} cm3/mol = 0.45*S2 - 12, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V2)))
                elif ZC == 2:
                    if cation_element in transition_metals:
                        V2 = 0.45*S2 - 12
                        eqs_used.append("V2 = {} cm3/mol = 0.45*S2 - 12, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V2)))
                    elif cation_element in alkaline_earths:
                        V2 = 0.16*S2 + 4.9
                        eqs_used.append("V2 = {} cm3/mol = 0.16*S2 + 4.9, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V2)))
                    else:
                        # cation element is not supported
                        V2 = float('NaN')
                elif ZC == 3:
                    if cation_element in first_row_transition_metals:
                        V2 = 0.45*S2 - 12
                        eqs_used.append("V2 = {} cm3/mol = 0.45*S2 - 12, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V2)))
                    else:
                        V2 = 0.16*S2 + 4.9
                        eqs_used.append("V2 = {} cm3/mol = 0.16*S2 + 4.9, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V2)))
                elif ZC == 4:
                    V2 = 0.16*S2 + 4.9
                    eqs_used.append("V2 = {} cm3/mol = 0.16*S2 + 4.9, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V2)))
                else:
                    # cations with 5+ charge are not supported
                    V2 = float("NaN")
                
                DELVR2 = V2 - VC - (2.0*VL) + V_water_298K
            
            else:
                DELVR2  = 0.11419*V1 + 8.9432 # Eq. 62 in Sverjensky et al 1997
                V2 = DELVR2 + DELVR1 + VC + (2.0*VL)
                eqs_used.append("DELVR2 = {} cm3/mol = 0.11419*V1 + 8.9432, Eq. 62 in Sverjensky et al 1997".format("{0:.5g}".format(DELVR2)))
                eqs_used.append("V2 = {} cm3/mol = DELVR2 + DELVR1 + VC + (2.0*VL), standard volume of the second complex".format("{0:.5g}".format(V2)))

    if not math.isnan(BETA3):
        ###### Calculations for the third complex
        eqs_used.append("Beginning calculations for the third complex...")
        Z = Z + ZL
        Z3 = Z
        
        if not math.isnan(G_3):
            G3 = G_3
            DELGR3 = -G3 + GC + (3.0*GL)
            eqs_used.append("G3 = {} cal/mol, standard Gibbs free energy of formation of the third complex".format("{0:.5g}".format(G3)))
            eqs_used.append("DELGR3 = {} cal/mol = -G3 + GC + (3.0*GL), Gibbs free energy of association of the third complex".format("{0:.5g}".format(DELGR3)))
        else:
            DELGR3 = 2.30259*1.98719*298.15*BETA3
            G3 = DELGR3 + GC + (3.0*GL)
            eqs_used.append("DELGR3 = {} cal/mol = 2.30259*1.98719*298.15*BETA3, Gibbs free energy of association of the third complex".format("{0:.5g}".format(DELGR3)))
            eqs_used.append("G3 = {} cal/mol = DELGR3 + GC + (3.0*GL), standard state partial molal Gibbs free energy of formation of the third complex".format("{0:.5g}".format(G3)))
    
        if ligand == "OH-":
            G3 = G3 - G_water_298K
            DELGR3 = G3 - GC - (3.0*GL) + G_water_298K
            # subtract the Gibbs free energy of 1 water
            eqs_used.append("G3 = {} cal/mol = G3 - G_water_298K, subtract the standard Gibbs free energy of H2O at 298.15K as per the convention outlined in Shock et al., 1997".format("{0:.5g}".format(G3)))
    
        # if sass3 is available, use it. If not, predict it.
        if not math.isnan(sass3) or not math.isnan(S_3):
            if math.isnan(S_3):
                DELSR3 = sass3
                DELS3 = DELSR3 - DELS1 - DELS2
                S3 = SC + 3.0*SL + DELS1 + DELS2 + DELS3
                eqs_used.append("Sass_23= {} cal/mol/K, entropy of association of the third complex".format(sass3))
                eqs_used.append("DELS3 = {} cal/mol/K = DELSR3 - DELS1 - DELS2, entropy of association of the third complex from the second complex".format("{0:.5g}".format(DELS3)))
                eqs_used.append("S3 = {} cal/mol/K = SC + 2.0*SL + DELS1 + DELS2, standard state partial molal third law entropy of the third complex".format("{0:.5g}".format(S3)))
            else:
                S3 = S_3
                DELSR3 = S3 - SC - (3.0*SL)
                DELS3 = DELSR3 - DELS1 - DELS2
                eqs_used.append("Third law entropy of the third complex S3 = {}, cation SC = {}, and ligand SL = {}, cal/mol/K".format(S_3, SC, SL))
                eqs_used.append("DELSR3 = {} cal/mol/K = S3 - SC - (3.0*SL), entropy of association of the third complex".format("{0:.5g}".format(DELSR3)))
                eqs_used.append("DELS3 = {} cal/mol/K = DELSR3 - DELS1 - DELS2, entropy of association of the third complex from the second complex".format("{0:.5g}".format(DELS3)))
        else:
            if ligand == "OH-":
                # correlation from Shock et al., 1997 for third hydroxide complexes (table 9)
                eqs_used.append("Beginning standard entropy estimation for the third hydroxide complex according to Shock et al., 1997...")
                if ZC == 1:
                    # no support for a third hydroxide complex for monovalent cations
                    msg = "The third hydroxide complex for a monovalent cation is not supported by existing estimation methods. Skipping third hydroxide complex of " + cation
                    warning_list.append(msg)
                    S3 = float('NaN')
                elif ZC == 2:
                    S3 = 1.52*SC + 15.5
                    eqs_used.append("S3 = {} cal/K/mol = 1.52*SC + 15.5, Eqn for divalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S3)))
                elif ZC == 3:
                    S3 = 1.52*SC + 123
                    eqs_used.append("S3 = {} cal/K/mol = 1.52*SC + 123, Eqn for trivalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S3)))
                elif ZC == 4:
                    S3 = 1.52*SC + 140
                    eqs_used.append("S3 = {} cal/K/mol = 1.52*SC + 140, Eqn for tetravalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S3)))
                else:
                    # cation charges of 5+ are not supported
                    msg = "Cations with 5+ charge are not supoported by existing estimation methods. Skipping " + cation
                    warning_list.append(msg)
                    S3 = float('NaN')
                
                DELSR3 = S3 - SC - (3.0*SL) + S_water_298K
                DELS3 = DELSR3 - DELSR2 - DELS1
                eqs_used.append("DELSR3 = {} cal/mol/K = S3 - SC - (3.0*SL) + S_water_298K, entropy of association of the third complex".format("{0:.5g}".format(DELSR3)))
                eqs_used.append("DELS3 = {} cal/mol/K = DELSR3 - DELSR2 - DELS1, entropy of association of the third complex from the second complex".format("{0:.5g}".format(DELS3)))
            else:
                # Eq 74-77 in Sverjensky et al 1997
                AZ= 0.016241*Z - 0.000479
                AZP= -0.36097*Z + 0.3209
                BZ= 0.32102*Z  - 0.05996
                BZP= 8.2198*Z - 1.557
                ALPHA = AZ*(SL+(-5.0*ZL)) + AZP
                BETA = BZ*(SL+(-5.0*ZL)) + BZP
                DELS3 = ALPHA*(S2+(-5.0*(ZC+(2*ZL)))) + BETA
                DELSR3= DELS1 + DELS2 + DELS3
                S3 = DELSR3 + SC + (3.0*SL)
                eqs_used.append("Calculating parameters to estimate the standard third law entropy of the third complex...")
                eqs_used.append("AZ = {} = 0.016241*Z - 0.000479, Eqn 74 in Sverjensky et al. 1997".format("{0:.5g}".format(AZ)))
                eqs_used.append("AZP = {} = -0.36097*Z + 0.3209, Eqn 75 in Sverjensky et al. 1997".format("{0:.5g}".format(AZP)))
                eqs_used.append("BZ = {} = 0.32102*Z  - 0.05996, Eqn 76 in Sverjensky et al. 1997".format("{0:.5g}".format(BZ)))
                eqs_used.append("BZP = {} = 8.2198*Z - 1.557, Eqn 77 in Sverjensky et al. 1997".format("{0:.5g}".format(BZP)))
                eqs_used.append("ALPHA = {} = AZ*(SL+(-5.0*ZL)) + AZP, Eqn pertaining to Sverjensky et al. 1997".format("{0:.5g}".format(ALPHA)))
                eqs_used.append("BETA = {} = BZ*(SL+(-5.0*ZL)) + BZP, Eqn pertaining to Sverjensky et al. 1997".format("{0:.5g}".format(BETA)))
                eqs_used.append("DELS3 = {} cal/mol/K = ALPHA*(S2+(-5.0*(ZC+(2*ZL)))) + BETA, Eqn pertaining to Sverjensky et al. 1997".format("{0:.5g}".format(DELS3)))
                eqs_used.append("DELSR3 = {} cal/mol/K = DELS1 + DELS2 + DELS3, entropy of association for the third complex (cal/mol/K)".format("{0:.5g}".format(DELSR3)))
                eqs_used.append("S3 = {} cal/mol/K = DELSR2 + SC + (2.0*SL), standard third law entropy of the third complex".format("{0:.5g}".format(S3)))

        if not math.isnan(H_3):
            H3 = H_3
            eqs_used.append("H3 = {} cal/mol, standard enthalpy of the third complex".format("{0:.5g}".format(H3)))
        else:
            DELHR3 = DELGR3 + (298.15*DELSR3)
            H3 = DELHR3 + HC + (3.0*HL)
            eqs_used.append("DELHR3 = {} cal/mol = DELGR3 + (298.15*DELSR3), enthalpy of the third association".format("{0:.5g}".format(DELHR3)))
            eqs_used.append("H3 = {} cal/mol = DELHR3 + HC + (3.0*HL), standard enthalpy of the third complex".format("{0:.5g}".format(H3)))
        
        if ligand == "OH-":
            H3 = H3 - H_water_298K
            eqs_used.append("H3 = {} cal/mol = H3 + H_water_298K, subtract the standard enthalpy of H2O at 298.15K as per the convention outlined in Shock et al., 1997".format("{0:.5g}".format(H3)))

        if not math.isnan(Cp_3):
            # if user provides a heat capacity
            CP3 = Cp_3
            DELCP3 = CP3 - CPC - CPL
            eqs_used.append("Isobaric heat capacity of the complex CP3 = {}, cation CPC = {}, and ligand CPL = {}, cal/mol/K".format(Cp_3, CPC, CPL))
            eqs_used.append("DELCP3 = {} cal/mol/K = CP3 - CPC - CPL, heat capacity of association of the third complex".format("{0:.5g}".format(DELCP3)))
        else:
            if ligand == "OH-":
                if ZC == 1:
                    # third hydroxide complex of a monovalent cation is not supported
                    CP3 = float('NaN')
                elif ZC in [2,3,4]:
                    CP3 = -2.28*S3 - 24.0
                    eqs_used.append("CP3 = {} cal/mol/K = -2.28*S3 - 24.0, Eqn for di-, tri-, tetravalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(CP2)))
                else:
                    # cations with >4 charge have been weeded out in entropy calculation already
                    CP3 = float('NaN')
                DELCPR3 = CP3 - CPC - (3.0*CPL)
                DELCP3 = DELCPR3 - DELCPR2 - DELCP1
            else:
                if SSH97_cp_eqns and ligand == "Cl-":
                    # GB added 2024
                    DELCP3 = (0.89*CPC - 4.9)*2 + DELCP1
                    CP3 = DELCP3 + CP2 + CPL
                    eqs_used.append("DELCP3 = {} cal/mol/K = (0.89*CPC - 4.9)*2 + DELCP1, Eqn 58 in Sverjensky et al. 1997".format("{0:.5g}".format(DELCP3)))
                    eqs_used.append("CP3 = {} cal/mol/K = DELCP1 + CPC + CPL, standard isobaric heat capacity of the third complex".format("{0:.5g}".format(CP3)))
                elif SSH97_cp_eqns and ligand == "acetate" and ZC == 2:
                    # GB modified 2024 to add equations as they appear in Sverjensky et al., 1997:
                    DELCP3 = (0.89*CPC + 20.6)*2 + DELCP1
                    CP3 = DELCP3 + CP2 + CPL
                    eqs_used.append("DELCP4 = {} cal/mol/K = (0.89*CPC + 20.6)*2 + DELCP1, Eqn 59 in Sverjensky et al. 1997".format("{0:.5g}".format(DELCP3)))
                    eqs_used.append("CP3 = {} cal/mol/K = DELCP1 + CPC + CPL, standard isobaric heat capacity of the third complex".format("{0:.5g}".format(CP3)))
                else:
                    # ES note from original fortran code: CHANGED 17 MARCH 1992!!!
                    gz = 0.89*CPC + 0.72*CPL + 16.3
                    DELCP3 = DELCP2 + gz
                    CP3 = DELCP3 + CP2 + CPL
                    eqs_used.append("gz = {} = 0.89*CPC + 0.72*CPL + 16.3, unpublished equation pertaining to Haas et al 1995, Shock et al. 1995, and Sverjensky et al. 1997".format("{0:.5g}".format(gz)))
                    eqs_used.append("DELCP3 = {} cal/mol/K = DELCP2 + gz, unpublished equation pertaining to Haas et al 1995, Shock et al. 1995, and Sverjensky et al. 1997".format("{0:.5g}".format(DELCP3)))
                    eqs_used.append("CP3 = {} cal/mol/K = DELCP3 + CP2 + CPL, standard isobaric heat capacity of the second complex".format("{0:.5g}".format(CP3)))

        if not math.isnan(V_3):
            # if user provides a volume
            V3 = V_3
            DELVR3 = V3 - VC - VL
            eqs_used.append("Volume of the complex V3 = {}, cation VC = {}, and ligand VL = {}, cal/mol/K".format(V3, VC, VL))
            eqs_used.append("DELVR3 = {} cm3/mol = V3 - VC - VL, volume of association of the third complex".format("{0:.5g}".format(DELVR3)))
        else:
            if ligand == "OH-":
                if ZC == 1:
                    # third complex of a monovalent cation is not supported
                    V3 = float("NaN")
                elif ZC == 2:
                    if cation_element in transition_metals:
                        V3 = 0.54*S3 - 4.8
                        eqs_used.append("V3 = {} cm3/mol = 0.54*S3 - 4.8, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V3)))
                    elif cation_element in alkaline_earths:
                        V3 = 0.25*S3 + 11.7
                        eqs_used.append("V3 = {} cm3/mol = 0.25*S3 + 11.7, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V3)))
                    else:
                        # cation element is not supported.
                        V3 = float('NaN')
                elif ZC == 3:
                    if cation_element in first_row_transition_metals:
                        V3 = 0.54*S3 - 4.8
                        eqs_used.append("V3 = {} cm3/mol = 0.54*S3 - 4.8, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V3)))
                    else:
                        V3 = 0.25*S3 + 11.7
                        eqs_used.append("V3 = {} cm3/mol = 0.25*S3 + 11.7, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V3)))
                elif ZC == 4:
                    V3 = 0.25*S3 + 11.7
                    eqs_used.append("V3 = {} cm3/mol = 0.25*S3 + 11.7, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V3)))
                else:
                    # cations with 5+ charge are not supported
                    V3 = float('NaN')
                    
                DELVR3 = V3 - VC - (3.0*VL) + V_water_298K
            else:
                DELVR3  = 0.11419*V2 + 8.9432
                V3 = DELVR3 + DELVR2 + DELVR1 + VC + (3.0*VL)
                eqs_used.append("DELVR3 = {} cm3/mol = 0.11419*V2+ 8.9432, Eq. 62 in Sverjensky et al 1997".format("{0:.5g}".format(DELVR3)))
                eqs_used.append("V3 = {} cm3/mol = DELVR3 + DELVR2 + DELVR1 + VC + (3.0*VL), standard volume of the third complex".format("{0:.5g}".format(V3)))

    if not math.isnan(BETA4):
        ##### Calculations for the fourth complex
        eqs_used.append("Beginning calculations for the fourth complex...")
        Z = Z + ZL
        Z4 = Z

        if not math.isnan(G_4):
            G4 = G_4
            DELGR4 = -G4 + GC + (4.0*GL)
            eqs_used.append("G4 = {} cal/mol, standard Gibbs free energy of formation of the fourth complex".format("{0:.5g}".format(G4)))
            eqs_used.append("DELGR4 = {} cal/mol = -G4 + GC + (4.0*GL), Gibbs free energy of association of the fourth complex".format("{0:.5g}".format(DELGR4)))
        else:
            DELGR4 = 2.30259*1.98719*298.15*BETA4
            G4 = DELGR4 + GC + (4.0*GL)
            eqs_used.append("DELGR4 = {} cal/mol = 2.30259*1.98719*298.15*BETA4, Gibbs free energy of association of the fourth complex".format("{0:.5g}".format(DELGR4)))
            eqs_used.append("G4 = {} cal/mol = DELGR4 + GC + (4.0*GL), standard state partial molal Gibbs free energy of formation of the fourth complex".format("{0:.5g}".format(G4)))
    
        if ligand == "OH-":
            G4 = G4 - 2*G_water_298K
            DELGR4 = G4 - GC - (4.0*GL) + 2*G_water_298K
            # subtract the Gibbs free energy of 2 waters
            eqs_used.append("G4 = {} cal/mol = G4 - 2*G_water_298K, subtract the standard Gibbs free energy of two H2O at 298.15K as per the convention outlined in Shock et al., 1997".format("{0:.5g}".format(G4)))
    
        
        # if sass4 is available, use it. If not, predict it.
        if not math.isnan(sass4) or not math.isnan(S_4):
            if math.isnan(S_4):
                DELSR4 = sass4
                DELS4 = DELSR4 - DELS1 - DELS2 - DELS3
                S4 = SC + 4.0*SL + DELS1 + DELS2 + DELS3 + DELS4
                eqs_used.append("Sass_4 = {} cal/mol/K, entropy of association of the third complex".format(sass4))
                eqs_used.append("DELS4 = {} cal/mol/K = DELSR3 - DELS1 - DELS2, entropy of association of the fourth complex from the third complex".format("{0:.5g}".format(DELS4)))
                eqs_used.append("S4 = {} cal/mol/K = SC + 4.0*SL + DELS1 + DELS2 + DELS3 + DELS4, standard state partial molal third law entropy of the third complex".format("{0:.5g}".format(S4)))
    
            else:
                S4 = S_4
                DELSR4 = S4 - SC - (4.0*SL)
                DELS4 = DELSR4 - DELS1 - DELS2 - DELS3
                eqs_used.append("Third law entropy of the fourth complex S4 = {}, cation SC = {}, and ligand SL = {}, cal/mol/K".format(S_4, SC, SL))
                eqs_used.append("DELSR4 = {} cal/mol/K = S4 - SC - (4.0*SL), entropy of association of the third complex".format("{0:.5g}".format(DELSR4)))
                eqs_used.append("DELS4 = {} cal/mol/K = DELSR3 - DELS1 - DELS2, entropy of association of the fourth complex from the third complex".format("{0:.5g}".format(DELS4)))
            
        else:
            if ligand == "OH-":
                # correlation from Shock et al., 1997 for fourth hydroxide complexes (table 9)
                eqs_used.append("Beginning standard entropy estimation for the fourth hydroxide complex according to Shock et al., 1997...")
                if ZC == 1:
                    # no support for a fourth complex for monovalent cations
                    # todo: exit out here
                    S4 = float('NaN')
                elif ZC == 2:
                    S4 = 1.62*SC + 11
                    eqs_used.append("S4 = {} cal/K/mol = 1.62*SC + 11, Eqn for divalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S4)))
                elif ZC == 3:
                    S4 = 1.62*SC + 118
                    eqs_used.append("S4 = {} cal/K/mol = 1.62*SC + 118, Eqn for trivalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S4)))
                elif ZC == 4:
                    S4 = 1.62*SC + 135
                    eqs_used.append("S4 = {} cal/K/mol = 1.62*SC + 135, Eqn for tetravalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(S4)))
                else:
                    # cation charges of 5+ are not supported but have already been weeded out...
                    S4 = float('NaN')
                
                DELSR4 = S4 - SC - (4.0*SL) + 2*S_water_298K
                DELS4 = DELSR4 - DELSR3 - DELSR2 - DELS1
                eqs_used.append("DELSR4 = {} cal/mol/K = S4 - SC - (4.0*SL) + 2*S_water_298K, entropy of association of the fourth complex".format("{0:.5g}".format(DELSR4)))
                eqs_used.append("DELS3 = {} cal/mol/K = DELSR4 - DELSR3 - DELSR2 - DELS1, entropy of association of the fourth complex from the third complex".format("{0:.5g}".format(DELS4)))
            else:
                # Eq 74-77 in Sverjensky et al 1997
                AZ= 0.016241*Z - 0.000479
                AZP= -0.36097*Z + 0.3209
                BZ= 0.32102*Z  - 0.05996
                BZP= 8.2198*Z - 1.557
                ALPHA = AZ*(SL+(-5.0*ZL)) + AZP
                BETA = BZ*(SL+(-5.0*ZL)) + BZP
                DELS4 = ALPHA*(S3+(-5.0*(ZC+(3*ZL)))) + BETA
                DELSR4= DELS1 + DELS2 + DELS3 + DELS4
                S4 = DELSR4 + SC + (4.0*SL)
                eqs_used.append("Calculating parameters to estimate the standard third law entropy of the fourth complex...")
                eqs_used.append("AZ = {} = 0.016241*Z - 0.000479, Eqn 74 in Sverjensky et al. 1997".format("{0:.5g}".format(AZ)))
                eqs_used.append("AZP = {} = -0.36097*Z + 0.3209, Eqn 75 in Sverjensky et al. 1997".format("{0:.5g}".format(AZP)))
                eqs_used.append("BZ = {} = 0.32102*Z  - 0.05996, Eqn 76 in Sverjensky et al. 1997".format("{0:.5g}".format(BZ)))
                eqs_used.append("BZP = {} = 8.2198*Z - 1.557, Eqn 77 in Sverjensky et al. 1997".format("{0:.5g}".format(BZP)))
                eqs_used.append("ALPHA = {} = AZ*(SL+(-5.0*ZL)) + AZP, Eqn pertaining to Sverjensky et al. 1997".format("{0:.5g}".format(ALPHA)))
                eqs_used.append("BETA = {} = BZ*(SL+(-5.0*ZL)) + BZP, Eqn pertaining to Sverjensky et al. 1997".format("{0:.5g}".format(BETA)))
                eqs_used.append("DELS4 = {} cal/mol/K = ALPHA*(S3+(-5.0*(ZC+(3*ZL)))) + BETA, Eqn pertaining to Sverjensky et al. 1997".format("{0:.5g}".format(DELS4)))
                eqs_used.append("DELSR4 = {} cal/mol/K = DELS1 + DELS2 + DELS3 + DELS4, entropy of association for the fourth complex (cal/mol/K)".format("{0:.5g}".format(DELSR4)))
                eqs_used.append("S2 = {} cal/mol/K = DELSR2 + SC + (2.0*SL), standard third law entropy of the fourth complex".format("{0:.5g}".format(S4)))

        if not math.isnan(H_4):
            H4 = H_4
            eqs_used.append("H4 = {} cal/mol, standard enthalpy of the fourth complex".format("{0:.5g}".format(H4)))
        else:
            DELHR4 = DELGR4 + (298.15*DELSR4)
            H4 = DELHR4 + HC + (4.0*HL)
            eqs_used.append("DELHR4 = {} cal/mol = DELGR4 + (298.15*DELSR4), enthalpy of the fourth association".format("{0:.5g}".format(DELHR4)))
            eqs_used.append("H4 = {} cal/mol = DELHR4 + HC + (4.0*HL), standard enthalpy of the fourth complex".format("{0:.5g}".format(H4)))
        
        if ligand == "OH-":
            H4 = H4 - 2*H_water_298K
            eqs_used.append("H4 = {} cal/mol = H4 - 2*H_water_298K, subtract the standard enthalpy of two H2O at 298.15K as per the convention outlined in Shock et al., 1997".format("{0:.5g}".format(H4)))

        if not math.isnan(Cp_4):
            # if user provides a heat capacity
            CP4 = Cp_4
            DELCP4 = CP4 - CPC - CPL
            eqs_used.append("Isobaric heat capacity of the complex CP4 = {}, cation CPC = {}, and ligand CPL = {}, cal/mol/K".format(Cp_4, CPC, CPL))
            eqs_used.append("DELCP4 = {} cal/mol/K = CP4 - CPC - CPL, heat capacity of association of the fourth complex".format("{0:.5g}".format(DELCP4)))
        else:
            if ligand == "OH-":
                if ZC == 1:
                    # fourth hydroxide complex of a monovalent cation is not supported
                    CP4 = float('NaN')
                elif ZC == 2:
                    CP4 = -2.28*S4 - 106.2
                    eqs_used.append("CP4 = {} cal/mol/K = 2.28*S4 - 106.2, Eqn for divalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(CP4)))
                elif ZC in [3,4]:
                    CP4 = -2.06*S4 - 34.5
                    eqs_used.append("CP4 = {} cal/mol/K = -2.06*S4 - 34.5, Eqn for tri- and tetravalent cations from Table 9 of Shock et al., 1997".format("{0:.5g}".format(CP4)))
                else:
                    # cations with >4 charge have been weeded out in entropy calculation already
                    CP4 = float('NaN')
                DELCPR4 = CP4 - CPC - (4.0*CPL)
                DELCP4 = DELCPR4 - DELCPR3 - DELCPR2 - DELCP1
            else:
                if SSH97_cp_eqns and ligand == "Cl-":
                    # GB added 2024
                    DELCP4 = (0.89*CPC - 4.9)*3 + DELCP1
                    CP4 = DELCP4 + CP3 + CPL
                    eqs_used.append("DELCP4 = {} cal/mol/K = (0.89*CPC - 4.9)*3 + DELCP1, Eqn 58 in Sverjensky et al. 1997".format("{0:.5g}".format(DELCP4)))
                    eqs_used.append("CP4 = {} cal/mol/K = DELCP1 + CPC + CPL, standard isobaric heat capacity of the fourth complex".format("{0:.5g}".format(CP4)))
                elif SSH97_cp_eqns and ligand == "acetate" and ZC == 2:
                    # GB modified 2024 to add equations as they appear in Sverjensky et al., 1997:
                    DELCP4 = (0.89*CPC + 20.6)*3 + DELCP1
                    CP4 = DELCP4 + CP3 + CPL
                    eqs_used.append("DELCP4 = (0.89*CPC + 20.6)*3 + DELCP1, Eqn 59 in Sverjensky et al. 1997".format("{0:.5g}".format(DELCP4)))
                    eqs_used.append("CP4 = {} cal/mol/K = DELCP4 + CP3 + CPL, standard isobaric heat capacity of the fourth complex".format("{0:.5g}".format(CP4)))
                else:
                    # ES note from original fortran code: CHANGED 17 MARCH 1992!!!
                    gz = 0.89*CPC + 0.72*CPL + 16.3
                    DELCP4 = DELCP3 + gz
                    CP4 = DELCP4 + CP3 + CPL
                    eqs_used.append("gz = {} = 0.89*CPC + 0.72*CPL + 16.3, unpublished equation pertaining to Haas et al 1995, Shock et al. 1995, and Sverjensky et al. 1997".format("{0:.5g}".format(gz)))
                    eqs_used.append("DELCP4 = {} cal/mol/K = DELCP3 + gz, unpublished equation pertaining to Haas et al 1995, Shock et al. 1995, and Sverjensky et al. 1997".format("{0:.5g}".format(DELCP4)))
                    eqs_used.append("CP4 = {} cal/mol/K = DELCP4 + CP3 + CPL, standard isobaric heat capacity of the fourth complex".format("{0:.5g}".format(CP4)))
        
        if not math.isnan(V_4):
            # if user provides a volume
            V4 = V_4
            DELVR4 = V4 - VC - VL
            eqs_used.append("Volume of the complex V4 = {}, cation VC = {}, and ligand VL = {}, cal/mol/K".format(V4, VC, VL))
            eqs_used.append("DELVR4 = {} cm3/mol = V4 - VC - VL, volume of association of the fourth complex".format("{0:.5g}".format(DELVR4)))
        else:
            if ligand == "OH-":
                if ZC == 1:
                    # fourth complex of a monovalent cation is not supported
                    V4 = float('NaN')
                elif ZC == 2:
                    if cation_element in transition_metals:
                        V4 = 0.54*S4 - 4.8
                        eqs_used.append("V4 = {} cm3/mol = 0.54*S4 - 4.8, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V4)))
                    elif cation_element in alkaline_earths:
                        V4 = 0.25*S4 + 11.7
                        eqs_used.append("V4 = {} cm3/mol = 0.25*S4 + 11.7, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V4)))
                    else:
                        # cation element is not supported.
                        V4 = float('NaN')
                elif ZC == 3:
                    if cation_element in first_row_transition_metals:
                        V4 = 0.54*S4 - 4.8
                        eqs_used.append("V4 = {} cm3/mol = 0.54*S4 - 4.8, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V4)))
                    else:
                        V4 = 0.25*S4 + 11.7
                        eqs_used.append("V4 = {} cm3/mol = 0.25*S4 + 11.7, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V4)))
                elif ZC == 4:
                    V4 = 0.25*S4 + 11.7
                    eqs_used.append("V4 = {} cm3/mol = 0.25*S4 + 11.7, see Table 9 in Shock et al., 1997".format("{0:.5g}".format(V4)))
                else:
                    # cations with 5+ charge are not supported
                    V4 = float('NaN')
                    
                DELVR4 = V4 - VC - (4.0*VL) + 2*V_water_298K
            else:
                DELVR4  = 0.11419*V3 + 8.9432 # Eq. 62 in Sverjensky et al 1997
                V4 = DELVR4 + DELVR3 + DELVR2 + DELVR1 + VC + (4.0*VL)
                eqs_used.append("DELVR4 = {} cm3/mol = 0.11419*V3 + 8.9432, Eq. 62 in Sverjensky et al 1997".format("{0:.5g}".format(DELVR4)))
                eqs_used.append("V4 = {} cm3/mol = DELVR4 + DELVR3 + DELVR2 + DELVR1 + VC + (4.0*VL), standard volume of the fourth complex".format("{0:.5g}".format(V4)))
        
        # override properties and parameters if they are supplied by the user
        if not math.isnan(G_1):
            G1 = G_1
        if not math.isnan(H_1):
            H1 = H_1
        if not math.isnan(S_1):
            S1 = S_1
        if not math.isnan(Cp_1):
            CP1 = Cp_1
        if not math.isnan(V_1):
            V1 = V_1
        
    if metal_organic_complex:
        organic=True
    else:
        organic=False
    
    if not math.isnan(BETA1):

        try:
            eqs_used.append("Estimating HKF parameters for the first complex...")
            G, H, S, a1, a2, a3, a4, c1, c2, wcon, hkf_eqns_1 = calc_params(Z1, G1, H1, S1, CP1, V1, a1_1, a2_1, a3_1, a4_1, c1_1, c2_1, omega_1, water_model=water_model, organic=organic)
            eqs_used.append("HKF parameters for the first complex were estimated...")
        except:
            eqs_used.append("HKF parameters for the first complex could not be estimated...")
            return pd.DataFrame(), eqs_used, warning_list, duplicate_list

        if not math.isnan(G):
            if sigfigs:
                G_out = sf(G, min(GC_sf, GL_sf))
                H_out = sf(H, min(HC_sf, HL_sf))
                S_out = sf(S, min(SC_sf, SL_sf))
                CP_out = sf(CP1, min(CPC_sf, CPL_sf))
                V_out = sf(V1, min(VC_sf, VL_sf))
                a1_out = sf(a1*10, min(VC_sf, VL_sf))
                a2_out = sf(a2/100, min(VC_sf, VL_sf))
                a3_out = sf(a3, min(VC_sf, VL_sf))
                a4_out = sf(a4/10000, min(VC_sf, VL_sf))
                c1_out = sf(c1, min(CPC_sf, CPL_sf))
                c2_out = sf(c2/10000, min(CPC_sf, CPL_sf))
                wcon_out = sf(wcon/100000, min(SC_sf, SL_sf))
            else:
                G_out = round(G, rt)
                H_out = round(H, rt)
                S_out = round(S, rt)
                CP_out = round(CP1, rt)
                V_out = round(V1, rt)
                a1_out = round(a1*10, rt)
                a2_out = round(a2/100, rt)
                a3_out = round(a3, rt)
                a4_out = round(a4/10000, rt)
                c1_out = round(c1, rt)
                c2_out = round(c2/10000, rt)
                wcon_out = round(wcon/100000, rt)
            Z_out = Z1

            df_out_1, duplicate_list_1 = __write_output(filename=out_name,
                         cation=cation, ligand=ligand, nth_complex=1, G=G_out, H=H_out, S=S_out,
                         CP=CP_out, V=V_out, a1=a1_out, a2=a2_out, a3=a3_out, a4=a4_out, c1=c1_out, c2=c2_out,
                         wcon=wcon_out, Z=Z_out, azero=azero,
                         cation_dissrxn_dict=cation_dissrxn_dict,
                         ligand_dissrxn_dict=ligand_dissrxn_dict,
                         thermo_data=thermo_data,
                         cation_formula_ox_dict=cation_formula_ox_dict,
                         ligand_formula_ox_dict=ligand_formula_ox_dict,
                         skip_duplicates=skip_duplicates,
                         ligand_name_abbrv_pairs=ligand_name_abbrv_pairs,
                         organic=organic)

        eqs_used = eqs_used + hkf_eqns_1
        duplicate_list = duplicate_list + duplicate_list_1
    
    # override properties and parameters if they are supplied by the user
    if not math.isnan(G_2):
        G2 = G_2
    if not math.isnan(H_2):
        H2 = H_2
    if not math.isnan(S_2):
        S2 = S_2
    if not math.isnan(Cp_2):
        CP2 = Cp_2
    if not math.isnan(V_2):
        V2 = V_2
    
    if not math.isnan(BETA2):
        try:
            eqs_used.append("Estimating HKF parameters for the second complex...")
            G, H, S, a1, a2, a3, a4, c1, c2, wcon, hkf_eqns_2 = calc_params(Z2, G2, H2, S2, CP2, V2, a1_2, a2_2, a3_2, a4_2, c1_2, c2_2, omega_2, water_model=water_model, organic=organic)
            eqs_used.append("HKF parameters for the second complex were estimated...")
        except:
            eqs_used.append("HKF parameters for the second complex could not be estimated...")
            return pd.DataFrame(), eqs_used, warning_list, duplicate_list

        if not math.isnan(G):
            if sigfigs:
                G_out = sf(G, min(GC_sf, GL_sf))
                H_out = sf(H, min(HC_sf, HL_sf))
                S_out = sf(S, min(SC_sf, SL_sf))
                CP_out = sf(CP1, min(CPC_sf, CPL_sf))
                V_out = sf(V1, min(VC_sf, VL_sf))
                a1_out = sf(a1*10, min(VC_sf, VL_sf))
                a2_out = sf(a2/100, min(VC_sf, VL_sf))
                a3_out = sf(a3, min(VC_sf, VL_sf))
                a4_out = sf(a4/10000, min(VC_sf, VL_sf))
                c1_out = sf(c1, min(CPC_sf, CPL_sf))
                c2_out = sf(c2/10000, min(CPC_sf, CPL_sf))
                wcon_out = sf(wcon/100000, min(SC_sf, SL_sf))
            else:
                G_out = round(G, rt)
                H_out = round(H, rt)
                S_out = round(S, rt)
                CP_out = round(CP2, rt)
                V_out = round(V2, rt)
                a1_out = round(a1*10, rt)
                a2_out = round(a2/100, rt)
                a3_out = round(a3, rt)
                a4_out = round(a4/10000, rt)
                c1_out = round(c1, rt)
                c2_out = round(c2/10000, rt)
                wcon_out = round(wcon/100000, rt)
            Z_out = Z2
        
            df_out_2, duplicate_list_2 = __write_output(filename=out_name,
                         cation=cation,
                         ligand=ligand, nth_complex=2,
                         G=G_out, H=H_out, S=S_out,
                         CP=CP_out, V=V_out, a1=a1_out,
                         a2=a2_out, a3=a3_out, a4=a4_out,
                         c1=c1_out, c2=c2_out,
                         wcon=wcon_out, Z=Z_out, azero=azero,
                         cation_dissrxn_dict=cation_dissrxn_dict,
                         ligand_dissrxn_dict=ligand_dissrxn_dict,
                         thermo_data=thermo_data,
                         cation_formula_ox_dict=cation_formula_ox_dict,
                         ligand_formula_ox_dict=ligand_formula_ox_dict,
                         skip_duplicates=skip_duplicates,
                         ligand_name_abbrv_pairs=ligand_name_abbrv_pairs,
                         organic=organic)

        eqs_used = eqs_used + hkf_eqns_2
        duplicate_list = duplicate_list + duplicate_list_2
    
    # override properties and parameters if they are supplied by the user
    if not math.isnan(G_3):
        G3 = G_3
    if not math.isnan(H_3):
        H3 = H_3
    if not math.isnan(S_3):
        S3 = S_3
    if not math.isnan(Cp_3):
        CP3 = Cp_3
    if not math.isnan(V_3):
        V3 = V_3
    
    if not math.isnan(BETA3):
        try:
            eqs_used.append("Estimating HKF parameters for the third complex...")
            G, H, S, a1, a2, a3, a4, c1, c2, wcon, hkf_eqns_3 = calc_params(Z3, G3, H3, S3, CP3, V3, a1_3, a2_3, a3_3, a4_3, c1_3, c2_3, omega_3, water_model=water_model, organic=organic)
            eqs_used.append("HKF parameters for the third complex were estimated...")
        except:
            eqs_used.append("HKF parameters for the third complex could not be estimated...")
            return pd.DataFrame(), eqs_used, warning_list, duplicate_list
        
        if not math.isnan(G):
            if sigfigs:
                G_out = sf(G, min(GC_sf, GL_sf))
                H_out = sf(H, min(HC_sf, HL_sf))
                S_out = sf(S, min(SC_sf, SL_sf))
                CP_out = sf(CP1, min(CPC_sf, CPL_sf))
                V_out = sf(V1, min(VC_sf, VL_sf))
                a1_out = sf(a1*10, min(VC_sf, VL_sf))
                a2_out = sf(a2/100, min(VC_sf, VL_sf))
                a3_out = sf(a3, min(VC_sf, VL_sf))
                a4_out = sf(a4/10000, min(VC_sf, VL_sf))
                c1_out = sf(c1, min(CPC_sf, CPL_sf))
                c2_out = sf(c2/10000, min(CPC_sf, CPL_sf))
                wcon_out = sf(wcon/100000, min(SC_sf, SL_sf))
            else:
                G_out = round(G, rt)
                H_out = round(H, rt)
                S_out = round(S, rt)
                CP_out = round(CP3, rt)
                V_out = round(V3, rt)
                a1_out = round(a1*10, rt)
                a2_out = round(a2/100, rt)
                a3_out = round(a3, rt)
                a4_out = round(a4/10000, rt)
                c1_out = round(c1, rt)
                c2_out = round(c2/10000, rt)
                wcon_out = round(wcon/100000, rt)
            Z_out = Z3

            df_out_3, duplicate_list_3 = __write_output(filename=out_name,
                         cation=cation,
                         ligand=ligand, nth_complex=3,
                         G=G_out, H=H_out, S=S_out,
                         CP=CP_out, V=V_out, a1=a1_out,
                         a2=a2_out, a3=a3_out, a4=a4_out,
                         c1=c1_out, c2=c2_out,
                         wcon=wcon_out, Z=Z_out, azero=azero,
                         cation_dissrxn_dict=cation_dissrxn_dict,
                         ligand_dissrxn_dict=ligand_dissrxn_dict,
                         thermo_data=thermo_data,
                         cation_formula_ox_dict=cation_formula_ox_dict,
                         ligand_formula_ox_dict=ligand_formula_ox_dict,
                         skip_duplicates=skip_duplicates,
                         ligand_name_abbrv_pairs=ligand_name_abbrv_pairs,
                         organic=organic)

        eqs_used = eqs_used + hkf_eqns_3
        duplicate_list = duplicate_list + duplicate_list_3
    
    # override properties and parameters if they are supplied by the user
    if not math.isnan(G_4):
        G4 = G_4
    if not math.isnan(H_4):
        H4 = H_4
    if not math.isnan(S_4):
        S4 = S_4
    if not math.isnan(Cp_4):
        CP4 = Cp_4
    if not math.isnan(V_4):
        V4 = V_4

    if not math.isnan(BETA4):
        try:
            eqs_used.append("Estimating HKF parameters for the fourth complex...")
            G, H, S, a1, a2, a3, a4, c1, c2, wcon, hkf_eqns_4 = calc_params(Z4, G4, H4, S4, CP4, V4, a1_4, a2_4, a3_4, a4_4, c1_4, c2_4, omega, water_model=water_model, organic=organic)
            eqs_used.append("HKF parameters for the fourth complex were estimated...")
        except:
            eqs_used.append("HKF parameters for the fourth complex could not be estimated...")
            return pd.DataFrame(), eqs_used, warning_list, duplicate_list

        if not math.isnan(G):
            if sigfigs:
                G_out = sf(G, min(GC_sf, GL_sf))
                H_out = sf(H, min(HC_sf, HL_sf))
                S_out = sf(S, min(SC_sf, SL_sf))
                CP_out = sf(CP1, min(CPC_sf, CPL_sf))
                V_out = sf(V1, min(VC_sf, VL_sf))
                a1_out = sf(a1*10, min(VC_sf, VL_sf))
                a2_out = sf(a2/100, min(VC_sf, VL_sf))
                a3_out = sf(a3, min(VC_sf, VL_sf))
                a4_out = sf(a4/10000, min(VC_sf, VL_sf))
                c1_out = sf(c1, min(CPC_sf, CPL_sf))
                c2_out = sf(c2/10000, min(CPC_sf, CPL_sf))
                wcon_out = sf(wcon/100000, min(SC_sf, SL_sf))
            else:
                G_out = round(G, rt)
                H_out = round(H, rt)
                S_out = round(S, rt)
                CP_out = round(CP4, rt)
                V_out = round(V4, rt)
                a1_out = round(a1*10, rt)
                a2_out = round(a2/100, rt)
                a3_out = round(a3, rt)
                a4_out = round(a4/10000, rt)
                c1_out = round(c1, rt)
                c2_out = round(c2/10000, rt)
                wcon_out = round(wcon/100000, rt)
            Z_out = Z4
            
            df_out_4, duplicate_list_4 = __write_output(filename=out_name,
                         cation=cation,
                         ligand=ligand, nth_complex=4,
                         G=G_out, H=H_out, S=S_out,
                         CP=CP_out, V=V_out, a1=a1_out,
                         a2=a2_out, a3=a3_out, a4=a4_out,
                         c1=c1_out, c2=c2_out,
                         wcon=wcon_out, Z=Z_out, azero=azero,
                         cation_dissrxn_dict=cation_dissrxn_dict,
                         ligand_dissrxn_dict=ligand_dissrxn_dict,
                         thermo_data=thermo_data,
                         cation_formula_ox_dict=cation_formula_ox_dict,
                         ligand_formula_ox_dict=ligand_formula_ox_dict,
                         skip_duplicates=skip_duplicates,
                         ligand_name_abbrv_pairs=ligand_name_abbrv_pairs,
                         organic=organic)


        eqs_used = eqs_used + hkf_eqns_4
        duplicate_list = duplicate_list + duplicate_list_4
    
    df_out = pd.concat([df_out_1, df_out_2, df_out_3, df_out_4], ignore_index=True)
    
    if len(duplicate_list) > 0:
        
        msg = ("Certain complexes that have just been estimated with the "
               "Complicator are likely duplicates of complexes that are already "
               "in the thermodynamic database. Set skip_duplicates=True "
               "and rerun to skip and exclude these complexes.")
        warning_list.append(msg)

        for dup_pair in duplicate_list:
            complex_name = dup_pair[0]
            complex_name_no_parentheses = dup_pair[1]
            
            warning_list.append(complex_name + " appears to be a duplicate of " + complex_name_no_parentheses + " already present in the thermodynamic database.")

    return df_out, eqs_used, warning_list, duplicate_list