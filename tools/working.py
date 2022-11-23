import numpy as np
from cbsyst import Csys

# variables used to calculate Q
Q_alpha = 1e-4
Q_beta = 7e-4
Q_k = 8.3e17

# salinity balance - the total moles of salt added or removed to the surface boxes
E = 1.52e14

def ocean_model_TS(lolat, hilat, deep, tmax, dt):
    """Run the ocean model for a given time period and return the results for each box.

    Parameters
    ----------
    lolat, hilat, deep : dict
        dictionaries containing the box properties
    tmax : int or float
        The maximum time to run the model for (yr)
    dt : float
        The time step of the model (yr)

    Returns
    -------
    tuple of (time, lolat, hilat, deep)
    """
    # set which variables will change with time
    model_vars = ['T', 'S']

    time = np.arange(0, tmax + dt, dt)  # the time axis for the model
    
    # create copies of the input dictionaries so we don't modify the originals
    lolat = lolat.copy()
    hilat = hilat.copy()
    deep = deep.copy()
    
    # turn all time-evolving variables into arrays containing the start values
    for box in [lolat, hilat, deep]:
        for k in model_vars:
            box[k] = np.full(time.shape, box[k])

    fluxes = {}

    for i in range(1, time.size):
        last = i - 1  # index of last model step
        
        # calculate circulation flux, Q
        dT = lolat['T'][last] - hilat['T'][last]
        dS = lolat['S'][last] - hilat['S'][last]
        Q_T = Q_k * (Q_alpha * dT - Q_beta * dS)
        
        # calculate mixing fluxes for T and S
        for var in model_vars:
            fluxes[f'Q_{var}_deep'] = Q_T * (hilat[var][last] - deep[var][last]) * dt  # mol dt-1
            fluxes[f'Q_{var}_hilat'] = Q_T * (lolat[var][last] - hilat[var][last]) * dt  # mol dt-1
            fluxes[f'Q_{var}_lolat'] = Q_T * (deep[var][last] - lolat[var][last]) * dt  # mol dt-1
            
            fluxes[f'vmix_{var}_hilat'] = hilat['V'] / hilat['tau_M'] * (hilat[var][last] - deep[var][last]) * dt  # mol dt-1
            fluxes[f'vmix_{var}_lolat'] = lolat['V'] / lolat['tau_M'] * (lolat[var][last] - deep[var][last]) * dt  # mol dt-1
        
        # calculate temperature exchange with each surface box
        for box in [hilat, lolat]:
            boxname = box['name']
            # temperature exchange with atmosphere
            fluxes[f'dT_{boxname}'] = box['V'] / box['tau_T'] * (box['T_atmos'] - box['T'][last]) * dt  # mol dt-1

        # update deep box
        for var in model_vars:
            deep[var][i] = deep[var][last] + (
                fluxes[f'Q_{var}_deep'] + fluxes[f'vmix_{var}_hilat'] + fluxes[f'vmix_{var}_lolat']
            ) / deep['V']
        
        # update surface boxes
        for box in [hilat, lolat]:
            boxname = box['name']
            box['S'][i] = box['S'][last] + (fluxes[f'Q_S_{boxname}'] - fluxes[f'vmix_S_{boxname}'] + box['E'] * dt) / box['V']
            box['T'][i] = box['T'][last] + (fluxes[f'Q_T_{boxname}'] - fluxes[f'vmix_T_{boxname}'] + fluxes[f'dT_{boxname}']) / box['V']
    
    return time, lolat, hilat, deep

def ocean_model_DIC(lolat, hilat, deep, atmos, tmax, dt):
    """Run the ocean model for a given time period and return the results for each box.

    Parameters
    ----------
    lolat, hilat, deep, atmos : dict
        dictionaries containing the box properties
    tmax : int or float
        The maximum time to run the model for (yr)
    dt : float
        The time step of the model (yr)

    Returns
    -------
    tuple of (time, lolat, hilat, deep)
    """

    # create the time scale for the model
    time = np.arange(0, tmax + dt, dt)
    
    # identify which variables will change with time
    model_vars = ['T', 'S', 'DIC', 'TA']
    atmos_model_vars = ['moles_CO2', 'pCO2']

    # create copies of the input dictionaries so we don't modify the originals
    lolat = lolat.copy()
    hilat = hilat.copy()
    deep = deep.copy()
    atmos = atmos.copy()
    
    # turn all time-evolving variables into arrays containing the start values
    for box in [lolat, hilat, deep]:
        for k in model_vars:
            box[k] = np.full(time.shape, box[k])
    for k in atmos_model_vars:
        atmos[k] = np.full(time.shape, atmos[k])
    if isinstance(atmos['GtC_emissions'], (int, float)):
        atmos['GtC_emissions'] = np.full(time.shape, atmos['GtC_emissions'])
    
    # calculate initial surface carbon chemistry in the surface boxes using Csys, and store a few key variables - CO2, pH, pCO2 and K0
    for box in [lolat, hilat]:
        csys = Csys(
            TA=1e3 * box['TA'],  # 1e3 is necessary because cbsyst expects units of umol L-1
            DIC=1e3 * box['DIC'],  # 1e3 is necessary because cbsyst expects units of umol L-1
            T_in=box['T'], S_in=box['S'],
            )
        box['CO2'] = csys.CO2 * 1e-3  # 1e3 converts back to mol m-3
        box['pH'] = csys.pHtot
        box['pCO2'] = csys.pCO2
        box['K0'] = csys.Ks.K0
    
    # Create a dictionary to keep track of the fluxes calculated at each step    
    fluxes = {}
    
    ### LOOP STARTS HERE ###
    for i in range(1, time.size):
        last = i - 1  # index of last model step
        
        # calculate circulation flux, Q
        dT = lolat['T'][last] - hilat['T'][last]
        dS = lolat['S'][last] - hilat['S'][last]
        Q_T = Q_k * (Q_alpha * dT - Q_beta * dS)
        
        # calculate mixing fluxes for model variables
        for var in model_vars:
            fluxes[f'Q_{var}_deep'] = Q_T * (hilat[var][last] - deep[var][last]) * dt  # mol dt-1
            fluxes[f'Q_{var}_hilat'] = Q_T * (lolat[var][last] - hilat[var][last]) * dt  # mol dt-1
            fluxes[f'Q_{var}_lolat'] = Q_T * (deep[var][last] - lolat[var][last]) * dt  # mol dt-1
            
            fluxes[f'vmix_{var}_hilat'] = hilat['V'] / hilat['tau_M'] * (hilat[var][last] - deep[var][last]) * dt  # mol dt-1
            fluxes[f'vmix_{var}_lolat'] = lolat['V'] / lolat['tau_M'] * (lolat[var][last] - deep[var][last]) * dt  # mol dt-1
        
        # calculate surface-specific fluxes
        for box in [hilat, lolat]:
            boxname = box['name']
            # temperature exchange with atmosphere
            fluxes[f'dT_{boxname}'] = box['V'] / box['tau_T'] * (box['T_atmos'] - box['T'][last]) * dt  # mol dt-1
            # CO2 exchange with atmosphere
            fluxes[f'dCO2_{boxname}'] = box['V'] / box['tau_CO2'] * (box['CO2'][last] - 1e-3 * atmos['pCO2'][last] * box['K0'][last]) * dt  # mol dt-1
        
        fluxes['dCO2_emissions'] = atmos['GtC_emissions'][last] * 1e15 / 12 * dt  # mol dt-1
        
        # update deep box
        for var in model_vars:
            deep[var][i] = deep[var][last] + (
                fluxes[f'Q_{var}_deep'] + fluxes[f'vmix_{var}_hilat'] + fluxes[f'vmix_{var}_lolat']
            ) / deep['V']
        
        # update surface boxes
        for box in [hilat, lolat]:
            boxname = box['name']
            box['S'][i] = box['S'][last] + (fluxes[f'Q_S_{boxname}'] - fluxes[f'vmix_S_{boxname}'] + box['E'] * dt) / box['V']
            box['T'][i] = box['T'][last] + (fluxes[f'Q_T_{boxname}'] - fluxes[f'vmix_T_{boxname}'] + fluxes[f'dT_{boxname}']) / box['V']
            box['DIC'][i] = box['DIC'][last] + (fluxes[f'Q_DIC_{boxname}'] - fluxes[f'vmix_DIC_{boxname}'] - fluxes[f'dCO2_{boxname}']) / box['V']
            box['TA'][i] = box['TA'][last] + (fluxes[f'Q_TA_{boxname}'] - fluxes[f'vmix_TA_{boxname}']) / box['V']
            
            # update carbon speciation
            csys = Csys(
                TA=1e3 * box['TA'][i],  # 1e3 is necessary because cbsyst expects units of umol L-1
                DIC=1e3 * box['DIC'][i],  # 1e3 is necessary because cbsyst expects units of umol L-1
                T_in=box['T'][i], S_in=box['S'][i],
                )
            box['CO2'][i] = csys.CO2 * 1e-3  # 1e3 converts back to mol m-3
            box['pH'][i] = csys.pHtot
            box['pCO2'][i] = csys.pCO2
            box['K0'][i] = csys.Ks.K0
        
        # update atmosphere
        atmos['moles_CO2'][i] = atmos['moles_CO2'][last] + fluxes['dCO2_hilat'] + fluxes['dCO2_lolat'] + fluxes['dCO2_emissions']
        atmos['pCO2'][i] = 1e6 * atmos['moles_CO2'][i] / atmos['moles_air']
        
    return time, lolat, hilat, deep, atmos


def ocean_model_bio(lolat, hilat, deep, atmos, tmax, dt):

    # create the time scale for the model
    time = np.arange(0, tmax + dt, dt)
    
    # identify which variables will change with time
    model_vars = ['T', 'S', 'DIC', 'TA', 'PO4']
    atmos_model_vars = ['moles_CO2', 'pCO2']

    # create copies of the input dictionaries so we don't modify the originals
    lolat = lolat.copy()
    hilat = hilat.copy()
    deep = deep.copy()
    atmos = atmos.copy()
    
    # turn all time-evolving variables into arrays containing the start values
    for box in [lolat, hilat, deep]:
        for k in model_vars:
            box[k] = np.full(time.shape, box[k])
    for k in atmos_model_vars:
        atmos[k] = np.full(time.shape, atmos[k])
    if isinstance(atmos['GtC_emissions'], (int, float)):
        atmos['GtC_emissions'] = np.full(time.shape, atmos['GtC_emissions'])
    
    # calculate initial surface carbon chemistry in the surface boxes using Csys, and store a few key variables - CO2, pH, pCO2 and K0
    for box in [lolat, hilat]:
        csys = Csys(
            TA=1e3 * box['TA'],  # 1e3 is necessary because cbsyst expects units of umol L-1
            DIC=1e3 * box['DIC'],  # 1e3 is necessary because cbsyst expects units of umol L-1
            T_in=box['T'], S_in=box['S'],
            )
        box['CO2'] = csys.CO2 * 1e-3  # 1e3 converts back to mol m-3
        box['pH'] = csys.pHtot
        box['pCO2'] = csys.pCO2
        box['K0'] = csys.Ks.K0
    
    # Create a dictionary to keep track of the fluxes calculated at each step    
    fluxes = {}

    for i in range(1, time.size):
        last = i - 1  # index of last model step
        
        # calculate circulation flux, Q
        dT = lolat['T'][last] - hilat['T'][last]
        dS = lolat['S'][last] - hilat['S'][last]
        Q = Q_k * (Q_alpha * dT - Q_beta * dS)
        
        # calculate mixing fluxes for model variables
        for var in model_vars:
            fluxes[f'Q_{var}_deep'] = Q * (hilat[var][last] - deep[var][last]) * dt  # mol dt-1
            fluxes[f'Q_{var}_hilat'] = Q * (lolat[var][last] - hilat[var][last]) * dt  # mol dt-1
            fluxes[f'Q_{var}_lolat'] = Q * (deep[var][last] - lolat[var][last]) * dt  # mol dt-1
            
            fluxes[f'vmix_{var}_hilat'] = hilat['V'] / hilat['tau_M'] * (hilat[var][last] - deep[var][last]) * dt  # mol dt-1
            fluxes[f'vmix_{var}_lolat'] = lolat['V'] / lolat['tau_M'] * (lolat[var][last] - deep[var][last]) * dt  # mol dt-1
        
        # calculate surface-specific fluxes
        for box in [hilat, lolat]:
            boxname = box['name']
            # temperature exchange with atmosphere
            fluxes[f'dT_{boxname}'] = box['V'] / box['tau_T'] * (box['T_atmos'] - box['T'][last]) * dt  # mol dt-1
            # CO2 exchange with atmosphere
            fluxes[f'dCO2_{boxname}'] = box['V'] / box['tau_CO2'] * (box['CO2'][last] - 1e-3 * atmos['pCO2'][last] * box['K0'][last]) * dt  # mol dt-1
            # organic matter production
            fluxes[f'export_PO4_{boxname}'] = box['PO4'][last] * box['V'] / box['tau_PO4'] * dt  # mol PO4 dt-1
            # DIC export by productivity :                                       redfield + calcification
            fluxes[f'export_DIC_{boxname}'] = fluxes[f'export_PO4_{boxname}'] * (106 + 106 * box['f_CaCO3'])  # mol DIC dt-1
            # TA export by productivity :                                       calcification
            fluxes[f'export_TA_{boxname}'] = fluxes[f'export_PO4_{boxname}'] * (2 * 106 * box['f_CaCO3'])  # mol TA dt-1
            
        fluxes['dCO2_emissions'] = atmos['GtC_emissions'][last] * 1e15 / 12 * dt  # mol dt-1

        # update deep box
        for var in model_vars:
            if var in ['T', 'S']:
                deep[var][i] = deep[var][last] + (
                    fluxes[f'Q_{var}_deep'] + fluxes[f'vmix_{var}_hilat'] + fluxes[f'vmix_{var}_lolat']
                ) / deep['V']
            else:
                deep[var][i] = deep[var][last] + (
                    fluxes[f'Q_{var}_deep'] + fluxes[f'vmix_{var}_hilat'] + fluxes[f'vmix_{var}_lolat'] + fluxes[f'export_{var}_hilat'] + fluxes[f'export_{var}_lolat']
                ) / deep['V']
        
        # update surface boxes
        for box in [hilat, lolat]:
            boxname = box['name']
            box['S'][i] = box['S'][last] + (fluxes[f'Q_S_{boxname}'] - fluxes[f'vmix_S_{boxname}'] + box['E'] * dt) / box['V']
            box['T'][i] = box['T'][last] + (fluxes[f'Q_T_{boxname}'] - fluxes[f'vmix_T_{boxname}'] + fluxes[f'dT_{boxname}']) / box['V']
            box['DIC'][i] = box['DIC'][last] + (fluxes[f'Q_DIC_{boxname}'] - fluxes[f'vmix_DIC_{boxname}'] - fluxes[f'dCO2_{boxname}'] - fluxes[f'export_DIC_{boxname}']) / box['V']
            box['TA'][i] = box['TA'][last] + (fluxes[f'Q_TA_{boxname}'] - fluxes[f'vmix_TA_{boxname}'] - fluxes[f'export_TA_{boxname}']) / box['V']
            box['PO4'][i] = box['PO4'][last] + (fluxes[f'Q_PO4_{boxname}'] - fluxes[f'vmix_PO4_{boxname}'] - fluxes[f'export_PO4_{boxname}']) / box['V']

            # update carbon speciation
            csys = Csys(
                TA=1e3 * box['TA'][i],  # 1e3 is necessary because cbsyst expects units of umol L-1
                DIC=1e3 * box['DIC'][i],  # 1e3 is necessary because cbsyst expects units of umol L-1
                T_in=box['T'][i], S_in=box['S'][i],
                )
            box['CO2'][i] = csys.CO2 * 1e-3  # 1e3 converts back to mol m-3
            box['pCO2'][i] = csys.pCO2
            box['pH'][i] = csys.pHtot
            box['K0'][i] = csys.Ks.K0
        
        # update atmosphere
        atmos['moles_CO2'][i] = atmos['moles_CO2'][last] + fluxes['dCO2_hilat'] + fluxes['dCO2_lolat'] + fluxes['dCO2_emissions']
        atmos['pCO2'][i] = 1e6 * atmos['moles_CO2'][i] / atmos['moles_air']
    
    return time, lolat, hilat, deep, atmos