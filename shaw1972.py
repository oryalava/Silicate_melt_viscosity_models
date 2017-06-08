import json
import math
import matplotlib.pyplot as plt
import numpy


def convert_to_mole(normalised_element_list):
    """ This function converts the normalised oxides weigh percent into
    mole percent that is used for the silicate melt viscosity model"""
    converted_element_list = []
    molar_mass = {
        'SiO2': 60.0843,
        'TiO2': 79.8658,
        'Al2O3': 101.961276,
        'Fe2O3': 159.69,
        'FeO': 71.8444,
        'MnO': 70.937449,
        'MgO': 40.3044,
        'CaO': 56.0774,
        'Na2O': 61.97894,
        'K2O': 94.2,
        'P2O5': 141.9446,
        'H2O': 18.01528}

    mole_SiO2 = normalised_element_list[0] / molar_mass['SiO2']
    converted_element_list.append(mole_SiO2)

    mole_TiO2 = normalised_element_list[1] / molar_mass['TiO2']
    converted_element_list.append(mole_TiO2)

    mole_Al2O3 = normalised_element_list[2] / molar_mass['Al2O3']
    converted_element_list.append(mole_Al2O3)

    mole_FeOtot = normalised_element_list[3]/ molar_mass['FeO']
    converted_element_list.append(mole_FeOtot)

    mole_MnO = normalised_element_list[4] / molar_mass['MnO']
    converted_element_list.append(mole_MnO)

    mole_MgO = normalised_element_list[5] / molar_mass['MgO']
    converted_element_list.append(mole_MgO)

    mole_CaO = normalised_element_list[6] / molar_mass['CaO']
    converted_element_list.append(mole_CaO)

    mole_Na2O = normalised_element_list[7] / molar_mass['Na2O']
    converted_element_list.append(mole_Na2O)

    mole_K2O = normalised_element_list[8] / molar_mass['K2O']
    converted_element_list.append(mole_K2O)

    # mole_P2O5 = normalised_element_list[9] / molar_mass['P2O5']
    # converted_element_list.append(mole_P2O5)

    mole_H2O =normalised_element_list[9] / molar_mass['H2O']
    converted_element_list.append(mole_H2O)

    #mole_F2O_1 = normalised_element_list[11] / molar_mass['F2O_1']
    #converted_element_list.append(mole_F2O_1)

    # return the created list with mole
    return converted_element_list


def compute_melt_viscosity_Shaw72(ax):
    """ This function calculates the parameters of the Arrheanian equation using the model
    of Shaw 1972 to estimate the viscosity of silicate liquid

    ln viscosity(Poise) = slope*(10000/T(K))-(1.5*slope)-6.4

    where slope is the intercept calculated from the chemical composotion of the silicate liquid

    Input data
    -----------
    json file containing the silicate melt chemical composition in oxide wt. %

    Returns
    ------------
    the slope and the viscosity in Pa.s

    Reference
    ---------
    Shaw (1972). Viscosity of magmatic silicate liquids: an empirical method of prediction.
    American Journal of science, Vol. 272, p. 870-893.

    """

    # TODO: here read the composition in oxide weight % from json file

    filename = 'Etna74.json'

    with open(filename) as data_file:
        data = json.load(data_file)
        temperature = float(data['temperature'])
        wt_SiO2 = float(data['melt_oxide_composition']['SiO2'])
        wt_TiO2 = float(data['melt_oxide_composition']['TiO2'])
        wt_Al2O3 = float(data['melt_oxide_composition']['Al2O3'])
        wt_FeO = float(data['melt_oxide_composition']['FeO'])
        wt_Fe2O3 = float(data['melt_oxide_composition']['Fe203'])
        wt_MnO = float(data['melt_oxide_composition']['MnO'])
        wt_MgO = float(data['melt_oxide_composition']['MgO'])
        wt_CaO = float(data['melt_oxide_composition']['CaO'])
        wt_Na2O = float(data['melt_oxide_composition']['Na2O'])
        wt_K2O = float(data['melt_oxide_composition']['K2O'])
        #wt_P2O5 = float(data['melt_oxide_composition']['P2O5'])
        wt_H2O = float(data['melt_oxide_composition']['H2O'])
        #wt_F2O_1 = float(data['melt_oxide_composition']['F2O_1'])

    # convert Fe2O3 in FeO and compute FeOtot
    wt_FeOtot = wt_Fe2O3 * ((2. * 71.844) / 159.688) + wt_FeO

    element_list = [wt_SiO2, wt_TiO2, wt_Al2O3, wt_FeOtot, wt_MgO, wt_CaO, wt_Na2O, wt_K2O, wt_H2O]

    # normalize values (weight normalised = wtn)
    sum_all_element = 0.
    for i in range(0, len(element_list)):
        sum_all_element += element_list[i]
    wtn_SiO2 = wt_SiO2 * 100. / sum_all_element
    wtn_TiO2 = wt_TiO2 * 100. / sum_all_element
    wtn_Al2O3 = wt_Al2O3 * 100. / sum_all_element
    wtn_FeOtot = wt_FeOtot * 100. / sum_all_element
    wtn_MnO = wt_MnO * 100. / sum_all_element
    wtn_MgO = wt_MgO * 100. / sum_all_element
    wtn_CaO = wt_CaO * 100. / sum_all_element
    wtn_Na2O = wt_Na2O * 100. / sum_all_element
    wtn_K2O = wt_K2O * 100. / sum_all_element
    #wtn_P2O5 = wt_P2O5 * 100. / sum_all_element
    wtn_H2O = wt_H2O * 100. / sum_all_element
    #wtn_F2O_1 = wt_F2O_1 * 100. / sum_all_element

    normalised_element_list = [wtn_SiO2, wtn_TiO2, wtn_Al2O3, wtn_FeOtot, wtn_MnO, wtn_MgO, wtn_CaO, wtn_Na2O, wtn_K2O, wtn_H2O]

    # convert in mole:
    element_list_in_mole = convert_to_mole(normalised_element_list)
    # normalize values and obtains mole fractions : mf_
    sum_all_element_in_mole = 0.
    for i in range(0, len(element_list_in_mole)):
        sum_all_element_in_mole += element_list_in_mole[i]

    mf_SiO2 = element_list_in_mole[0] * 100. / sum_all_element_in_mole
    mf_TiO2 = element_list_in_mole[1] * 100. / sum_all_element_in_mole
    mf_Al2O3 = element_list_in_mole[2] * 100. / sum_all_element_in_mole
    mf_FeOtot = element_list_in_mole[3] * 100. / sum_all_element_in_mole
    mf_MnO = element_list_in_mole[4] * 100. / sum_all_element_in_mole
    mf_MgO = element_list_in_mole[5] * 100. / sum_all_element_in_mole
    mf_CaO = element_list_in_mole[6] * 100. / sum_all_element_in_mole
    mf_Na2O = element_list_in_mole[7] * 100. / sum_all_element_in_mole
    mf_K2O = element_list_in_mole[8] * 100. / sum_all_element_in_mole
    #mf_P2O5 = element_list_in_mole[9] * 100. / sum_all_element_in_mole
    mf_H2O = element_list_in_mole[9] * 100. / sum_all_element_in_mole
   #mf_F2O_1= element_list_in_mole[11] * 100. / sum_all_element_in_mole

    xsi = mf_SiO2 / 100.0
    xal = mf_Al2O3 / 100.0
    xfm = (mf_FeOtot + mf_MnO + mf_MgO) / 100.0
    xct = (mf_CaO + mf_TiO2) / 100.0
    xnk = (mf_Na2O + mf_K2O) / 100.0
    xh = mf_H2O / 100.0
    sum = xsi + xal + xfm + xct + xnk + xh
    if sum != 1:
        print('something is wrong in Shaw model')

    p1 = xal * xsi * 6.7
    p2 = xfm * xsi * 3.4
    p3 = xct * xsi * 4.5
    p4 = xnk * xsi * 2.8
    p5 = xh * xsi * 2.0

    total_p_values = p1+p2+p3+p4+p5
    shaw_slope = total_p_values /(1- xsi)

    lnvisco= shaw_slope * (10000/(temperature + 273))-(1.5 * shaw_slope) - 6.4 #in poise
    logvisco = lnvisco/2.303 - 1 # in Pa.s
    viscosity = 10 ** logvisco
    print('Slope =', shaw_slope)
    print('viscosity =', viscosity, 'Pa.s (log(viscosity)=', logvisco,') at T =',temperature,'(°C)')

    inverse_temperature = []
    logviscosity = []

    for temp in range(700, 1400):
        logvisco = (shaw_slope * (10000/(temp + 273))-(1.5 * shaw_slope) - 6.4)/2.303 - 1
        inverse_temperature.append(10000.0 / (temp + 273.15))
        logviscosity.append(logvisco)

    ax.plot(inverse_temperature, logviscosity, '-', color='b', label='Shaw72 model')
    ax.legend(loc=1, prop={'size': 8})
    ax.set_xlabel("10000 / T (K)")
    ax.set_ylabel("log viscosity (Pa.s)")
    plt.show()
    return shaw_slope


if __name__ == "__main__":
    fig1 = plt.figure()

    ax = fig1.add_subplot(111)
    liquid_viscosity_model = compute_melt_viscosity_Shaw72(ax)
