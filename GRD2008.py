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
        'H2O': 18.01528,
        'F2O_1': 37.9968}

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

    mole_P2O5 = normalised_element_list[9] / molar_mass['P2O5']
    converted_element_list.append(mole_P2O5)

    mole_H2O =normalised_element_list[10] / molar_mass['H2O']
    converted_element_list.append(mole_H2O)

    mole_F2O_1 = normalised_element_list[11] / molar_mass['F2O_1']
    converted_element_list.append(mole_F2O_1)

    # return the created list with mole
    return converted_element_list


def compute_melt_viscosity_GRD08(ax):

    """ This function calculates the parameters of the non-Arrheanian VFT equation using the model
    of Giordano et al. 2008:

    log viscosity(Pa.s) = A + B / (T(K) - C),

    where A is a constant independent of composition and B and C are adjustable parameters

    Input data
    -----------
    json file containing the silicate melt chemical composition in oxide wt. %

    Returns
    ------------
    B and C of the VFT equation
    Tg: the glass transition temperature at viscosity = e12 Pa.s

    Reference
    ---------
    Giordano, D., Russell, J. K., & Dingwell, D. B. (2008). Viscosity of magmatic liquids: a model.
    Earth and Planetary Science Letters, 271(1), 123-134.

    """

    # TODO: here enter the path to the json file containing the composition in oxide weight %

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
        wt_P2O5 = float(data['melt_oxide_composition']['P2O5'])
        wt_H2O = float(data['melt_oxide_composition']['H2O'])
        wt_F2O_1 = float(data['melt_oxide_composition']['F2O_1'])

    # convert Fe2O3 in FeO and compute FeOtot
    wt_FeOtot = wt_Fe2O3 * ((2. * 71.844) / 159.688) + wt_FeO

    element_list = [wt_SiO2, wt_TiO2, wt_Al2O3, wt_FeOtot, wt_MnO, wt_MgO, wt_CaO, wt_Na2O, wt_K2O, wt_P2O5, wt_H2O,
                    wt_F2O_1]

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
    wtn_P2O5 = wt_P2O5 * 100. / sum_all_element
    wtn_H2O = wt_H2O * 100. / sum_all_element
    wtn_F2O_1 = wt_F2O_1 * 100. / sum_all_element

    normalised_element_list = [wtn_SiO2, wtn_TiO2, wtn_Al2O3, wtn_FeOtot, wtn_MnO, wtn_MgO, wtn_CaO, wtn_Na2O, wtn_K2O,
                               wtn_P2O5, wtn_H2O, wtn_F2O_1]

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
    mf_P2O5 = element_list_in_mole[9] * 100. / sum_all_element_in_mole
    mf_H2O = element_list_in_mole[10] * 100. / sum_all_element_in_mole

    # fitting parameters:
    B1 = 159.60
    B2 = -173.30
    B3 = 72.10
    B4 = 75.70
    B5 = -39.00
    B6 = -84.10
    B7 = 141.50
    B11 = -2.43
    B12 = -0.91
    B13 = 17.60

    C1 = 2.75
    C2 = 15.70
    C3 = 8.30
    C4 = 10.20
    C5 = -12.30
    C6 = -99.50
    C11 = 0.30

    # calculate VFT paramters
    A = -4.55
    B = (mf_SiO2 + mf_TiO2) * B1 + \
        mf_Al2O3 * B2 + \
        (mf_FeOtot + mf_MnO + mf_P2O5) * B3 + \
        mf_MgO * B4 + \
        mf_CaO * B5 + \
        (mf_Na2O + mf_H2O) * B6 + \
        ((mf_H2O) + math.log(1. + mf_H2O)) * B7 + \
        (mf_SiO2 + mf_TiO2) * (mf_FeOtot + mf_MnO + mf_MgO) * B11 + \
        (mf_SiO2 + mf_TiO2 + mf_Al2O3 + mf_P2O5) * (mf_Na2O + mf_K2O + mf_H2O) * B12 + \
        mf_Al2O3 * (mf_Na2O + mf_K2O) * B13

    C = mf_SiO2 * C1 + \
        (mf_TiO2 + mf_Al2O3) * C2 + \
        (mf_FeOtot + mf_MnO + mf_MgO) * C3 + \
        mf_CaO * C4 + \
        (mf_Na2O + mf_K2O) * C5 + \
        math.log(1. + mf_H2O) * C6 + \
        (mf_Al2O3 + mf_FeOtot + mf_MnO + mf_MgO + mf_CaO - mf_P2O5) * (mf_Na2O + mf_K2O + mf_H2O) * C11

    Tg = B / (12 - A) + C
    logviscosity_liquid_GRD08 = A + B / (temperature + 273.15 - C)
    print('A =', A)
    print('B =', B)
    print('C =', C)
    print('Tg =', Tg)
    viscosity_liquid_GRD08 = 10**logviscosity_liquid_GRD08
    print('viscosity =', viscosity_liquid_GRD08, 'Pa.s (log(viscosity)=', logviscosity_liquid_GRD08,') at T =',temperature,'(°C)')

    inverse_temperature = []
    logviscosity_GRD = []

    for temp in range(700, 1400):
        logvisco_GRD = A + B/ (temp + 273.15 - C)
        inverse_temperature.append(10000.0 / (temp + 273.15))
        logviscosity_GRD.append(logvisco_GRD)

    #fig1 = plt.figure()
    #ax = fig1.add_subplot(111)
    #plt.title('GRD model')
    ax.plot(inverse_temperature, logviscosity_GRD, '-', color='r', label='GRD model')
    ax.legend(loc=1, prop={'size': 8})
    ax.set_xlabel("10000 / T (K)")
    ax.set_ylabel("log viscosity (Pa.s)")
    plt.show()

    return viscosity_liquid_GRD08, A, B, C


if __name__ == "__main__":
    fig1 = plt.figure()

    ax = fig1.add_subplot(111)

    liquid_viscosity_model_GRD = compute_melt_viscosity_GRD08(ax)


