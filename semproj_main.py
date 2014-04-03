from numpy import *


# Creates the ID array used in all the other scripts
def ID_array():
    ID = array(['00358', '01590', '04054', '04078', '04913', '10394', '12242',
                '12291', '12713', '12862', '12863', '12864', '12865', '12866',
                '12867', '12868', '12869', '12870', '12871', '12872', '12873',
                '13896', '13898', '13899'])
    RA = array([213.311845, 253.245244, 204.565931, 331.758169, 209.011858,
                32.4103586, 198.821929, 186.444822, 253.245172, 6.38405016,
                174.122207, 90.5434278, 136.153654, 38.6575224, 100.048517,
                31.5664714, 116.037870, 37.0600602, 171.400592, 10.7198346,
                176.418473, 171.702015, 15.1453511, 18.4625008])
    DEC = array([-3.20753360, 2.40109491, 4.54279789, 10.2336514, 18.3719250,
                 -10.1466462, 44.4073438, 12.6622092, 2.40101865, 68.3626943,
                 21.5962976, 28.4721918, 55.6008422, -8.78774935, -25.8951959,
                 -0.291330007, 29.2476806, 31.3116767, 54.3828728, -23.5409277,
                 -18.4540326, 35.2508226, -47.8675907, 13.2717757])
    z = array([0.0062, 0.024, 0.023, 0.027, 0.050, 0.013, 0.037, 0.0084,
               0.024, 0.012, 0.030, 0.033, 0.037, 0.043, 0.025, 0.042,
               0.016, 0.017, 0.021, 0.022, 0.033, 0.032, 0.048, 0.050])
    name = array(['NGC 5506', 'NGC 6240', 'NGC 5252', 'NGC 7212', 'Mrk 463',
                  'NGC 0838', 'UGC 08327', 'NGC 4388', 'NGC 6240',
                  '2MASX J00253292+6821442', 'Mrk 739', 'IRAS 05589+2828',
                  '2MASX J09043699+5536025', 'NGC 0985', 'ESO 490-IG026',
                  'Mrk 1018', 'UGC 03995 NOTES01', 'NGC 0931', 'ARP 151',
                  'NGC 0235A', '2MASX J11454045-1827149', 'Mrk 423',
                  'ESO 195-IG 021 NED03', 'Mrk 975'])
    dist = array([27.2, 1.5, 54.8, 8.9, 3.9, 7.7, 26.1, 49.2, 1.5, 1.2, 3.4, 8,
                  9, 44.7, 16.4, 44.4, 10.3, 6, 8.2, 9, 14, 5.9, 17.9, 13.1])
    area_ext = array([33.62351955, 38.8780399, 30.73492265, 38.91548009,
                      34.30817824, 38.91937507, 38.9014255, 38.9311909,
                      38.8780399, 38.92611961, 30.88875137, 32.10554299,
                      30.10141226, 26.25091718, 30.23754831, 28.49346697,
                      38.94176865, 33.9022708, 32.2574816, 38.89774288,
                      30.02613024, 31.04877957, 38.94595966, 38.87812167])
    ret = vstack((ID, RA, DEC, z, name, dist, area_ext))
    return ret


# Creates the list of double sources
def double_list():
    double = ['04913', '12863']
    return double


# Creates the list of piledup sources
def piledup():
    piled = ['04054', '12242', '12291', '12863', '12864', '12865', '12866',
             '12867', '12868', '12869', '12870', '12871', '12872', '12873',
             '13896', '13899']
    return piled
