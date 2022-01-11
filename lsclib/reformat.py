import pandas as pd
import os
import csv

# GLOBAL VARIABLES
# extrinsic params
light_forms = ['direct', 'diffuse', 'ground']
output_light_formats = ['beam_circ', 'iso', 'grnd']  # output directory names
phis = [0.001, 15.0, 30.0, 45.0, 60.0, 75.0]
thetas = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.999, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0]
tilts = [5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]
n_sample = 100


# convert list of results to phi & theta matrix
def listToMatrix(l1):
    result = []
    row = []
    for i in range(0, len(l1)):
        row.append(float(l1[i]))
        if (i + 1) % 13 == 0:
            result.append(row)
            row = []
    # insert phi labels (row label)
    for i in range(len(result)):
        result[i].insert(0, phis[i])
    return result


# write to csv file
def writeCSV(fname, fields, rows):
    with open(fname, 'w') as f:
        write = csv.writer(f)
        write.writerow(fields)
        write.writerows(rows)
    return


if __name__ == "__main__":
    # Generate data directory for this light form
    if not os.path.isdir("./reformatted_data"):
        os.mkdir("./reformatted_data")
    for light_f in output_light_formats:
        if not os.path.isdir("./reformatted_data/{}".format(light_f)):
            os.mkdir("./reformatted_data/{}".format(light_f))
            os.mkdir("./reformatted_data/{}/opt_eff".format(light_f))
            os.mkdir("./reformatted_data/{}/spect_mismatch".format(light_f))

    parameters = pd.DataFrame()
    # Length	Short Side Position	Film Concentration (kg/m^3)
    csv_data_params = pd.read_csv('./sim_data/direct/{}_{}_output.csv'.format(phis[0], thetas[0]))
    parameters['Length'] = csv_data_params['length']
    parameters['Short Side Position'] = csv_data_params['short_side_height']
    parameters['Film Concentration (kg/m^3)'] = csv_data_params['film_concentration']
    parameters.to_csv('./reformatted_data/parameters.csv')

    # read direct light csv files
    direct_dicts_list = []  # list that store all 78 dicts(storing information from csv)
    for phi in phis:
        for theta in thetas:
            csv_data = pd.read_csv('./sim_data/direct/{}_{}_output.csv'.format(phi, theta))  # read output data
            dict = csv_data.to_dict()  # dict that store a single input csv file
            direct_dicts_list.append(dict)
    # print(direct_dicts_list[0])

    l_opt_effs_direct = []  # output files data list               visual: [file1, file2, file3, ....]
    l_spect_mismatchs_direct = []  # output files data list
    for i in range(n_sample):
        l_opt_effs_direct.append([])
        l_spect_mismatchs_direct.append([])

    for i in range(n_sample):
        for dict in direct_dicts_list:
            l_opt_effs_direct[i].append(dict['opt_eff'][i])
            l_spect_mismatchs_direct[i].append(dict['spect_mismatch'][i])

    # Do writing
    fields_direct = thetas[:]
    fields_direct.insert(0, '')
    for i in range(len(l_opt_effs_direct)):
        fname_opt_eff = './reformatted_data/beam_circ/opt_eff/opt_eff_{}.csv'.format(i)
        fname_spect_mismatch = './reformatted_data/beam_circ/spect_mismatch/spect_mismatch_{}.csv'.format(i)
        opt_eff_rows = listToMatrix(l_opt_effs_direct[i])
        spect_mismatch_rows = listToMatrix(l_spect_mismatchs_direct[i])
        writeCSV(fname_opt_eff,fields_direct,opt_eff_rows)
        writeCSV(fname_spect_mismatch, fields_direct, spect_mismatch_rows)
    print("Finish writing csv files for direct lights")


    # do the same with diffuse & ground
    diffuse_dicts_list = []
    ground_dicts_list = []
    for tilt in tilts:
        csv_data_diffuse = pd.read_csv('./sim_data/diffuse/{}_output.csv'.format(tilt))  # read output data
        csv_data_ground = pd.read_csv('./sim_data/ground/{}_output.csv'.format(tilt))  # read output data
        dict_diffuse = csv_data_diffuse.to_dict()  # dict that store a single input csv file
        dict_ground = csv_data_ground.to_dict()  # dict that store a single input csv file
        diffuse_dicts_list.append(dict_diffuse)
        ground_dicts_list.append(dict_ground)

    l_opt_effs_diffuse = []  # output files data list               visual: [file1, file2, file3, ....]
    l_spect_mismatchs_diffuse = []  # output files data list
    l_opt_effs_ground = []  # output files data list               visual: [file1, file2, file3, ....]
    l_spect_mismatchs_ground = []  # output files data list
    for i in range(n_sample):
        l_opt_effs_diffuse.append([])
        l_opt_effs_ground.append([])
        l_spect_mismatchs_diffuse.append([])
        l_spect_mismatchs_ground.append([])

    for i in range(n_sample):
        for dict in diffuse_dicts_list:
            l_opt_effs_diffuse[i].append([dict['tilt'][i], dict['opt_eff'][i]])
            l_spect_mismatchs_diffuse[i].append([dict['tilt'][i], dict['spect_mismatch'][i]])
        for dict in ground_dicts_list:
            l_opt_effs_ground[i].append([dict['tilt'][i], dict['opt_eff'][i]])
            l_spect_mismatchs_ground[i].append([dict['tilt'][i], dict['spect_mismatch'][i]])


    # Do writing
    fields_diffuse_ground = ['tilt', 'lsc_objects']
    for i in range(len(l_opt_effs_diffuse)):
        fname_opt_eff_diffuse = './reformatted_data/iso/opt_eff/opt_eff_{}.csv'.format(i)
        fname_opt_eff_ground = './reformatted_data/grnd/opt_eff/opt_eff_{}.csv'.format(i)
        fname_spect_mismatch_diffuse = './reformatted_data/iso/spect_mismatch/spect_mismatch_{}.csv'.format(i)
        fname_spect_mismatch_ground = './reformatted_data/grnd/spect_mismatch/spect_mismatch_{}.csv'.format(i)

        writeCSV(fname_opt_eff_diffuse, fields_diffuse_ground, l_opt_effs_diffuse[i])
        writeCSV(fname_opt_eff_ground, fields_diffuse_ground, l_opt_effs_ground[i])
        writeCSV(fname_spect_mismatch_diffuse, fields_diffuse_ground, l_spect_mismatchs_diffuse[i])
        writeCSV(fname_spect_mismatch_ground, fields_diffuse_ground, l_spect_mismatchs_ground[i])
    print("Finish writing csv files for diffuse & ground lights")

