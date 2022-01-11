from math import floor
import numpy as np
import os
import pandas as pd
from pyDOE import lhs
import threading
import xmltodict

import dask.delayed
from dask.distributed import Client
from dask_jobqueue import SLURMCluster

from run import wedge, wedge2


"""
    sim_driver.py
    Read parameter sweep ranges and cluster parameters
    Run simulations to record neff, spectral mismatch for each parameter combo and incidence type
    Save results to csv
"""

def LHS_sample(lsc_params, lsc_consts, sim_params):
    # Reduced set of params with LHS
    combined_keys = list(lsc_params.keys()) + list(lsc_consts.keys())
    samples = pd.DataFrame(lhs(len(combined_keys), int(sim_params["nsamples"])),
                           columns=combined_keys)
    for key in lsc_params.keys():
        samples[key] = samples[key] * \
                       (float(lsc_params[key]["ubound"]) - float(lsc_params[key]["lbound"])) + \
                       float(lsc_params[key]["lbound"])
    for key in lsc_consts.keys():
        samples[key] = float(lsc_consts[key]["const"])
    return samples


def parse_range(range_str, floatcast=False):
    if floatcast:
        return [float(val) for val in range_str.split(',')]
    return range_str.split(',')


def phi_theta_sim(phis, thetas, light_f, samples_dict, n_bundles, lsc_func):
    results = []
    labels = []
    for p in phis:
        for t in thetas:
            for sample in samples_dict:
                labels.append({'phi': p, 'theta': t})
                lsc_outputs = dask.delayed(lsc_func)(n_bundles,
                                                  light_form=light_f,
                                                  theta_o=t,
                                                  phi_o=p,
                                                  input_params=sample)
                results.append(lsc_outputs)
    futures = dask.persist(*results)
    results = dask.compute(*futures)
    results_to_csv(results, light_f, labels)
    return results


def tilt_sim(l_tilts, light_f, samples_dict, n_bundles, lsc_func):
    results = []
    labels = []
    for t in l_tilts:
        for sample in samples_dict:
            labels.append({'tilt': t})
            lsc_outputs = dask.delayed(lsc_func)(n_bundles,
                                              light_form=light_f,
                                              tilt=t,
                                              input_params=sample)
            results.append(lsc_outputs)
    futures = dask.persist(*results)
    results = dask.compute(*futures)
    results_to_csv(results, light_f, labels)
    return


def results_to_csv(results, light_f, labels):
    # Generate data directory for this light form
    if not os.path.isdir("./sim_data"):
        os.mkdir("./sim_data")
    if not os.path.isdir("./sim_data/{}".format(light_f)):
        os.mkdir("./sim_data/{}".format(light_f))

    # These shouldn't be hardcoded, but it should be find as long as simulation return values don't change
    columns = list(results[0]['parameters'].keys())
    columns = columns + ['opt_eff', 'spect_mismatch']
    results_df = pd.DataFrame(columns=columns)

    # Assemble full results df
    label_count = 0
    for result in results:
        row = result['parameters']
        row['opt_eff'] = result['opt_eff']
        row['spect_mismatch'] = result['spect_mismatch']
        row.update(labels[label_count])
        print(row)
        results_df = results_df.append(row, ignore_index=True)
        label_count += 1

    # For each insolation label, filter results df and save to csv
    for label in labels:
        label_df = results_df
        label_str = ""
        for label_key in list(label.keys()):
            label_df = label_df.loc[(label_df[label_key] == label[label_key])]
            label_str = label_str + str(label[label_key]) + "_"
        label_df.to_csv('./sim_data/{}/{}output.csv'.format(light_f, label_str), index=False)
    return


if __name__ == "__main__":
    # Read simulation configs
    fconfig = input("Enter simulation config (IE: params_small.xml):")
    conf_s = "".join([line.rstrip('\n') for line in open(fconfig)])
    sim_params = xmltodict.parse(conf_s)["params"]["sim_params"]
    insolation_params = xmltodict.parse(conf_s)["params"]["insolation_params"]
    lsc_params = xmltodict.parse(conf_s)["params"]["lsc_params"]
    lsc_consts = xmltodict.parse(conf_s)["params"]["lsc_consts"]

    # Read the desired LSC geometry
    lsc_func = None
    if sim_params["geometry"] == "wedge_mirror":
        lsc_func = wedge
    elif sim_params["geometry"] == "wedge_no_mirror":
        lsc_func = wedge2
    else:
        print("Unknown geometry. Exiting")
        exit(1)


    # Unpack light forms
    light_forms = parse_range(sim_params["light_forms"])

    # Unpack phi/theta/tilt ranges
    phi_v = parse_range(insolation_params["direct_ranges"]["phi"], floatcast=True)
    theta_v = parse_range(insolation_params["direct_ranges"]["theta"], floatcast=True)
    tilt_v = parse_range(insolation_params["diffuse_ranges"]["tilt"], floatcast=True)

    # Set up Slurm cluster context
    cluster = None
    client = None
    if sim_params['hostcluster'] == 'slurm':
        cluster = SLURMCluster(cores=int(sim_params["cores"]),
                               processes=int(sim_params["processes"]),
                               memory="512GiB",
                               project="SLSC",
                               walltime="06:00:00",
                               job_extra=['--gres=gpu:4'])
        cluster.scale(int(sim_params["clusterscale"]))
        client = Client(cluster)
    else:
        # Or a local cluster
        client = Client()

    # Reduced set of param combos with LHS
    samples = LHS_sample(lsc_params, lsc_consts, sim_params)
    samples_dict = samples.to_dict(orient='records')

    threads = []
    # We don't really need to thread here and it seems like an antipattern with parallelism, but hear me out:
    #   Threading allows different light forms to produce outputs when they're ready, instead of waiting for the slowest
    #       light form to finish.
    #       Diffuse and Ground sims require far fewer incidence trials since they depend on tilt only, so their data
    #       comes out sooner.
    #   Threading saves us from having to filter an aggregate dataframe for light form.
    #   Threading also happens to make phi/theta/tilt labeling a little easier.
    if "direct" in light_forms:
        t = threading.Thread(target=phi_theta_sim,
                             args=[phi_v, theta_v, "direct", samples_dict, int(sim_params["nbundles"]), lsc_func])
        threads.append(t)
    if "diffuse" in light_forms:
        t = threading.Thread(target=tilt_sim,
                             args=[tilt_v, "diffuse", samples_dict, int(sim_params["nbundles"]), lsc_func])
        threads.append(t)
    if "ground" in light_forms:
        t = threading.Thread(target=tilt_sim,
                             args=[tilt_v, "ground", samples_dict, int(sim_params["nbundles"]), lsc_func])
        threads.append(t)

    # Kick off, then destroy threads once they finish
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
