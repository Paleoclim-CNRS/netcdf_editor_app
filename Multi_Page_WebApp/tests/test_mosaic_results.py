# coding: utf-8

import os
import xarray as xr
import numpy

import sys, getopt


def compare_xarrays(ds1, ds2, name):
    print(ds1.equals(ds2), "\t", name)
    if not ds1.equals(ds2):
        print("testing all close")
        nb_diff = 0
        for var in ds1.data_vars:
            try:
                if not numpy.allclose(ds1[var], ds2[var], 1e-15):
                    nb_diff += 1
            except ValueError as e:
                print(e)
                nb_diff += 1
        print(nb_diff < 1, "\t", name)
    print("\n")


def test_files(base_data_dir, docker_data_dir):

    base_data_name = os.listdir(os.path.join(base_data_dir, "ATM", "START"))[-1]
    docker_data_name = os.listdir(os.path.join(docker_data_dir, "ATM", "START"))[-1]

    filenames = os.listdir(
        os.path.join(base_data_dir, "CPL", "IPSLCM5A2", base_data_name)
    )

    for name in filenames:
        if name.split(".")[-1] != "nc":
            continue
        base = os.path.join(base_data_dir, "CPL", "IPSLCM5A2", base_data_name, name)
        docker = os.path.join(
            docker_data_dir, "CPL", "IPSLCM5A2", docker_data_name, name
        )
        ds_base = xr.open_dataset(base)
        ds_docker = xr.open_dataset(docker)
        compare_xarrays(ds_docker, ds_base, name)

    filenames = os.listdir(os.path.join(base_data_dir, "ATM", "START", base_data_name))

    for name in filenames:
        if name.split(".")[-1] != "nc":
            continue
        base = os.path.join(base_data_dir, "ATM", "START", base_data_name, name)
        docker = os.path.join(docker_data_dir, "ATM", "START", docker_data_name, name)
        ds_base = xr.open_dataset(base)
        ds_docker = xr.open_dataset(docker)
        compare_xarrays(ds_docker, ds_base, name)


def main(argv):
    base_data_dir = ""
    docker_data_dir = ""
    # We have 2 options
    print(len(argv))
    if len(argv) == 2:
        # Make sure Folder exists and it is a Folder
        for i in range(2):
            assert os.path.exists(argv[i]), f"{argv[i]} Does not exist"
            assert os.path.isdir(argv[i]), f"{argv[i]} is not a Folder"
        base_data_dir = argv[0]
        docker_data_dir = argv[1]
    try:
        opts, args = getopt.getopt(argv, "ho:n:", ["origdir=", "newdir="])
    except getopt.GetoptError:
        print("test_mosaic_results.py -o <originaldir> -n <newdir>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("test_mosaic_results.py -o <originaldir> -n <newdir>")
            sys.exit()
        elif opt in ("-o", "--origdir"):
            base_data_dir = arg
        elif opt in ("-n", "--newdir"):
            docker_data_dir = arg
    print('base_data_dir is "', base_data_dir)
    print('docker_data_dir file is "', docker_data_dir)
    if len(base_data_dir) and len(docker_data_dir):
        test_files(base_data_dir, docker_data_dir)
    else:
        print("No folders specified")


if __name__ == "__main__":
    main(sys.argv[1:])
