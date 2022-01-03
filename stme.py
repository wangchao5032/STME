import argparse
import os
import pandas as pd
from stme_module import STME

parser = argparse.ArgumentParser(description="stme.py")

parser.add_argument("-k", type=int, default=20,
    help="Number of neighbors searched, determine the size of the spatiotemporal neighbourhood of each event")

parser.add_argument("-tWindow", type=int, default=5,
    help="Time window, determine the size of the spatiotemporal neighbourhood of each event")

parser.add_argument("-kt", type=int, default=10,
    help="Number of shared neighbors,determine whether Oj is a refined neighbour of Oi")

parser.add_argument("-MinPts", type=int, default=10,
    help="Number of shared neighbors,used to identify the core event")

parser.add_argument("-distK_sigma_multi", type=int, default=100,
    help="点数最多的类型可以放宽k距离的阈值")

parser.add_argument("-distK_sigma_times", type=int, default=3,
    help="")

parser.add_argument("-data", default="data",
    help="Path to clustering data")

parser.add_argument("-sheet", default="Sheet1",
    help="sheet name")

parser.add_argument("-name",  default="data",
    help="file name")



args = parser.parse_args()


if __name__ == "__main__":
    # 读取数据
    data_path = "{0}/{1}.xlsx".format(args.data, args.name)
    if not os.path.exists(data_path):
        print("can not find data")
        exit()

    data = pd.read_excel(data_path, engine='openpyxl', sheet_name=args.sheet)

    # 执行stme算法
    data, clusterDensityList_nor = STME(data, args.k, args.kt, args.tWindow, args.MinPts, args.distK_sigma_multi, args.distK_sigma_times)

    # 保存数据
    out_path = "{0}/output.xlsx".format(args.data)
    data.to_excel(out_path, header=False, index=False)
    exit()
