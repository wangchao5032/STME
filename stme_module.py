import math
import sys

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, convex_hull_plot_2d

def Dist(data, i, j, tWindow):
    """
    计算两点之间的距离（暂时没用）
    :param data: 数据集
    :param i: 数据集中的第i个对象
    :param j: 数据集中的第j个对象
    :param tWindow: 时间窗口
    :return: 数据集中的第i个对象和第j个对象之间的距离
    """
    dist_i_j = float("inf")
    if i == j:
        return dist_i_j
    # 提取第i和j个对象的时间、经度、维度
    t1, t2 = data.values[i][0], data.values[j][0]
    x1, x2 = data.values[i][1], data.values[j][1]
    y1, y2 = data.values[i][2], data.values[j][2]

    if abs(t1 - t2) <= tWindow:
        dist_i_j = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

    return dist_i_j


def GetKN(data, k, tWindow):
    """
    计算所有点的k近邻距离，以及K近邻
    :param data: 数据集
    :param k: k
    :param tWindow: 时间窗口
    :return:
    """
    rows = data.shape[0]  # 行数
    types = max(data['type']) + 1  # 数据集中点的类型数

    dist = np.zeros((rows, rows))  # 距离矩阵
    dist_topK = np.zeros((rows, types+1))  # 所有点的k近邻矩阵：第0列是不考虑类型的k距离，第1列是第1种类型的k距离，以此类推
    topK_Neighbor = np.zeros((rows, k))  # 所有点的k近邻编号

    # 计算距离矩阵 (Method 1 ---too slow)
    # dist2 = np.zeros((rows, rows))  # 距离矩阵
    # for i in range(rows):
    #     if i % 100 == 0:
    #         print("正在处理第{0}个轨迹点".format(i))
    #     for j in range(i, rows):
    #         dist2[i][j] = Dist(data, i, j, tWindow)
    #         dist2[j][i] = dist2[i][j]

    # 计算距离矩阵 (Method 2)
    for i in range(rows):
        if i % 100 == 0:
            print("正在处理第{0}个轨迹点".format(i))
        validData = (np.abs(data.values[:, 0] - data.values[i, 0]) <= tWindow)
        invalidData = np.logical_not(validData)
        dist[i, validData] = np.linalg.norm((data.values[validData, 1:3] - data.values[i, 1:3]), axis=1)
        dist[i, invalidData] = float("inf")
        dist[i, i] = float("inf")

    # 查找每个点的k近邻及距离
    for i in range(rows):
        dist_i_sort = sorted(dist[i])  # 排序后的第i个对象的k近邻距离
        dist_i_sort_id = np.argsort(dist[i])  # 排序后的第i个对象的k近邻距离 对应的 数据编号
        # 包含所有类型的，第k个最近的距离及编号
        dist_topK[i][0] = dist_i_sort[k-1]  # 包含所有类型的第k个最近的距离
        topK_Neighbor[i] = dist_i_sort_id[:k]  # 包含所有类型的第k个最近的距离 对应的 数据编号
        # 查找每个类型的，第k个最近的距离
        for ty in range(types):  # ty类型
            # 计算多类型k距离(Method 1 ---too slow)
            # temp = []  # 存放所有ty类型的数据和第i个对象的距离
            # for j in range(rows):
            #     if data.iloc[j].at['type'] == ty:
            #         temp.append(dist[i][j])
            # 计算多类型k距离(Method 2)
            selectedData = (data.values[:, 3] == ty)
            temp2 = dist[i, selectedData].tolist()

            dist_topK[i][ty+1] = sorted(temp2)[k-1]  # 对所有ty类型的数据和第i个对象的距离排序，查找第k个
    topK_Neighbor = topK_Neighbor.astype(int)  # 把编号转为int类型
    return dist_topK, topK_Neighbor


def GetRN(topK_Neighbor, kt):
    """
    计算每个点的共享近邻
    :param topK_Neighbor:所有点的k近邻编号
    :param kt:kt
    :return: 每个点的共享近邻列表
    """
    rows = topK_Neighbor.shape[0]  # 数据集中点的数目
    rstn = []  # 每个点的共享近邻列表   [[] for _ in range(rows)]
    for i in range(rows):
        rstn_i = []  # 第i个点的共享近邻列表
        neighbor_i = topK_Neighbor[i]  # 第i个点的k近邻
        for j in neighbor_i:  # 第i个点的j邻居
            neighbor_j = topK_Neighbor[int(j)]  # 第j个点的k近邻
            neighbor_i_j = list(set(neighbor_i).intersection(set(neighbor_j)))  # 第i个点和第j个点的k近邻交集
            # 如果i和j的共享近邻个数大于kt，且i也是j的k近邻
            if len(neighbor_i_j) >= kt and (i in neighbor_j):
                rstn_i.append(j)
        rstn.append(rstn_i)
    return rstn


def STME(data, k, kt, t_window, min_pts, distK_sigma_multi, distK_sigma_times, spatiotemoral=False):
    rows = data.shape[0]  # 行数
    types = max(data['type']) + 1  # 点的类型数
    print("共{0}个轨迹点,{1}种类型".format(rows, types))

    # 1. 计算每个点的k近邻距离，以及最近的K个邻居
    print("正在计算k近邻距离...")
    dist_topK, topK_Neighbor = GetKN(data, k, t_window)


    # 2. 计算聚类优先级（每个点第k个邻居的升序排序后的距离）
    print("正在计算聚类优先级...")
    dist_k = dist_topK[:, 0]  # 所有数据的k近邻距离
    dist_k_sort = sorted(dist_k)
    priority = np.argsort(dist_k)


    # 3. 计算每个点的共享近邻
    print("正在计算共享近邻...")
    RSTN = GetRN(topK_Neighbor, kt)

    # 4. 聚类
    clusterLabel = 0  # 初始化簇标号
    clusterDensityList = [[] for _ in range(types+1)]  # 簇有效密度列表：第0列是不考虑类型的有效密度，第1列是第1种类型的有效密度，以此类推
    clusterNumList = [[] for _ in range(types+1)]  # 簇中点数：第0列存的是簇中点的总数，第1列是簇中第1种类型点的数目，以此类推
    clusterCoreList = []  # 簇的第一个核心点列表
    labelKList = []  # 按照标记的顺序记录被标记点的k距离
    clusterDensityList_nor = []  # 簇的归一化混合密度
    print("开始聚类...")
    for p in priority:  # 按照第k个邻居距离的升序进行枚举，也就是第i个点一定是目前所有点里k邻域范围最小的
        # 判断第p个点是否已经被处理过了
        if data.iloc[p].at['cluster'] != 0:
            continue

        # 判断第p个点是否是噪声
        r_neighbor = RSTN[p]  # 第p个点的共享邻居 / 直接可达点
        if len(r_neighbor) < min_pts:  # 如果第p个点是噪声，不是核心点
            data.iloc[p, 4] = -1
            continue

        # 判断共享近邻中free_pts的数目
        RSTN_free_pts = data.iloc[r_neighbor].query('cluster==0 | cluster==-1')  # 查找第p个点的共享近邻中标签为0或-1的点
        tmp = data.iloc[r_neighbor]['cluster'].values.tolist()
        RSTN_most_label = max(tmp, key=tmp.count)  # 第p个点的共享近邻中出现最多的类型
        if len(RSTN_free_pts) < min_pts and RSTN_most_label != 0:  # 第p个点的共享近邻中标签为0或-1的点的数目小于MinPts，且出现最多的类型不是0
            data.iloc[p, 4] = RSTN_most_label
            continue

        # 第p个点是核心点，开始一个新的簇！
        clusterLabel += 1
        print("cluster:{0}".format(clusterLabel))
        clusterValidList = []  # 参与计算当前簇有效体积的集合，存放当前簇的所有核心点。和type_buf_core的作用一样，其实可以删除的。
        clusterCoreList.append(p)  # 当前簇的第一个核心点
        # 构造初始队列
        queue = [p] + r_neighbor  # 将第p个点及其RSTN加入队列
        distK_buf = dist_topK[queue]  # 队列中所有点的考虑类型的k距离
        type_buf_init = data.iloc[queue]["type"].values.tolist()  # 簇的有效点1:簇的初始队列中的所有点
        type_buf_core = []  # 簇的有 效点2: 后续加入队列的点中的核心点

        # 计算初始队列中所有点的均值和标准差
        distK_mu = np.average(distK_buf[:, 1:], axis=0)  # 计算“队列中所有点的考虑类型的k距离”中，点的 均值
        distK_sigma_init = np.std(distK_buf[:, 1:], axis=0) # 计算“队列中所有点的考虑类型的k距离”中，点的 标准差
        # 遍历队列
        while len(queue)>0:
            #  (1) 统计当前对列中各类型点的数目
            type_buf = type_buf_init + type_buf_core  # 簇的有效点：初始队列中所有点+后续加入队列的点中的核心点
            type_hist = []  # 统计当前队列中有效点的所有类型
            for ty in range(types):
                type_hist.append(type_buf.count(ty))
            # 当前队列的有效点中，如果某种类型的点数目大于MinPts，则认为这些类型是有影响力的
            type_valid = [0] * types
            for ty in range(types):
                if type_hist[ty] >= min_pts:
                    type_valid[ty] = 1
            type_valid = [True if type_valid[i] == 1 else False for i in range(len(type_valid))]
            # 找出当前队列的有效点中，影响力最大的类型
            dominate_type = type_hist.index(max(type_hist))
            distK_sigma = distK_sigma_init.copy()
            distK_sigma[dominate_type] = distK_sigma[dominate_type] * distK_sigma_multi  # 影响力最大的类型可以放宽k距离的阈值（为了使模拟数据边界密度较低的点包含进来）
            # (2) 取出当前点并标记
            ptCurrent = queue[0]  # 取出队列中的第一个点作为当前点
            queue.remove(queue[0])  # 把当前点移出队列
            if data.iloc[ptCurrent].at['cluster'] > 0:  #如果当前点被标记过，这个点只可能是作为共享近邻被加入队列了
                continue
            # 队列中的点有以下情况:
            # - 初始簇形成时依赖的第1个点，一定是没有被标记也不是噪声的核心点。
            # - 初始簇形成时依赖的第1个点的共享近邻，不一定是核心点，可能被标记过。会在后文进一步筛选是否是核心点（if length(X2) >= MinPts），如果是再将其k邻居加入队列
            # - 后期加入队列中的点，一定是没有被标记过的点或噪声，不一定是核心点。会在后文进一步筛选是否是核心点（if length(X2) >= MinPts），如果是再将其k邻居加入队列
            data.iloc[ptCurrent, 4] = clusterLabel
            labelKList.append(dist_topK[ptCurrent][0])

            # (3) 遍历当前点的k邻居，加入队列
            topK_Neighbor_pt = topK_Neighbor[ptCurrent]  # 队列中当前点的k个邻居
            RSTN_pt = RSTN[ptCurrent]  # 队列中当前点的共享邻居/直接可达点
            if len(RSTN_pt) >= min_pts:  # 如果当前点是核心点
                type_buf_core.append(data.iloc[ptCurrent].at['type'])
                clusterValidList.append(ptCurrent)
                for neigh in topK_Neighbor_pt:
                    # 计算当前点第j个邻居的考虑类型的k距离 和 初始簇考虑类型的k距离均值 的差值
                    distK_diff = abs(dist_topK[:, 1:][neigh] - distK_mu)
                    # 如果当前点的第j个邻居:a) 是噪声或还没有被标记过; b)不在对列中;
                    # c) 在只考虑有影响力的类型下，当前点第j个邻居的考虑类型的k距离（1*2） 和 初始簇考虑类型的k距离均值（1*2）的差值distK_diff(type_valid)（1*2），若差值在初始簇的3倍标准差以内
                    if (data.iloc[neigh].at['cluster'] == -1 or data.iloc[neigh].at['cluster'] == 0) \
                            and (not (neigh in queue)) \
                            and all(distK_diff[type_valid] < distK_sigma[type_valid] * distK_sigma_times):
                        queue.append(neigh)
        # (4) 计算当前簇的点数
        cluster = data[data['cluster'] == clusterLabel]  # 当前簇
        clusterNumList[0].append(len(cluster))  # 当前簇的点的总数
        for i in range(types):
            clusterNumList[i+1].append(len(cluster[cluster['type'] == i]))

        # (5) 计算当前簇有效密度
        if spatiotemoral == False or all(cluster['t'].values == np.zeros(len(cluster))):  # 如果计算的是二维凸包或时间维度一列都是0
            hull = ConvexHull(cluster[['x', 'y']].values)  # 当前簇的最小凸包
        else:
            hull = ConvexHull(cluster[['t', 'x', 'y']].values)
        V = hull.volume
        clusterDensityList[0].append(clusterNumList[0][0] / V)  # 当前簇的密度=当前簇的数量/当前簇的体积
        for i in range(types):  # 当前簇的第i类点的密度 = 当前簇中该类点的数量/当前簇的体积
            clusterDensityList[i+1].append(clusterNumList[i+1][0] / V)

    # 5. 计算簇的归一化混合密度
    print("计算簇的归一化混合密度")
    nClusters = len(clusterDensityList[0])
    if nClusters > 1:
        for i in range(nClusters):
            tmp = (clusterDensityList[0][i] - min(clusterDensityList[0])) / (max(clusterDensityList[0]) - min(clusterDensityList[0]))
            clusterDensityList_nor.append(tmp)
    else:
        clusterDensityList_nor.append(1)

    return data, clusterDensityList_nor
