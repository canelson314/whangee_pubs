#!/bin/python
import os, sys
import argparse
from Bio import SeqIO
import cPickle as pic
import gzip
import numpy as np
from scipy import stats
#import seqmatch
import math
import re
import matplotlib.pyplot as plt
import pandas as pd
#import seaborn as sns
#from sklearn.neighbors.kde import KernelDensity
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
from scipy import cluster
import operator


class MidpointNormalize(Normalize):
    # Defines the midpoint of diverging colormap
    # Adjusts the midpoint of a color bar

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # Ignore masked values and edge cases to make a simple example...
        # Find the midpoint
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def make_indices_dict():
    indices_dict ={}
    indexFile = '/data/whangee/data/index.txt'
    #makes dict of indecies pointing to timepoints 
    with open(indexFile,"rU") as indices_handle:
        indices_handle.readline()
        for line in indices_handle:
            linelist = (line.strip()).split("\t")
            timepoint = linelist[0]
            seq_index = linelist[6]
            generations =linelist[9]
            indices_dict[seq_index] = [int(linelist[2].strip()),int(linelist[3].strip()),timepoint, generations]
            #need to add timepoint to index file- maybe they should be the ODs
            #tp_dict[timepoint] = [seq_index, int(linelist[2].strip()),int(linelist[3].strip())]
            #time_dict[(linelist[2].strip()+linelist[3].strip())] = [int(linelist[7].strip()), int(linelist[8].strip())]
    pic.dump(indices_dict, open('indices_dict.pkl', 'w'))
    return indices_dict

def get_aa_to_codon_dict(translate):
    aa_to_codon_dict = {}
    for codon, aa in translate.items():
        if aa not in aa_to_codon_dict:
            aa_to_codon_dict[aa] = []
        aa_to_codon_dict[aa].append(codon)
    return aa_to_codon_dict

def get_timepoint_dict(indices_dict):
    tp_dict = {}
    # group num-> [ 0.  ,  1.85,  4.1 ,  0.  ,  1.97,  3.45]
    for index, indexList in indices_dict.items():
        group = indices_dict[index][0]
        experiment = indices_dict[index][1]
        tp = indices_dict[index][3]
        if group not in tp_dict:
            tp_dict[group] = np.zeros(6)
        tp_dict[group][experiment] = float(tp)
    return tp_dict

def get_wt_bc_and_counts_lists(translate, mut_to_BC_dict, wtCodonList):
    wtBarcodeList = []
    wtCountList = []
    for mut, bc_list in mut_to_BC_dict.items():
        pos= int("".join(re.findall('[0-9]*' ,mut)))
        codon = "".join(re.findall('[ACGU]{3}|WT+' ,mut))
        for bc in bc_list:
            if (codon == 'WT') or (translate[codon] == translate[wtCodonList[pos-1]]):
                wtBarcodeList.append(bc)
                sum_bc = 0
                for group in xrange(0,6):
                    sum_bc+=sum(barcode_to_count_dict[bc][3][group])
                wtCountList.append(sum_bc)
    #plt.hist((sorted(np.array(wtCountList), reverse=True))[30:], bins= 30, color='#3399FF', edgecolor = '#2E8AE6')
    #plt.show()
    return wtBarcodeList, wtCountList

def get_wt_matrix(barcode_to_count_dict, wtBarcodeList):
    wt_count_matrix = np.zeros((7,6))
    for group in xrange(0,7):
        group_tmp_wt_lists = [[],[],[],[],[],[]]
        for wt_bc in wtBarcodeList:
            for tp in xrange(0,6):
                if barcode_to_count_dict[wt_bc][3][group][tp] > 0:
                    group_tmp_wt_lists[tp].append(barcode_to_count_dict[wt_bc][3][group][tp])
        for tp in xrange(0,6):
            if group == 6 and tp>2:
                continue
            else:
                wt_count_matrix[group][tp] = sum(group_tmp_wt_lists[tp])/len(group_tmp_wt_lists[tp])
    return wt_count_matrix

def get_slope_dict(barcode_to_count_dict, wt_count_matrix, tp_dict):
    barcode_to_slope_dict ={}
    #bc_list = [[],[],[],[],[],[]]
    #difference_list = [[],[],[],[],[],[]]
    count = 0
    std_err_dict = {0:[],1:[]}
    for bc_key, bc_val in barcode_to_count_dict.items():
        barcode = bc_key
        barcode_to_slope_dict[barcode] = [bc_val[0], bc_val[1], bc_val[2], np.zeros((7,2))]
        for group in xrange(0, 7):
            day1_norm_count = []
            day2_norm_count = []
            day1_counts = bc_val[3][group][:3]
            day2_counts = bc_val[3][group][3:]
            day1_wt = wt_count_matrix[group][:3]
            day2_wt = wt_count_matrix[group][3:]
            time1 = []
            time2 = []
            for tp in xrange(0,3):
                if day1_counts[tp] != 0:
                    day1_norm_count.append(math.log((day1_counts[tp]/day1_wt[tp]),2))
                    time1.append(tp_dict[group][tp])
                if day2_counts[tp] != 0:
                    day2_norm_count.append(math.log((day2_counts[tp]/day2_wt[tp]),2))
                    time2.append(tp_dict[group][tp+3])
            if group < 6:
                if (len(day1_norm_count) > 1) and (len(day2_norm_count) > 1):
                    slopeOverall_1, interceptOverall_1, r_valOverall_1, p_valOverall_1, std_errOverall_1 = stats.linregress(time1, day1_norm_count)
                    slopeOverall_2, interceptOverall_2, r_valOverall_2, p_valOverall_2, std_errOverall_2 = stats.linregress(time2, day2_norm_count)
                    std1 = "".join(re.findall('inf|nan', str(std_errOverall_1)))
                    std2 = "".join(re.findall('inf|nan', str(std_errOverall_2)))
                    """
                    if std1 != 'inf' and std1 !='nan':
                        std_err_dict[0].append(std_errOverall_1)
                    if std2 != 'inf' and std2 != 'nan':
                        std_err_dict[1].append(std_errOverall_2)
                    if group < 6:
                        bc_list[group].append(barcode)
                        difference_list[group].append(((slopeOverall_1 - slopeOverall_2)**2)/slopeOverall_2)
                """
                barcode_to_slope_dict[barcode][3][group] = [slopeOverall_1, slopeOverall_2]
                """
                if (std_errOverall_1 < 0.4) and (std_errOverall_2 < 0.4):
                    count+=1
                    barcode_to_slope_dict[barcode][3][group] = [slopeOverall_1, slopeOverall_2]
                else:
                    barcode_to_slope_dict[barcode][3][group] = [float('NaN'), float('NaN')]
                """
            else:
                if len(day1_norm_count) > 1:
                    barcode_to_slope_dict[barcode][3][group] = [slopeOverall_1, slopeOverall_2]


    pic.dump(barcode_to_slope_dict, open('slope_no_filter.pkl', 'w'))
    return barcode_to_slope_dict

def get_codon_to_slope_dict(pos_codon_to_bar_dict, barcode_to_slope_dict):
    pos_codon_to_slope_dict ={}
    for pos_codon, barcodeList in pos_codon_to_bar_dict.items():
        slopeMatrix1 = np.zeros((7,len(barcodeList)))
        slopeMatrix2 = np.zeros((7,len(barcodeList)))
        count_list = np.zeros(7)
        avgSlopeMatrix = np.zeros((7,2))
        slopeList1 = []
        slopeList2 = []
        notInDICT = []
        for n in xrange(0, len(barcodeList)):
            for group in xrange(0,7):
                barcode = barcodeList[n]
                if barcode in barcode_to_slope_dict:
                    slope1 = barcode_to_slope_dict[barcode][3][group][0]
                    slope2 = barcode_to_slope_dict[barcode][3][group][1]
                    inf_nan_1 = "".join(re.findall('inf|nan', str(slope1)))
                    inf_nan_2 = "".join(re.findall('inf|nan', str(slope2)))
                    #if str(slope1)!=inf_nan_1 and str(slope2)!=inf_nan_2:
                    slopeMatrix1[group][n] = barcode_to_slope_dict[barcode][3][group][0]
                    slopeMatrix2[group][n] = barcode_to_slope_dict[barcode][3][group][1]
                    count_list[group]+=1
                else:
                    notInDICT.append(barcode)
        for group in xrange(0,7):
            avgSlopeMatrix[group][0] = np.average(np.ma.array(slopeMatrix1[group]))
            avgSlopeMatrix[group][1] = np.average(np.ma.array(slopeMatrix2[group]))
        posCodon = "".join(re.findall('[0-9]*' ,pos_codon) ) + "_"+"".join(re.findall('[ACGU]{3}|WT+' ,pos_codon))
        pos_codon_to_slope_dict[posCodon] = avgSlopeMatrix
    pic.dump(pos_codon_to_slope_dict, open('codon_pos_to_fit_no_filter.pkl', 'w'))
    return pos_codon_to_slope_dict

def get_aa_to_slope_dict(pos_codon_to_slope_dict, translate):
    aa_to_fit_dict = {}
    aa_count_dict = {}
    for pos_codon in pos_codon_to_slope_dict.items():
        (pos, codon) = pos_codon[0].split("_")
        #print pos, codon
        aa = ''
        if codon == 'WT':
            aa = 'WT'
        else:
            aa = translate[codon]
        pos_aa = pos + '_'+ aa
        if pos_aa not in aa_to_fit_dict:
            aa_to_fit_dict[pos_aa] = np.zeros((7,2))
            aa_count_dict[pos_aa] = np.zeros(7)
        for group in xrange(0,7):
            #print aa_to_fit_dict[pos_aa][group][0], pos_codon_to_slope_dict[pos_codon][group][0]
            slope1 = pos_codon_to_slope_dict[pos_codon[0]][group][0]
            slope2 = pos_codon_to_slope_dict[pos_codon[0]][group][1]
            inf_nan_1 = "".join(re.findall('inf|nan', str(slope1)))
            inf_nan_2 = "".join(re.findall('inf|nan', str(slope2)))
            #if str(slope1)!=inf_nan_1 and str(slope2)!=inf_nan_2:
            aa_to_fit_dict[pos_aa][group][0]+= pos_codon_to_slope_dict[pos_codon[0]][group][0]
            aa_to_fit_dict[pos_aa][group][1]+= pos_codon_to_slope_dict[pos_codon[0]][group][1]
            aa_count_dict[pos_aa][group]+=1
    for aa , aa_matrix in aa_to_fit_dict.items():
        for group in xrange(0,7):
            aa_to_fit_dict[aa][group][0]= np.average(np.ma.array(aa_to_fit_dict[aa][group][0]))
            aa_to_fit_dict[aa][group][1]= np.average(np.ma.array(aa_to_fit_dict[aa][group][1]))

    pic.dump(aa_to_fit_dict, open('aa_pos_to_fit_no_filter.pkl', 'w'))
    return aa_to_fit_dict

def avg_aa_day1_2(pos_aa_to_slope_dict):
    avg_aa_dict = {}
    for aa, aa_matrix in pos_aa_to_slope_dict.items():
        avg_aa_dict[aa] = np.zeros(7)
        for group in xrange(0,7):
            if group == 6:
                avg_aa_dict[aa][group]=np.average(aa_matrix[group][0])
            else:
                avg_aa_dict[aa][group]=np.average(np.ma.array(aa_matrix[group]))
    for aa, aa_val in avg_aa_dict.items():
        print aa, aa_val[3]
    pic.dump(avg_aa_dict, open('avg_aa_fit_no_filter.pkl', 'w'))
    return avg_aa_dict
    


def print_bc_scatter(bc_to_fit):
    for group in xrange(0,6):
        fit1 = []
        fit2 = []
        for bc in bc_to_fit.items():
            fit1.append(bc[1][3][group][0])
            fit2.append(bc[1][3][group][1])
        plt.scatter(fit1, fit2)
        plt.show()  

def print_codon_or_aa_scatter(residue_to_fit):
    for group in xrange(0,6):
        fit1 = []
        fit2 = []
        for residue in residue_to_fit.items():
            fit1.append(residue[1][group][0])
            fit2.append(residue[1][group][1])
        plt.scatter(fit1, fit2)
        plt.show()  

def residue_dict_to_df(residue_slope_dict):
    df1 = pd.DataFrame()
    temp_residue_dict1 = {}
    temp_residue_dict2 = {}
    residue_list = []
    for group in xrange(0,6):
        temp_residue_dict1[group] = []
        temp_residue_dict2[group] = []
    for group in xrange(0,6):
        res_list = []
        for residue, residue_value in residue_slope_dict.items():
            temp_residue_dict1[group].append(residue_value[0])
            temp_residue_dict2[group].append(residue_value[1])
            res_list.append(residue)
        residue_list = res_list
        tempDF2 = pd.DataFrame(temp_residue_dict1, index=[str(group)+'_day1'], columns = residue_list)
        df1 = pd.concat((df1, tempDF2))
        tempDF2 = pd.DataFrame(temp_residue_dict2, index=[str(group)+'_day2'], columns = residue_list)
        df1 = pd.concat((df1, tempDF2))
    print df1
    return df1

def print_distribution_of_all(dictionary):
    groupList = ['onion', 'apc', 'shmoo', 'whangee', 'pynd', 'etoh', 'control']
    dict_by_group = {}
    f = plt.figure()
    for group in xrange(0,7):
        dict_by_group[group] = []
    for dict_key, dict_value in dictionary.items():
        for group in xrange(0,7):
            if str(dict_value[group]) == 'nan' or str(dict_value[group]) == 'inf':
                dict_by_group[group].append(0)
            else:
                dict_by_group[group].append(dict_value[group])
    for group in xrange(0,7):
        sns.distplot(dict_by_group[group], hist=False, label = groupList[group])
        plt.legend(bbox_to_anchor=(1.3, 0.5), loc="center right")
    sns.despine()
    sns.axlabel("Amino Acid Fitness", "Density")
    plt.tight_layout()
    aa_string = 'all_group_aa_fit_dist.png'
    plt.savefig(aa_string, dpi=300)
    plt.show()
    

def pert_vs_mut(aa_dict, group):
    groupList = ['Onion', 'Apc', 'Shmoo', 'Whangee', 'Pynd', 'Etoh', 'Control']
    wt_aa = aa_dict['0_WT'][group]
    wt_mut_diff_list = []
    pert_control_diff_list = []
    print wt_aa
    for aa, aa_fit in aa_dict.items():
        pert_fitness = aa_fit[group]
        control_fitness = aa_fit[6]
        if str(pert_fitness)!='nan' and str(control_fitness)!='nan':
            wt_mut_diff_list.append(pert_fitness-wt_aa)
            pert_control_diff_list.append(pert_fitness-control_fitness)
    plt.axis([min(wt_mut_diff_list)-1, max(wt_mut_diff_list)+1, min(pert_control_diff_list)-1, max(pert_control_diff_list)+1])
    plt.plot([0,0],[min(pert_control_diff_list)-1,max(pert_control_diff_list)+1], color='#C0C0C0')
    plt.plot([min(wt_mut_diff_list)-1, max(wt_mut_diff_list)+1], [0,0], color='#C0C0C0')
    plt.title(groupList[group] + ' Impact of Mutation vs Perturbation')
    sns.axlabel("Impact of Mutation", "Impact of Perturbation")
    plt.scatter(wt_mut_diff_list, pert_control_diff_list, color='#48D1CC')
    name_string = groupList[group] + 'MutationvsPerturbation.png'
    plt.savefig(name_string, dpi=300)
    plt.show()

def aa_by_pos_matrix(codon_fit, avg_aa_dict, number):
    #make aa by pos matrix
    avg_aa_matrix = np.zeros((21, 78))
    vari_matrix = np.zeros((21, 78))
    clean_aa_dict = {}
    clean_codon_dict = {}
    norm = MidpointNormalize(midpoint=0)
    for pos_aa, aa_fit in avg_aa_dict.items():
        position = int(pos_aa.split('_')[0])
        aa = pos_aa.split('_')[1]
        #vari = np.var(np.ma.array(aa_fit))
        all_list = []

        """
        if pos_aa == '0_STOP':
            print pos_aa, aa_fit
        """
        
        if aa != 'WT':
            avg_aa_matrix[number[aa]][position] = aa_fit[3]

            codon_list_fit = []
            new_fit_list = []
            for codon in aa_to_codon_dict[aa]:
                pos_codon = str(position) +'_'+ codon
                if pos_codon in codon_fit:
                    codon_list_fit.append(codon_fit[pos_codon][3][0])
                    codon_list_fit.append(codon_fit[pos_codon][3][1])
            if len(codon_list_fit) > 2:
                mean = np.mean(np.ma.masked_array(codon_list_fit))
                var = np.var(np.ma.masked_array(codon_list_fit))
                negSigma = mean - (var**2)
                posSigma = mean + (var**2)
                med = np.median(np.ma.masked_array(codon_list_fit))
                q75, q25 = np.percentile(np.ma.masked_array(codon_list_fit),[75, 25])
                for fit in codon_list_fit:
                    if fit > q25 and fit < q75:
                        new_fit_list.append(fit)    
            else:
                new_fit_list = codon_list_fit
            vari_matrix[number[aa]][position] = np.mean(new_fit_list)
            clean_aa_dict[pos_aa] = np.mean(new_fit_list)
            clean_aa_dict[pos_aa] = np.mean(new_fit_list)
            if len(new_fit_list) == 0:
                vari_matrix[number[aa]][position] = 0


        """
        else:
            new_fit_list = []
            avg_aa_matrix[0][position] = aa_fit[3]
            med = np.median(np.ma.masked_array(codon_fit['0_WT']))
            q75, q25 = np.percentile(np.ma.masked_array(codon_fit['0_WT']),[75, 25])
            for fit in codon_list_fit:
                if fit > q25 and fit < q75:
                    new_fit_list.append(fit)
            vari_matrix[0][position] = np.mean(new_fit_list)
        """
    #pic.dump(clean_codon_to_fit_dict, open('aa_pos_to_fit_no_filter.pkl', 'w'))
    return avg_aa_matrix, vari_matrix, clean_aa_dict

#make aa by aa matrix

def make_aa_by_aa_matrix(avg_aa_dict, number):
    aa_by_aa_matrix = np.zeros((21,21))
    aa_dict = {}
    ubseq = list('MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG*')
    row_name_aa = []
    for aa, keys in sorted(number.items(), key=operator.itemgetter(1)):
        row_name_aa.append(aa)

    for pos_aa, aa_fit in avg_aa_dict.items():
        position = int(pos_aa.split('_')[0])
        aa = pos_aa.split('_')[1]
        if aa != 'WT' and aa != ubseq[position-1]:
            if aa not in aa_dict:
                aa_dict[aa] = np.zeros(78)
            if type(aa_fit) != np.float64:
                aa_dict[aa][position] = aa_fit[3]
            else:
                aa_dict[aa][position] = aa_fit
    print 'AA,'+','.join(row_name_aa)
    row_index = 0
    for aa1, aa_num1 in number.items():
        aa_row_list = []
        for aa2, aa_num2 in number.items():
            aa_list1 = aa_dict[aa1]
            aa_list2 = aa_dict[aa2]
            new_aa_list1 = []
            new_aa_list2 = []
            for pos in xrange(0, len(aa_list1)):
                if aa_list1[pos]!=0 and aa_list2[pos]!=0:
                    new_aa_list1.append(aa_list1[pos])
                    new_aa_list2.append(aa_list2[pos])
            slopeOverall_1, interceptOverall_1, r_valOverall_1, p_valOverall_1, std_errOverall_1 = stats.linregress(new_aa_list1, new_aa_list2)
            aa_by_aa_matrix[number[aa1]][number[aa2]]=r_valOverall_1
            aa_row_list.append(str(r_valOverall_1))
        print ','.join(aa_row_list)


    plt.matshow(aa_by_aa_matrix)
    plt.show()
    return aa_by_aa_matrix

def print_matrix(row, col, row_label, col_label, matrix_1, matrix_2, number, same_row_col):
    f, ax = plt.subplots(2, sharex = True)
    ax[0] = plt.subplot2grid((2,2), (0,0), colspan = 2)
    ax[1] = plt.subplot2grid((2,2), (1,0), colspan = 2)
        
    # Assigns color according to values in fitness_heat arrays
    norm = MidpointNormalize(midpoint=0)
    #im1 = ax[0].imshow(matrix_1, cmap=plt.cm.seismic, interpolation='none')
    im1 = ax[0].imshow(matrix_1, norm=norm, cmap=plt.cm.seismic, interpolation='none')
    divider1 = make_axes_locatable(ax[0])
    cax1 = divider1.append_axes("right", size="2%", pad=0.05)

    im2 = ax[1].imshow(matrix_2, cmap=plt.cm.seismic, interpolation='none')
    #im2 = ax[1].imshow(matrix_2, norm=norm, cmap=plt.cm.seismic, interpolation='none')
    divider2 = make_axes_locatable(ax[1])
    cax2 = divider2.append_axes("right", size="2%", pad=0.05)

    # Sorting out labels for y-axis in correct order and placing in new list 'yax'
    yax1 = sorted(number.items(), key=lambda (k, v): v)
    yax = []
    for item in yax1:
        if item[0] =='STOP':
            yax.append ('*')
        else:
            yax.append (item[0])
    
    if same_row_col == True:
        xax = yax
        x = np.arange(0, row+0.5, 1)
        ax[0].set_xticks(x)
        ax[0].set_xticklabels(xax, rotation='horizontal', fontsize = 10)
        ax[0].set_xlabel(row_label)
        ax[1].set_xticks(x)
        ax[1].set_xticklabels(xax, rotation='horizontal', fontsize = 10)
        ax[1].set_xlabel(row_label)
    else:
        x = np.arange(0.5, col+0.5, 5)
        xlab = np.arange(1, col, 5)
        plt.xticks(x, xlab, rotation='horizontal')
        ax[1].set_xlabel(col_label)
        
    # Setting up y-axis ticks...want AA labels, centered
    y = np.arange(0, row, 1)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels(yax, rotation='horizontal', fontsize = 10)
    ax[0].set_ylabel(row_label)
    ax[1].set_yticks(y)
    ax[1].set_yticklabels(yax, rotation='horizontal', fontsize = 10)
    ax[1].set_ylabel(row_label)

    # Setting up x-axis ticks...want residue number, centered
    

    # Limit axes to appropriate values
    plt.axis([0, col, 0, row])

    plt.colorbar(im1, cax = cax1)
    plt.colorbar(im2, cax = cax2)

    plt.show()

def print_csv(aa_dict, translate, number, aa_type_dict):
    groupList = ['Onion', 'Apc', 'Shmoo', 'Whangee', 'Pynd', 'Etoh', 'Control']
    """
    print ','.join(['pos_aa', 'Group','AA', 'AA Type', 'Position',  'Fitness'])
    indv_group_dict = {}
    for pos_aa, aa_fit in aa_dict.items():
        pos, aa = pos_aa.split('_')
        for group in xrange(0, 7):
            g_list = [pos_aa, groupList[group], aa, aa_type_dict[aa],  pos,  str(aa_fit[group])]
            print ','.join(g_list)
   
    print ','.join(['pos_aa', 'Position', 'Codon', 'AA', 'Day 1 Fitness', 'Day 2 Fitness'])
    for pos_aa, aa_fit in aa_dict.items():
        pos, codon = pos_aa.split('_')
        aa = ''
        if codon != 'WT':
            aa = translate[codon]
        else:
            aa = 'WT'

        g_list = [pos_aa, pos, codon, aa, str(aa_fit[3][0]), str(aa_fit[3][0])]
        print ','.join(g_list)
    """
    print ','.join(['pos_aa', 'AA', 'AA Type', 'Position',  'Fitness'])
    for pos_aa, aa_fit in aa_dict.items():
        pos, aa = pos_aa.split('_')
        g_list = [pos_aa, aa, aa_type_dict[aa],  pos,  str(aa_fit[3])]
        print ','.join(g_list)

def clust():
    df= pd.DataFrame.from_csv('clean_groups.csv')
    clusters = sns.clustermap(df, columns= groupList)
    plt.show()

    

def make_aa_type_dict():
    aa_type_dict = {}
    with open('aa_type.txt', 'r') as aa_type_file:
        for lines in aa_type_file:
            aa, aa_type = (lines.strip()).split('\t')
            aa_type_dict[aa] = aa_type
    return aa_type_dict





wtseq = 'ATGCAGATTTTCGTCAAGACTTTGACCGGTAAAACCATAACATTGGAAGTTGAATCTTCCGATACCATCGACAACGTTAAGTCGAAAATTCAAGACAAGGAAGGTATCCCTCCAGATCAACAAAGATTGATCTTTGCCGGTAAGCAGCTAGAAGACGGTAGAACGCTGTCTGATTACAACATTCAGAAGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGTTGA'
wtRNAseq = wtseq.replace('T', 'U')
#77 wt codons
wtCodonList = [ wtRNAseq[i:i+3] for i in range(0, len(wtRNAseq), 3) ]
stopBarcodeList = []
translate = pic.load(open("translate.pkl", "rb"))
aa_to_codon_dict = get_aa_to_codon_dict(translate)
aa_type_dict = make_aa_type_dict()
#load other dicts
"""
allele = pic.load(open("/data/whangee/data/allele_dic_with_WT.pkl", "rb"))



barcode_to_count_dict = pic.load(open("/data/whangee/data/clean_all_count_dict.pkl.pkl", "rb"))
mut_to_BC_dict = pic.load(open('/data/whangee/data/clean_all_mut_BC_dict.pkl', 'rb'))
indices_dict = pic.load(open("/data/whangee/data/indices_dict.pkl", "rb"))


tp_dict = get_timepoint_dict(indices_dict)
wtBarcodeList, wtCountList = get_wt_bc_and_counts_lists(translate, mut_to_BC_dict, wtCodonList)
wt_count_matrix = get_wt_matrix(barcode_to_count_dict, wtBarcodeList)
"""
number = pic.load(open("aminotonumber.pkl", "rb"))
#barcode_to_slope_dict = get_slope_dict(barcode_to_count_dict, wt_count_matrix, tp_dict)
#barcode_to_slope_dict = pic.load(open('slope_no_filter.pkl', 'rb'))
#print_bc_scatter(barcode_to_slope_dict)

#pos_codon_to_slope_dict = get_codon_to_slope_dict(mut_to_BC_dict, barcode_to_slope_dict)
pos_codon_to_slope_dict = pic.load(open('codon_pos_to_fit_no_filter.pkl', 'rb'))
#print_codon_or_aa_scatter(pos_codon_to_slope_dict)

#pos_aa_to_slope_dict = get_aa_to_slope_dict(pos_codon_to_slope_dict, translate)
pos_aa_to_slope_dict = pic.load(open('aa_pos_to_fit_no_filter.pkl', 'rb'))
#print_codon_or_aa_scatter(pos_aa_to_slope_dict)

#avg_aa_dict = avg_aa_day1_2(pos_aa_to_slope_dict)
avg_aa_dict = pic.load(open('avg_aa_fit_no_filter.pkl', 'rb'))
#avg_aa_matrix, vari_matrix, clean_aa_dict = aa_by_pos_matrix(pos_codon_to_slope_dict, avg_aa_dict, number)


"""
pert_vs_mut(avg_aa_dict, 3)
#print_distribution_of_all(avg_aa_dict)

aa_to_codon_dict = get_aa_to_codon_dict(translate)

aa_by_aa_matrix = make_aa_by_aa_matrix(avg_aa_dict, number)
clean_aa_by_aa_matrix = make_aa_by_aa_matrix(clean_aa_dict, number)
"""

avg_aa_matrix, vari_matrix, clean_aa_dict = aa_by_pos_matrix(pos_codon_to_slope_dict, avg_aa_dict, number)

row_name_aa = []
for aa, keys in sorted(number.items(), key=operator.itemgetter(1)):
    row_name_aa.append(aa)
pos_list = xrange(1,79)
#print 'AA,'+','.join(pos_list)

rosetta_dict = {}
with open('uby_1ubq.csv.txt', 'r') as rosetta:
    for lines in rosetta:
        lines = lines.strip()
        linelist = lines.split(',')
        aa_pos = linelist[1][-3:-1] + '_' + linelist[1][-1]

        ddg = linelist[2]
        rank = linelist[3]
        rosetta_dict[aa_pos] = [ddg, rank]

print ','.join(['Position', 'AA', 'Fitness', 'DDG', 'Rank'])
for pos_aa, fit in clean_aa_dict.items():
    if pos_aa in rosetta_dict:
        pos, aa = pos_aa.split('_')
        row_list = [pos, aa, str(fit), rosetta_dict[pos_aa][0], rosetta_dict[pos_aa][1]]
        print ','.join(row_list)




"""
for row in xrange(0,21):
    aa_row_list = [row_name_aa[row]]
    for col in xrange(0,78):
        aa_row_list.append(avg_aa_dict[row][col])
"""






"""
row_name_aa = []
for aa, keys in sorted(number.items(), key=operator.itemgetter(1)):
    row_name_aa.append(aa)
aa_by_aa_matrix = make_aa_by_aa_matrix(avg_aa_dict, number)
print 'AA,' + ','.join(row_name_aa)
for row_num in xrange(0, len(row_name_aa)):
    row_list = []
    for col in xrange(0, len(row_name_aa)):
        row_list.append(str(aa_by_aa_matrix[row_num][col]))
    print row_name_aa[row_num] +','+ ','.join(row_list)

print_matrix(21, 21, "Amino Acid", "Amino Acid", aa_by_aa_matrix, aa_by_aa_matrix, number, True)
#print_csv(avg_aa_dict, translate, number, aa_type_dict)
"""
