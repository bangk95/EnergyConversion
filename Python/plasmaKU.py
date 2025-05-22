import numpy as np
import pandas as pd
import os
os.environ["CDF_LIB"] = "/home/local/cdf/lib/"
from spacepy import pycdf
from datetime import datetime
import glob


def clweb(file, mode, var):
    file = str(file)
    mode = str(mode)
    cdffile0 = pycdf.CDF(file)
    epoch = cdffile0['epoch'][...]

    if var=='E':
        ex = cdffile0['edp_dce_gse_ex_'+mode+'_l2'][...]
        ey = cdffile0['edp_dce_gse_ey_'+mode+'_l2'][...]
        ez = cdffile0['edp_dce_gse_ez_'+mode+'_l2'][...]
        return ex,ey,ez,epoch

    if var=='B':
        bx = cdffile0['fgm_bx_gse_brst_l2'][...]
        by = cdffile0['fgm_by_gse_brst_l2'][...]
        bz = cdffile0['fgm_bz_gse_brst_l2'][...]
        bt = cdffile0['fgm_bt_gse_brst_l2'][...]
        return bx,by,bz,bt,epoch
    
    if var=='ui':
        uix = cdffile0['dis_bulkvx_gse_'+mode][...]
        uiy = cdffile0['dis_bulkvy_gse_'+mode][...]
        uiz = cdffile0['dis_bulkvz_gse_'+mode][...]
        return uix,uiy,uiz,epoch
    
    if var=='pi':
        pixx = cdffile0['dis_prestensorxx_gse_'+mode][...]
        pixy = cdffile0['dis_prestensorxy_gse_'+mode][...]
        pixz = cdffile0['dis_prestensorxz_gse_'+mode][...]
        piyx = cdffile0['dis_prestensoryx_gse_'+mode][...]
        piyy = cdffile0['dis_prestensoryy_gse_'+mode][...]
        piyz = cdffile0['dis_prestensoryz_gse_'+mode][...]
        pizx = cdffile0['dis_prestensorzx_gse_'+mode][...]
        pizy = cdffile0['dis_prestensorzy_gse_'+mode][...]
        pizz = cdffile0['dis_prestensorzz_gse_'+mode][...]
        return pixx,pixy,pixz,piyx,piyy,piyz,pizx,pizy,pizz,epoch
        
    if var=='ni':
        ni = cdffile0['dis_numberdensity_'+mode][...]
        return ni,epoch
        
    if var=='ti':
        tipara = cdffile0['dis_temppara_'+mode][...]
        tiperp = cdffile0['dis_tempperp_'+mode][...]
        return tipara,tiperp,epoch

    if var=='ue':
        uex = cdffile0['des_bulkvx_gse_'+mode][...]
        uey = cdffile0['des_bulkvy_gse_'+mode][...]
        uez = cdffile0['des_bulkvz_gse_'+mode][...]
        return uex,uey,uez,epoch
        
    if var=='pe':
        pexx = cdffile0['des_prestensorxx_gse_'+mode][...]
        pexy = cdffile0['des_prestensorxy_gse_'+mode][...]
        pexz = cdffile0['des_prestensorxz_gse_'+mode][...]
        peyx = cdffile0['des_prestensoryx_gse_'+mode][...]
        peyy = cdffile0['des_prestensoryy_gse_'+mode][...]
        peyz = cdffile0['des_prestensoryz_gse_'+mode][...]
        pezx = cdffile0['des_prestensorzx_gse_'+mode][...]
        pezy = cdffile0['des_prestensorzy_gse_'+mode][...]
        pezz = cdffile0['des_prestensorzz_gse_'+mode][...]
        return pexx,pexy,pexz,peyx,peyy,peyz,pezx,pezy,pezz,epoch

    if var=='te':
        tepara = cdffile0['des_temppara_'+mode][...]
        teperp = cdffile0['des_tempperp_'+mode][...]
        return tepara,teperp,epoch

    if var=='ne':
        ne = cdffile0['des_numberdensity_'+mode][...]
        return ne,epoch
        
    if var=='j':
        jx = cdffile0['jx_gse_A_m2'][...]
        jy = cdffile0['jy_gse_A_m2'][...]
        jz = cdffile0['jz_gse_A_m2'][...]
        jt = cdffile0['j_A_m2'][...]
        return jx,jy,jz,jt,epoch
        
    if var=='jcurl':
        jcurlx = cdffile0['jx_gse_A_m2'][...]
        jcurly = cdffile0['jy_gse_A_m2'][...]
        jcurlz = cdffile0['jz_gse_A_m2'][...]
        jcurlt = cdffile0['j_A_m2'][...]
        return jcurlx,jcurly,jcurlz,jcurlt,epoch



def lasp(mode, n, tstr, tend, manual=False, yok=False):
    tstr = str(tstr)
    tend = str(tend)
    
    start_Y = tstr[:4]
    start_M = tstr[4:6]
    start_D = tstr[6:8]
    start_H = tstr[9:11]
    start_m = tstr[11:13]
    start_s = tstr[13:]
    
    end_Y = tend[:4]
    end_M = tend[4:6]
    end_D = tend[6:8]
    end_H = tend[9:11]
    end_m = tend[11:13]
    end_s = tend[13:]
    
    print('Start : ', start_M+'/'+start_D+'/'+start_Y, start_H+':'+start_m+':'+start_s)
    print('Stop  : ', end_M+'/'+end_D+'/'+end_Y, end_H+':'+end_m+':'+end_s)
    
    mode = str(mode)

    n = str(n)
    mms = 'mms'+n
    E_nan_path = ''

    if manual:
        start_Y = input('Input start year (YYYY) : ')
        start_M = input('Input start month(MM) : ')
        start_D = input('Input start date (DD) : ')
        start_H = input('Input start hour (HH) : ')
        start_m = input('Input start minute (MM) : ')
        start_s = input('Input start second (SS) : ')

        end_Y = start_Y
        end_M = start_M
        end_D = start_D
        end_H = input('Input end hour (HH) : ')
        end_m = input('Input end minute (MM) : ')
        end_s = input('Input end second (SS) : ')

    start_time1 = int(start_Y+start_M+start_D+start_H+start_m+start_s)
    end_time1 = int(end_Y+end_M+end_D+end_H+end_m+end_s)
    start_time2 = int(start_Y+start_M+start_D)
    end_time2 = int(end_Y+end_M+end_D)
    
    interval_file1 = range(start_time1, end_time1+1)
    interval_file1 = np.array(interval_file1)
    interval_file1 = list(map(str,interval_file1))
    
    interval_file2 = range(start_time2, end_time2+1)
    interval_file2 = np.array(interval_file2)
    interval_file2 = list(map(str,interval_file2))
    
    #print(interval_file1)
    #print(interval_file2)

    fileE = [] #EDP
    fileB = [] #FGM
    filePi = [] #FPI-i
    filePe = [] #FPI-e
    
    if yok:
        path_add = '../../yok/MMS/'
    else:
        path_add = "/Users/thanaponaiamsai/Desktop/Plasma/"
    
    if mode in ['fast'] or mode in ['srvy']:
        modeb = 'srvy'
        for elem2 in interval_file2:
            for file1 in glob.glob(path_add + "Data/LASP/mms/data/" + mms + "/edp/fast/l2/dce/" + start_Y + "/" + start_M + "/" + mms + "_edp_fast_l2_dce_" + elem2 + "*.cdf"):
                #print(file1)
                fileE.append(file1)
            for file2 in glob.glob(path_add + "Data/LASP/mms/data/" + mms + "/fgm/srvy/l2/" + start_Y + "/" + start_M + "/" + mms + "_fgm_srvy_l2_" + elem2 + "*.cdf"):
                #print(file2)
                fileB.append(file2)
        for elem1 in interval_file1:
            for file3 in glob.glob(path_add + "Data/LASP/mms/data/" + mms + "/fpi/fast/l2/" + "dis-moms" + "/" + start_Y + "/" + start_M  + "/" + mms + "_fpi_fast_l2_" + "dis-moms" + "_" + elem1 + "*.cdf"):
                #print(file3)
                filePi.append(file3)
            for file4 in glob.glob(path_add + "Data/LASP/mms/data/" + mms + "/fpi/fast/l2/" + "des-moms" + "/" + start_Y + "/" + start_M  + "/" + mms + "_fpi_fast_l2_" + "des-moms" + "_" + elem1 + "*.cdf"):
                #print(file4)
                filePe.append(file4)
                
    if mode in ['brst']:
        modeb = 'brst'
        for elem1 in interval_file1:
            for file1 in glob.glob(path_add + "Data/LASP/mms/data/" + mms + "/edp/brst/l2/dce/" + start_Y + "/" + start_M + "/" + start_D + "/" + mms + "_edp_brst_l2_dce_" + elem1 + "*.cdf"):
                #print(file1)
                fileE.append(file1)
            for file2 in glob.glob(path_add + "Data/LASP/mms/data/" + mms + "/fgm/brst/l2/" + start_Y + "/" + start_M + "/" + start_D + "/" + mms + "_fgm_brst_l2_" + elem1 + "*.cdf"):
                #print(file2)
                fileB.append(file2)
            for file3 in glob.glob(path_add + "Data/LASP/mms/data/" + mms + "/fpi/brst/l2/" + "dis-moms" + "/" + start_Y + "/" + start_M  + "/" + start_D  + "/" +  mms + "_fpi_brst_l2_" + "dis-moms" + "_" + elem1 + "*.cdf"):
                #print(file3)
                filePi.append(file3)
            for file4 in glob.glob(path_add + "Data/LASP/mms/data/" + mms + "/fpi/brst/l2/" + "des-moms" + "/" + start_Y + "/" + start_M  + "/" + start_D  + "/" +  mms + "_fpi_brst_l2_" + "des-moms" + "_" + elem1 + "*.cdf"):
                #print(file4)
                filePe.append(file4)

    print('EDP files:', len(fileE))
    print('FGM files:', len(fileB))
    print('FPI-i files:', len(filePi))
    print('FPI-e files:', len(filePe))
        
    cdffileE = []
    cdffileB = []
    cdffilePi = []
    cdffilePe = []
    for i in range(len(fileE)):
        cdffileE.append(pycdf.CDF(fileE[i]))
    for i in range(len(fileB)):
        cdffileB.append(pycdf.CDF(fileB[i]))
    for i in range(len(filePi)):
        cdffilePi.append(pycdf.CDF(filePi[i]))
        cdffilePe.append(pycdf.CDF(filePe[i]))

#     print(cdffileE[0], '\n')
#     print(cdffileB[0], '\n')
#     print(cdffilePi[0], '\n')
#     print(cdffilePe[0], '\n')

    timeE = []
    timeb = []
    timer = []
    timei = []
    timee = []
    
    x = []
    y = []
    z = []
    bx = []
    by = []
    bz = []
    
    ex = []
    ey = []
    ez = []
    
    edp_qual = []
    fgm_flag = []
    
    uix = []
    uiy = []
    uiz = []
    uex = []
    uey = []
    uez = []
    uix_err = []
    uiy_err = []
    uiz_err = []
    uex_err = []
    uey_err = []
    uez_err = []
    uix_spt = []
    uiy_spt = []
    uiz_spt = []
    uex_spt = []
    uey_spt = []
    uez_spt = []

    pixx = []
    pixy = []
    pixz = []
    piyx = []
    piyy = []
    piyz = []
    pizx = []
    pizy = []
    pizz = []
    pexx = []
    pexy = []
    pexz = []
    peyx = []
    peyy = []
    peyz = []
    pezx = []
    pezy = []
    pezz = []
    
    energy_i01 = []
    energy_i02 = []
    energy_i03 = []
    energy_i04 = []
    energy_i05 = []
    energy_i06 = []
    energy_i07 = []
    energy_i08 = []
    energy_i09 = []
    energy_i10 = []
    energy_i11 = []
    energy_i12 = []
    energy_i13 = []
    energy_i14 = []
    energy_i15 = []
    energy_i16 = []
    energy_i17 = []
    energy_i18 = []
    energy_i19 = []
    energy_i20 = []
    energy_i21 = []
    energy_i22 = []
    energy_i23 = []
    energy_i24 = []
    energy_i25 = []
    energy_i26 = []
    energy_i27 = []
    energy_i28 = []
    energy_i29 = []
    energy_i30 = []
    energy_i31 = []
    energy_i32 = []
    
    energy_e01 = []
    energy_e02 = []
    energy_e03 = []
    energy_e04 = []
    energy_e05 = []
    energy_e06 = []
    energy_e07 = []
    energy_e08 = []
    energy_e09 = []
    energy_e10 = []
    energy_e11 = []
    energy_e12 = []
    energy_e13 = []
    energy_e14 = []
    energy_e15 = []
    energy_e16 = []
    energy_e17 = []
    energy_e18 = []
    energy_e19 = []
    energy_e20 = []
    energy_e21 = []
    energy_e22 = []
    energy_e23 = []
    energy_e24 = []
    energy_e25 = []
    energy_e26 = []
    energy_e27 = []
    energy_e28 = []
    energy_e29 = []
    energy_e30 = []
    energy_e31 = []
    energy_e32 = []
    
    pixx_err = []
    pixy_err = []
    pixz_err = []
    piyx_err = []
    piyy_err = []
    piyz_err = []
    pizx_err = []
    pizy_err = []
    pizz_err = []
    pexx_err = []
    pexy_err = []
    pexz_err = []
    peyx_err = []
    peyy_err = []
    peyz_err = []
    pezx_err = []
    pezy_err = []
    pezz_err = []
    
    err_flagi = []
    err_flage = []

    numdeni = []
    numdene = []
    numdeni_err = []
    numdene_err = []

    tempi_para = []
    tempi_perp = []
    tempe_para = []
    tempe_perp = []
    
    tixx = []
    tixy = []
    tixz = []
    tiyx = []
    tiyy = []
    tiyz = []
    tizx = []
    tizy = []
    tizz = []
    texx = []
    texy = []
    texz = []
    teyx = []
    teyy = []
    teyz = []
    tezx = []
    tezy = []
    tezz = []
    
    tixx_err = []
    tixy_err = []
    tixz_err = []
    tiyx_err = []
    tiyy_err = []
    tiyz_err = []
    tizx_err = []
    tizy_err = []
    tizz_err = []
    texx_err = []
    texy_err = []
    texz_err = []
    teyx_err = []
    teyy_err = []
    teyz_err = []
    tezx_err = []
    tezy_err = []
    tezz_err = []

    for i in range(len(filePi)):
        timei = np.append(timei, cdffilePi[i]["Epoch"][...])
        timee = np.append(timee, cdffilePe[i]["Epoch"][...])

        uix = np.append(uix, cdffilePi[i][mms + "_" + "dis" + "_bulkv_gse_" + mode][:,0])
        uiy = np.append(uiy, cdffilePi[i][mms + "_" + "dis" + "_bulkv_gse_" + mode][:,1])
        uiz = np.append(uiz, cdffilePi[i][mms + "_" + "dis" + "_bulkv_gse_" + mode][:,2])
        uex = np.append(uex, cdffilePe[i][mms + "_" + "des" + "_bulkv_gse_" + mode][:,0])
        uey = np.append(uey, cdffilePe[i][mms + "_" + "des" + "_bulkv_gse_" + mode][:,1])
        uez = np.append(uez, cdffilePe[i][mms + "_" + "des" + "_bulkv_gse_" + mode][:,2])
        
        uix_err = np.append(uix_err, cdffilePi[i][mms + "_" + "dis" + "_bulkv_err_" + mode][:,0])
        uiy_err = np.append(uiy_err, cdffilePi[i][mms + "_" + "dis" + "_bulkv_err_" + mode][:,1])
        uiz_err = np.append(uiz_err, cdffilePi[i][mms + "_" + "dis" + "_bulkv_err_" + mode][:,2])
        uex_err = np.append(uex_err, cdffilePe[i][mms + "_" + "des" + "_bulkv_err_" + mode][:,0])
        uey_err = np.append(uey_err, cdffilePe[i][mms + "_" + "des" + "_bulkv_err_" + mode][:,1])
        uez_err = np.append(uez_err, cdffilePe[i][mms + "_" + "des" + "_bulkv_err_" + mode][:,2])
        
        uix_spt = np.append(uix_spt, cdffilePi[i][mms + "_" + "dis" + "_bulkv_spintone_gse_" + mode][:,0])
        uiy_spt = np.append(uiy_spt, cdffilePi[i][mms + "_" + "dis" + "_bulkv_spintone_gse_" + mode][:,1])
        uiz_spt = np.append(uiz_spt, cdffilePi[i][mms + "_" + "dis" + "_bulkv_spintone_gse_" + mode][:,2])
        uex_spt = np.append(uex_spt, cdffilePe[i][mms + "_" + "des" + "_bulkv_spintone_gse_" + mode][:,0])
        uey_spt = np.append(uey_spt, cdffilePe[i][mms + "_" + "des" + "_bulkv_spintone_gse_" + mode][:,1])
        uez_spt = np.append(uez_spt, cdffilePe[i][mms + "_" + "des" + "_bulkv_spintone_gse_" + mode][:,2])
        
        err_flagi = np.append(err_flagi, cdffilePi[i][mms + "_" + "dis" + "_errorflags_" + mode][:])
        err_flage = np.append(err_flage, cdffilePe[i][mms + "_" + "des" + "_errorflags_" + mode][:])
        
        energy_i01 = np.append(energy_i01, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,0])
        energy_i02 = np.append(energy_i02, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,1])
        energy_i03 = np.append(energy_i03, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,2])
        energy_i04 = np.append(energy_i04, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,3])
        energy_i05 = np.append(energy_i05, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,4])
        energy_i06 = np.append(energy_i06, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,5])
        energy_i07 = np.append(energy_i07, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,6])
        energy_i08 = np.append(energy_i08, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,7])
        energy_i09 = np.append(energy_i09, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,8])
        energy_i10 = np.append(energy_i10, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,9])
        energy_i11 = np.append(energy_i11, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,10])
        energy_i12 = np.append(energy_i12, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,11])
        energy_i13 = np.append(energy_i13, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,12])
        energy_i14 = np.append(energy_i14, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,13])
        energy_i15 = np.append(energy_i15, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,14])
        energy_i16 = np.append(energy_i16, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,15])
        energy_i17 = np.append(energy_i17, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,16])
        energy_i18 = np.append(energy_i18, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,17])
        energy_i19 = np.append(energy_i19, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,18])
        energy_i20 = np.append(energy_i20, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,19])
        energy_i21 = np.append(energy_i21, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,20])
        energy_i22 = np.append(energy_i22, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,21])
        energy_i23 = np.append(energy_i23, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,22])
        energy_i24 = np.append(energy_i24, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,23])
        energy_i25 = np.append(energy_i25, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,24])
        energy_i26 = np.append(energy_i26, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,25])
        energy_i27 = np.append(energy_i27, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,26])
        energy_i28 = np.append(energy_i28, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,27])
        energy_i29 = np.append(energy_i29, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,28])
        energy_i30 = np.append(energy_i30, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,29])
        energy_i31 = np.append(energy_i31, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,30])
        energy_i32 = np.append(energy_i32, cdffilePi[i][mms + "_" + "dis" + "_energy_" + mode][:,31])
        
        energy_e01 = np.append(energy_e01, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,0])
        energy_e02 = np.append(energy_e02, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,1])
        energy_e03 = np.append(energy_e03, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,2])
        energy_e04 = np.append(energy_e04, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,3])
        energy_e05 = np.append(energy_e05, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,4])
        energy_e06 = np.append(energy_e06, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,5])
        energy_e07 = np.append(energy_e07, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,6])
        energy_e08 = np.append(energy_e08, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,7])
        energy_e09 = np.append(energy_e09, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,8])
        energy_e10 = np.append(energy_e10, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,9])
        energy_e11 = np.append(energy_e11, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,10])
        energy_e12 = np.append(energy_e12, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,11])
        energy_e13 = np.append(energy_e13, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,12])
        energy_e14 = np.append(energy_e14, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,13])
        energy_e15 = np.append(energy_e15, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,14])
        energy_e16 = np.append(energy_e16, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,15])
        energy_e17 = np.append(energy_e17, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,16])
        energy_e18 = np.append(energy_e18, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,17])
        energy_e19 = np.append(energy_e19, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,18])
        energy_e20 = np.append(energy_e20, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,19])
        energy_e21 = np.append(energy_e21, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,20])
        energy_e22 = np.append(energy_e22, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,21])
        energy_e23 = np.append(energy_e23, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,22])
        energy_e24 = np.append(energy_e24, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,23])
        energy_e25 = np.append(energy_e25, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,24])
        energy_e26 = np.append(energy_e26, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,25])
        energy_e27 = np.append(energy_e27, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,26])
        energy_e28 = np.append(energy_e28, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,27])
        energy_e29 = np.append(energy_e29, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,28])
        energy_e30 = np.append(energy_e30, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,29])
        energy_e31 = np.append(energy_e31, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,30])
        energy_e32 = np.append(energy_e32, cdffilePe[i][mms + "_" + "des" + "_energy_" + mode][:,31])
         
        pixx = np.append(pixx, cdffilePi[i][mms + "_" + "dis" + "_prestensor_gse_" + mode][:,0,0])
        pixy = np.append(pixy, cdffilePi[i][mms + "_" + "dis" + "_prestensor_gse_" + mode][:,0,1])
        pixz = np.append(pixz, cdffilePi[i][mms + "_" + "dis" + "_prestensor_gse_" + mode][:,0,2])
        piyx = np.append(piyx, cdffilePi[i][mms + "_" + "dis" + "_prestensor_gse_" + mode][:,1,0])
        piyy = np.append(piyy, cdffilePi[i][mms + "_" + "dis" + "_prestensor_gse_" + mode][:,1,1])
        piyz = np.append(piyz, cdffilePi[i][mms + "_" + "dis" + "_prestensor_gse_" + mode][:,1,2])
        pizx = np.append(pizx, cdffilePi[i][mms + "_" + "dis" + "_prestensor_gse_" + mode][:,2,0])
        pizy = np.append(pizy, cdffilePi[i][mms + "_" + "dis" + "_prestensor_gse_" + mode][:,2,1])
        pizz = np.append(pizz, cdffilePi[i][mms + "_" + "dis" + "_prestensor_gse_" + mode][:,2,2])
        pexx = np.append(pexx, cdffilePe[i][mms + "_" + "des" + "_prestensor_gse_" + mode][:,0,0])
        pexy = np.append(pexy, cdffilePe[i][mms + "_" + "des" + "_prestensor_gse_" + mode][:,0,1])
        pexz = np.append(pexz, cdffilePe[i][mms + "_" + "des" + "_prestensor_gse_" + mode][:,0,2])
        peyx = np.append(peyx, cdffilePe[i][mms + "_" + "des" + "_prestensor_gse_" + mode][:,1,0])
        peyy = np.append(peyy, cdffilePe[i][mms + "_" + "des" + "_prestensor_gse_" + mode][:,1,1])
        peyz = np.append(peyz, cdffilePe[i][mms + "_" + "des" + "_prestensor_gse_" + mode][:,1,2])
        pezx = np.append(pezx, cdffilePe[i][mms + "_" + "des" + "_prestensor_gse_" + mode][:,2,0])
        pezy = np.append(pezy, cdffilePe[i][mms + "_" + "des" + "_prestensor_gse_" + mode][:,2,1])
        pezz = np.append(pezz, cdffilePe[i][mms + "_" + "des" + "_prestensor_gse_" + mode][:,2,2])
        
        pixx_err = np.append(pixx_err, cdffilePi[i][mms + "_" + "dis" + "_prestensor_err_" + mode][:,0,0])
        pixy_err = np.append(pixy_err, cdffilePi[i][mms + "_" + "dis" + "_prestensor_err_" + mode][:,0,1])
        pixz_err = np.append(pixz_err, cdffilePi[i][mms + "_" + "dis" + "_prestensor_err_" + mode][:,0,2])
        piyx_err = np.append(piyx_err, cdffilePi[i][mms + "_" + "dis" + "_prestensor_err_" + mode][:,1,0])
        piyy_err = np.append(piyy_err, cdffilePi[i][mms + "_" + "dis" + "_prestensor_err_" + mode][:,1,1])
        piyz_err = np.append(piyz_err, cdffilePi[i][mms + "_" + "dis" + "_prestensor_err_" + mode][:,1,2])
        pizx_err = np.append(pizx_err, cdffilePi[i][mms + "_" + "dis" + "_prestensor_err_" + mode][:,2,0])
        pizy_err = np.append(pizy_err, cdffilePi[i][mms + "_" + "dis" + "_prestensor_err_" + mode][:,2,1])
        pizz_err = np.append(pizz_err, cdffilePi[i][mms + "_" + "dis" + "_prestensor_err_" + mode][:,2,2])
        pexx_err = np.append(pexx_err, cdffilePe[i][mms + "_" + "des" + "_prestensor_err_" + mode][:,0,0])
        pexy_err = np.append(pexy_err, cdffilePe[i][mms + "_" + "des" + "_prestensor_err_" + mode][:,0,1])
        pexz_err = np.append(pexz_err, cdffilePe[i][mms + "_" + "des" + "_prestensor_err_" + mode][:,0,2])
        peyx_err = np.append(peyx_err, cdffilePe[i][mms + "_" + "des" + "_prestensor_err_" + mode][:,1,0])
        peyy_err = np.append(peyy_err, cdffilePe[i][mms + "_" + "des" + "_prestensor_err_" + mode][:,1,1])
        peyz_err = np.append(peyz_err, cdffilePe[i][mms + "_" + "des" + "_prestensor_err_" + mode][:,1,2])
        pezx_err = np.append(pezx_err, cdffilePe[i][mms + "_" + "des" + "_prestensor_err_" + mode][:,2,0])
        pezy_err = np.append(pezy_err, cdffilePe[i][mms + "_" + "des" + "_prestensor_err_" + mode][:,2,1])
        pezz_err = np.append(pezz_err, cdffilePe[i][mms + "_" + "des" + "_prestensor_err_" + mode][:,2,2])

        numdeni = np.append(numdeni, cdffilePi[i][mms + "_" + "dis" + "_numberdensity_" + mode][...])
        numdene = np.append(numdene, cdffilePe[i][mms + "_" + "des" + "_numberdensity_" + mode][...])
        numdeni_err = np.append(numdeni_err, cdffilePi[i][mms + "_" + "dis" + "_numberdensity_err_" + mode][...])
        numdene_err = np.append(numdene_err, cdffilePe[i][mms + "_" + "des" + "_numberdensity_err_" + mode][...])
        
        tempi_para = np.append(tempi_para, cdffilePi[i][mms + "_" + "dis" + "_temppara_" + mode][...])
        tempi_perp = np.append(tempi_perp, cdffilePi[i][mms + "_" + "dis" + "_tempperp_" + mode][...])
        tempe_para = np.append(tempe_para, cdffilePe[i][mms + "_" + "des" + "_temppara_" + mode][...])
        tempe_perp = np.append(tempe_perp, cdffilePe[i][mms + "_" + "des" + "_tempperp_" + mode][...])
        
        tixx = np.append(tixx, cdffilePi[i][mms + "_" + "dis" + "_temptensor_gse_" + mode][:,0,0])
        tixy = np.append(tixy, cdffilePi[i][mms + "_" + "dis" + "_temptensor_gse_" + mode][:,0,1])
        tixz = np.append(tixz, cdffilePi[i][mms + "_" + "dis" + "_temptensor_gse_" + mode][:,0,2])
        tiyx = np.append(tiyx, cdffilePi[i][mms + "_" + "dis" + "_temptensor_gse_" + mode][:,1,0])
        tiyy = np.append(tiyy, cdffilePi[i][mms + "_" + "dis" + "_temptensor_gse_" + mode][:,1,1])
        tiyz = np.append(tiyz, cdffilePi[i][mms + "_" + "dis" + "_temptensor_gse_" + mode][:,1,2])
        tizx = np.append(tizx, cdffilePi[i][mms + "_" + "dis" + "_temptensor_gse_" + mode][:,2,0])
        tizy = np.append(tizy, cdffilePi[i][mms + "_" + "dis" + "_temptensor_gse_" + mode][:,2,1])
        tizz = np.append(tizz, cdffilePi[i][mms + "_" + "dis" + "_temptensor_gse_" + mode][:,2,2])
        texx = np.append(texx, cdffilePe[i][mms + "_" + "des" + "_temptensor_gse_" + mode][:,0,0])
        texy = np.append(texy, cdffilePe[i][mms + "_" + "des" + "_temptensor_gse_" + mode][:,0,1])
        texz = np.append(texz, cdffilePe[i][mms + "_" + "des" + "_temptensor_gse_" + mode][:,0,2])
        teyx = np.append(teyx, cdffilePe[i][mms + "_" + "des" + "_temptensor_gse_" + mode][:,1,0])
        teyy = np.append(teyy, cdffilePe[i][mms + "_" + "des" + "_temptensor_gse_" + mode][:,1,1])
        teyz = np.append(teyz, cdffilePe[i][mms + "_" + "des" + "_temptensor_gse_" + mode][:,1,2])
        tezx = np.append(tezx, cdffilePe[i][mms + "_" + "des" + "_temptensor_gse_" + mode][:,2,0])
        tezy = np.append(tezy, cdffilePe[i][mms + "_" + "des" + "_temptensor_gse_" + mode][:,2,1])
        tezz = np.append(tezz, cdffilePe[i][mms + "_" + "des" + "_temptensor_gse_" + mode][:,2,2])
        
        tixx_err = np.append(tixx_err, cdffilePi[i][mms + "_" + "dis" + "_temptensor_err_" + mode][:,0,0])
        tixy_err = np.append(tixy_err, cdffilePi[i][mms + "_" + "dis" + "_temptensor_err_" + mode][:,0,1])
        tixz_err = np.append(tixz_err, cdffilePi[i][mms + "_" + "dis" + "_temptensor_err_" + mode][:,0,2])
        tiyx_err = np.append(tiyx_err, cdffilePi[i][mms + "_" + "dis" + "_temptensor_err_" + mode][:,1,0])
        tiyy_err = np.append(tiyy_err, cdffilePi[i][mms + "_" + "dis" + "_temptensor_err_" + mode][:,1,1])
        tiyz_err = np.append(tiyz_err, cdffilePi[i][mms + "_" + "dis" + "_temptensor_err_" + mode][:,1,2])
        tizx_err = np.append(tizx_err, cdffilePi[i][mms + "_" + "dis" + "_temptensor_err_" + mode][:,2,0])
        tizy_err = np.append(tizy_err, cdffilePi[i][mms + "_" + "dis" + "_temptensor_err_" + mode][:,2,1])
        tizz_err = np.append(tizz_err, cdffilePi[i][mms + "_" + "dis" + "_temptensor_err_" + mode][:,2,2])
        texx_err = np.append(texx_err, cdffilePe[i][mms + "_" + "des" + "_temptensor_err_" + mode][:,0,0])
        texy_err = np.append(texy_err, cdffilePe[i][mms + "_" + "des" + "_temptensor_err_" + mode][:,0,1])
        texz_err = np.append(texz_err, cdffilePe[i][mms + "_" + "des" + "_temptensor_err_" + mode][:,0,2])
        teyx_err = np.append(teyx_err, cdffilePe[i][mms + "_" + "des" + "_temptensor_err_" + mode][:,1,0])
        teyy_err = np.append(teyy_err, cdffilePe[i][mms + "_" + "des" + "_temptensor_err_" + mode][:,1,1])
        teyz_err = np.append(teyz_err, cdffilePe[i][mms + "_" + "des" + "_temptensor_err_" + mode][:,1,2])
        tezx_err = np.append(tezx_err, cdffilePe[i][mms + "_" + "des" + "_temptensor_err_" + mode][:,2,0])
        tezy_err = np.append(tezy_err, cdffilePe[i][mms + "_" + "des" + "_temptensor_err_" + mode][:,2,1])
        tezz_err = np.append(tezz_err, cdffilePe[i][mms + "_" + "des" + "_temptensor_err_" + mode][:,2,2])
        
    for i in range(len(fileE)):
        timeE = np.append(timeE, cdffileE[i][mms + "_edp_epoch_" + mode + "_l2"][...])
        ex = np.append(ex, cdffileE[i][mms + "_edp_dce_gse_" + mode + "_l2"][:,0])
        ey = np.append(ey, cdffileE[i][mms + "_edp_dce_gse_" + mode + "_l2"][:,1])
        ez = np.append(ez, cdffileE[i][mms + "_edp_dce_gse_" + mode + "_l2"][:,2])
        edp_qual = np.append(edp_qual, cdffileE[i][mms + "_edp_quality_" + mode + "_l2"][:])
        
    for i in range(len(fileB)):
        timeb = np.append(timeb, cdffileB[i]["Epoch"][...])
        bx = np.append(bx, cdffileB[i][mms + "_fgm_b_gse_" + modeb + "_l2"][:,0])
        by = np.append(by, cdffileB[i][mms + "_fgm_b_gse_" + modeb + "_l2"][:,1])
        bz = np.append(bz, cdffileB[i][mms + "_fgm_b_gse_" + modeb + "_l2"][:,2])
        timer = np.append(timer, cdffileB[i]["Epoch_state"][...])
        x = np.append(x, cdffileB[i][mms + "_fgm_r_gse_" + modeb + "_l2"][:,0])
        y = np.append(y, cdffileB[i][mms + "_fgm_r_gse_" + modeb + "_l2"][:,1])
        z = np.append(z, cdffileB[i][mms + "_fgm_r_gse_" + modeb + "_l2"][:,2])
        fgm_flag = np.append(fgm_flag, cdffileB[i][mms + "_fgm_flag_" + modeb + "_l2"][:])


    
    dfenergyi = pd.concat([pd.DataFrame(timei), 
                          pd.DataFrame(energy_i01), pd.DataFrame(energy_i02), pd.DataFrame(energy_i03), pd.DataFrame(energy_i04), pd.DataFrame(energy_i05), 
                          pd.DataFrame(energy_i06), pd.DataFrame(energy_i07), pd.DataFrame(energy_i08), pd.DataFrame(energy_i09), pd.DataFrame(energy_i10), 
                          pd.DataFrame(energy_i11), pd.DataFrame(energy_i12), pd.DataFrame(energy_i13), pd.DataFrame(energy_i14), pd.DataFrame(energy_i15), 
                          pd.DataFrame(energy_i16), pd.DataFrame(energy_i17), pd.DataFrame(energy_i18), pd.DataFrame(energy_i19), pd.DataFrame(energy_i20), 
                          pd.DataFrame(energy_i21), pd.DataFrame(energy_i22), pd.DataFrame(energy_i23), pd.DataFrame(energy_i24), pd.DataFrame(energy_i25), 
                          pd.DataFrame(energy_i26), pd.DataFrame(energy_i27), pd.DataFrame(energy_i28), pd.DataFrame(energy_i29), pd.DataFrame(energy_i30), 
                          pd.DataFrame(energy_i31), pd.DataFrame(energy_i32)], axis=1)
    dfenergyi.columns = ['epoch','energy01i'+n, 'energy02i'+n, 'energy03i'+n, 'energy04i'+n, 'energy05i'+n, 
                                 'energy06i'+n, 'energy07i'+n, 'energy08i'+n, 'energy09i'+n, 'energy10i'+n, 
                                 'energy11i'+n, 'energy12i'+n, 'energy13i'+n, 'energy14i'+n, 'energy15i'+n, 
                                 'energy16i'+n, 'energy17i'+n, 'energy18i'+n, 'energy19i'+n, 'energy20i'+n, 
                                 'energy21i'+n, 'energy22i'+n, 'energy23i'+n, 'energy24i'+n, 'energy25i'+n, 
                                 'energy26i'+n, 'energy27i'+n, 'energy28i'+n, 'energy29i'+n, 'energy30i'+n, 
                                 'energy31i'+n, 'energy32i'+n]
    
    dfenergye = pd.concat([pd.DataFrame(timee), 
                          pd.DataFrame(energy_e01), pd.DataFrame(energy_e02), pd.DataFrame(energy_e03), pd.DataFrame(energy_e04), pd.DataFrame(energy_e05), 
                          pd.DataFrame(energy_e06), pd.DataFrame(energy_e07), pd.DataFrame(energy_e08), pd.DataFrame(energy_e09), pd.DataFrame(energy_e10), 
                          pd.DataFrame(energy_e11), pd.DataFrame(energy_e12), pd.DataFrame(energy_e13), pd.DataFrame(energy_e14), pd.DataFrame(energy_e15), 
                          pd.DataFrame(energy_e16), pd.DataFrame(energy_e17), pd.DataFrame(energy_e18), pd.DataFrame(energy_e19), pd.DataFrame(energy_e20), 
                          pd.DataFrame(energy_e21), pd.DataFrame(energy_e22), pd.DataFrame(energy_e23), pd.DataFrame(energy_e24), pd.DataFrame(energy_e25), 
                          pd.DataFrame(energy_e26), pd.DataFrame(energy_e27), pd.DataFrame(energy_e28), pd.DataFrame(energy_e29), pd.DataFrame(energy_e30), 
                          pd.DataFrame(energy_e31), pd.DataFrame(energy_e32)], axis=1)
    dfenergye.columns = ['epoch','energy01e'+n, 'energy02e'+n, 'energy03e'+n, 'energy04e'+n, 'energy05e'+n, 
                                 'energy06e'+n, 'energy07e'+n, 'energy08e'+n, 'energy09e'+n, 'energy10e'+n, 
                                 'energy11e'+n, 'energy12e'+n, 'energy13e'+n, 'energy14e'+n, 'energy15e'+n, 
                                 'energy16e'+n, 'energy17e'+n, 'energy18e'+n, 'energy19e'+n, 'energy20e'+n, 
                                 'energy21e'+n, 'energy22e'+n, 'energy23e'+n, 'energy24e'+n, 'energy25e'+n, 
                                 'energy26e'+n, 'energy27e'+n, 'energy28e'+n, 'energy29e'+n, 'energy30e'+n, 
                                 'energy31e'+n, 'energy32e'+n]
    
    dfe = pd.concat([pd.DataFrame(timeE), pd.DataFrame(ex), pd.DataFrame(ey), pd.DataFrame(ez), pd.DataFrame(edp_qual)], axis=1)
    dfe.columns = ['epoch','ex'+n,'ey'+n,'ez'+n,'qual'+n]
    dfb = pd.concat([pd.DataFrame(timeb), pd.DataFrame(bx), pd.DataFrame(by), pd.DataFrame(bz), pd.DataFrame(fgm_flag)], axis=1)
    dfb.columns = ['epoch','bx'+n,'by'+n,'bz'+n,'flag'+n]
    dfr = pd.concat([pd.DataFrame(timer), pd.DataFrame(x), pd.DataFrame(y), pd.DataFrame(z)], axis=1)
    dfr.columns = ['epoch','x'+n,'y'+n,'z'+n]
    dfui = pd.concat([pd.DataFrame(timei), pd.DataFrame(uix), pd.DataFrame(uiy), pd.DataFrame(uiz), pd.DataFrame(uix_err), pd.DataFrame(uiy_err), pd.DataFrame(uiz_err), pd.DataFrame(uix_spt), pd.DataFrame(uiy_spt), pd.DataFrame(uiz_spt), pd.DataFrame(err_flagi)], axis=1)
    dfui.columns = ['epoch','uix'+n,'uiy'+n,'uiz'+n,'uix_err'+n,'uiy_err'+n,'uiz_err'+n,'uix_spt'+n,'uiy_spt'+n,'uiz_spt'+n,'flagi'+n]
    dfpi = pd.concat([pd.DataFrame(timei), pd.DataFrame(pixx), pd.DataFrame(pixy), pd.DataFrame(pixz), pd.DataFrame(piyx), pd.DataFrame(piyy), pd.DataFrame(piyz), pd.DataFrame(pizx), pd.DataFrame(pizy), pd.DataFrame(pizz),
                      pd.DataFrame(pixx_err), pd.DataFrame(pixy_err), pd.DataFrame(pixz_err), pd.DataFrame(piyx_err), pd.DataFrame(piyy_err), pd.DataFrame(piyz_err), pd.DataFrame(pizx_err), pd.DataFrame(pizy_err), pd.DataFrame(pizz_err)], axis=1)
    dfpi.columns = ['epoch','pixx'+n,'pixy'+n,'pixz'+n,'piyx'+n,'piyy'+n,'piyz'+n,'pizx'+n,'pizy'+n,'pizz'+n,
                    'pixx_err'+n,'pixy_err'+n,'pixz_err'+n,'piyx_err'+n,'piyy_err'+n,'piyz_err'+n,'pizx_err'+n,'pizy_err'+n,'pizz_err'+n]
    dfni = pd.concat([pd.DataFrame(timei), pd.DataFrame(numdeni), pd.DataFrame(numdeni_err)], axis=1)
    dfni.columns = ['epoch','ni'+n,'ni_err'+n]
    dfti = pd.concat([pd.DataFrame(timei), pd.DataFrame(tempi_para), pd.DataFrame(tempi_perp), 
                      pd.DataFrame(tixx),pd.DataFrame(tixy),pd.DataFrame(tixz),pd.DataFrame(tiyx),pd.DataFrame(tiyy),pd.DataFrame(tiyz),pd.DataFrame(tizx),pd.DataFrame(tizy),pd.DataFrame(tizz),
                      pd.DataFrame(tixx_err),pd.DataFrame(tixy_err),pd.DataFrame(tixz_err),pd.DataFrame(tiyx_err),pd.DataFrame(tiyy_err),pd.DataFrame(tiyz_err),pd.DataFrame(tizx_err),pd.DataFrame(tizy_err),pd.DataFrame(tizz_err)], axis=1)
    dfti.columns = ['epoch','Ti_para'+n,'Ti_perp'+n,'Tixx'+n,'Tixy'+n,'Tixz'+n,'Tiyx'+n,'Tiyy'+n,'Tiyz'+n,'Tizx'+n,'Tizy'+n,'Tizz'+n,
                    'Tixx_err'+n,'Tixy_err'+n,'Tixz_err'+n,'Tiyx_err'+n,'Tiyy_err'+n,'Tiyz_err'+n,'Tizx_err'+n,'Tizy_err'+n,'Tizz_err'+n]
    
    dfue = pd.concat([pd.DataFrame(timee), pd.DataFrame(uex), pd.DataFrame(uey), pd.DataFrame(uez), pd.DataFrame(uex_err), pd.DataFrame(uey_err), pd.DataFrame(uez_err), pd.DataFrame(uex_spt), pd.DataFrame(uey_spt), pd.DataFrame(uez_spt), pd.DataFrame(err_flage)], axis=1)
    dfue.columns = ['epoch','uex'+n,'uey'+n,'uez'+n,'uex_err'+n,'uey_err'+n,'uez_err'+n,'uex_spt'+n,'uey_spt'+n,'uez_spt'+n,'flage'+n]
    dfpe = pd.concat([pd.DataFrame(timee), pd.DataFrame(pexx), pd.DataFrame(pexy), pd.DataFrame(pexz), pd.DataFrame(peyx), pd.DataFrame(peyy), pd.DataFrame(peyz), pd.DataFrame(pezx), pd.DataFrame(pezy), pd.DataFrame(pezz),
                      pd.DataFrame(pexx_err), pd.DataFrame(pexy_err), pd.DataFrame(pexz_err), pd.DataFrame(peyx_err), pd.DataFrame(peyy_err), pd.DataFrame(peyz_err), pd.DataFrame(pezx_err), pd.DataFrame(pezy_err), pd.DataFrame(pezz_err)], axis=1)
    dfpe.columns = ['epoch','pexx'+n,'pexy'+n,'pexz'+n,'peyx'+n,'peyy'+n,'peyz'+n,'pezx'+n,'pezy'+n,'pezz'+n,
                    'pexx_err'+n,'pexy_err'+n,'pexz_err'+n,'peyx_err'+n,'peyy_err'+n,'peyz_err'+n,'pezx_err'+n,'pezy_err'+n,'pezz_err'+n]
    dfne = pd.concat([pd.DataFrame(timee), pd.DataFrame(numdene), pd.DataFrame(numdene_err)], axis=1)
    dfne.columns = ['epoch','ne'+n,'ne_err'+n]
    dfte = pd.concat([pd.DataFrame(timee), pd.DataFrame(tempe_para), pd.DataFrame(tempe_perp),
                      pd.DataFrame(texx),pd.DataFrame(texy),pd.DataFrame(texz),pd.DataFrame(teyx),pd.DataFrame(teyy),pd.DataFrame(teyz),pd.DataFrame(tezx),pd.DataFrame(tezy),pd.DataFrame(tezz),
                      pd.DataFrame(texx_err),pd.DataFrame(texy_err),pd.DataFrame(texz_err),pd.DataFrame(teyx_err),pd.DataFrame(teyy_err),pd.DataFrame(teyz_err),pd.DataFrame(tezx_err),pd.DataFrame(tezy_err),pd.DataFrame(tezz_err)], axis=1)
    dfte.columns = ['epoch','Te_para'+n,'Te_perp'+n,'Texx'+n,'Texy'+n,'Texz'+n,'Teyx'+n,'Teyy'+n,'Teyz'+n,'Tezx'+n,'Tezy'+n,'Tezz'+n,
                    'Texx_err'+n,'Texy_err'+n,'Texz_err'+n,'Teyx_err'+n,'Teyy_err'+n,'Teyz_err'+n,'Tezx_err'+n,'Tezy_err'+n,'Tezz_err'+n]
    
    dfe_notnan = dfe
    dfe.loc[dfe['qual'+n] == 0, ['ex'+n, 'ey'+n, 'ez'+n]] = np.nan
    dfe.loc[dfe['qual'+n] == 1, ['ex'+n, 'ey'+n, 'ez'+n]] = np.nan
    E_nan_path = '_withnan'
    
    month_str = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    if mode in ['brst']:
        mode_dir = 'brst/'
    else:
        mode_dir = 'srvy_fast/'
    date_dir = start_D + month_str[int(start_M)-1] + start_Y + '/'
    path_dir = "/Users/thanaponaiamsai/Desktop/Plasma/"+'DataFrame/' + mode_dir + date_dir


    dfe_notnan.to_csv(path_dir + 'e' + n + '.csv', index=False)
    dfe.to_csv(path_dir + 'e' + n + E_nan_path + '.csv', index=False)
    dfb.to_csv(path_dir + 'b' + n + '.csv', index=False)
    dfr.to_csv(path_dir + 'r' + n + '.csv', index=False)
    dfui.to_csv(path_dir + 'ui' + n + '.csv', index=False)
    dfpi.to_csv(path_dir + 'pi' + n + '.csv', index=False)
    dfni.to_csv(path_dir + 'ni' + n + '.csv', index=False)
    dfti.to_csv(path_dir + 'ti' + n + '.csv', index=False)
    dfue.to_csv(path_dir + 'ue' + n + '.csv', index=False)
    dfpe.to_csv(path_dir + 'pe' + n + '.csv', index=False)
    dfne.to_csv(path_dir + 'ne' + n + '.csv', index=False)
    dfte.to_csv(path_dir + 'te' + n + '.csv', index=False)
    dfenergyi.to_csv(path_dir + 'energyi' + n + '.csv', index=False)
    dfenergye.to_csv(path_dir + 'energye' + n + '.csv', index=False)
    print('Success for MMS'+n)
    return 'Success!'



### Cumulative

def acc(q, average=False):
    cumu = np.zeros(len(q)+1)
    latest = 0
    for i in range(1,len(q)+1):
        if np.isnan(q[i-1]):
            cumu[i] = np.nan
        else:
            if np.isnan(cumu[i-1]):
                cumu[i] = latest + q[i-1]
                latest = cumu[i]
            else:
                cumu[i] = cumu[i-1] + q[i-1]
                latest = cumu[i]
    return cumu

def acc_time(time):
    dt = time[2] - time[1]
    t0 = time[0] - dt
    t0 = np.array([t0])
    time_new = np.append(t0, time)
    return time_new


### Convert to Datetime

def readtime1(Ctime):
    return [datetime.strptime(x, '%Y-%m-%d %H:%M:%S.%f') for x in Ctime]
def readtime2(Ctime):
    return [datetime.strptime(x, '%Y-%m-%d %H:%M:%S') for x in Ctime]



### Time Window

def window(df, t1, t2, epoch='epoch', reset_index=True, dropna=False):
    epoch = epoch
    if type(t1)!=datetime:
        t1 = pd.to_datetime(t1, format='%Y-%m-%d %H:%M:%S.%f')
        t2 = pd.to_datetime(t2, format='%Y-%m-%d %H:%M:%S.%f')
    window = (df[epoch] > t1) & (df[epoch] <= t2)
    if reset_index:
        if dropna:
            df1 = df.loc[window].dropna().reset_index().drop('index', axis=1)
        elif not dropna:
            df1 = df.loc[window].reset_index().drop('index', axis=1)
    elif not reset_index:
        if dropna:
            df1 = df.loc[window].dropna()
        elif not dropna:
            df1 = df.loc[window]
    return df1



### Average Data

def average_data(date, n, mode='brst', Enan=True):
    if mode == 'srvy':
        mode = 'srvy_fast'
    path = 'DataFrame/'+mode+'/'+str(date)+'/'
    n = str(n)
    dfui = pd.read_csv(str(path) + 'ui'+n+'.csv')
    dfni = pd.read_csv(str(path) + 'ni'+n+'.csv')
    dfpi = pd.read_csv(str(path) + 'pi'+n+'.csv')
    dfti = pd.read_csv(str(path) + 'ti'+n+'.csv')
    dfenergyi = pd.read_csv(str(path) + 'energyi'+n+'.csv')
    dfenergye = pd.read_csv(str(path) + 'energye'+n+'.csv')
    dfue = pd.read_csv(str(path) + 'ue'+n+'.csv')
    dfne = pd.read_csv(str(path) + 'ne'+n+'.csv')
    dfpe = pd.read_csv(str(path) + 'pe'+n+'.csv')
    dfte = pd.read_csv(str(path) + 'te'+n+'.csv')
    if Enan:
        dfe = pd.read_csv(str(path) + 'e'+n+'_withnan.csv')
    else:
        dfe = pd.read_csv(str(path) + 'e'+n+'.csv')
    dfb = pd.read_csv(str(path) + 'b'+n+'.csv')
    dfr = pd.read_csv(str(path) + 'r'+n+'.csv')
    
    dfui['epoch'] = readtime1(dfui['epoch'])
    dfni['epoch'] = readtime1(dfni['epoch'])
    dfpi['epoch'] = readtime1(dfpi['epoch'])
    dfti['epoch'] = readtime1(dfti['epoch'])
    dfue['epoch'] = readtime1(dfue['epoch'])
    dfne['epoch'] = readtime1(dfne['epoch'])
    dfpe['epoch'] = readtime1(dfpe['epoch'])
    dfte['epoch'] = readtime1(dfte['epoch'])
    dfe['epoch'] = readtime1(dfe['epoch'])
    dfb['epoch'] = readtime1(dfb['epoch'])
    dfr['epoch'] = readtime1(dfr['epoch'])
    dfenergyi['epoch'] = readtime1(dfenergyi['epoch'])
    dfenergye['epoch'] = readtime1(dfenergye['epoch'])
    dfr.set_index('epoch', inplace=True)
    df1r = dfr.resample('150ms').mean().interpolate(method='linear')
    dfr = df1r.reset_index()
    time = dfui['epoch'].values
    dt = (time[1]-time[0])/2
    pos = 0
    ta = time[pos]-dt
    tb = time[pos]+dt
    dfue_window = window(dfue,ta,tb)
    dfe_window = window(dfe,ta,tb)
    dfb_window = window(dfb,ta,tb)
    dfr_window = window(dfr,ta,tb)
    while (len(dfb_window) == 0) or (len(dfe_window) == 0) or (len(dfue_window) == 0) or (len(dfr_window) == 0) :
        ta = time[pos]-dt
        tb = time[pos]+dt
        dfue_window = window(dfue,ta,tb)
        dfe_window = window(dfe,ta,tb)
        dfb_window = window(dfb,ta,tb)
        dfr_window = window(dfr,ta,tb)
        pos +=1
    ta = time[pos]-dt
    tb = time[pos]+dt
    dfue_window = window(dfue,ta,tb)
    dfe_window = window(dfe,ta,tb)
    dfb_window = window(dfb,ta,tb)
    dfr_window = window(dfr,ta,tb)
    print(pos)
    l1 = pos
    l2 = len(dfui)-pos-1
    print('l1:\t',l1)
    print('l2:\t',l2)
    
    uex = []
    uey = []
    uez = []
    x  = []
    y  = []
    z  = []
    ex = []
    ey = []
    ez = []
    bx = []
    by = []
    bz = []
    ne = []
    tepara = []
    teperp = []
    texx = []
    texy = []
    texz = []
    teyx = []
    teyy = []
    teyz = []
    tezx = []
    tezy = []
    tezz = []

    pexx = []
    pexy = []
    pexz = []
    peyx = []
    peyy = []
    peyz = []
    pezx = []
    pezy = []
    pezz = []
    
    energy_e01 = []
    energy_e02 = []
    energy_e03 = []
    energy_e04 = []
    energy_e05 = []
    energy_e06 = []
    energy_e07 = []
    energy_e08 = []
    energy_e09 = []
    energy_e10 = []
    energy_e11 = []
    energy_e12 = []
    energy_e13 = []
    energy_e14 = []
    energy_e15 = []
    energy_e16 = []
    energy_e17 = []
    energy_e18 = []
    energy_e19 = []
    energy_e20 = []
    energy_e21 = []
    energy_e22 = []
    energy_e23 = []
    energy_e24 = []
    energy_e25 = []
    energy_e26 = []
    energy_e27 = []
    energy_e28 = []
    energy_e29 = []
    energy_e30 = []
    energy_e31 = []
    energy_e32 = []
    
    for i in range(l1, l2):
        ta = time[i]-dt
        tb = time[i]+dt

        dfue_window = window(dfue,ta,tb)
        dfne_window = window(dfne,ta,tb)
        dfte_window = window(dfte,ta,tb)
        dfpe_window = window(dfpe,ta,tb)
        dfenergye_window = window(dfenergye,ta,tb)
        dfe_window = window(dfe,ta,tb)
        dfb_window = window(dfb,ta,tb)
        dfr_window = window(dfr,ta,tb)
        
        uex.append(np.mean(dfue_window['uex'+n].values))
        uey.append(np.mean(dfue_window['uey'+n].values))
        uez.append(np.mean(dfue_window['uez'+n].values))

        ne.append(np.mean(dfne_window['ne'+n].values))

        tepara.append(np.mean(dfte_window['Te_para'+n].values))
        teperp.append(np.mean(dfte_window['Te_perp'+n].values))
        texx.append(np.mean(dfte_window['Texx'+n].values))
        texy.append(np.mean(dfte_window['Texy'+n].values))
        texz.append(np.mean(dfte_window['Texz'+n].values))
        teyx.append(np.mean(dfte_window['Teyx'+n].values))
        teyy.append(np.mean(dfte_window['Teyy'+n].values))
        teyz.append(np.mean(dfte_window['Teyz'+n].values))
        tezx.append(np.mean(dfte_window['Tezx'+n].values))
        tezy.append(np.mean(dfte_window['Tezy'+n].values))
        tezz.append(np.mean(dfte_window['Tezz'+n].values))

        pexx.append(np.mean(dfpe_window['pexx'+n].values))
        pexy.append(np.mean(dfpe_window['pexy'+n].values))
        pexz.append(np.mean(dfpe_window['pexz'+n].values))
        peyx.append(np.mean(dfpe_window['peyx'+n].values))
        peyy.append(np.mean(dfpe_window['peyy'+n].values))
        peyz.append(np.mean(dfpe_window['peyz'+n].values))
        pezx.append(np.mean(dfpe_window['pezx'+n].values))
        pezy.append(np.mean(dfpe_window['pezy'+n].values))
        pezz.append(np.mean(dfpe_window['pezz'+n].values))
        
        energy_e01.append(np.mean(dfenergye_window['energy01e'+n].values))
        energy_e02.append(np.mean(dfenergye_window['energy02e'+n].values))
        energy_e03.append(np.mean(dfenergye_window['energy03e'+n].values))
        energy_e04.append(np.mean(dfenergye_window['energy04e'+n].values))
        energy_e05.append(np.mean(dfenergye_window['energy05e'+n].values))
        energy_e06.append(np.mean(dfenergye_window['energy06e'+n].values))
        energy_e07.append(np.mean(dfenergye_window['energy07e'+n].values))
        energy_e08.append(np.mean(dfenergye_window['energy08e'+n].values))
        energy_e09.append(np.mean(dfenergye_window['energy09e'+n].values))
        energy_e10.append(np.mean(dfenergye_window['energy10e'+n].values))
        energy_e11.append(np.mean(dfenergye_window['energy11e'+n].values))
        energy_e12.append(np.mean(dfenergye_window['energy12e'+n].values))
        energy_e13.append(np.mean(dfenergye_window['energy13e'+n].values))
        energy_e14.append(np.mean(dfenergye_window['energy14e'+n].values))
        energy_e15.append(np.mean(dfenergye_window['energy15e'+n].values))
        energy_e16.append(np.mean(dfenergye_window['energy16e'+n].values))
        energy_e17.append(np.mean(dfenergye_window['energy17e'+n].values))
        energy_e18.append(np.mean(dfenergye_window['energy18e'+n].values))
        energy_e19.append(np.mean(dfenergye_window['energy19e'+n].values))
        energy_e20.append(np.mean(dfenergye_window['energy20e'+n].values))
        energy_e21.append(np.mean(dfenergye_window['energy21e'+n].values))
        energy_e22.append(np.mean(dfenergye_window['energy22e'+n].values))
        energy_e23.append(np.mean(dfenergye_window['energy23e'+n].values))
        energy_e24.append(np.mean(dfenergye_window['energy24e'+n].values))
        energy_e25.append(np.mean(dfenergye_window['energy25e'+n].values))
        energy_e26.append(np.mean(dfenergye_window['energy26e'+n].values))
        energy_e27.append(np.mean(dfenergye_window['energy27e'+n].values))
        energy_e28.append(np.mean(dfenergye_window['energy28e'+n].values))
        energy_e29.append(np.mean(dfenergye_window['energy29e'+n].values))
        energy_e30.append(np.mean(dfenergye_window['energy30e'+n].values))
        energy_e31.append(np.mean(dfenergye_window['energy31e'+n].values))
        energy_e32.append(np.mean(dfenergye_window['energy32e'+n].values))
        
        ex.append(np.mean(dfe_window['ex'+n].values))
        ey.append(np.mean(dfe_window['ey'+n].values))
        ez.append(np.mean(dfe_window['ez'+n].values))

        bx.append(np.mean(dfb_window['bx'+n].values))
        by.append(np.mean(dfb_window['by'+n].values))
        bz.append(np.mean(dfb_window['bz'+n].values))

        x.append(np.mean(dfr_window['x'+n].values))
        y.append(np.mean(dfr_window['y'+n].values))
        z.append(np.mean(dfr_window['z'+n].values))
        
    x, y, z = np.array(x), np.array(y), np.array(z)
    bx, by, bz = np.array(bx), np.array(by), np.array(bz)
    ex, ey, ez = np.array(ex), np.array(ey), np.array(ez)
    pexx, pexy, pexz = np.array(pexx), np.array(pexy), np.array(pexz)
    peyx, peyy, peyz = np.array(peyx), np.array(peyy), np.array(peyz)
    pezx, pezy, pezz = np.array(pezx), np.array(pezy), np.array(pezz)
    uex, uey, uez = np.array(uex), np.array(uey), np.array(uez)
    tepara, teperp = np.array(tepara), np.array(teperp)
    texx, texy, texz = np.array(texx), np.array(texy), np.array(texz)
    teyx, teyy, teyz = np.array(teyx), np.array(teyy), np.array(teyz)
    tezx, tezy, tezz = np.array(tezx), np.array(tezy), np.array(tezz)
    ne = np.array(ne)
    energy_e01, energy_e02, energy_e03, energy_e04, energy_e05 = np.array(energy_e01), np.array(energy_e02), np.array(energy_e03), np.array(energy_e04), np.array(energy_e05)
    energy_e06, energy_e07, energy_e08, energy_e09, energy_e10 = np.array(energy_e06), np.array(energy_e07), np.array(energy_e08), np.array(energy_e09), np.array(energy_e10)
    energy_e11, energy_e12, energy_e13, energy_e14, energy_e15 = np.array(energy_e11), np.array(energy_e12), np.array(energy_e13), np.array(energy_e14), np.array(energy_e15)
    energy_e16, energy_e17, energy_e18, energy_e19, energy_e20 = np.array(energy_e16), np.array(energy_e17), np.array(energy_e18), np.array(energy_e19), np.array(energy_e20)
    energy_e21, energy_e22, energy_e23, energy_e24, energy_e25 = np.array(energy_e21), np.array(energy_e22), np.array(energy_e23), np.array(energy_e24), np.array(energy_e25)
    energy_e26, energy_e27, energy_e28, energy_e29, energy_e30 = np.array(energy_e26), np.array(energy_e27), np.array(energy_e28), np.array(energy_e29), np.array(energy_e30)
    energy_e31, energy_e32 = np.array(energy_e31), np.array(energy_e32)
    
    print('Len B : ',len(bx), len(by), len(bz))
    print('Len E : ',len(ex), len(ey), len(ez))
    print('Len R : ',len(x), len(y), len(z))
    print('Len U : ',len(uex), len(uey), len(uez))
    print('Len P : ',len(pexx), len(pexy), len(pexz), len(peyx), len(peyy), len(peyz), len(pezx), len(pezy), len(pezz))
    print('Len N : ',len(ne))
    print('Len T : ',len(tepara), len(teperp))
    print('Len T : ',len(texx), len(texy), len(texz), len(teyx), len(teyy), len(teyz), len(tezx), len(tezy), len(tezz))
    print('Len En: ',len(energy_e04), len(energy_e08), len(energy_e12), len(energy_e16), len(energy_e20), len(energy_e24), len(energy_e28), len(energy_e32))
    print()
    
    timestamp = dfui['epoch'].values[l1:l2]
    uix = dfui['uix'+n].values[l1:l2]
    uiy = dfui['uiy'+n].values[l1:l2]
    uiz = dfui['uiz'+n].values[l1:l2]
    
    pixx = dfpi['pixx'+n].values[l1:l2]
    pixy = dfpi['pixy'+n].values[l1:l2]
    pixz = dfpi['pixz'+n].values[l1:l2]
    piyx = dfpi['piyx'+n].values[l1:l2]
    piyy = dfpi['piyy'+n].values[l1:l2]
    piyz = dfpi['piyz'+n].values[l1:l2]
    pizx = dfpi['pizx'+n].values[l1:l2]
    pizy = dfpi['pizy'+n].values[l1:l2]
    pizz = dfpi['pizz'+n].values[l1:l2]
    
    tipara = dfti['Ti_para'+n].values[l1:l2]
    tiperp = dfti['Ti_perp'+n].values[l1:l2]
    tixx = dfti['Tixx'+n].values[l1:l2]
    tixy = dfti['Tixy'+n].values[l1:l2]
    tixz = dfti['Tixz'+n].values[l1:l2]
    tiyx = dfti['Tiyx'+n].values[l1:l2]
    tiyy = dfti['Tiyy'+n].values[l1:l2]
    tiyz = dfti['Tiyz'+n].values[l1:l2]
    tizx = dfti['Tizx'+n].values[l1:l2]
    tizy = dfti['Tizy'+n].values[l1:l2]
    tizz = dfti['Tizz'+n].values[l1:l2]
    
    energy_i01 = dfenergyi['energy01i'+n].values[l1:l2]
    energy_i02 = dfenergyi['energy02i'+n].values[l1:l2]
    energy_i03 = dfenergyi['energy03i'+n].values[l1:l2]
    energy_i04 = dfenergyi['energy04i'+n].values[l1:l2]
    energy_i05 = dfenergyi['energy05i'+n].values[l1:l2]
    energy_i06 = dfenergyi['energy06i'+n].values[l1:l2]
    energy_i07 = dfenergyi['energy07i'+n].values[l1:l2]
    energy_i08 = dfenergyi['energy08i'+n].values[l1:l2]
    energy_i09 = dfenergyi['energy09i'+n].values[l1:l2]
    energy_i10 = dfenergyi['energy10i'+n].values[l1:l2]
    energy_i11 = dfenergyi['energy11i'+n].values[l1:l2]
    energy_i12 = dfenergyi['energy12i'+n].values[l1:l2]
    energy_i13 = dfenergyi['energy13i'+n].values[l1:l2]
    energy_i14 = dfenergyi['energy14i'+n].values[l1:l2]
    energy_i15 = dfenergyi['energy15i'+n].values[l1:l2]
    energy_i16 = dfenergyi['energy16i'+n].values[l1:l2]
    energy_i17 = dfenergyi['energy17i'+n].values[l1:l2]
    energy_i18 = dfenergyi['energy18i'+n].values[l1:l2]
    energy_i19 = dfenergyi['energy19i'+n].values[l1:l2]
    energy_i20 = dfenergyi['energy20i'+n].values[l1:l2]
    energy_i21 = dfenergyi['energy21i'+n].values[l1:l2]
    energy_i22 = dfenergyi['energy22i'+n].values[l1:l2]
    energy_i23 = dfenergyi['energy23i'+n].values[l1:l2]
    energy_i24 = dfenergyi['energy24i'+n].values[l1:l2]
    energy_i25 = dfenergyi['energy25i'+n].values[l1:l2]
    energy_i26 = dfenergyi['energy26i'+n].values[l1:l2]
    energy_i27 = dfenergyi['energy27i'+n].values[l1:l2]
    energy_i28 = dfenergyi['energy28i'+n].values[l1:l2]
    energy_i29 = dfenergyi['energy29i'+n].values[l1:l2]
    energy_i30 = dfenergyi['energy30i'+n].values[l1:l2]
    energy_i31 = dfenergyi['energy31i'+n].values[l1:l2]
    energy_i32 = dfenergyi['energy32i'+n].values[l1:l2]
    
    df = pd.DataFrame()
    df['epoch'] = timestamp
    df['uix'+n] = uix
    df['uiy'+n] = uiy
    df['uiz'+n] = uiz
    df['uex'+n] = uex
    df['uey'+n] = uey
    df['uez'+n] = uez
    
    df['pixx'+n] = pixx
    df['pixy'+n] = pixy
    df['pixz'+n] = pixz
    df['piyx'+n] = piyx
    df['piyy'+n] = piyy
    df['piyz'+n] = piyz
    df['pizx'+n] = pizx
    df['pizy'+n] = pizy
    df['pizz'+n] = pizz
    df['pexx'+n] = pexx
    df['pexy'+n] = pexy
    df['pexz'+n] = pexz
    df['peyx'+n] = peyx
    df['peyy'+n] = peyy
    df['peyz'+n] = peyz
    df['pezx'+n] = pezx
    df['pezy'+n] = pezy
    df['pezz'+n] = pezz
    
    
    df['ex'+n] = ex
    df['ey'+n] = ey
    df['ez'+n] = ez
    df['bx'+n] = bx
    df['by'+n] = by
    df['bz'+n] = bz
    df['x'+n] = x
    df['y'+n] = y
    df['z'+n] = z
    
    df['ne'+n] = ne
        
    df['Te_para'+n] = tepara
    df['Te_perp'+n] = teperp
    df['Ti_para'+n] = tipara
    df['Ti_perp'+n] = tiperp
    
    df['Tixx'+n] = tixx
    df['Tixy'+n] = tixy
    df['Tixz'+n] = tixz
    df['Tiyx'+n] = tiyx
    df['Tiyy'+n] = tiyy
    df['Tiyz'+n] = tiyz
    df['Tizx'+n] = tizx
    df['Tizy'+n] = tizy
    df['Tizz'+n] = tizz
    df['Texx'+n] = texx
    df['Texy'+n] = texy
    df['Texz'+n] = texz
    df['Teyx'+n] = teyx
    df['Teyy'+n] = teyy
    df['Teyz'+n] = teyz
    df['Tezx'+n] = tezx
    df['Tezy'+n] = tezy
    df['Tezz'+n] = tezz
    
    df['energy01i'+n] = energy_i01
    df['energy02i'+n] = energy_i02
    df['energy03i'+n] = energy_i03
    df['energy04i'+n] = energy_i04
    df['energy05i'+n] = energy_i05
    df['energy06i'+n] = energy_i06
    df['energy07i'+n] = energy_i07
    df['energy08i'+n] = energy_i08
    df['energy09i'+n] = energy_i09
    df['energy10i'+n] = energy_i10
    df['energy11i'+n] = energy_i11
    df['energy12i'+n] = energy_i12
    df['energy13i'+n] = energy_i13
    df['energy14i'+n] = energy_i14
    df['energy15i'+n] = energy_i15
    df['energy16i'+n] = energy_i16
    df['energy17i'+n] = energy_i17
    df['energy18i'+n] = energy_i18
    df['energy19i'+n] = energy_i19
    df['energy20i'+n] = energy_i20
    df['energy21i'+n] = energy_i21
    df['energy22i'+n] = energy_i22
    df['energy23i'+n] = energy_i23
    df['energy24i'+n] = energy_i24
    df['energy25i'+n] = energy_i25
    df['energy26i'+n] = energy_i26
    df['energy27i'+n] = energy_i27
    df['energy28i'+n] = energy_i28
    df['energy29i'+n] = energy_i29
    df['energy30i'+n] = energy_i30
    df['energy31i'+n] = energy_i31
    df['energy32i'+n] = energy_i32
    
    df['energy01e'+n] = energy_e01
    df['energy02e'+n] = energy_e02
    df['energy03e'+n] = energy_e03
    df['energy04e'+n] = energy_e04
    df['energy05e'+n] = energy_e05
    df['energy06e'+n] = energy_e06
    df['energy07e'+n] = energy_e07
    df['energy08e'+n] = energy_e08
    df['energy09e'+n] = energy_e09
    df['energy10e'+n] = energy_e10
    df['energy11e'+n] = energy_e11
    df['energy12e'+n] = energy_e12
    df['energy13e'+n] = energy_e13
    df['energy14e'+n] = energy_e14
    df['energy15e'+n] = energy_e15
    df['energy16e'+n] = energy_e16
    df['energy17e'+n] = energy_e17
    df['energy18e'+n] = energy_e18
    df['energy19e'+n] = energy_e19
    df['energy20e'+n] = energy_e20
    df['energy21e'+n] = energy_e21
    df['energy22e'+n] = energy_e22
    df['energy23e'+n] = energy_e23
    df['energy24e'+n] = energy_e24
    df['energy25e'+n] = energy_e25
    df['energy26e'+n] = energy_e26
    df['energy27e'+n] = energy_e27
    df['energy28e'+n] = energy_e28
    df['energy29e'+n] = energy_e29
    df['energy30e'+n] = energy_e30
    df['energy31e'+n] = energy_e31
    df['energy32e'+n] = energy_e32
    df.to_csv(str(path)+'df'+n+'_raw.csv', index=False)
    
    print('success MMS'+n)



### Import DataFrame

def import_df0(path = 'DataFrame/brst/08Sep2015/'):
    files = glob.glob(path+"df*.csv")
    file = []
    mms = []
    for i in range(len(files)):
        for j in range(1,5):
            if files[i] == str(path)+"df"+str(j)+".csv":
                mms.append(j)
                file.append(files[i])
    df0 = []
    for i in range(len(file)):
        df0.append(pd.read_csv(file[i]))
    for i in range(len(df0)):
        df0[i]['epoch'+str(mms[i])] = readtime1(df0[i]['epoch'+str(mms[i])])
    return df0, mms
def data(df0, mms, settime=False, time1=None, time2=None):
    t = []
    b = []
    e = []
    ui = []
    ue = []
    r = []
    ti = []
    te = []
    n = []
    energy_i = []
    energy_e = []
    ji = []
    je = []
    ji_para = []
    je_para = []
    e_para = []
    ji_perp = []
    je_perp = []
    e_perp = []
    jiE_para = []
    jeE_para = []
    jiE_perp = []
    jeE_perp = []
    jE_para = []
    jE_perp = []
    jE_tot = []
    file = []
    
    sorted_mms = sorted(mms)
    index = []
    for i in sorted_mms:
        index.append(mms.index(i))
    dfs = []
    for i in index:
        dfs.append(df0[i])
    df0 = dfs
    mms = sorted_mms
    
    if settime:
        time1 = time1
        time2 = time2
        for i in range(len(df0)):
            df0[i] = window(df0[i], time1[i], time2[i], epoch='epoch'+str(mms[i]))
    
    for i in range(len(df0)):
        t.append(df0[i]['epoch'+str(mms[i])].values)
        b.append(df0[i][['bx'+str(mms[i]),'by'+str(mms[i]),'bz'+str(mms[i])]].values)
        e.append(df0[i][['ex'+str(mms[i]),'ey'+str(mms[i]),'ez'+str(mms[i])]].values)
        ui.append(df0[i][['uix'+str(mms[i]),'uiy'+str(mms[i]),'uiz'+str(mms[i])]].values)
        ue.append(df0[i][['uex'+str(mms[i]),'uey'+str(mms[i]),'uez'+str(mms[i])]].values)
        r.append(df0[i][['x'+str(mms[i]),'y'+str(mms[i]),'z'+str(mms[i])]].values)
        ti.append(df0[i][['Ti_para'+str(mms[i]),'Ti_perp'+str(mms[i])]].values)
        te.append(df0[i][['Te_para'+str(mms[i]),'Te_perp'+str(mms[i])]].values)
        n.append(df0[i][['ne'+str(mms[i])]].values)
        energy_i.append(df0[i][['energy01i'+str(mms[i]), 'energy02i'+str(mms[i]), 'energy03i'+str(mms[i]), 'energy04i'+str(mms[i]), 'energy05i'+str(mms[i]),
                                'energy06i'+str(mms[i]), 'energy07i'+str(mms[i]), 'energy08i'+str(mms[i]), 'energy09i'+str(mms[i]), 'energy10i'+str(mms[i]),
                                'energy11i'+str(mms[i]), 'energy12i'+str(mms[i]), 'energy13i'+str(mms[i]), 'energy14i'+str(mms[i]), 'energy15i'+str(mms[i]),
                                'energy16i'+str(mms[i]), 'energy17i'+str(mms[i]), 'energy18i'+str(mms[i]), 'energy19i'+str(mms[i]), 'energy20i'+str(mms[i]),
                                'energy21i'+str(mms[i]), 'energy22i'+str(mms[i]), 'energy23i'+str(mms[i]), 'energy24i'+str(mms[i]), 'energy25i'+str(mms[i]),
                                'energy26i'+str(mms[i]), 'energy27i'+str(mms[i]), 'energy28i'+str(mms[i]), 'energy29i'+str(mms[i]), 'energy30i'+str(mms[i]),
                                'energy31i'+str(mms[i]), 'energy32i'+str(mms[i])]].values)
        energy_e.append(df0[i][['energy01e'+str(mms[i]), 'energy02e'+str(mms[i]), 'energy03e'+str(mms[i]), 'energy04e'+str(mms[i]), 'energy05e'+str(mms[i]),
                                'energy06e'+str(mms[i]), 'energy07e'+str(mms[i]), 'energy08e'+str(mms[i]), 'energy09e'+str(mms[i]), 'energy10e'+str(mms[i]),
                                'energy11e'+str(mms[i]), 'energy12e'+str(mms[i]), 'energy13e'+str(mms[i]), 'energy14e'+str(mms[i]), 'energy15e'+str(mms[i]),
                                'energy16e'+str(mms[i]), 'energy17e'+str(mms[i]), 'energy18e'+str(mms[i]), 'energy19e'+str(mms[i]), 'energy20e'+str(mms[i]),
                                'energy21e'+str(mms[i]), 'energy22e'+str(mms[i]), 'energy23e'+str(mms[i]), 'energy24e'+str(mms[i]), 'energy25e'+str(mms[i]),
                                'energy26e'+str(mms[i]), 'energy27e'+str(mms[i]), 'energy28e'+str(mms[i]), 'energy29e'+str(mms[i]), 'energy30e'+str(mms[i]),
                                'energy31e'+str(mms[i]), 'energy32e'+str(mms[i])]].values)
        ji.append(df0[i][['jix'+str(mms[i]),'jiy'+str(mms[i]),'jiz'+str(mms[i])]].values)
        je.append(df0[i][['jex'+str(mms[i]),'jey'+str(mms[i]),'jez'+str(mms[i])]].values)
        ji_para.append(df0[i][['ji_parax'+str(mms[i]),'ji_paray'+str(mms[i]),'ji_paraz'+str(mms[i])]].values)
        je_para.append(df0[i][['je_parax'+str(mms[i]),'je_paray'+str(mms[i]),'je_paraz'+str(mms[i])]].values)
        ji_perp.append(df0[i][['ji_perpx'+str(mms[i]),'ji_perpy'+str(mms[i]),'ji_perpz'+str(mms[i])]].values)
        je_perp.append(df0[i][['je_perpx'+str(mms[i]),'je_perpy'+str(mms[i]),'je_perpz'+str(mms[i])]].values)
        e_para.append(df0[i][['e_parax'+str(mms[i]),'e_paray'+str(mms[i]),'e_paraz'+str(mms[i])]].values)
        e_perp.append(df0[i][['e_perpx'+str(mms[i]),'e_perpy'+str(mms[i]),'e_perpz'+str(mms[i])]].values)

        jiE_para.append(df0[i]['jiE_para'+str(mms[i])].values)
        jeE_para.append(df0[i]['jeE_para'+str(mms[i])].values)
        jiE_perp.append(df0[i]['jiE_perp'+str(mms[i])].values)
        jeE_perp.append(df0[i]['jeE_perp'+str(mms[i])].values)
        jE_para.append(df0[i]['jE_para'+str(mms[i])].values)
        jE_perp.append(df0[i]['jE_perp'+str(mms[i])].values)
        jE_tot.append(df0[i]['jE_tot'+str(mms[i])].values)
    
    jE_net = []
    for i in range(len(df0)):
        jE_net.append(np.sum(jE_tot[i][~np.isnan(jE_tot[i])]))
    
    return df0, t, b, e, r, energy_i, energy_e, ui, ue, ti, te, n, ji, je, ji_para, je_para, ji_perp, je_perp, e_para, e_perp, jiE_para, jiE_perp, jeE_para, jeE_perp, jE_para, jE_perp, jE_tot, jE_net

def import_df(file, path = 'DataFrame/brst/08Sep2015/', ConvertTime=False):
    df = pd.read_csv(str(path) + str(file) + '.csv')
    if ConvertTime:
        df['epoch'] = readtime1(df['epoch'])
    return df

