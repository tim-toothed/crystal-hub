import numpy as np

#######READ FROM FILE########

def cfp_from_file(conf):

    prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

    if conf[0]=='d':
        file = open('tables/cfp_d_conf.txt', 'r').readlines()
    elif conf[0]=='f':
        file = open('tables/cfp_f_conf.txt', 'r').readlines()[1:] #skip the first line
    else:
        raise ValueError('conf must be dn or fn')
    
    check = False
    cfp = []
    for i,line in enumerate(file):
        if conf in line:
            check = True
        if check == True:
            if r'#' not in line:
                cfp.append(line)
            elif r'#' in line:
                break

    cfp_dic = {}
    for i,line in enumerate(cfp):
        if i>0:
            splitline = line.split('\t')
            factor = int(splitline[2])
            number = 1
            if len(splitline)>3:
                splitline[-1] = splitline[-1].split()
                for j in range(len(splitline[-1])):
                    number *= prime[j]**int(splitline[-1][j])
            number = np.sqrt(number)
            number *= factor
            try:
                cfp_dic[splitline[0]].update({splitline[1]:number})
                #cfp_dic[splitline[0]] += [splitline[1],number]
            except:
                cfp_dic[splitline[0]] = {splitline[1]: number}

    return cfp_dic

def read_matrix_from_file(conf_print, closed_shell=False):

    file = open('tables/tables_'+conf_print[0]+'conf.txt').readlines()
    
    dizionario = {}
    conf = None
    for i,line in enumerate(file):
        if 'CONFIGURATION' in line:
            line = line.strip()
            conf = line.split()[-1]
            dizionario[conf] = {}
        elif 'MATRIX' in line:
            line = line.strip()
            key1 = line.split()[-1]
            dizionario[conf][key1] = {}
        elif '#' in line:
            pass
        else:
            splitline = line.split()
            if splitline[0] not in dizionario[conf][key1].keys():
                dizionario[conf][key1][splitline[0]] = {}
            dizionario[conf][key1][splitline[0]][splitline[1]] = float(splitline[-1])

            L = state_legend(splitline[0][1])
            S = int(splitline[0][0])
            L1 = state_legend(splitline[1][1])
            S1 = int(splitline[1][0])

            if key1=='V11':
                if splitline[1] not in dizionario[conf][key1].keys():
                    dizionario[conf][key1][splitline[1]] = {}
                dizionario[conf][key1][splitline[1]][splitline[0]] = float(splitline[-1])*(-1)**(L-L1-S/2+S1/2)
            else:
                if splitline[1] not in dizionario[conf][key1].keys():
                    dizionario[conf][key1][splitline[1]] = {}
                dizionario[conf][key1][splitline[1]][splitline[0]] = float(splitline[-1])*(-1)**(L-L1)

    if closed_shell==False:
        conf_str = conf_print
    else:
        conf_n = almost_closed_shells(conf_print)
        conf_str = conf[0]+str(conf_n)

    return dizionario[conf_str]

def read_ee_int(conf, closed_shell):

    if closed_shell==False:
        conf_str = conf
    else:
        conf_n = almost_closed_shells(conf)
        conf_str = conf[0]+str(conf_n)

    file = open('tables/dic_ee_values.txt', 'r').readlines()
    dic_ee_loaded = {}
    conf = None
    for line in file:
        if line=='#':
            conf = None
        else:
            splitline = line.split()
            if len(splitline)==1:
                conf = splitline[0].strip()
                dic_ee_loaded[conf] = {}
            else:
                if splitline[0] not in dic_ee_loaded[conf].keys():
                    dic_ee_loaded[conf][splitline[0]] = {}
               
                if splitline[0]!=splitline[1]:
                    if splitline[1] not in dic_ee_loaded[conf].keys():
                        dic_ee_loaded[conf][splitline[1]] = {}
                    dic_ee_loaded[conf][splitline[1]].update({splitline[0]:np.array(splitline[2:], dtype=float)})
                dic_ee_loaded[conf][splitline[0]].update({splitline[1]:np.array(splitline[2:], dtype=float)})

   
    return dic_ee_loaded[conf_str]

def calc_ee_int(conf, closed_shell=False):

    if closed_shell==False:
        conf_str = conf
    else:
        conf_n = almost_closed_shells(conf)
        conf_str = conf[0]+str(conf_n)

    import sympy

    def calc_dic_ek(conf, dic_ek_out):

        def omegaUL(U,L):
            gU = {99:0,
                10:6,
                11:12,
                20:14,
                21:21,
                30:24,
                22:30,
                31:32,
                40:36}
            return 1/2*state_legend(L)*(state_legend(L)+1) - gU[U]

        S_list = [7/2, 3, 5/2, 5/2, 2, 2, 3/2, 3/2, 1, 1/2]
        v_list = [7, 6, 5, 7, 4, 6, 5, 7, 6, 7]
        Sv_list = [(S_list[i],v_list[i]) for i in range(len(S_list))]

        AB_label = {'f5': {'2F':[6,7], '2H':[6,7], '2I':[4,5], '2K':[4,5]},
                    'f6': {'3F':[8,9], '3H':[8,9], '3I':[5,6], '3K':[5,6], '1G':[7,8], '1I':[6,7], '1L':[3,4]},
                    'f7': {'2F':[6,7], '2G':[9,0], '2H':[6,7], '2I':[4,5,8,9], '2K':[4,5], '2L':[4,5]}}

        #W, U, U1
        x_g = {200:
                {20:{20:2}},  #questo è back calcolato da 1D1:1D1 sulla base delle tabelle di ryley
               210:
                {11:
                {21:12*np.sqrt(455)},
                20:
                {20:-6/7, 21:6*np.sqrt(66)/7},
                21:
                {11:12*np.sqrt(455), 20:6*np.sqrt(66)/7, 21:[3/7, 0]}},
               211:
                {10:
                {30:-20*np.sqrt(143)},
                11:
                {21:10*np.sqrt(182), 30:10},
                20:
                {20:-8/7, 21:4*np.sqrt(33)/7, 30:4*np.sqrt(3)},
                21:
                {11:10*np.sqrt(182), 20:4*np.sqrt(33)/7, 21:[4/7, 3], 30:2},
                30:
                {10:-20*np.sqrt(143), 11:10, 20:4*np.sqrt(3), 21:2, 30:2}},
               220:
                {20:
                {20:3/14, 21:3*np.sqrt(55)/7, 22:-3*np.sqrt(5/28)},
                21:
                {20:3*np.sqrt(55)/7, 21:[-6/7,-3], 22:3/np.sqrt(7)},
                22:
                {20:-3*np.sqrt(5/28), 21:3/np.sqrt(7), 22:3/2}},
               221:
                {10:
                {30:5*np.sqrt(143), 31:-15*np.sqrt(429)},
                11:
                {21:14*np.sqrt(910/11), 30:2*np.sqrt(10), 31:2*np.sqrt(39)/11},
                20:
                {20:2/7, 21:-10*np.sqrt(6)/7, 30:np.sqrt(3), 31:9*np.sqrt(3/7)},
                21:
                {11:14*np.sqrt(910/11), 20:-10*np.sqrt(6)/7, 21:[-1/7,12/11], 30: 5*np.sqrt(2/11), 31:3*np.sqrt(2)/11},
                30:
                {10: 5*np.sqrt(143), 11:2*np.sqrt(10), 20: np.sqrt(3), 21:5*np.sqrt(2/11), 30:-1/2, 31:3/(2*np.sqrt(11))},
                31:
                {10:-15*np.sqrt(429), 11:2*np.sqrt(39)/11, 20:9*np.sqrt(3/7), 21:3*np.sqrt(2)/11, 30:3/(2*np.sqrt(11)), 31:1/22}},
               222:
                {99:
                {40:-30*np.sqrt(143)},
                10:
                {30:-3*np.sqrt(1430), 40:9*np.sqrt(1430)},
                20:
                {20:6/11, 30:-3*np.sqrt(42/11), 40:9*np.sqrt(2)/11},
                30:
                {10:-3*np.sqrt(1430), 20:-3*np.sqrt(42/11), 30:-3, 40:1/np.sqrt(11)},
                40:
                {99:-30*np.sqrt(143), 10:9*np.sqrt(1430), 20:9*np.sqrt(2)/11, 30:1/np.sqrt(11), 40:3/11}}}
        #U1 L U
        chi_L = {20:
                {"D":{20:143},
                "G":{20:-130},
                "I":{20:35}},
                21:
                {"H":{11:1, 20:0, 21:[49,-75]},
                "D":{20:-39*np.sqrt(2), 21:[377, 13]},
                "F":{21:[455, -65]},
                "G":{20:4*np.sqrt(65), 21:[-561, 55]},
                "K":{21:[-315, 133]},
                "L":{21:[245, -75]}},
                30:
                {"P":{11:-13*np.sqrt(11), 30:-52},
                "F":{10:1, 21:12*np.sqrt(195), 30:38},
                "G":{20:-13*np.sqrt(5), 21:8*np.sqrt(143), 30:-52},
                "H":{11:np.sqrt(39), 21:11*np.sqrt(42), 30:88},
                "I":{20:30, 30:25},
                "K":{21:-4*np.sqrt(17), 30:-94},
                "M":{30:25}},
                22:
                {"S":{22:260},
                "D":{20:3*np.sqrt(429), 21:45*np.sqrt(78), 22:-25},
                "G":{20:-38*np.sqrt(65), 21:12*np.sqrt(11), 22:94},
                "H":{21:-12*np.sqrt(546), 22:104},
                "I":{20:21*np.sqrt(85), 22:-181},
                "L":{21:-8*np.sqrt(665), 22:-36},
                "N":{22:40}},
                31:
                {"P":{11:11*np.sqrt(330), 30:76*np.sqrt(143), 31:-6644},
                "D":{20:-8*np.sqrt(78), 21:-60*np.sqrt(39/7), 31:4792},
                "F":{10:[0,1], 21:[-312*np.sqrt(5), 12*np.sqrt(715)], 30:[-48*np.sqrt(39), -98*np.sqrt(33)], 31:[4420, -902, 336*np.sqrt(143)]},
                "G":{20:5*np.sqrt(65), 21:2024/np.sqrt(7), 30:20*np.sqrt(1001), 31:-2684},
                "H":{11:[11*np.sqrt(85), -25*np.sqrt(77)], 21:[31*np.sqrt(1309/3), 103*np.sqrt(5/3)], 30:[-20*np.sqrt(374), -44*np.sqrt(70)], 31:[-2024,2680,-48*np.sqrt(6545)]},
                "I":{20:[10*np.sqrt(21),0], 30:[-57*np.sqrt(33), 18*np.sqrt(1122)], 31:[-12661/5,17336/5,-3366*np.sqrt(34)/5]},
                "K":{21:[-52*np.sqrt(323/23), -336*np.sqrt(66/23)], 30:[-494*np.sqrt(19/23), 73*np.sqrt(1122/23)], 31:[123506/23, -85096/23, 144*np.sqrt(21318)/23]},
                "L":{21:-24*np.sqrt(190), 31:-4712},
                "M":{30:-21*np.sqrt(385), 31:-473},
                "N":{31:1672},
                "O":{31:220}},
                40:
                {"S":{99:1, 40:-1408},
                "D":{20:-88*np.sqrt(13), 40:-44},
                "F":{10:1, 30:90*np.sqrt(11), 40:1078},
                "G":{20:[53*np.sqrt(715/27), 7*np.sqrt(15470/27)], 30:[-16*np.sqrt(1001), 64*np.sqrt(442)], 40:[-16720/9, 10942/9, -34*np.sqrt(2618)/9]},
                "H":{30:-72*np.sqrt(462), 40:-704},
                "I":{20:[34*np.sqrt(1045/31), -12*np.sqrt(1785/31)], 30:[-9*np.sqrt(21945/31), 756*np.sqrt(85/31)], 40:[-2453/31, 36088/31, 60*np.sqrt(74613)/31]},
                "K":{30:-84*np.sqrt(33), 40:-132},
                "L":{40:[-4268/31, 11770/31, 924*np.sqrt(1995)/31]},
                "M":{30:-99*np.sqrt(15), 40:-1067},
                "N":{40:528},
                "Q":{40:22}}}
        #U1 L U
        phi_L = {11:
                {"P":{11:-11},
                "H":{11:3}},
                20:
                {"D":{20:-11},
                "G":{20:-4},
                "I":{20:7}},
                21:
                {"D":{20:6*np.sqrt(2), 21:-57},
                "F":{10:1, 21:63},
                "G":{20:np.sqrt(65), 21:55},
                "H":{21:-105},
                "K":{21:-14},
                "L":{21:42}},
                30:
                {"P":{11:np.sqrt(11), 30:83},
                "F":{21:np.sqrt(195), 30:-72},
                "G":{20:2*np.sqrt(5), 21:-np.sqrt(143), 30:20},
                "H":{11:np.sqrt(39), 21:-2*np.sqrt(42), 30:-15},
                "I":{20:3, 30:42},
                "K":{21:-4*np.sqrt(17), 30:-28},
                "M":{30:6}},
                22:
                {"S":{99:1, 22:144},
                "D":{20:3*np.sqrt(429), 22:69},
                "G":{20:4*np.sqrt(65), 22:-148},
                "H":{22:72},
                "I":{20:3*np.sqrt(85), 22:39},
                "L":{22:-96},
                "N":{22:56}},
                31:
                {"P":{11:np.sqrt(330), 30:17*np.sqrt(143), 31:209},
                "D":{21:12*np.sqrt(273), 31:-200},
                "F":{10:[1,0], 21:[-36*np.sqrt(5), -3*np.sqrt(715)], 30:[-16*np.sqrt(39), 24*np.sqrt(33)], 31:[624, -616, -80*np.sqrt(143)]},
                "G":{21:11*np.sqrt(7), 30:4*np.sqrt(1001), 31:836},
                "H":{11:[np.sqrt(85), np.sqrt(77)], 21:[-2*np.sqrt(1309/3), -74*np.sqrt(5/3)], 30:[np.sqrt(187/2), 31*np.sqrt(35/2)], 31:[-1353/2, 703/2, -5*np.sqrt(6545)/2]},
                "I":{30:[30*np.sqrt(33), 0], 31:[-2662/5, -88/5, 528*np.sqrt(34)/5]},
                "K":{21:[-28*np.sqrt(323/23), 42*np.sqrt(66/23)], 30:[4*np.sqrt(437),0], 31:[6652/23, -5456/23, 96*np.sqrt(21318)/23]},
                "L":{21:-6*np.sqrt(190), 31:-464},
                "M":{30:-6*np.sqrt(385), 31:814},
                "N":{31:-616},
                "O":{31:352}},
                40:
                {"S":{22:2*np.sqrt(2145)},
                "D":{20:11*np.sqrt(13), 21:-6*np.sqrt(26), 22:9*np.sqrt(33)},
                "F":{21:3*np.sqrt(455)},
                "G":{20:[-4*np.sqrt(715/27),np.sqrt(15470/27)], 21:[-131*np.sqrt(11/27), 17*np.sqrt(238/27)], 22:[-4*np.sqrt(11/27), -17*np.sqrt(238/27)]},
                "H":{21:-12*np.sqrt(21), 22:3*np.sqrt(286)},
                "I":{20:[7*np.sqrt(1045/31),3*np.sqrt(1785/31)], 22:[3*np.sqrt(3553/31),75*np.sqrt(21/31)]},
                "K":{21:-2*np.sqrt(119)},
                "L":{21:[22*np.sqrt(105/31), -84*np.sqrt(19/31)], 22:[4*np.sqrt(627/31), 12*np.sqrt(385/31)]},
                "N":{22:-np.sqrt(2530)}}}
        #conf v:(2*S+1):U v1:(2*S1+1):U1
        y_g = {"f3":
                {"1:2:10":{"3:2:21":-6*np.sqrt(22)},
                "3:2:11":{"3:2:11":2},
                "3:2:20":{"3:2:20":10/7, "3:2:21":2*np.sqrt(66)/7},
                "3:2:21":{"3:2:20":2*np.sqrt(66)/7, "3:2:21":2/7}},
                "f4":
                {"2:3:10":{"4:3:21":-12*np.sqrt(33/5)},
                "2:3:11":{"4:3:11":6/5, "4:3:30":6},
                "4:3:10":{"4:3:21":8*np.sqrt(11/15)},
                "4:3:11":{"4:3:11":29/15,"4:3:30":-1/3},
                "4:3:20":{"4:3:20":6/7, "4:3:21":-8*np.sqrt(11/147), "4:3:30":4/np.sqrt(3)},
                "4:3:21":{"4:3:10":8*np.sqrt(11/15), "4:3:20":-8*np.sqrt(11/147), "4:3:21":-2/21, "4:3:30":-4/3},
                "4:3:30":{"4:3:11":-1/3, "4:3:20":4/np.sqrt(3), "4:3:21":-4/3, "4:3:30":1/3},
                "0:1:99":{"4:1:22":-12*np.sqrt(22)},
                "2:1:20":{"4:1:20":3*np.sqrt(3/175), "4:1:21":-4*np.sqrt(33/35), "4:1:22":-np.sqrt(3/5)},
                "4:1:20":{"4:1:20":221/140, "4:1:21":8*np.sqrt(11/245), "4:1:22":-np.sqrt(7/80)},
                "4:1:21":{"4:1:20":8*np.sqrt(11/245), "4:1:21":2/7},
                "4:1:22":{"4:1:20":-np.sqrt(7/80), "4:1:22":1/4}},
                "f5":
                {"3:4:10":{"5:4:21":9*np.sqrt(11)},
                "3:4:20":{"5:4:20":3/np.sqrt(7), "5:4:21":np.sqrt(33/7), "5:4:30":-2*np.sqrt(21)},
                "5:4:10":{"5:4:21":-np.sqrt(55/3)},
                "5:4:11":{"5:4:11":-1/3, "5:4:30":-5/3},
                "5:4:20":{"5:4:20":5/7, "5:4:21":5*np.sqrt(11/147), "5:4:30":2/np.sqrt(3)},
                "5:4:21":{"5:4:10":-np.sqrt(55/3), "5:4:20":5*np.sqrt(11/147), "5:4:21":-4/21, "5:4:30":-2/3},
                "5:4:30":{"5:4:11":-5/3, "5:4:20":2/np.sqrt(3), "5:4:21":-2/3, "5:4:30":-1/3},
                "1:2:10":{"5:2:21":36/np.sqrt(5), "5:2:31":-36*np.sqrt(2)},
                "3:2:11":{"5:2:11":3/np.sqrt(2), "5:2:30":3*np.sqrt(5)/2, "5:2:31":-np.sqrt(39/8)},
                "3:2:20":{"5:2:20":3/7, "5:2:21":-11*np.sqrt(6)/7, "5:2:30":-4*np.sqrt(3)},
                "3:2:21":{"5:2:10":3*np.sqrt(33/10), "5:2:20":-3*np.sqrt(33/98), "5:2:21":3/(7*np.sqrt(11)), "5:2:30":-3/(2*np.sqrt(2)), "5:2:31":3/(2*np.sqrt(22))},
                "5:2:10":{"5:2:21":43/np.sqrt(30), "5:2:31":4*np.sqrt(3)},
                "5:2:11":{"5:2:11":-5/6, "5:2:30":-5*np.sqrt(5/72), "5:2:31":-np.sqrt(13/48)},
                "5:2:20":{"5:2:20":11/7, "5:2:21":-11/(7*np.sqrt(6)), "5:2:30":4/np.sqrt(3)},
                "5:2:21":{"5:2:10":43/np.sqrt(30), "5:2:20":-11/(7*np.sqrt(6)), "5:2:21":25/231, "5:2:30":29/(6*np.sqrt(22)), "5:2:31":1/(22*np.sqrt(2))},
                "5:2:30":{"5:2:11":-5*np.sqrt(5/72), "5:2:20":4/np.sqrt(3), "5:2:21":29/(6*np.sqrt(22)), "5:2:30":-1/12, "5:2:31":1/(4*np.sqrt(11))},
                "5:2:31":{"5:2:10":4*np.sqrt(3), "5:2:11":-np.sqrt(13/48), "5:2:21":1/(22*np.sqrt(2)), "5:2:30":1/(4*np.sqrt(11)), "5:2:31":1/44}},
                "f6":
                {"4:5:10":{"6:5:21":-6*np.sqrt(11)},
                "4:5:20":{"6:5:20":-2*np.sqrt(2/7), "6:5:21":2*np.sqrt(33/7)},
                "2:3:10":{"6:3:21":-48*np.sqrt(2/5), "6:3:31":-36},
                "2:3:11":{"6:3:11":np.sqrt(6/5), "6:3:30":np.sqrt(3), "6:3:31":3*np.sqrt(13/10)},
                "4:3:10":{"6:3:21":46/np.sqrt(15), "6:3:31":-8*np.sqrt(6)},
                "4:3:11":{"6:3:11":11/(3*np.sqrt(5)), "6:3:30":-19/(3*np.sqrt(2)), "6:3:31":np.sqrt(13/60)},
                "4:3:20":{"6:3:20":-6*np.sqrt(2)/7, "6:3:21":-22/(7*np.sqrt(3)), "6:3:30":8*np.sqrt(2/3)},
                "4:3:21":{"6:3:10":-np.sqrt(110/3), "6:3:20":np.sqrt(22/147), "6:3:21":-16/(21*np.sqrt(11)), "6:3:30":5/(3*np.sqrt(2)), "6:3:31":1/np.sqrt(22)},
                "4:3:30":{"6:3:11":-np.sqrt(5)/3, "6:3:20":4*np.sqrt(2/3), "6:3:21":4/(3*np.sqrt(11)), "6:3:30":1/(3*np.sqrt(2)), "6:3:31":-1/np.sqrt(22)},
                "2:1:20":{"6:1:20":6/np.sqrt(55), "6:1:30":2*np.sqrt(42/5), "6:1:40":6*np.sqrt(2/55)},
                "4:1:20":{"6:1:20":-61/np.sqrt(770), "6:1:30":8*np.sqrt(3/5), "6:1:40":-6/np.sqrt(385)},
                "4:1:21":{"6:1:10":3*np.sqrt(22), "6:1:20":np.sqrt(2/7), "6:1:30":-np.sqrt(3), "6:1:40":1/np.sqrt(7)},
                "4:1:22":{"6:1:99":-4*np.sqrt(33/5), "6:1:20":-1/np.sqrt(22), "6:1:40":2/np.sqrt(11)}},
                "f7":
                {"3:4:99":{"7:4:22":-12*np.sqrt(11)},
                "3:4:10":{"7:4:21":6*np.sqrt(33)},
                "3:4:20":{"7:4:20":-np.sqrt(5/7), "7:4:21":2*np.sqrt(11/7), "7:4:22":-1},
                "3:2:11":{"7:2:30":2*np.sqrt(10)},
                "3:2:20":{"7:2:20":-16/np.sqrt(77), "7:2:30":-2*np.sqrt(6), "7:2:40":6*np.sqrt(2/77)},
                "3:2:21":{"7:2:10":-np.sqrt(66), "7:2:20":np.sqrt(6/7), "7:2:30":1, "7:2:40":np.sqrt(3/7)}}}

        def calc_ek(conf, label1, label2, S, L, dic_ek): 

            #dic_ek [label1:n:v:U:2S+1:L][label2:n:v1:U1:2S+1:L] 
            v,W,U = terms_labels_symm(conf)[label1]
            v1,W1,U1 = terms_labels_symm(conf)[label2]

            ek_coeff = np.zeros(4)
            n = int(conf[1:]) 
            if label1==label2:
                ek_coeff[0] = n*(n-1)/2 
                ek_coeff[1] = 9*(n - v)/2 + 1/4*v*(v+2) - S*(S+1)

            if v==v1:

                if v!=2*S and W==W1 and int(str(W)[0])==2:
                    factor1 = x_g.get(W, {}).get(U, {}).get(U1)
                    if factor1 is None: 
                        factor1 = x_g.get(W, {}).get(U1, {}).get(U)
                    factor2 = chi_L.get(U1, {}).get(L, {}).get(U)
                    if factor2 is None:  
                        factor2 = chi_L.get(U, {}).get(L, {}).get(U1)
                    if factor1 is not None and factor2 is not None:
                        if (isinstance(factor1, float) or isinstance(factor1, int)) and (isinstance(factor2, float) or isinstance(factor2, int)):
                            ek_coeff[2] = factor1*factor2
                        elif (isinstance(factor1, float) or isinstance(factor1, int)) and isinstance(factor2, list):
                            idx = 1
                            try:       
                                if AB_label[conf][label1[:2]].index(int(label1[-1]))%2==0:
                                    idx = 0
                            except KeyError:
                                if AB_label[conf][label2[:2]].index(int(label2[-1]))%2==0:
                                    idx = 0
                            if label1!=label2 and v==v1 and U==U1 and W==W1:
                                idx = -1
                            ek_coeff[2] = factor2[idx]*factor1
                        else:
                            ek_coeff[2] = np.sum(np.array(factor1)*np.array(factor2))
                        
                    if (S,v) in Sv_list:
                        ek_coeff[2] *= -1
                    
                if n==2*S and U1==U:
                    ek_coeff[3] = -3*omegaUL(U,L)+omegaUL(U,L)  #I'll remove it later
                if v==n and (v==6 or v==7):
                    ek_coeff[3] = 0
                key1 = ':'.join([label1,str(v),str(v),str(U),str(int(2*S+1)),L])
                key2 = ':'.join([label2,str(v),str(v),str(U1),str(int(2*S+1)),L])
                key1_cut = ':'.join([label1[:2],str(v),str(v),str(U),str(int(2*S+1)),L])
                key2_cut = ':'.join([label2[:2],str(v),str(v),str(U1),str(int(2*S+1)),L])

                if dic_ek['f'+str(v)].get(key1, {}).get(key2) is not None:
                    prev_int = dic_ek['f'+str(v)][key1][key2][3].copy()
                    if U1==U and label1==label2:
                        prev_int += omegaUL(U,L)
                    if n==v+2:
                        ek_coeff[3] = prev_int*(1-v)/(7-v)
                    elif n==v+4:
                        ek_coeff[3] = prev_int*(-4)/(7-v)
                elif dic_ek['f'+str(v)].get(key1_cut, {}).get(key2_cut) is not None:
                    prev_int = dic_ek['f'+str(v)][key1_cut][key2_cut][3].copy()
                    if U1==U and label1==label2:
                        prev_int += omegaUL(U,L)
                    if n==v+2:
                        ek_coeff[3] = prev_int*(1-v)/(7-v)
                    elif n==v+4:
                        ek_coeff[3] = prev_int*(-4)/(7-v)
                        
            else:

                if n==5 and ((v,v1)==(1,3) or (v1,v)==(1,3)) and (2*S+1)==2:
                    key1 = ':'.join([label1,'3',str(v),str(U),str(int(2*S+1)),L])
                    key2 = ':'.join([label2,'3',str(v1),str(U1),str(int(2*S+1)),L])
                    if dic_ek['f3'].get(key1, {}).get(key2) is not None:
                        ek_coeff[3] = dic_ek['f3'][key1][key2][3]*np.sqrt(2/5)
                    elif dic_ek['f3'].get(key2, {}).get(key1) is not None:
                        ek_coeff[3] = dic_ek['f3'][key2][key1][3]*np.sqrt(2/5)
                if n==6 and ((v,v1)==(0,4) or (v,v1)==(4,0)) and (2*S+1)==1:
                    key1 = ':'.join([label1,'4',str(v),str(U),str(int(2*S+1)),L])
                    key2 = ':'.join([label2,'4',str(v1),str(U1),str(int(2*S+1)),L])
                    if dic_ek['f4'].get(key1, {}).get(key2) is not None:
                        ek_coeff[3] = dic_ek['f4'][key1][key2][3]*np.sqrt(9/5)
                    elif dic_ek['f4'].get(key2, {}).get(key1) is not None:
                        ek_coeff[3] = dic_ek['f4'][key2][key1][3]*np.sqrt(9/5)
                if n==6 and ((v,v1)==(2,4) or (v,v1)==(4,2)):
                    key1 = ':'.join([label1,'4',str(v),str(U),str(int(2*S+1)),L])
                    key2 = ':'.join([label2,'4',str(v1),str(U1),str(int(2*S+1)),L])
                    if dic_ek['f4'].get(key1, {}).get(key2) is not None:
                        ek_coeff[3] = dic_ek['f4'][key1][key2][3]*np.sqrt(1/6)
                    elif dic_ek['f4'].get(key2, {}).get(key1) is not None:
                        ek_coeff[3] = dic_ek['f4'][key2][key1][3]*np.sqrt(1/6)
                if n==7 and ((v,v1)==(1,5) or (v,v1)==(5,1)) and (2*S+1)==2:
                    key1 = ':'.join([label1,'5',str(v),str(U),str(int(2*S+1)),L])
                    key2 = ':'.join([label2,'5',str(v1),str(U1),str(int(2*S+1)),L])
                    if dic_ek['f5'].get(key1, {}).get(key2) is not None:
                        ek_coeff[3] = dic_ek['f5'][key1][key2][3]*np.sqrt(3/2)
                    elif dic_ek['f5'].get(key2, {}).get(key1) is not None:
                        ek_coeff[3] = dic_ek['f5'][key2][key1][3]*np.sqrt(3/2)

            key1 = ':'.join([str(v),str(int(2*S+1)),str(U)])
            key2 = ':'.join([str(v1),str(int(2*S+1)),str(U1)])
            factor1 = y_g.get(conf, {}).get(key1, {}).get(key2)
            if factor1 is None: 
                factor1 = y_g.get(conf, {}).get(key2, {}).get(key1)
            factor2 = phi_L.get(U1, {}).get(L, {}).get(U)
            if factor2 is None:  
                factor2 = phi_L.get(U, {}).get(L, {}).get(U1)
            if factor1 is not None and factor2 is not None:
                if (isinstance(factor1, float) or isinstance(factor1, int)) and (isinstance(factor2, float) or isinstance(factor2, int)):
                    ek_coeff[3] = factor1*factor2
                else:
                    idx = 1
                    try:       
                        if AB_label[conf][label1[:2]].index(int(label1[-1]))%2==0:
                            idx = 0
                    except KeyError:
                        if AB_label[conf][label2[:2]].index(int(label2[-1]))%2==0:
                            idx = 0
                    if label1!=label2 and v==v1 and U==U1 and W==W1:
                        idx = -1

                    ek_coeff[3] = factor2[idx]*factor1

            if U1==U and v1==v and label1==label2:
                ek_coeff[3] -= omegaUL(U,L)

            key1 = ':'.join([label1,str(n),str(v),str(U),str(int(2*S+1)),L])
            key2 = ':'.join([label2,str(n),str(v1),str(U1),str(int(2*S+1)),L])
            if key1 not in dic_ek[conf].keys():
                dic_ek[conf][key1] = {}
            if key2 not in dic_ek[conf].keys():
                dic_ek[conf][key2] = {}
            else:
                if key2 not in dic_ek[conf][key1].keys():
                    dic_ek[conf][key1][key2] = ek_coeff
                    dic_ek[conf][key2][key1] = ek_coeff

            return dic_ek 

        calc = calculation(conf, TAB=True, wordy=False)
        basis = calc.basis
        dic_ek_out[conf] = {}
        labels_list = []
        for i in range(basis.shape[0]):
            statei = basis[i]
            Si = statei[0]/2.
            Li = statei[1]
            Ji = statei[2]/2.
            MJi = statei[3]/2.
            labeli = calc.dic_LS[':'.join([f'{qq}' for qq in statei])]

            for j in range(0,i+1):
                statej = basis[j]
                Sj = statej[0]/2.
                Lj = statej[1]
                Jj = statej[2]/2.
                MJj = statej[3]/2.
                labelj = calc.dic_LS[':'.join([f'{qq}' for qq in statej])]

                if Ji==Jj and MJi==MJj and Li == Lj and Si == Sj:
                    if labeli+':'+labelj not in labels_list and labelj+':'+labeli not in labels_list:
                        labels_list.append(labeli+':'+labelj)
                        dic_ek_out = calc_ek(conf, labeli, labelj, Si, state_legend(str(Li), inv=True), dic_ek_out)
                

        return dic_ek_out

    def calc_ee_int(conf, label, label1, S, L, dic_ek):

        n = int(conf[1:])
        v,W,U = terms_labels_symm(conf)[label]
        v1,W1,U1 = terms_labels_symm(conf)[label1]
        key1 = ':'.join([label,str(n),str(v),str(U),str(int(2*S+1)),L])
        key2 = ':'.join([label1,str(n),str(v1),str(U1),str(int(2*S+1)),L])
        ek_coeff = dic_ek[key1][key2]

        F0, F2, F4, F6 = sympy.Symbol("F0, F2, F4, F6")
        coeff = [F0, F2, F4, F6]  

        Ek_coeff = []
        Ek_coeff.append(coeff[0] - 10/225*coeff[1] - 33/1089*coeff[2] - 286*25/184041*coeff[3])
        Ek_coeff.append(1/9*(70/225*coeff[1] + 231/1089*coeff[2] + 2002*25/184041*coeff[3]))
        Ek_coeff.append(1/9*(1/225*coeff[1] - 3/1089*coeff[2] + 7*25/184041*coeff[3]))
        Ek_coeff.append(1/3*(5/225*coeff[1] + 6/1089*coeff[2] - 91*25/184041*coeff[3]))

        ee_int = sympy.simplify(Ek_coeff[0]*ek_coeff[0] + Ek_coeff[1]*ek_coeff[1] + Ek_coeff[2]*ek_coeff[2] + Ek_coeff[3]*ek_coeff[3])

        return ee_int, ek_coeff

    #create the dictionry with e-e interaction
    conf_list = ['f3', 'f4', 'f5', 'f6', 'f7']
    dic_ek_conf = {'f0':{},
                   'f1':{},
                   'f2':
                   {'3P:2:2:11:3:P':{'3P:2:2:11:3:P':np.array([0,0,0,22+11])},
                   '3F:2:2:10:3:F':{'3F:2:2:10:3:F':np.array([0,0,0,0])},
                   '3H:2:2:11:3:H':{'3H:2:2:11:3:H':np.array([0,0,0,-6-3])},
                   '1S:2:0:99:1:S':{'1S:2:0:99:1:S':np.array([0,0,0,0])},
                   '1D:2:2:20:1:D':{'1D:2:2:20:1:D':np.array([0,0,0,-22+11])},
                   '1G:2:2:20:1:G':{'1G:2:2:20:1:G':np.array([0,0,0,-8+4])},
                   '1I:2:2:20:1:I':{'1I:2:2:20:1:I':np.array([0,0,0,14-7])}}}
    
    dic_ee_expr = {}
    dic_ee_values = {}

    for conf in conf_list:

        dic_ek_conf = calc_dic_ek(conf, dic_ek_conf)
        dic_ee_expr[conf] = {}
        dic_ee_values[conf] = {}

        calc = calculation(conf, TAB=True, wordy=False)
        basis = calc.basis
        for i in range(basis.shape[0]):
            statei = basis[i]
            Si = statei[0]/2.
            Li = statei[1]
            Ji = statei[2]/2.
            MJi = statei[3]/2.
            labeli = calc.dic_LS[':'.join([f'{qq}' for qq in statei])]

            for j in range(0,i+1):
                statej = basis[j]
                Sj = statej[0]/2.
                Lj = statej[1]
                Jj = statej[2]/2.
                MJj = statej[3]/2.
                labelj = calc.dic_LS[':'.join([f'{qq}' for qq in statej])]

                if Ji==Jj and MJi==MJj and Li == Lj and Si == Sj:
                    ee_value, ee_c = calc_ee_int(conf, labeli, labelj, Si, state_legend(str(Li), inv=True), dic_ek_conf[conf])
                    dic_ee_expr[conf][labeli+':'+labelj] = ee_value
                    dic_ee_values[conf][labeli+':'+labelj] = ee_c
                    if labeli!=labelj:
                        dic_ee_expr[conf][labelj+':'+labeli] = ee_value
                        dic_ee_values[conf][labelj+':'+labeli] = ee_c

    return dic_ee_expr[conf_str], dic_ee_values[conf_str]
