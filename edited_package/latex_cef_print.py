import numpy as np

def printLaTexCEFparams(Bs):
    precision = 5
    '''prints CEF parameters in the output that Latex can read'''
    print('\\begin{table}\n\\caption{Fitted vs. Calculated CEF parameters for ?}')
    print('\\begin{ruledtabular}')
    print('\\begin{tabular}{c|'+'c'*len(Bs)+'}')
    # Create header
    print('$B_n^m$ (meV)' +' & Label'*len(Bs)
        +' \\tabularnewline\n \\hline ')
    for i, (n,m) in enumerate([[n,m] for n in range(2,8,2) for m in range(0,n+1, 3)]):
        print('$ B_'+str(n)+'^'+str(m)+'$ &', 
              ' & '.join([str(np.around(bb[i],decimals=precision)) for bb in Bs]),
              '\\tabularnewline')
    print('\\end{tabular}\\end{ruledtabular}')
    print('\\label{flo:CEF_params}\n\\end{table}')
