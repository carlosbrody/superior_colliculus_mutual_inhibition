import numpy as np
import scipy.linalg as s 
from tqdm import tqdm
import h5py
import matplotlib.pyplot as plt

# study = run_study(10000)
# plot_study(study, modes=True)
# plot_study(study, modes=False)

def load_mat(filename):
    f = h5py.File(filename)
    d = dict()
    d['weights'] = np.array(f['paramvals'])   
    d['arg_names']=[]
    ref = f['argnames']
    for i in range(0,26):
        obj = f[ref[i]]
        d['arg_names'].append(''.join(chr(i) for i in obj[:]))
    return d

def get_solution_w(d, dex):
    # Need to make sure I have cwa1/cwa2 mapped right. 
    # should probably make this algorithmic
    # Do I have the mapping correct?
    dwpa1 = d['weights'][9,dex]
    dwpa2 = d['weights'][21,dex]
    dwa1p = d['weights'][11,dex]
    dwa2p = d['weights'][19,dex]
    vwpa1 = d['weights'][6,dex]
    vwpa2 = d['weights'][20,dex]
    vwa1p = d['weights'][1,dex]
    vwa2p = d['weights'][22,dex]
    swp = d['weights'][15,dex]
    swa1 = d['weights'][13,dex]
    swa2 = d['weights'][25,dex]
    hwp = d['weights'][4,dex] 
    hwa1a2 = d['weights'][5,dex]
    hwa2a2 = d['weights'][8,dex]
    hwa1a1 = d['weights'][23,dex]
    hwa2a1 = d['weights'][2,dex]
    cwa1 = d['weights'][17,dex]
    cwa2 = d['weights'][3,dex]

    W = np.array([
        [swp,   vwpa1,  vwpa2,  dwpa2,  dwpa1,  hwp],
        [vwa1p, swa1,   cwa1,   hwa1a2, hwa1a1, dwa1p],
        [vwa2p, cwa2,   swa2,   hwa2a2, hwa2a1, dwa2p],
        [dwa2p, hwa2a1, hwa2a2, swa2,   cwa2,   vwa2p],
        [dwa1p, hwa1a1, hwa1a2, cwa1,   swa1,   vwa1p],
        [hwp,   dwpa1,  dwpa2,  vwpa2,  vwpa1,  swp]
        ])
    return W

def run_study(num2test):
    random = test_random(num2test)
    sym    = test_random(num2test,force_symmetry=True)
    bound  = test_random(num2test,force_symmetry=True, force_bounds=True)
    solutions = test_random(num2test, solutions=True)
    return random, sym, bound, solutions

def plot_study(study, modes=True):
    plt.figure()
    w = 1
    dw = .25
    count = 0
    n=66
    if 'Other' in study[2][0]:
        study[2][0].pop('Other')
    if 'Other' in study[3][0]:
        study[3][0].pop('Other')

    if modes:
        setdex = 0
        ystr = 'Mode Frequency (% of solutions)'
        filename = 'modes'
        chance = 6/8*100
    else:
        setdex = 1
        ystr = 'Mode Stability (% of modes)'
        filename = 'stability'
        chance = 0.5*100

    plt.axhline(chance, color='k',linestyle='--', alpha=.5)
    for key in study[2][setdex].keys():
        if np.mod(count,2)==0:
            plt.axvspan(count*w-dw,(count+1)*w-dw,color='k', alpha=0.1)
        if count ==0:
            plt.plot([count*w, count*w+dw], [study[2][setdex][key]*100, study[2][setdex][key]*100],'k',linewidth=3,label='Random')
            plt.plot([count*w+dw,count*w+dw*2], [study[3][setdex][key]*100, study[3][setdex][key]*100],'r',linewidth=3,label='Solutions')
        else:
            plt.plot([count*w, count*w+dw], [study[2][setdex][key]*100, study[2][setdex][key]*100],'k',linewidth=3)
            plt.plot([count*w+dw,count*w+dw*2], [study[3][setdex][key]*100, study[3][setdex][key]*100],'r',linewidth=3)   
        p = study[3][setdex][key]
        err = 1.98*np.sqrt((p*(1-p))/n)
        plt.plot([count*w+dw*1.5,count*w+dw*1.5], [study[3][setdex][key]*100-err*100, study[3][setdex][key]*100+err*100],'r',linewidth=1)
        count+=1
    plt.ylabel(ystr,fontsize=14)
    plt.xticks(np.array(range(0,len(study[2][setdex].keys())))+dw, study[2][setdex].keys(), rotation=60,fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(0,100)
    plt.xlim(left=-.25)
    if modes:
        plt.xlim(right=7.75)
    else:
        plt.xlim(right=7.75)
    plt.legend(loc='lower left',fontsize=14)
    plt.tight_layout()
    plt.savefig(filename+".svg")
    plt.savefig(filename+".png")

def get_w(symmetry, bounds):
    if bounds: # Rand within bounds
        # I double checked this. 
        # all are -3,3 ,except self weights, which are 0 to 3
        dwpa1 = np.random.rand()*6-3
        dwpa2 = np.random.rand()*6-3
        dwa1p = np.random.rand()*6-3
        dwa2p = np.random.rand()*6-3
        vwpa1 = np.random.rand()*6-3
        vwpa2 = np.random.rand()*6-3
        vwa1p = np.random.rand()*6-3
        vwa2p = np.random.rand()*6-3
        swp = np.random.rand()*3
        swa1 = np.random.rand()*3
        swa2 = np.random.rand()*3
        hwp = np.random.rand()*6-3
        hwa1a2 = np.random.rand()*6-3
        hwa2a2 = np.random.rand()*6-3
        hwa1a1 = np.random.rand()*6-3
        hwa2a1 = np.random.rand()*6-3
        cwa1 = np.random.rand()*6-3
        cwa2 = np.random.rand()*6-3

        W = np.array([
            [swp,   vwpa1,  vwpa2,  dwpa2,  dwpa1,  hwp],
            [vwa1p, swa1,   cwa1,   hwa1a2, hwa1a1, dwa1p],
            [vwa2p, cwa2,   swa2,   hwa2a2, hwa2a1, dwa2p],
            [dwa2p, hwa2a1, hwa2a2, swa2,   cwa2,   vwa2p],
            [dwa1p, hwa1a1, hwa1a2, cwa1,   swa1,   vwa1p],
            [hwp,   dwpa1,  dwpa2,  vwpa2,  vwpa1,  swp]
            ])
    elif symmetry: # Randn with symmetry
        dwpa1 = np.random.randn()
        dwpa2 = np.random.randn()
        dwa1p = np.random.randn()
        dwa2p = np.random.randn()
        vwpa1 = np.random.randn()
        vwpa2 = np.random.randn()
        vwa1p = np.random.randn()
        vwa2p = np.random.randn()
        swp = np.random.randn()
        swa1 = np.random.randn()
        swa2 = np.random.randn()
        hwp = np.random.randn()
        hwa1a2 = np.random.randn()
        hwa2a2 = np.random.randn()
        hwa1a1 = np.random.randn()
        hwa2a1 = np.random.randn()
        cwa1 = np.random.randn()
        cwa2 = np.random.randn()

        W = np.array([
            [swp,   vwpa1,  vwpa2,  dwpa2,  dwpa1,  hwp],
            [vwa1p, swa1,   cwa1,   hwa1a2, hwa1a1, dwa1p],
            [vwa2p, cwa2,   swa2,   hwa2a2, hwa2a1, dwa2p],
            [dwa2p, hwa2a1, hwa2a2, swa2,   cwa2,   vwa2p],
            [dwa1p, hwa1a1, hwa1a2, cwa1,   swa1,   vwa1p],
            [hwp,   dwpa1,  dwpa2,  vwpa2,  vwpa1,  swp]
            ])
    else: # Pure Random
        W = np.random.randn(6,6)
    return W

def test_random(num2test, force_symmetry=False, force_bounds=False, solutions=False):

    # Set up dictionaries of mode and stability counts
    modes = {'All':0, 'Side':0,'Task':0,'Diagonal':0,'A1':0,'A2':0,'Mixed_1':0,'Mixed_2':0,'Other':0}
    stability = {'All':0, 'Side':0,'Task':0,'Diagonal':0,'A1':0,'A2':0,'Mixed_1':0,'Mixed_2':0,'Other':0}

    # Define modes 
    # Units are          [Pro_L  Anti_1L  Anti_2L  Anti_2R  Anti_1R  Pro_R
    side_mask = np.array([True,  True,    True,    False,   False,   False])
    task_mask = np.array([True,  False,   False,   False,   False,   True])
    diag_mask = np.array([True,  False,   False,   True,    True,    False])
    a1_mask =   np.array([True,  False,   True,    True,    False,   True])
    a2_mask =   np.array([True,  True,    False,   False,   True,    True])
    mxa1_mask = np.array([True,  False,   True,    False,   True,    False])
    mxa2_mask = np.array([True,  True,    False,   True,    False,   False])
   
    if solutions:
        d = load_mat('solutions6.mat')
        num2test = np.shape(d['weights'])[1]
 
    # Iterate random samples 
    for count in tqdm(range(0, num2test)):

        # Generate network and do schur decomposition
        if solutions:
            W = get_solution_w(d, count)
        else:
            W = get_w(force_symmetry, force_bounds)
        vals, vecs = s.schur(W)
        
        # Iterate modes
        for dex in range(0,6):
            # Extract relavant mode and eigenvalue
            vec = vecs[:,dex]
            val = vals[dex,dex]
    
            # Classify mode
            if np.all(vec > 0) or np.all(vec < 0):
                mode = 'All'
            elif np.all([x==y for (x,y) in zip(vec < 0, side_mask)]) or np.all([x==y for (x,y) in zip(vec < 0, ~side_mask)]):
                mode = 'Side'
            elif np.all([x==y for (x,y) in zip(vec < 0, task_mask)]) or np.all([x==y for (x,y) in zip(vec < 0, ~task_mask)]):
                mode = 'Task'
            elif np.all([x==y for (x,y) in zip(vec < 0, diag_mask)]) or np.all([x==y for (x,y) in zip(vec < 0, ~diag_mask)]):
                mode = 'Diagonal'
            elif np.all([x==y for (x,y) in zip(vec < 0, a1_mask)]) or np.all([x==y for (x,y) in zip(vec < 0, ~a1_mask)]):
                mode = 'A1' # "Its an all vs a1"
            elif np.all([x==y for (x,y) in zip(vec < 0, a2_mask)]) or np.all([x==y for (x,y) in zip(vec < 0, ~a2_mask)]):
                mode = 'A2' # "Its an all vs a2"
            elif np.all([x==y for (x,y) in zip(vec < 0, mxa1_mask)]) or np.all([x==y for (x,y) in zip(vec < 0, ~mxa1_mask)]):
                mode = 'Mixed_1' # "Its a side mode, but a1 is mixed"
            elif np.all([x==y for (x,y) in zip(vec < 0, mxa2_mask)]) or np.all([x==y for (x,y) in zip(vec < 0, ~mxa2_mask)]):
                mode = 'Mixed_2' # "Its a side mode, but a2 is mixed"
            else:
                mode = 'Other' 
        
            # Tabulate mode and stability
            modes[mode] +=1
            if val > 0:
                stability[mode] +=1

    
    # Finished sampling networks, post-processing
    for key in modes.keys():
        # Convert to fraction stable    
        if modes[key] > 0:
            stability[key] = np.round(stability[key]/modes[key],2)
        else:
            # None of this mode, stability is not defined
            stability.pop(key)  
            
        # Convert to fraction with mode
        if modes[key] > 0:
            modes[key] = np.round(modes[key]/num2test,2)

    return modes, stability           

