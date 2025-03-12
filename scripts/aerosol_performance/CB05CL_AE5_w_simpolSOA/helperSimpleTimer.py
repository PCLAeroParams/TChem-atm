import numpy as np
import json
# read a database produced by CSPlib
def readDataBase(file_names, first_name='' ):
    data=[]
    for f in file_names:
        f1 = open(first_name+f,)
        try:
            data += [json.load(f1)]
#             break
        except :
            print (first_name+f)

        f1.close()
    return data
# read a database produced by the kokko profiling tools: simple timer json
def readDataBaseJson(file_names,first_name):
    data=[]
    read_files=[]
    for i, f in enumerate(file_names):
        try:
            f1 = open(first_name+f,)
            data += [json.load(f1)]
            f1.close()
            read_files += [f]
        except :
            print (first_name+f)
    return data, read_files
# get data in array form from data that is produced from readDataBaseJson
def getTimeTask(data, selected_tasks):
    Nfiles = len(data)
    Ntask = len(selected_tasks)
    task_times = np.empty([Nfiles, Ntask])

    # loop over files
    for i in range(Nfiles):
        # get a dic with task index
        # the index of task can change between files,
        #so I cannot assume that this index is equal for all files
        task_index = getDataJson(data[i]['kokkos-kernel-data'])
        for j,itask in enumerate(selected_tasks):
            # find index for selected task
            ## get time of select task
            task_times[i,j] = data[i]['kokkos-kernel-data']["kernel-perf-info"][task_index[itask]]['time-per-call']
    return task_times
# make a list of file name e.g., simple_timer10_csp_times.dat
def makeFileNames(increment,first_name, middle, last_name):
    files= []
    for i in increment:
        files += [first_name +str(i)+middle+last_name]
    return files

# get data in array form from data that is produce from readDataBase: format from CSPlib
def getTimeTask2(data, list_tasks, class_name):
    Nfiles = len(data)
    Ntask =len(list_tasks)
    task_times = np.empty([Nfiles,Ntask])
    for j,itask in enumerate(list_tasks):
        for i in range(Nfiles):
            task_times[i,j] = data[i][class_name][itask]/data[i]['Number of state vectors']
    return task_times
#get a dic with kernel names and index
def getDataJson(data):
    # get data
    # other line
    kernel_info = data["kernel-perf-info"]
    nk=len(kernel_info)
    # loop over kernel name, and save index for CSPlib, TChem, and Tines, do not include kokkos
    task_index = {}
    for i in range(nk):
        kernel_name = data["kernel-perf-info"][i]['kernel-name']
        library=kernel_name.split('::')[0]
        print_it=False
        if library=='CSPlib':
            print_it=True
        elif  library=='TChem':
            print_it=True
        elif library=='Tines':
            print_it=True
        if print_it:
            task_index.update({kernel_name:i})
    return task_index
