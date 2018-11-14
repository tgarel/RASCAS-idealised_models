
# ----------------------------------------------------------------------
# functions to write RASCAS parameter files from dictionaries.
# ----------------------------------------------------------------------

def write_param_section(f,sectionName,paramsDict):
    f.write('[%s] \n'%sectionName)
    for k in paramsDict.keys():
        f.write('  %s = %s \n'%(k,paramsDict[k]))

def write_parameter_file(filename,paramdict):
    f = open(filename,'w')
    for s in paramdict.keys():
        write_param_section(f,s,paramdict[s])
    f.close()

# ----------------------------------------------------------------------

