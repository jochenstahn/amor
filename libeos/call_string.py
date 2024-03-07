# module to assemble an orso 'call' string form command line arguments and some default values.

def call_string():
    base = 'python eos.py '

    inpt = ''
    if clas.year:
        inpt += f' --year {clas.year}'
    else:
        inpt += f' --year {datetime.now().year}'
    if clas.dataPath:
        inpt += f' --dataPath {clas.dataPath}'
    if clas.subtract:
        inpt += f' -subtract {clas.subtract}'
    if clas.normalisationFileIdentifier:
        inpt += f' -r {clas.normalisationFileIdentifier}'
    # get file list somehow
    if ...:
        inpt += f' -n {flie_list}'
    else:
        inpt += f' -n {clas.fileIdentifier}'

    otpt = ''
    if outputFormats != 'Rqz.ort':
        otpt =  f" -of  '{outputFormats}'"
    if clas.outputName:
        otpt += f' -o {clas.outputName}'
    else:
        pass
        # default name
        
    mask = ''    
    if clas.yRange:
        mask += f' -y {clas.yRange}'
    if clas.lambdaRange:
        mask += f' -l {clas.lambdaRange}'
    if clas.thetaRange:
        mask += f' -- thetaRange {clas.thetaRange}'
    elif clas.thetaRangeR:
        mask += f' -t {clas.thetaRangeR}'
    if clas.qzRange:
        mask += f' -q {clas.qzRange}'
    if clas.qResolution:
        mask += f' -a {clas.qResolution}'

    para = ''
    if clas.chopperPhase:
        para += f' --chopperPhase {clas.chopperPhase}'
    if clas.chopperPhaseOffset:
        para += f' --chopperPhaseOffset {clas.chopperPhaseOffset}'
    if clas.mu:
        para += f' --mu {clas.mu}'
    elif clas.muOffset:
        para += f' --muOffset {clas.muOffset}'
    if clas.nu:
        para += f' --nu {clas.nu}'

    if clas.sampleModel:
        modl =  f" --sampleModel '{clas.sampleModel}'"

    acts = ''
    if clas.autoscale:
        acts += f' --autoscale {clas.autoscale}'
    if clas.scale:
        acts += f' --scale {clas.scale}'
    if clas.timeSlize:
        acts += f' --timeSlize {clas.timeSlize}'

    mlst = base + '\n' + inpt + '\n' + outp 
    if mask:
        mlst += '\n' + mask
    if para:
        mlst += '\n' + para
    if acts:
        mlst += '\n' + acts
    if mask:
        modl += '\n' + modl

    return  mlst
