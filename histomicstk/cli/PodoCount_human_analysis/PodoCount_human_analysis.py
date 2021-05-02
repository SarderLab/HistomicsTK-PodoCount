import os
from histomicstk.cli.utils import CLIArgumentParser

def main(args2):     
    
    cmd = "python3 ../BAS_folder_human/pod_quant_main_serv.py -A '{}' -B '{}' -D {} -P '{}' -L '{}' -M '{}'".format(args2.inputImageFile, args2.inputAnnotationFile1, args2.slider, args2.outputAnnotationFile1, args2.outputAnnotationFile2 , args2.csvFile)

    os.system(cmd)    

if __name__ == "__main__":
    main(CLIArgumentParser().parse_args())
