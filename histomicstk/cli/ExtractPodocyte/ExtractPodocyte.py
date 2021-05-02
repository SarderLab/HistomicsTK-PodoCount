import os
from histomicstk.cli.utils import CLIArgumentParser

def main(args2):     
    
    cmd = "python3 ../BAS_folder/pod_quantification_main.py -A '{}' -B '{}' -C '{}' -D {} -E {} -F {} -P '{}' -L '{}' -M '{}'".format(args2.inputImageFile, args2.inputAnnotationFile1,args2.inputAnnotationFile2,args2.slider, args2.ihc_gauss_sd, args2.num_sections, args2.outputAnnotationFile1 ,args2.outputAnnotationFile2 , args2.csvFile)
    os.system(cmd)    

if __name__ == "__main__":
    main(CLIArgumentParser().parse_args())
