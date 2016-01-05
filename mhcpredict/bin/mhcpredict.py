#!/usr/bin/env python

import mhcpredict.predict
import mhcpredict.util

if __name__ == "__main__":
    result = mhcpredict.predict.predictPeptides(
        "VIFRLMRTNFL", 
        alleles="HLA-DRB1*0101",
        methods=["NetMHCpan"],
        config=mhcpredict.util.DictConfig(dict(NetMHCpan=dict(
            executable="/Users/didionjp/Downloads/netMHCIIpan-3.1/netMHCIIpan",
            tempdir=".")))
    )
    result.to_csv("output.txt", sep="\t")
    
import mhcpredict.predict
import mhcpredict.util    
result = mhcpredict.predict.predictPeptides(
    "VIFRLMRTNFL",
    alleles="HLA-A01:01",
    methods=["NetMHCpan"],
    config=mhcpredict.util.DictConfig(dict(NetMHCpan=dict(
        executable="/home/ubuntu/software/netMHCpan-2.8/Linux_x86_64/bin/netMHCpan",
        tempdir="/home/ubuntu/software/Epitopes_from_TCRs/mhcpredict")))
)