#!/usr/bin/env bash
if [[ ! -e init.1.done ]]; then
    Rscript load_jun2018_data.R
    Rscript load_and_save_PD2182_data.R
    Rscript load_and_save_PD2182sequel_data.R
    Rscript load_jun2018_data_abcd_kl.R
    Rscirpt call_extreme_ipd.R
    Rscript call_extreme_ipd.PD2182.R
    Rscript call_modification.R
    date > init.1.done
fi
if [[ ! -e init.2.done ]]; then
    Rscript load_and_save_greer_data.R
    date > init.2.done
fi
