init_dir () {
    dir=$1
    fixed_position=$2
    motif=$3
    mkdir -p $dir
    ln -sf ../subdir.makefile $dir/Makefile
    ln -sf $fixed_position $dir/fixed_position.gff
    echo $motif > $dir/PATTERN
}

init_dir greer_pacbio_m6A_cov50 ../../../greer_pacbio_no_chunk/ipd_summary.m6A_cov50.gff A
init_dir greer_pacbio_orig_m6A_cov50 ../../../greer_pacbio_orig_no_chunk/ipd_summary.m6A_cov50.gff A
