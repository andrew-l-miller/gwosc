from pathlib import Path

def make_input_file(ifo,ffl_fname,sciseg_fname, channel_name,input_file_name,subsamp_factor='1',nsamps='4194304',interlace_ffts='2',window_type='5'):
    # ifo: H1, L1
    # ffl_fname:'H1_O3a.ffl'; list of frames with paths to them; create with make_ffl.py
    # sciseg_fname: 'H1_O3a_sciseg_for_sfdb.txt' contains seg number, start time, end time, duration; create with convert_sciseg_file 
    # channel_name: 'H1:GWOSC-4KHZ_R1_STRAIN'; channel in GWOSC to use
    # subsamp_fact: downsampling factor [1 = no downsampling]; for 4 KHz frames, no downsampling needed
    # nsamps: number of samples per FFT; def = 4194304, which means TFFT = 4194304/4096 = 1024 s
    # interlace_ffts: 2 for interlacing by 50%, 1 for no interlacing; 0 for ??
    # window_type: 0=no,1=Hann,2=Hamm,3=MAP, 4=Blackmann flatcos; 5=flat top,cosine edge. Sugg. 5) 
    flag_sfdb = '2'
    fact_evf = '2'
    verb_lvl = '1'
    max_num_FFTs_total = '1000000' ## maximum number of FFTs to do total; if large, does full data set (be careful for small TFFT)
    subsamp_fact_ar_spec = '128' ## factor by which AR specrum is downsammpled w.r.t. FFTs
    veto_freq = '100' ## frequency of veto
    subsamp_fact_veto = '-1' ##always keep this
    max_num_ffts_per_file = '100' ## max number of FFTs per SFDB

    lines = [
        str(ifo),                   # 1) detector
        str(flag_sfdb),             # 2) 2 or 3
        str(fact_evf),              # 3) "2" in template (keep as-is)
        str(ffl_fname),             # 4) .ffl file
        str(sciseg_fname),          # 5) science segments
        str(channel_name),          # 6) channel
        str(subsamp_factor),        # 7) subsampling factor
        str(verb_lvl),              # 8) extra 1 line in template (keep as-is)
        str(nsamps),                # 9) FFT length in samples
        str(interlace_ffts),        # 10) overlap/interlace
        str(max_num_FFTs_total),    # 11) max total FFTs
        str(subsamp_fact_ar_spec),  # 12) header spectrum subsampling
        str(window_type),           # 13) window type
        str(veto_freq),             # 14) veto frequency
        str(subsamp_fact_veto),     # 15) always -1
        str(max_num_ffts_per_file), # 16) FFTs per SFDB file
    ]

    outpath = Path(input_file_name)
    outpath.write_text("\n".join(lines) + "\n")

